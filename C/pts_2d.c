#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include "simple_wavwrite.h"

// --- audio file constants
#define SAMPLE_RATE 48000
#define INV_SAMPLE_RATE 0.00002083333333

#define MAX_SPEAKERS 512

typedef struct
{
    int channel;
    double angle;
    int range_min;
    int range_max;
} t_speaker_2D;

// void pts_2d(int offset, double *tables, int channels, int bits, double *gains);
void pts_2d(int offset, double *tables, int channels, int bits, double *gains, int *active_lspk);

void build_panning_tables(double* tables, int bits, int channels, t_speaker_2D *speakers, int *active_lspk);
int spk_cmp(const void *a, const void *b);
double vbap_2d_equation(double phi, double phi_0);

// --- main
int main(int argc, char **argv)
{
    // synthesis parameters
    int channels = 0;
    double spk_angles[MAX_SPEAKERS]; // max number of speakers
    double duration = 1.;
    double input_f = 100.;
    double horizontal_f = 0.;
    double horizontal_p = 0.;
    int bits = 513;
    int trials = 10;
    char *filename = NULL;
    int b_verbose = 0;

    // parse commandline input
    int opt;
    while((opt = getopt(argc, argv, "hc:t:l:i:f:p:b:d:o:v")) != -1)
    {
        switch(opt)
        {
            case 'h':
                printf("\nHere is how you use this magnificent piece of software:\
                    \n\n\t-c NUM: the number of channels in a ring (default: 3). Provide several arguments of -c and each value becomes the angle in channel order.\
                    \n\t-t NUM: the number of trials of which a median execution time is determined (default: 10)\
                    \n\t-l SECONDS: defines the sample length in seconds that is computed (default: 1.)\
                    \n\t-i FREQ: input sound frequency (default: 100.)\
                    \n\t-f FREQ: rotation freuqency (default: 0.)\
                    \n\t-p PHASE: rotaton phase in angles. Useful if rotation frequency is 0 (default: 0.)\
                    \n\t-b BITS: the panning-table size in bits (default: 513)\
                    \n\t-o FILENAME: if you provide a filename, then the result for each channel with be written as a sequence of mono wav files. (default: none)\
                    \n\t-v: activates verbose output\
                    \n");
                return 0;
            case 'c':
                if(channels < MAX_SPEAKERS) {
                    spk_angles[channels] = atof(optarg);
                    channels += 1;
                } else {
                    printf("Warning: max channel count (MAX_SPEAKERS) reached, skipping...");
                }
                break;
            case 't':
                trials = atoi(optarg);
                break;
            case 'l':
                duration = atof(optarg);
                break;
            case 'i':
                input_f = atof(optarg);
                break;
            case 'f':
                horizontal_f = atof(optarg);
                break;
            case 'p':
                horizontal_p = atof(optarg);
                break;
            case 'b':
                bits = atoi(optarg);
                break;
            case 'o':
                filename = optarg;
                break;
            case 'v':
                b_verbose = 1;
                break;
            case ':':
                printf("Argument required!\n");
                break;
            case '?':
                printf("Unknown option: %c\n", optopt);
                break;
        }
    }

    if(channels == 0) // default to 3 channels, but simulate user input first...
    {
        channels = 1;
        spk_angles[0] = 3;
    }
   
    if(channels == 1) // single argument passed means number of speakers, not speaker angle
    {
        channels = (int) spk_angles[0];
        for(int ch=0; ch<channels; ch++)
        {
            spk_angles[ch] = 360. / channels * ch; // we simulate user input, thus we use degrees
        }
    }

    // prepare speaker array
    t_speaker_2D speakers[channels];
    for(int ch=0; ch<channels; ch++)
    {
        speakers[ch].channel = ch;
        speakers[ch].angle = fmod(fmod(spk_angles[ch], 360.) + 360., 360.) / 180. * M_PI;
    }

    // prepare output array
    int samples = (int) (SAMPLE_RATE * duration);
    double *output = malloc(channels * samples * sizeof(double));

    // prepare panning tables per channel
    double *tables = malloc(channels * bits * sizeof(double));
    int *active_lspk = malloc(2 * bits * sizeof(int));
    build_panning_tables(tables, bits, channels, (t_speaker_2D *)speakers, active_lspk);

    horizontal_p = fmod(fmod(horizontal_p, 360.) + 360., 360.) / 360. * bits;

    // calculate gains
    int offset;
    double step, input, gains[channels];
    double input_per_second = 2*M_PI * input_f;
    double bits_per_second = bits * horizontal_f;
    clock_t start, end;
    double cpu_time_used, time_sum = 0;
    for(int j=0; j<trials; j++)
    {
        for(int i=0; i<samples; i++)
        {
            // time measurement
            start = clock();

            step = i * INV_SAMPLE_RATE;
            offset = (int) fmod(horizontal_p + bits_per_second * step, bits); // NOTE: option (1) using modulo
            // offset = atan(tan(M_PI_2 - M_PI * step)) * M_1_PI + 0.5; // NOTE: option (2) using sawtooth wave

            // pts_2d(offset, tables, channels, bits, gains);
            pts_2d(offset, tables, channels, bits, gains, active_lspk);

            // stop clock
            end = clock();
            cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
            if(b_verbose) printf("Time used: %0.10f\n", cpu_time_used);
            time_sum += cpu_time_used; // sum for average

            // output...
            input = cos(input_per_second * step * 2*M_PI);
            for(int ch=0; ch<channels; ch++)
            {
                output[ch * samples + i] = input * gains[ch];
            }
        }
    }
    printf("== Average time: %0.10f\n", time_sum/(trials*samples));

    // output results, if necessary
    if(filename != NULL) {
        for(int ch=0; ch<channels; ch++)
        {
            char outname[strlen(filename) + 6];
            sprintf(outname, "%s%02d.wav", filename, ch);
            write_mono_wav_file_d(outname,
                                  output+ch*samples,
                                  samples,
                                  SAMPLE_RATE);
        }
    }

    free(output); // output array
    free(tables); // panning tables
    free(active_lspk); // active loudspeaker pair

    return 0;
}



// void pts_2d(int offset, double *tables, int channels, int bits, double *gains)
void pts_2d(int offset, double *tables, int channels, int bits, double *gains, int *active_lspk)
{
    // // constant complexity version using active_lspk
    // // NOTE: only works if it's known that only two loudspeakers are active at any time!
    // int index = 2*offset;
    // int ch = active_lspk[index];
    // gains[ch] = tables[ch * bits + offset];
    // ch = active_lspk[index+1];
    // gains[ch] = tables[ch * bits + offset];

    // linear complexity version ignoring active_lspk
    // NOTE: needed, if panning function is more complex and more than 2 loudspeakers could sound at some given position
    for(int ch=0; ch<channels; ++ch)
    {
        gains[ch] = tables[ch * bits + offset];
    }
}



// --- function to build many individual panning tables
// void build_panning_tables(double* tables, int bits, int channels, t_speaker_2D *speakers)
void build_panning_tables(double* tables, int bits, int channels, t_speaker_2D *speakers, int *active_lspk)
{
    int neighbors[channels][2];
    int i, ch;

    qsort(speakers, channels, sizeof(speakers[0]), spk_cmp);

    // NOTE: used to denote that an active speaker for this sample has already been found,
    // thus we can increase the index
    int active_spk_offset[bits];
    for(i=0; i<bits; i++)
        active_spk_offset[i] = 0;

    // calculate phi_0 in each direction and populate buffer per speaker
    int ch_left, ch_right, table_offset;
    double phi_0_left, phi_0_right, step, step_centered, x, y;
    for(ch=0; ch<channels; ch++)
    {
        ch_left = (ch+channels-1) % channels;
        ch_right = (ch+1) % channels;

        // phi_0 to the left side of the buffer for this speaker
        phi_0_left = (speakers[ch].angle - speakers[ch_left].angle)/2;
        if(phi_0_left<0)
            phi_0_left += M_PI;

        // phi_0 to the right side of the buffer for this speaker
        phi_0_right = (speakers[ch_right].angle - speakers[ch].angle)/2;
        if(phi_0_right<0)
            phi_0_right += M_PI;

        // populate table
        for(i=0; i<bits; i++)
        {
            step = 2*M_PI/bits * i;
            step_centered = fmod(step - speakers[ch].angle + 3*M_PI, 2*M_PI) - M_PI;
            table_offset = speakers[ch].channel * bits + i;
            if(step_centered == 0)
                tables[table_offset] = 1;
            else if(step_centered < 0 && step_centered > phi_0_left*-2)
                tables[table_offset] = vbap_2d_equation(step_centered + phi_0_left, phi_0_left);
            else if(step_centered > 0 && step_centered < phi_0_right*2)
                tables[table_offset] = vbap_2d_equation((step_centered - phi_0_right)*-1, phi_0_right);
            else
                tables[table_offset] = 0;

            // NOTE: watch out, only works with correct array size of active_lspk (ex. 2 in this case)!
            if(tables[table_offset] > 0)
            {
                int index = 2*i + active_spk_offset[i];
                active_lspk[index] = ch;
                active_spk_offset[i]++;
            }
        }
    }
}

int spk_cmp(const void *a, const void *b)
{
    t_speaker_2D *a1 = (t_speaker_2D *)a;
    t_speaker_2D *a2 = (t_speaker_2D *)b;
    if ((*a1).angle > (*a2).angle)
        return 1;
    else if ((*a1).angle < (*a2).angle)
        return -1;
    else
        return 0;
}

double vbap_2d_equation(double phi, double phi_0)
{
    return sqrt(( pow(tan(phi_0),2) + 2*tan(phi_0)*tan(phi) + pow(tan(phi),2) ) / (2*( pow(tan(phi_0),2) + pow(tan(phi),2) )));
}

