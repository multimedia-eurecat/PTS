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

#define MAX_SPEAKERS 100

static float vbap_atorad = (2.0f * M_PI) / 360.0f ;

typedef struct t_speaker_3D
{
    int channel;
    double azi;
    double ele;
    float x;
    float y;
    float z;
} t_speaker_3D;

typedef struct t_speaker_3D_set
{
  int channels[3];
  float inv_matrix[9];
  struct t_speaker_3D_set *next;
} t_speaker_3D_set;

// void pts_3d(double *output, int samples, double *tables, int bits, int channels, double horizontal_f, double vertical_f, t_speaker_3D *speakers);
void pts_3d(int offset, double *tables, int channels, int page_size, double *gains);
// void pts_3d(int offset, double *tables, int channels, int page_size, double *gains, int *active_lspk);
void vogel_distribution(double *azi, double *ele, int channels);
double m_fmod(double dividend, double divisor);
void build_panning_tables(double* tables, int h_bits, int v_bits, int channels, t_speaker_3D *speakers);

t_speaker_3D_set* vbap_choose_ls_triplets(t_speaker_3D *speakers, int channels);
void vbap_angle_to_cart(double azi, double ele, float res[3]);
float vbap_vec_length(t_speaker_3D v1);
float vbap_dot_product(t_speaker_3D v1, t_speaker_3D v2);
float vbap_vec_angle(t_speaker_3D v1, t_speaker_3D v2);
void vbap_cross_prod(t_speaker_3D v1, t_speaker_3D v2, t_speaker_3D *res);
int vbap_lines_intersect(int i, int j, int k, int l, t_speaker_3D *speakers);
int vbap_any_ls_inside_triplet(int a, int b, int c, t_speaker_3D *speakers, int channels);
void vbap_calculate_3x3_matrixes(t_speaker_3D_set *speaker_set, t_speaker_3D *speakers, int channels);
void vbap(float source_cart[3], t_speaker_3D_set *speaker_set, float g[3], int ls[3]);

// --- main
int main(int argc, char **argv)
{
    // synthesis parameters
    int channels = 0; // will default to 4
    int trials = 10;
    int c_param_count = 0, f_param_count = 0, p_param_count = 0;
    double azi[MAX_SPEAKERS];
    double ele[MAX_SPEAKERS];
    double duration = 1.;
    double input_f = 100.;
    double horizontal_f = 0.;
    double vertical_f = 0.;
    double horizontal_p = 0.;
    double vertical_p = 0.;
    int h_bits = 0; // will be set to 513 later
    int v_bits = 0;
    char round_to_integer = 0;
    char *filename = NULL;
    int b_verbose = 0;

    // parse commandline input
    int opt, tmparg;
    while((opt = getopt(argc, argv, "hc:t:l:i:f:p:b:mo:v")) != -1)
    {
        switch(opt)
        {
            case 'h':
                printf("\nHere is how you use this magnificent piece of software:\
                    \n\n\t-c NUM: the number of channels (default: 4). Provide several arguments of -c and each value becomes the azimuth or elevation consecutively in channel order.\
                    \n\t-t NUM: the number of trials of which a median execution time is determined (default: 10)\
                    \n\t-l SECONDS: defines the sample length in seconds that is computed (default: 1.)\
                    \n\t-i FREQ: input sound frequency (default: 100.)\
                    \n\t-f FREQ: rotation freuqency. Provide it a second time for elevation! (default: 0.)\
                    \n\t-p PHASE: rotaton phase in angles. Useful if rotation frequency is 0. Provide it a second time for elevation! (default: 0.)\
                    \n\t-b BITS: the horizontal panning-table size in bits. Give a second -b parameter for the vertical bit count (both default: 513 horizontal and 257 vertical)\
                    \n\t-m: convert spherical loudspeaker positions to integers to match MaxMSP external behaviour.\
                    \n\t-o FILENAME: if you provide a filename, then the result for each channel with be written as a sequence of mono wav files. (default: none)\
                    \n\t-v: activates verbose output\
                    \n");
                return 0;
            case 'c':
                if(channels < MAX_SPEAKERS)
                {
                    if(c_param_count % 2 == 0)
                        azi[channels] = atof(optarg);
                    else
                    {
                        ele[channels] = atof(optarg);
                        channels++;
                    }
                    c_param_count++;
                } else
                    printf("Warning: max channel count (MAX_SPEAKERS) reached, skipping...");
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
                if(f_param_count > 0)
                    vertical_f = atof(optarg);
                else
                {
                    horizontal_f = atof(optarg);
                    f_param_count++;
                }
                break;
            case 'p':
                if(p_param_count > 0)
                    vertical_p = atof(optarg);
                else
                {
                    horizontal_p = atof(optarg);
                    p_param_count++;
                }
                break;
            case 'm':
                round_to_integer = 1;
                break;
            case 'b':
                tmparg = atoi(optarg);
                if(tmparg < 1)
                {
                    printf("Error table size must be bigger 0!");
                    return 0;
                }
                if(h_bits < 1)
                    h_bits = tmparg;
                else
                    v_bits = tmparg;
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

    // no argument supplied, default to 4 speakers
    if(c_param_count == 0)
    {
        azi[0] = 4; // this simulates user input as 4 channels
        c_param_count = 1;
    }

    // single argument passed means azi[0] is number of speakers!
    if(c_param_count == 1)
    {
        channels = (int) azi[0];
        vogel_distribution((double *)azi, (double *)ele, channels);
    }
    else
    {
        if(c_param_count % 2 != 0)
        {
            printf("Warning: missing last elevation value; defaulting to 0°!\n");
            ele[channels] = 0.;
            channels++;
        }
    }

    // round loudspeaker positions to integers
    if(round_to_integer)
    {
        for(int ch=0; ch<channels; ch++)
        {
            azi[ch] = round(azi[ch]);
            ele[ch] = round(ele[ch]);
        }
    }

    // output channels if verbose
    if(b_verbose)
    {
        printf("== Channel positions: ");
        for(int i=0; i<channels; i++)
        {
            if(round_to_integer)
                printf("%d %d ", (int)azi[i], (int)ele[i]);
            else
                printf("\n%0.2f %0.2f ", azi[i], ele[i]);
        }
        printf("\n\n");
    }

    // check buffer sizes
    if(h_bits < 1)
        h_bits = 513;
    if(v_bits < 1)
        v_bits = 257;

    // normalise phase offset
    double h_bits_m_1 = h_bits-1;
    double v_bits_m_1 = v_bits-1;
    horizontal_p = fmod(fmod(horizontal_p, 360.) + 360., 360.) / 360. * h_bits_m_1;
    vertical_p = fmod(fmod((vertical_p - 90.) * -1, 360.) + 360., 360.) / 180. * v_bits_m_1;

    // prepare speaker array
    t_speaker_3D speakers[channels];
    for(int ch=0; ch<channels; ch++)
    {
        speakers[ch].channel = ch;
        speakers[ch].azi = fmod(fmod(azi[ch], 360.) + 360., 360.) / 180. * M_PI;
        speakers[ch].ele = fmod(fmod((ele[ch] - 90.) * -1, 360.) + 360., 360.) / 180. * M_PI;
        speakers[ch].x = cos(speakers[ch].azi) * sin(speakers[ch].ele);
        speakers[ch].y = sin(speakers[ch].azi) * sin(speakers[ch].ele);
        speakers[ch].z = cos(speakers[ch].ele);
    }

    // // prepare output array
    int samples = (int) (SAMPLE_RATE * duration);
    double *output = malloc(channels * samples * sizeof(double));

    // // prepare panning tables per channel
    double *tables = malloc(channels * v_bits * h_bits * sizeof(double));
    build_panning_tables(tables, h_bits, v_bits, channels, (t_speaker_3D *)speakers);

    // calculate gains
    int offset, v_offset, h_offset, page_size = h_bits*v_bits;
    double step, input, gains[channels];
    double input_per_second = 2*M_PI * input_f;
    double h_bits_per_second = h_bits * horizontal_f;
    double v_bits_per_second = 2*v_bits_m_1 * vertical_f;
    clock_t start, end;
    double cpu_time_used, time_sum = 0;
    for(int j=0; j<trials; j++)
    {
        for(int i=0; i<samples; i++)
        {
            // time measurement
            start = clock();

            step = i * INV_SAMPLE_RATE;
            v_offset = fmod(vertical_p + v_bits_per_second * step, 2*v_bits_m_1);
            h_offset = 0;
            if(v_offset < 0)
                v_offset += 2*v_bits_m_1;
            if(v_offset > v_bits_m_1)
            {
                v_offset = 2*v_bits_m_1 - v_offset;
                h_offset = h_bits/2;
            }
            offset = (int) fmod(horizontal_p + h_offset + h_bits_per_second * step, h_bits) + v_offset * h_bits;

            pts_3d(offset, tables, channels, page_size, gains);
            // pts_3d(offset, tables, channels, page_size, gains, active_lspk);

            // stop clock
            end = clock();
            cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
            if(b_verbose) { printf("Time used: %0.10f\n", cpu_time_used); }
            time_sum += cpu_time_used;

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

    return 0;
}

void vogel_distribution(double *azi, double *ele, int channels)
{
    float golden_angle = M_PI * (3. - sqrt(5.));
    float z_max = 1-1./ channels;
    float z_step = 2. / channels;

    for(int i=0; i<channels; i++)
    {
        float t = golden_angle * i;
        float z = z_max - z_step * i;
        float s = sqrt(1 - z * z);
        float x = s * cos(t);
        float y = s * sin(t);

        // convert to spherical
        float r = sqrt(x*x + y*y + z*z);
        float phi = atan2(y, x);
        float theta = acos(z/r);

        azi[i] = phi * 180. / M_PI; // we need to simulate user input, thus angles...
        ele[i] = (theta * 180. / M_PI - 90.) * -1;
    }
}

void pts_3d(int offset, double *tables, int channels, int page_size, double *gains)
// void pts_3d(int offset, double *tables, int channels, int page_size, double *gains, int *active_lspk)
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
        gains[ch] = tables[ch * page_size + offset];
    }
}

// this seems to work more efficiently for our purposes than fmod alone
double m_fmod(double dividend, double divisor)
{
    double ratio = dividend / divisor;
    double decimal = fmod(ratio, 1.);
    return decimal * divisor;
}

// --- function to build many individual panning tables
void build_panning_tables(double* tables,
                        int h_bits,
                        int v_bits,
                        int channels,
                        t_speaker_3D *speakers)
{
    t_speaker_3D_set *speaker_set = vbap_choose_ls_triplets(speakers, channels);
    if(speaker_set == NULL)
    {
        printf("Error: Creating speaker set failed.");
        return;
    }
    vbap_calculate_3x3_matrixes(speaker_set, speakers, channels);

    int i, ch, table_index, ele_step, azi_step, ls[3];
    float source_cart[3], g[3];
    double azi, ele;
    for(i=0; ch < channels*v_bits*h_bits; ch++)
    {
        tables[i] = 0;
    }
    for(ele_step=0; ele_step<v_bits; ele_step++)
    {
        ele = (180. * ele_step/v_bits - 90.) * -1;
        for(azi_step=0; azi_step<h_bits; azi_step++)
        {
            azi = 360. * azi_step/h_bits;
            vbap_angle_to_cart(azi, -ele, source_cart);
            vbap(source_cart, speaker_set, g, ls);
            for(ch=0; ch<3; ch++)
            {
                table_index = ls[ch]*v_bits*h_bits + ele_step*h_bits + azi_step;
                tables[table_index] = g[ch];
            }
        }
    }

        // char filename[] = "testdoubleout";
        // for(int ch=0; ch<channels; ch++)
        // {
        //     char outname[strlen(filename) + 6];
        //     sprintf(outname, "%s%02d.wav", filename, ch);
        //     write_mono_wav_file_d(outname,
        //                           tables+ch*v_bits*h_bits,
        //                           v_bits*h_bits,
        //                           SAMPLE_RATE);
        // }


    t_speaker_3D_set *prev_set;
    while(speaker_set != NULL)
    {
        prev_set = speaker_set;
        speaker_set = speaker_set->next;
        free(prev_set);
    }
}




t_speaker_3D_set* vbap_choose_ls_triplets(t_speaker_3D *speakers, int channels)
{
    int i, j, k, l, table_size;
    int connections[MAX_SPEAKERS][MAX_SPEAKERS];
    float distance_table[((MAX_SPEAKERS * (MAX_SPEAKERS - 1)) / 2)];
    int distance_table_i[((MAX_SPEAKERS * (MAX_SPEAKERS - 1)) / 2)];
    int distance_table_j[((MAX_SPEAKERS * (MAX_SPEAKERS - 1)) / 2)];
    float distance;
    t_speaker_3D_set *speaker_set, *prev_speaker_set = NULL;

    if (channels == 0) {
        printf("ERROR: define-loudspeakers: Number of loudspeakers is zero\n");
        return NULL;
    }

    for(i=0; i<channels; i++)
    {
        for(j=i+1; j<channels; j++)
        {
            for(k=j+1; k<channels; k++)
            { // NOTE: removed check for tiny triangle
                connections[i][j]=1;
                connections[j][i]=1;
                connections[i][k]=1;
                connections[k][i]=1;
                connections[j][k]=1;
                connections[k][j]=1;

                t_speaker_3D_set *new_speaker_set = malloc(sizeof(t_speaker_3D_set));
                new_speaker_set->next = NULL;
                new_speaker_set->channels[0] = i;
                new_speaker_set->channels[1] = j;
                new_speaker_set->channels[2] = k;
                if(prev_speaker_set == NULL)
                    speaker_set = new_speaker_set;
                else
                    prev_speaker_set->next = new_speaker_set;
                prev_speaker_set = new_speaker_set;
            }
        }
    }

    // calculate distancies between all lss and sorting them
    table_size = ((channels * (channels - 1)) / 2);
    for(i=0; i<table_size; i++)
        distance_table[i] = 100000.0;

    for(i=0; i<channels; i++)
    {
        for(j=i+1; j<channels; j++)
        {
            if(connections[i][j] == 1)
            {
                distance = fabs(vbap_vec_angle(speakers[i], speakers[j]));

                k = 0;
                while(distance_table[k] < distance)
                    k++;

                for(l = (table_size - 1); l>k; l--)
                {
                    distance_table[l] = distance_table[l-1];
                    distance_table_i[l] = distance_table_i[l-1];
                    distance_table_j[l] = distance_table_j[l-1];
                }

                distance_table[k] = distance;
                distance_table_i[k] = i;
                distance_table_j[k] = j;
            }
            else
                table_size--;
        }
    }

    /* disconnecting connections which are crossing shorter ones,
       starting from shortest one and removing all that cross it,
       and proceeding to next shortest */
    for(i=0; i<(table_size); i++)
    {
        int fst_ls = distance_table_i[i];
        int sec_ls = distance_table_j[i];
        if(connections[fst_ls][sec_ls] == 1)
        {
            for(j=0; j<channels ; j++)
            {
                for(k=j+1; k<channels; k++)
                {
                    if( (j!=fst_ls) && (k != sec_ls) && (k!=fst_ls) && (j != sec_ls))
                    {
                        if(vbap_lines_intersect(fst_ls, sec_ls, j, k, speakers) == 1)
                        {
                            connections[j][k] = 0;
                            connections[k][j] = 0;
                        }
                    }
                }
            }
        }
    }

    /* remove triangles which had crossing sides
     with smaller triangles or include loudspeakers*/
    t_speaker_3D_set *tmp_ptr, *speaker_set_ptr = speaker_set;
    prev_speaker_set = NULL;
    while (speaker_set_ptr != NULL)
    {
        i = speaker_set_ptr->channels[0];
        j = speaker_set_ptr->channels[1];
        k = speaker_set_ptr->channels[2];
        if(connections[i][j] == 0 ||
           connections[i][k] == 0 ||
           connections[j][k] == 0 ||
           vbap_any_ls_inside_triplet(i, j, k, speakers, channels) == 1 )
        {
            if(prev_speaker_set != NULL)
            {
                prev_speaker_set->next = speaker_set_ptr->next;
                tmp_ptr = speaker_set_ptr;
                speaker_set_ptr = speaker_set_ptr->next;
                free(tmp_ptr);
            }
            else
            {
                speaker_set = speaker_set_ptr->next;
                tmp_ptr = speaker_set_ptr;
                speaker_set_ptr = speaker_set_ptr->next;
                free(tmp_ptr);
                // free(speaker_set_ptr);
                // speaker_set_ptr = speaker_set;
            }
        }
        else
        {
            prev_speaker_set = speaker_set_ptr;
            speaker_set_ptr = speaker_set_ptr->next;
        }
    }

    return speaker_set;
}
void vbap_angle_to_cart(double azi, double ele, float res[3])
{
    res[0] = cos(azi * vbap_atorad) * cos(ele * vbap_atorad);
    res[1] = sin(azi * vbap_atorad) * cos(ele * vbap_atorad);
    res[2] = sin(ele * vbap_atorad);
}
float vbap_vec_angle(t_speaker_3D v1, t_speaker_3D v2)
{
    float inner = vbap_dot_product(v1, v2) / (vbap_vec_length(v1) * vbap_vec_length(v2));
    if(inner > 1.0)
        inner= 1.0;
    if (inner < -1.0)
        inner = -1.0;
    return fabs(acos(inner));
}

float vbap_dot_product(t_speaker_3D v1, t_speaker_3D v2)
{
    return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

float vbap_vec_length(t_speaker_3D v1)
{
  return (sqrt(vbap_dot_product(v1, v1)));
}

void vbap_cross_prod(t_speaker_3D v1, t_speaker_3D v2, t_speaker_3D *res)
{
  res->x = (v1.y * v2.z ) - (v1.z * v2.y);
  res->y = (v1.z * v2.x ) - (v1.x * v2.z);
  res->z = (v1.x * v2.y ) - (v1.y * v2.x);

  float length = vbap_vec_length(*res);
  res->x /= length;
  res->y /= length;
  res->z /= length;
}

int vbap_lines_intersect(int i, int j, int k, int l, t_speaker_3D *speakers)
{
    t_speaker_3D v1;
    t_speaker_3D v2;
    t_speaker_3D v3, neg_v3;
    //float angle;
    float dist_ij, dist_kl, dist_iv3, dist_jv3, dist_inv3, dist_jnv3;
    float dist_kv3, dist_lv3, dist_knv3, dist_lnv3;

    vbap_cross_prod(speakers[i],speakers[j], &v1);
    vbap_cross_prod(speakers[k],speakers[l], &v2);
    vbap_cross_prod(v1, v2, &v3);

    neg_v3.x = 0.0 - v3.x;
    neg_v3.y = 0.0 - v3.y;
    neg_v3.z = 0.0 - v3.z;

    dist_ij = (vbap_vec_angle(speakers[i],speakers[j]));
    dist_kl = (vbap_vec_angle(speakers[k],speakers[l]));
    dist_iv3 = (vbap_vec_angle(speakers[i],v3));
    dist_jv3 = (vbap_vec_angle(v3,speakers[j]));
    dist_inv3 = (vbap_vec_angle(speakers[i],neg_v3));
    dist_jnv3 = (vbap_vec_angle(neg_v3,speakers[j]));
    dist_kv3 = (vbap_vec_angle(speakers[k],v3));
    dist_lv3 = (vbap_vec_angle(v3,speakers[l]));
    dist_knv3 = (vbap_vec_angle(speakers[k],neg_v3));
    dist_lnv3 = (vbap_vec_angle(neg_v3,speakers[l]));

    // if one of loudspeakers is close to crossing point, don't do anything
    if(fabsf(dist_iv3) <= 0.01 || fabsf(dist_jv3) <= 0.01 ||
       fabsf(dist_kv3) <= 0.01 || fabsf(dist_lv3) <= 0.01 ||
       fabsf(dist_inv3) <= 0.01 || fabsf(dist_jnv3) <= 0.01 ||
       fabsf(dist_knv3) <= 0.01 || fabsf(dist_lnv3) <= 0.01 )
        return 0;

    // if crossing point is on line between both loudspeakers return 1
    if (((fabsf(dist_ij - (dist_iv3 + dist_jv3)) <= 0.01 ) &&
         (fabsf(dist_kl - (dist_kv3 + dist_lv3))  <= 0.01)) ||
        ((fabsf(dist_ij - (dist_inv3 + dist_jnv3)) <= 0.01)  &&
         (fabsf(dist_kl - (dist_knv3 + dist_lnv3)) <= 0.01 )))
        return 1;
    else
        return 0;
}

/*--------------------------------------------------------------------------*/
// returns 1 if there is loudspeaker(s) inside given ls triplet
int vbap_any_ls_inside_triplet(int a, int b, int c, t_speaker_3D *speakers, int channels)
{
    float invdet;
    t_speaker_3D *lp1, *lp2, *lp3;
    float invmx[9];
    int i,j;
    float tmp;
    int any_ls_inside, this_inside;

    lp1 =  &(speakers[a]);
    lp2 =  &(speakers[b]);
    lp3 =  &(speakers[c]);

    /* matrix inversion */
    invdet = 1.0 / (  lp1->x * ((lp2->y * lp3->z) - (lp2->z * lp3->y))
                    - lp1->y * ((lp2->x * lp3->z) - (lp2->z * lp3->x))
                    + lp1->z * ((lp2->x * lp3->y) - (lp2->y * lp3->x)));
 
    invmx[0] = ((lp2->y * lp3->z) - (lp2->z * lp3->y)) * invdet;
    invmx[3] = ((lp1->y * lp3->z) - (lp1->z * lp3->y)) * -invdet;
    invmx[6] = ((lp1->y * lp2->z) - (lp1->z * lp2->y)) * invdet;
    invmx[1] = ((lp2->x * lp3->z) - (lp2->z * lp3->x)) * -invdet;
    invmx[4] = ((lp1->x * lp3->z) - (lp1->z * lp3->x)) * invdet;
    invmx[7] = ((lp1->x * lp2->z) - (lp1->z * lp2->x)) * -invdet;
    invmx[2] = ((lp2->x * lp3->y) - (lp2->y * lp3->x)) * invdet;
    invmx[5] = ((lp1->x * lp3->y) - (lp1->y * lp3->x)) * -invdet;
    invmx[8] = ((lp1->x * lp2->y) - (lp1->y * lp2->x)) * invdet;

    any_ls_inside = 0;
    for(i=0; i< channels; i++)
    {
        if (i != a && i!=b && i != c)
        {
            this_inside = 1;
            for(j=0; j< 3; j++)
            {
                tmp = speakers[i].x * invmx[0 + j*3];
                tmp += speakers[i].y * invmx[1 + j*3];
                tmp += speakers[i].z * invmx[2 + j*3];
                if(tmp < -0.001)
                    this_inside = 0;
            }
            if(this_inside == 1)
                any_ls_inside = 1;
        }
    }
    return any_ls_inside;
}

// Calculates the inverse matrices for 3D
void vbap_calculate_3x3_matrixes(t_speaker_3D_set *speaker_set, t_speaker_3D *speakers, int channels)
{ 
    float invdet;
    float *invmx;//, *at;
    t_speaker_3D_set *tr_ptr = speaker_set;
    int i;//, speaker_sets = 0;//, pointer, list_length=0;
    t_speaker_3D *lp1, *lp2, *lp3;

    if (tr_ptr == NULL)
    {
        printf("ERROR: define-loudspeakers: Not valid 3-D configuration\n");
        return;
    }
   
    // while(tr_ptr != NULL)
    // {
    //     speaker_sets++;
    //     tr_ptr = tr_ptr->next;
    // }
    // tr_ptr = speaker_set;
    // at = (float *) malloc((speaker_sets * 21 + 2) * sizeof(float));

    // at[0] = 3.0; // always 3rd dimension;
    // at[1] = (float) channels;
    // pointer = 2;

    while(tr_ptr != NULL)
    {
        lp1 =  &(speakers[tr_ptr->channels[0]]);
        lp2 =  &(speakers[tr_ptr->channels[1]]);
        lp3 =  &(speakers[tr_ptr->channels[2]]);

        /* matrix inversion */
        invmx = tr_ptr->inv_matrix;
        invdet = 1.0 / (  lp1->x * ((lp2->y * lp3->z) - (lp2->z * lp3->y))
                        - lp1->y * ((lp2->x * lp3->z) - (lp2->z * lp3->x))
                        + lp1->z * ((lp2->x * lp3->y) - (lp2->y * lp3->x)));

        invmx[0] = ((lp2->y * lp3->z) - (lp2->z * lp3->y)) * invdet;
        invmx[3] = ((lp1->y * lp3->z) - (lp1->z * lp3->y)) * -invdet;
        invmx[6] = ((lp1->y * lp2->z) - (lp1->z * lp2->y)) * invdet;
        invmx[1] = ((lp2->x * lp3->z) - (lp2->z * lp3->x)) * -invdet;
        invmx[4] = ((lp1->x * lp3->z) - (lp1->z * lp3->x)) * invdet;
        invmx[7] = ((lp1->x * lp2->z) - (lp1->z * lp2->x)) * -invdet;
        invmx[2] = ((lp2->x * lp3->y) - (lp2->y * lp3->x)) * invdet;
        invmx[5] = ((lp1->x * lp3->y) - (lp1->y * lp3->x)) * -invdet;
        invmx[8] = ((lp1->x * lp2->y) - (lp1->y * lp2->x)) * invdet;

        // for(i=0;i<3;i++)
        // {
        //     at[pointer] = (float) (tr_ptr->channels[i]+1);
        //     pointer++;
        // }

        // for(i=0;i<9;i++)
        // {
        //     at[pointer] = invmx[i];
        //     pointer++;
        // }

        // at[pointer] = lp1->x; pointer++;
        // at[pointer] = lp2->x; pointer++;
        // at[pointer] = lp3->x; pointer++;
        // at[pointer] = lp1->y; pointer++;
        // at[pointer] = lp2->y; pointer++;
        // at[pointer] = lp3->y; pointer++;
        // at[pointer] = lp1->z; pointer++;
        // at[pointer] = lp2->z; pointer++;
        // at[pointer] = lp3->z; pointer++;
    
        tr_ptr = tr_ptr->next;
    }

    // free(at);
}


void vbap(float source_cart[3], t_speaker_3D_set *speaker_set, float g[3], int ls[3])
{
    /* calculates gain factors using loudspeaker setup and given direction */
    float power;
    int i,j,k, gains_modified;
    float small_g;
    float big_sm_g, gtmp[3];
    long winner_set = 0;
    // float cartdir[3];
    float new_cartdir[3];
    float new_angle_dir[3];
    // long dim = 3;
    long neg_g_am, best_neg_g_am;
    t_speaker_3D_set *current_speaker_set, *winner_speaker_set;

    big_sm_g = -100000.0;   // initial value for largest minimum gain value
    best_neg_g_am=3;        // how many negative values in this set

    // for(i=0; i<speaker_sets; i++)
    // {
    //     small_g = 10000000.0;
    //     neg_g_am = 3;
    //     for(j=0;j<3;j++)
    //     {
    //         gtmp[j]=0.0;
    //         for(k=0;k<3;k++)
    //             gtmp[j] += source_cart[k] * x_set_inv_matx[i][k+j*3];
    //         if(gtmp[j] < small_g)
    //             small_g = gtmp[j];
    //         if(gtmp[j]>= -0.01)
    //             neg_g_am--;
    //     }

    //     if(small_g > big_sm_g && neg_g_am <= best_neg_g_am)
    //     {
    //         big_sm_g = small_g;
    //         best_neg_g_am = neg_g_am;
    //         winner_set=i;
    //         g[0]=gtmp[0]; g[1]=gtmp[1];
    //         ls[0]= x_lsset[i][0]; ls[1]= x_lsset[i][1];
    //         if(3==3)
    //         {
    //             g[2]=gtmp[2];
    //             ls[2]= x_lsset[i][2];
    //         }
    //         else
    //         {
    //             g[2]=0.0;
    //             ls[2]=0;
    //         }
    //     }
    // }

    current_speaker_set = speaker_set;
    while(current_speaker_set != NULL)
    {
        small_g = 10000000.0;
        neg_g_am = 3;

        for(j=0; j<3; j++)
        {
            gtmp[j]=0.0;
            for(k=0; k<3; k++)
                gtmp[j] += source_cart[k] * current_speaker_set->inv_matrix[k+j*3];

            if(gtmp[j] < small_g)
                small_g = gtmp[j];
            if(gtmp[j] >= -0.01)
                neg_g_am--;
        }

        if(small_g > big_sm_g && neg_g_am <= best_neg_g_am)
        {
            big_sm_g = small_g;
            best_neg_g_am = neg_g_am;
            winner_speaker_set = current_speaker_set;

            g[0] = gtmp[0];
            g[1] = gtmp[1];
            g[2] = gtmp[2];

            ls[0]= current_speaker_set->channels[0];
            ls[1]= current_speaker_set->channels[1];
            ls[2]= current_speaker_set->channels[2];
        }

        current_speaker_set = current_speaker_set->next;
    }
 
    // If chosen set produced a negative value, make it zero and
    // calculate direction that corresponds to these new
    // gain values. This happens when the virtual source is outside of
    // all loudspeaker sets.
    gains_modified=0;
    for(i=0; i<3; i++)
        if(g[i]<-0.01)
        {
            g[i]=0.0001;
            gains_modified=1;
        }

    if(gains_modified==1)
    {
        new_cartdir[0] =  winner_speaker_set->inv_matrix[0] * g[0]
                        + winner_speaker_set->inv_matrix[1] * g[1]
                        + winner_speaker_set->inv_matrix[2] * g[2];
        new_cartdir[1] =  winner_speaker_set->inv_matrix[3] * g[0]
                        + winner_speaker_set->inv_matrix[4] * g[1]
                        + winner_speaker_set->inv_matrix[5] * g[2];
        new_cartdir[2] =  winner_speaker_set->inv_matrix[6] * g[0]
                        + winner_speaker_set->inv_matrix[7] * g[1]
                        + winner_speaker_set->inv_matrix[8] * g[2];
        // cart_to_angle(new_cartdir,new_angle_dir);
        // x_azi = (long) (new_angle_dir[0] + 0.5);
        // x_ele = (long) (new_angle_dir[1] + 0.5);
    }

    power = sqrt(g[0]*g[0] + g[1]*g[1] + g[2]*g[2]);
    g[0] /= power;
    g[1] /= power;
    g[2] /= power;
}

