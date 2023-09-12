#define SUBCHUNK1SIZE   (16) // For PCM, see http://soundfile.sapp.org/doc/WaveFormat/ for reference
#define AUDIO_FORMAT    (1) // For PCM



// -----------------------------------------
// WAV header
// -----------------------------------------

// --- The header of a wav file modified from: https://ccrma.stanford.edu/courses/422/projects/WaveFormat/
typedef struct wavfile_header_s
{
    char    ChunkID[4];     // 4
    int32_t ChunkSize;      // 4
    char    Format[4];      // 4

    char    Subchunk1ID[4]; // 4
    int32_t Subchunk1Size;  // 4
    int16_t AudioFormat;    // 2
    int16_t NumChannels;    // 2
    int32_t SampleRate;     // 4
    int32_t ByteRate;       // 4
    int16_t BlockAlign;     // 2
    int16_t BitsPerSample;  // 2

    char    Subchunk2ID[4];
    int32_t Subchunk2Size;
} wavfile_header_t;


// --- Helper method to write the WAV header
int write_PCM_header(FILE*   file_p,
                     int32_t sample_rate,
                     int32_t bit_depth,
                     int32_t num_channels,
                     int32_t frame_count)
{
    int ret;

    wavfile_header_t wav_header;
    int32_t subchunk2_size;
    int32_t chunk_size;
    int32_t byte_rate;
    int32_t block_align;

    size_t write_count;

    byte_rate = sample_rate * num_channels * bit_depth/8;
    block_align = num_channels * bit_depth/8;
    subchunk2_size  = frame_count * num_channels * bit_depth/8;
    chunk_size = 4 + (8 + SUBCHUNK1SIZE) + (8 + subchunk2_size);

    wav_header.ChunkID[0] = 'R';
    wav_header.ChunkID[1] = 'I';
    wav_header.ChunkID[2] = 'F';
    wav_header.ChunkID[3] = 'F';

    wav_header.ChunkSize = chunk_size;

    wav_header.Format[0] = 'W';
    wav_header.Format[1] = 'A';
    wav_header.Format[2] = 'V';
    wav_header.Format[3] = 'E';

    wav_header.Subchunk1ID[0] = 'f';
    wav_header.Subchunk1ID[1] = 'm';
    wav_header.Subchunk1ID[2] = 't';
    wav_header.Subchunk1ID[3] = ' ';

    wav_header.Subchunk1Size = SUBCHUNK1SIZE;
    wav_header.AudioFormat = AUDIO_FORMAT;
    wav_header.NumChannels = num_channels;
    wav_header.SampleRate = sample_rate;
    wav_header.ByteRate = byte_rate;
    wav_header.BlockAlign = block_align;
    wav_header.BitsPerSample = bit_depth;

    wav_header.Subchunk2ID[0] = 'd';
    wav_header.Subchunk2ID[1] = 'a';
    wav_header.Subchunk2ID[2] = 't';
    wav_header.Subchunk2ID[3] = 'a';
    wav_header.Subchunk2Size = subchunk2_size;

    write_count = fwrite(&wav_header, 
                         sizeof(wavfile_header_t),
                         1,
                         file_p);

    ret = (1 != write_count)? -1 : 0;

    return ret;
}



// -----------------------------------------
// WAV Data
// -----------------------------------------

// --- Data structure to hold a single frame (sample) with one channel
typedef struct PCM16_mono_s
{
    int16_t channel;
} PCM16_mono_t;

// --- Helper method to allocate the necessary space in memory for a mono signal
PCM16_mono_t *allocate_PCM16_mono_buffer(int32_t FrameCount)
{
    return (PCM16_mono_t *)malloc(sizeof(PCM16_mono_t) * FrameCount);
}

// --- Data structure to hold a single frame with two channels (two samples, on for each channel)
typedef struct PCM16_stereo_s
{
    int16_t left;
    int16_t right;
} PCM16_stereo_t;

// --- Helper method to allocate the necessary space in memory for a stereo signal
PCM16_stereo_t *allocate_PCM16_stereo_buffer(int32_t FrameCount)
{
    return (PCM16_stereo_t *)malloc(sizeof(PCM16_stereo_t) * FrameCount);
}

// --- Convert any buffer of doubles into a PCM16_mono_t (16 bit) buffer, in order to write it to file later
int convert_double_to_int16_buffer(double *buffer,
                                   int samples,
                                   PCM16_mono_t *buffer_pcm)
{
    for(int k = 0; k<samples; ++k)
    {
        buffer_pcm[k].channel = (int16_t) (double)SHRT_MAX * *(buffer+k);
    }

    return 0;
}

// --- Returns the number of audio frames sucessfully written
size_t  write_PCM_mono_data(FILE* file_p,
                            int32_t FrameCount,
                            PCM16_mono_t *buffer_p)
{
    size_t ret;

    ret = fwrite(buffer_p, 
                 sizeof(PCM16_mono_t),
                 FrameCount,
                 file_p);

    return ret;
}

// --- Returns the number of audio frames sucessfully written
size_t  write_PCM_stereo_data(FILE* file_p,
                              int32_t FrameCount,
                              PCM16_stereo_t *buffer_p)
{
    size_t ret;

    ret = fwrite(buffer_p, 
                 sizeof(PCM16_stereo_t),
                 FrameCount,
                 file_p);

    return ret;
}



// -----------------------------------------
// Write file helper function
// -----------------------------------------

int write_mono_wav_file_d(char* filepath, double* buffer, int FrameCount, int SampleRate)
{
    int ret;

    FILE *file_p;
    size_t written;

    // PCM16_stereo_t *buffer_p = NULL;
    PCM16_mono_t *buffer_p = NULL;

    /*Open the wav file*/
    file_p = fopen(filepath, "w");
    if(NULL == file_p)
    {
        perror("fopen failed in main");
        ret = -1;
        goto error0;
    }

    // Allocate the data buffer
    buffer_p = allocate_PCM16_mono_buffer(FrameCount);
    if(NULL == buffer_p)
    {
        perror("fopen failed in main");
        ret = -1;
        goto error1;        
    }

    // fill PCM 16 bit buffer
    ret = convert_double_to_int16_buffer(buffer,
                                         FrameCount,
                                         buffer_p);
    if(ret < 0)
    {
        fprintf(stderr, "generate_dual_sawtooth failed in main\n");
        ret = -1;
        goto error2;
    }

    // Write the wav file header
    ret = write_PCM_header(file_p,
                           SampleRate,
                           16, // bit depth
                           1, // channels
                           FrameCount);
    if(ret < 0)
    {
        perror("write_PCM16_stereo_header failed in main");
        ret = -1;
        goto error2;
    }

    // Write the data out to file
    written = write_PCM_mono_data(file_p,
                                  FrameCount,
                                  buffer_p);
    if(written < FrameCount)
    {
        perror("write_PCM16wav_data failed in main");
        ret = -1;
        goto error2;
    }

error2:
    free(buffer_p);
error1:
    fclose(file_p);
error0:
    return ret; 
}



// -----------------------------------------
// Some simple signal generator examples
// -----------------------------------------

// Generate two saw-tooth signals at two frequencies and amplitudes
int generate_sawtooth( double frequency,
                            double amplitude,
                            int32_t SampleRate,
                            int32_t FrameCount,
                            PCM16_mono_t *buffer_p)
{
    int ret = 0;
    double SampleRate_d = (double)SampleRate;
    double SamplePeriod = 1.0/SampleRate_d;

    double Period;
    double phase;
    double Slope;

    int32_t k;

    // Check for the violation of the Nyquist limit
    if( (frequency*2 >= SampleRate_d) )
    {
        ret = -1;
        goto error0;
    }

    // Compute the period
    Period = 1.0/frequency;

    // Compute the slope
    Slope = (amplitude * (double)SHRT_MAX)/Period;

    for(k = 0, phase = 0.0; 
        k < FrameCount; 
        ++k)
    {
        phase += SamplePeriod;
        phase = (phase > Period)? (phase - Period) : phase;

        buffer_p[k].channel = (int16_t)(phase * Slope);
    }

error0:
    return ret;
}

// Generate two saw-tooth signals at two frequencies and amplitudes
int generate_dual_sawtooth( double frequency1,
                            double amplitude1,
                            double frequency2,
                            double amplitude2,
                            int32_t SampleRate,
                            int32_t FrameCount,
                            PCM16_stereo_t  *buffer_p)
{
    int ret = 0;
    double SampleRate_d = (double)SampleRate;
    double SamplePeriod = 1.0/SampleRate_d;

    double Period1, Period2;
    double phase1, phase2;
    double Slope1, Slope2;

    int32_t k;

    // Check for the violation of the Nyquist limit
    if( (frequency1*2 >= SampleRate_d) || (frequency2*2 >= SampleRate_d) )
    {
        ret = -1;
        goto error0;
    }

    // Compute the period
    Period1 = 1.0/frequency1;
    Period2 = 1.0/frequency2;

    // Compute the slope
    Slope1  = (amplitude1 * (double)SHRT_MAX)/Period1;
    Slope2  = (amplitude2 * (double)SHRT_MAX)/Period2;

    for(k = 0, phase1 = 0.0, phase2 = 0.0; 
        k < FrameCount; 
        ++k)
    {
        phase1 += SamplePeriod;
        phase1 = (phase1 > Period1)? (phase1 - Period1) : phase1;

        phase2 += SamplePeriod;
        phase2 = (phase2 > Period2)? (phase2 - Period2) : phase2;

        buffer_p[k].left    = (int16_t)(phase1 * Slope1);
        buffer_p[k].right   = (int16_t)(phase2 * Slope2);
    }

error0:
    return ret;
}