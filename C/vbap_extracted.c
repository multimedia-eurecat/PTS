/*
 * This code is taken and modified from the original Pure Data or Max MSP code written in C by Ville Pulkki
 * updated by Nathan Wolek for 64bit MaxMSP from https://github.com/nwolek/vbap/
 * The intention is to only remove those parts specific to the Pure Data or Max MSP framework
 * which would inhibit execution of this code as a standalone script,
 * but the code added allows execution as a script and in a repetitive manner over x trials to determine a median
 * using a predefined loudspeaker distribution
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include "simple_wavwrite.h"



/*--------------------------------------------------------------------------*/
// preprocessor constants (or not...)
#define SAMPLE_RATE 48000
#define INV_SAMPLE_RATE 0.00002083333333

#define MAX_LS_AMOUNT 100 // maximum amount of loudspeakers, can be increased, but see comments next to MAX_LS_SETS above
#define MAX_LS_SETS 1000 // maximum number of loudspeaker sets (triplets or pairs) allowed - This can crash when too many speakers are defined

#define MIN_VOL_P_SIDE_LGTH 0.01

static float rad2ang = 360.0 / ( 2.0f * M_PI );
static float atorad = (2.0f * M_PI) / 360.0f ;



/*--------------------------------------------------------------------------*/
// A struct for a loudspeaker instance (taken from VBAP object)
typedef struct 
{                   // distance value is 1.0 == unit vectors
  float x;  // cartesian coordinates
  float y;
  float z;
  float azi;  // polar coordinates
  float ele;
  int channel_nbr;  // which speaker channel number 
} t_ls;

/* A struct for all loudspeaker sets */
typedef struct t_ls_set 
{
  int ls_nos[3];  // channel numbers
  float inv_mx[9]; // inverse 3x3 or 2x2 matrix
  struct t_ls_set *next;  // next set (triplet or pair)
} t_ls_set;



/*--------------------------------------------------------------------------*/
// Function defs
void vogel_distribution(double *azi, double *ele, int channels);

// VBAP functions extracted from external (no object *x used anywhere)
void vbap_set_azimuth(double n);
void vbap_set_elevation(double n);
void angle_to_cart(double azi, double ele, float res[3]);
void cart_to_angle(float cvec[3], float avec[3]);
void vect_cross_prod(float v1[3], float v2[3], float v3[3]);
void ls_cross_prod(t_ls v1,t_ls v2, t_ls *res);
void ls_angles_to_cart(t_ls *ls);
float vec_angle(t_ls v1, t_ls v2);
float vec_length(t_ls v1);
float vec_prod(t_ls v1, t_ls v2);
int lines_intersect(int i,int j,int k,int l,t_ls  lss[MAX_LS_AMOUNT]);

void vbap_bang(float *final_gs, int x_ls_amount, int x_dimension);
void vbap(float g[3], long ls[3]);

void initContent_ls_directions(int ac, float *av);

void choose_ls_tuplets();
void sort_2D_lss(t_ls lss[MAX_LS_AMOUNT], int sorted_lss[MAX_LS_AMOUNT], int ls_amount);
int calc_2D_inv_tmatrix(float azi1,float azi2, float inv_mat[4],float mat[4]);

void choose_ls_triplets();
float vol_p_side_lgth(int i, int j,int k, t_ls lss[MAX_LS_AMOUNT]);
void add_ldsp_triplet(int i, int j, int k);
int any_ls_inside_triplet(int a, int b, int c,t_ls lss[MAX_LS_AMOUNT],int ls_amount);
void calculate_3x3_matrixes();

void vbap_matrix(int ac, float *av);



/*--------------------------------------------------------------------------*/
// VBAP OBJECT, Max port to script format from original struct
double x_azi = 0;                       // panning direction azimuth
double x_ele = 0;                       // panning direction elevation               
float x_set_inv_matx[MAX_LS_SETS][9];   // inverse matrice for each loudspeaker set
float x_set_matx[MAX_LS_SETS][9];       // matrice for each loudspeaker set
long x_lsset[MAX_LS_SETS][3];           // channel numbers of loudspeakers in each LS set 
long x_lsset_available;                 // have loudspeaker sets been defined with define_loudspeakers
long x_lsset_amount=0;                  // amount of loudspeaker sets
long x_ls_amount;                       // amount of loudspeakers
long x_dimension=2;                     // 2 or 3
double x_spread;                        // speading amount of virtual source (0-100)
double x_gain;                          // general gain control (0-2)
float x_spread_base[3];                 // used to create uniform spreading

// define_loudspeaker data
long x_ls_read;                         // 1 if loudspeaker directions have been read
long x_triplets_specified=0;            // 1 if loudspeaker triplets have been chosen
t_ls x_ls[MAX_LS_AMOUNT];               // loudspeakers
t_ls_set *x_ls_set;                     // loudspeaker sets
long x_def_ls_amount;                   // number of loudspeakers
long x_def_ls_dimension;                // 2 (horizontal arrays) or 3 (3d setups)



/*--------------------------------------------------------------------------*/
// MAIN
int main(int argc, char **argv)
{
    // synthesis parameters
    int channels = 0; // will default to 3 or 4
    int trials = 1;
    int c_param_count = 0, f_param_count = 0, p_param_count = 0;
    double azi[MAX_LS_AMOUNT];
    double ele[MAX_LS_AMOUNT];
    double duration = 1.;
    double input_f = 100.;
    double horizontal_f = 0.;
    double vertical_f = 0.;
    double horizontal_p = 0.;
    double vertical_p = 0.;
    int dimension = 2;
    char round_to_integer = 0;
    char *filename = NULL;
    char b_verbose = 0;

    // parse commandline input
    int opt;
    int next_f = 0;
    while((opt = getopt(argc, argv, "hc:t:l:i:f:p:d:mo:v")) != -1) 
    { 
        switch(opt) 
        { 
            case 'h':
                printf("\nHere is how you use this magnificent piece of software:\
                    \n\n\t-c NUM: the number of channels in the array (default: 3)\
                    \n\t-t NUM: the number of trials of which a median execution time is determined (default: 1)\
                    \n\t-l SECONDS: defines the sample length in seconds (default: 1.)\
                    \n\t-i FREQ: input sound freuqency (default: 100.).\
                    \n\t-f FREQ: rotation freuqency (default: 1.). Supply -f a second time to define vertical rotation frequency (with DIM=3)\
                    \n\t-p PHASE: rotaton phase in angles. Useful if rotation frequency is 0 (default: 0.)\
                    \n\t-d DIM: dimension, can be either 2 or 3 (default: 2)\
                    \n\t-m: convert spherical loudspeaker positions to integers to match MaxMSP external behaviour.\
                    \n\t-o FILENAME: if you provide a filename, then the output of the first channel will be written as a WAV file to that name. (default: none)\
                    \n\t-v: activates verbose output\
                    \n");
                return 0;
            case 'c':
                if(channels < MAX_LS_AMOUNT)
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
            case 'd':
                dimension = atoi(optarg);
                if(dimension < 2) {
                    printf("ERROR: -d dimension cannot be smaller than 2!\n");
                    return 1;
                }
                if(dimension > 3) {
                    printf("ERROR: -d dimension cannot be greater than 3!\n");
                    return 1;
                }
                break;
            case 'm':
                round_to_integer = 1;
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

    int samples = (int) (SAMPLE_RATE * duration);

    // --- Loudspeakers setup -------------------------------------------------
    long ac;
    float *av;
    if(dimension==2) // 2D, spread evenly over circle
    {
        if(c_param_count == 0)
        {
            c_param_count = 1;
            channels = 3;
        }
        else if (c_param_count == 1)
            channels = azi[0];
        else
            channels = c_param_count;

        ac = 1+channels;
        av = malloc(ac * sizeof(float));
        av[0] = 2.0;
        
        if(c_param_count == 1)
        {
            float angular_distance = 360./channels;
            for(int i=0; i<channels; i++)
            {
                av[i+1] = i * angular_distance;
            }
        }
        else
        {
            for(int i=0; i<channels; i++)
            {
                if(i%2==0)
                    av[i+1] = azi[i - i/2];
                else
                    av[i+1] = ele[i - (i/2+1)];
            }
        }
    }
    else if(dimension==3) // 3D, distribution on sphere using "Vogel's" method
    {
        if(c_param_count == 0)
        {
            c_param_count = 1;
            channels = 4;
        }
        
        if (c_param_count == 1)
        {
            channels = azi[0];
            vogel_distribution((double *)azi, (double *)ele, channels);
        }
        else if(c_param_count % 2 != 0)
        {
            printf("Warning: missing last elevation value, setting to 0.\n");
            ele[channels] = 0.;
            channels++;
        }

        ac = 1+channels*2;
        av = malloc(ac * sizeof(float));
        av[0] = 3.0;
        for (int i=0; i<channels; ++i)
        {
            av[i*2+1] = azi[i];
            av[i*2+2] = ele[i];
        }
    }

    if(round_to_integer)
    {
        for(int i=1; i<ac; i++)
            av[i] = round(av[i]);
    }

    if(b_verbose)
    {
        printf("== Channel positions: ");
        for(int i=0; i<channels; i++)
        {
            if(round_to_integer)
                printf("%d %d ", (int)av[i*2+1], (int)av[i*2+2]);
            else
                printf("\n%0.2f %0.2f ", av[i*2+1], av[i*2+2]);
        }
        printf("\n\n");
    }

    // horizontal_p = fmod(fmod(horizontal_p, 360.) + 360., 360.);
    // vertical_p = fmod(fmod(vertical_p, 360.) + 360., 360.);

    /*-----------------------------------------------*/
    // ==> code from vbap object starts here
    /*-----------------------------------------------*/

    initContent_ls_directions(ac, av);

    x_lsset_available = 0;
    x_gain = 1.0;
    x_spread = 0.0; // TODO: spread currently not regarded
    x_spread_base[0] = 0.0; // TODO: spread currently not regarded
    x_spread_base[1] = 1.0; // TODO: spread currently not regarded
    x_spread_base[2] = 0.0; // TODO: spread currently not regarded

    /*-----------------------------------------------*/
    // from def_ls_bang(t_def_ls *x)
    if(x_def_ls_dimension == 3)
    {
        if(x_triplets_specified==0)
            choose_ls_triplets();
        calculate_3x3_matrixes();
    } 
    else if(x_def_ls_dimension == 2)
    {
        choose_ls_tuplets();
    }

    double *output = malloc(samples * x_ls_amount * sizeof(double));//[channels][samples];
    double step, input;
    double input_per_second = 2*M_PI * input_f;
    double h_angles_per_second = 360. * horizontal_f;
    double v_angles_per_second = 360. * vertical_f;
    double time_sum = 0;
    for(int j=0;j<trials;j++)
    {

        /*-----------------------------------------------*/
        // from vbap_bang
        // t_atom at[MAX_LS_AMOUNT]; 
        float final_gs[x_ls_amount];
        long i;
        // float *final_gs = (float *) sysmem_newptr (x_ls_amount * sizeof(float));

        clock_t start, end;
        double cpu_time_used;
        if(x_lsset_available ==1)
        {
            for(int i=0; i<samples; ++i)
            {
                // time measurement
                start = clock();

                step = i * INV_SAMPLE_RATE;

                if(x_def_ls_dimension == 3)
                    x_ele = fmod(v_angles_per_second * step + vertical_p, 360.);
                else if(x_def_ls_dimension == 2)
                    x_ele = 0; 

                x_azi = fmod(h_angles_per_second * step + horizontal_p, 360.);

                vbap_bang(final_gs, x_ls_amount, dimension);

                // stop clock
                end = clock();
                cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
                if(b_verbose) printf("-- Time used: %0.10f\n", cpu_time_used); 
                time_sum += cpu_time_used;

                input = cos(input_per_second * step * 2*M_PI);
                for(int ch=0; ch<x_ls_amount; ch++)
                {
                    output[ch * samples + i] = input * final_gs[ch];
                }

                // This is MaxMSP specific output, we replace it with print statements
                // for(i=0;i<x_ls_amount;i++) 
                // {
                //         // atom_setlong(&at[0], i);
                //         // atom_setfloat(&at[1], final_gs[i]*x_gain); // gain is applied here
                //         // outlet_list(x_outlet0, 0L, 2, at);
                // }
                // outlet_float(x_outlet1, x_azi);
                // outlet_float(x_outlet2, x_ele);
                // outlet_float(x_outlet3, x_spread);
                // outlet_float(x_outlet4, x_gain);
            }
        }
        else
            // object_error((t_object*)x, "vbap: Configure loudspeakers first!");
            printf("ERROR: vbap: Configure loudspeakers first!\n");

        // sysmem_freeptr(final_gs); // bug fix added 9/00
        // free(final_gs);
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

    free(av); // loudspeaker array
    free(output); // sample buffer

    // free memory allocated by void add_ldsp_triplet(int i, int j, int k)
    t_ls_set *trip_ptr,  *tmp_ptr, *prev;
    prev = NULL;
    trip_ptr = x_ls_set;
    while (trip_ptr != NULL)
    {
        tmp_ptr = trip_ptr;
        trip_ptr = trip_ptr->next;
        free(tmp_ptr);
    }
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



/*--------------------------------------------------------------------------*/
// HELPER: panning angle azimuth
// void vbap_set_azimuth(t_vbap *x, double n) { x_azi = n; } // original
void vbap_set_azimuth(double n) { x_azi = n; }
/*--------------------------------------------------------------------------*/
// HELPER: panning angle elevation
// void vbap_set_elevation(t_vbap *x, double n) { x_ele = n; } // original
void vbap_set_elevation(double n) { x_ele = n; }
/*--------------------------------------------------------------------------*/
// HELPER: converts angular coordinates to cartesian
void angle_to_cart(double azi, double ele, float res[3])
{
    res[0] = cos(azi * atorad) * cos(ele * atorad);
    res[1] = sin(azi * atorad) * cos(ele * atorad);
    res[2] = sin(ele * atorad);
}
/*--------------------------------------------------------------------------*/
// HELPER: converts cartesian coordinates to angular
void cart_to_angle(float cvec[3], float avec[3])
{
  //float tmp, tmp2, tmp3, tmp4;
  //float power;
  float dist, atan_y_per_x, atan_x_pl_y_per_z;
  float azi, ele;
  
  if(cvec[0]==0.0)
    atan_y_per_x = M_PI / 2;
  else
    atan_y_per_x = atan(cvec[1] / cvec[0]);
  azi = atan_y_per_x / atorad;
  if(cvec[0]<0.0)
    azi +=180.0;
  dist = sqrt(cvec[0]*cvec[0] + cvec[1]*cvec[1]);
  if(cvec[2]==0.0)
    atan_x_pl_y_per_z = 0.0;
  else
    atan_x_pl_y_per_z = atan(cvec[2] / dist);
  if(dist == 0.0)
    {
    if(cvec[2]<0.0)
      atan_x_pl_y_per_z = -M_PI/2.0;
    else
      atan_x_pl_y_per_z = M_PI/2.0;
    }
  ele = atan_x_pl_y_per_z / atorad;
  dist = sqrtf(cvec[0] * cvec[0] +cvec[1] * cvec[1] +cvec[2]*cvec[2]);
  avec[0]=azi;
  avec[1]=ele;
  avec[2]=dist;
}
/*--------------------------------------------------------------------------*/
// HELPER: vector cross product 
void vect_cross_prod(float v1[3], float v2[3], float v3[3]) 
{
  float length;
  v3[0] = (v1[1] * v2[2] ) - (v1[2] * v2[1]);
  v3[1] = (v1[2] * v2[0] ) - (v1[0] * v2[2]);
  v3[2] = (v1[0] * v2[1] ) - (v1[1] * v2[0]);

  length= sqrt(v3[0]*v3[0] + v3[1]*v3[1] + v3[2]*v3[2]);
  v3[0] /= length;
  v3[1] /= length;
  v3[2] /= length;
}
/*--------------------------------------------------------------------------*/
// HELPER: vector cross product for loudspeakers
void ls_cross_prod(t_ls v1,t_ls v2, t_ls *res)
{
  float length;
  res->x = (v1.y * v2.z ) - (v1.z * v2.y);
  res->y = (v1.z * v2.x ) - (v1.x * v2.z);
  res->z = (v1.x * v2.y ) - (v1.y * v2.x);

  length= vec_length(*res);
  res->x /= length;
  res->y /= length;
  res->z /= length;
}
/*--------------------------------------------------------------------------*/
// HELPER: spherical to cartesian conversion for a set of loudspeakers
void ls_angles_to_cart(t_ls *ls)
// convert angular direction to cartesian
{
  float azi = ls->azi;
  float ele = ls->ele;
  ls->x = cos((float) azi * atorad) * cos((float) ele * atorad);
  ls->y = sin((float) azi * atorad) * cos((float) ele * atorad);
  ls->z = sin((float) ele * atorad);
}
/*--------------------------------------------------------------------------*/
// HELPER: angle between two loudspeakers
float vec_angle(t_ls v1, t_ls v2)
{
  float inner= ((v1.x*v2.x + v1.y*v2.y + v1.z*v2.z)/
              (vec_length(v1) * vec_length(v2)));
  if(inner > 1.0)
    inner= 1.0;
  if (inner < -1.0)
    inner = -1.0;
  return fabs( acos( inner));
}
/*--------------------------------------------------------------------------*/
// HELPER: length of a vector
float vec_length(t_ls v1)
{
  return (sqrt(v1.x*v1.x + v1.y*v1.y + v1.z*v1.z));
}
/*--------------------------------------------------------------------------*/
// HELPER: vector dot product
float vec_prod(t_ls v1, t_ls v2)
{
  return (v1.x*v2.x + v1.y*v2.y + v1.z*v2.z);
}
/*--------------------------------------------------------------------------*/
// HELPER: checks if two lines intersect on 3D sphere 
int lines_intersect(int i,int j,int k,int l,t_ls  lss[MAX_LS_AMOUNT])
{
  t_ls v1;
  t_ls v2;
  t_ls v3, neg_v3;
  //float angle;
  float dist_ij,dist_kl,dist_iv3,dist_jv3,dist_inv3,dist_jnv3;
  float dist_kv3,dist_lv3,dist_knv3,dist_lnv3;

  ls_cross_prod(lss[i],lss[j],&v1);
  ls_cross_prod(lss[k],lss[l],&v2);
  ls_cross_prod(v1,v2,&v3);

  neg_v3.x= 0.0 - v3.x;
  neg_v3.y= 0.0 - v3.y;
  neg_v3.z= 0.0 - v3.z;

  dist_ij = (vec_angle(lss[i],lss[j]));
  dist_kl = (vec_angle(lss[k],lss[l]));
  dist_iv3 = (vec_angle(lss[i],v3));
  dist_jv3 = (vec_angle(v3,lss[j]));
  dist_inv3 = (vec_angle(lss[i],neg_v3));
  dist_jnv3 = (vec_angle(neg_v3,lss[j]));
  dist_kv3 = (vec_angle(lss[k],v3));
  dist_lv3 = (vec_angle(v3,lss[l]));
  dist_knv3 = (vec_angle(lss[k],neg_v3));
  dist_lnv3 = (vec_angle(neg_v3,lss[l]));

  /* if one of loudspeakers is close to crossing point, don't do anything*/
  if(fabsf(dist_iv3) <= 0.01 || fabsf(dist_jv3) <= 0.01 || 
         fabsf(dist_kv3) <= 0.01 || fabsf(dist_lv3) <= 0.01 ||
     fabsf(dist_inv3) <= 0.01 || fabsf(dist_jnv3) <= 0.01 || 
     fabsf(dist_knv3) <= 0.01 || fabsf(dist_lnv3) <= 0.01 )
    return(0);

  // if crossing point is on line between both loudspeakers return 1
  if (((fabsf(dist_ij - (dist_iv3 + dist_jv3)) <= 0.01 ) &&
       (fabsf(dist_kl - (dist_kv3 + dist_lv3))  <= 0.01)) ||
      ((fabsf(dist_ij - (dist_inv3 + dist_jnv3)) <= 0.01)  &&
       (fabsf(dist_kl - (dist_knv3 + dist_lnv3)) <= 0.01 ))) {
    return (1);
  } else {
    return (0);
  }
}



void vbap_bang(float *final_gs, int x_ls_amount, int x_dimension)
{
  float g[3];
  long ls[3];
  long i;

  vbap(g,ls);

  for(i=0;i<x_ls_amount;i++)
      final_gs[i]=0.0;

  for(i=0;i<x_dimension;i++)
  {
      final_gs[ls[i]-1]=g[i];
  }

  // TODO: no spreading for now...
  // if(x_spread != 0.0)
  // {
  //     spread_it(x,final_gs);
  // }
}



/*--------------------------------------------------------------------------*/
// Main VBAP function
// void vbap(float g[3], long ls[3], t_vbap *x) // original
void vbap(float g[3], long ls[3])
{
  /* calculates gain factors using loudspeaker setup and given direction */
  float power;
  int i,j,k, gains_modified;
  float small_g;
  float big_sm_g, gtmp[3];
  long winner_set = 0;
  float cartdir[3];
  float new_cartdir[3];
  float new_angle_dir[3];
  long dim = x_dimension;
  long neg_g_am, best_neg_g_am;
  
  // transfering the azimuth angle to a decent value
  while(x_azi > 180.0)
    x_azi -= 360.0;
  while(x_azi <= -180.0)
    x_azi += 360.0;
    
  // transferring the elevation to a decent value
  if(dim == 3){
    while(x_ele > 180.0)
        x_ele -= 360.0;
    while(x_ele <= -180.0)
        x_ele += 360.0;
  } else
    x_ele = 0.0;
  
  
  // go through all defined loudspeaker sets and find the set which
  // has all positive values. If such is not found, set with largest
  // minimum value is chosen. If at least one of gain factors of one LS set is negative
  // it means that the virtual source does not lie in that LS set. 
  
  angle_to_cart(x_azi,x_ele,cartdir);
  big_sm_g = -100000.0;   // initial value for largest minimum gain value
  best_neg_g_am=3;        // how many negative values in this set
  
  for(i=0;i<x_lsset_amount;i++)
    {
    small_g = 10000000.0;
    neg_g_am = 3;
    for(j=0;j<dim;j++)
        {
      gtmp[j]=0.0;
      for(k=0;k<dim;k++)
        gtmp[j]+=cartdir[k]* x_set_inv_matx[i][k+j*dim];
      if(gtmp[j] < small_g)
        small_g = gtmp[j];
      if(gtmp[j]>= -0.01)
        neg_g_am--;
    }
    if(small_g > big_sm_g && neg_g_am <= best_neg_g_am)
        {
      big_sm_g = small_g;
      best_neg_g_am = neg_g_am; 
      winner_set=i;
      g[0]=gtmp[0]; g[1]=gtmp[1];
      ls[0]= x_lsset[i][0]; ls[1]= x_lsset[i][1];
      if(dim==3)
            {
        g[2]=gtmp[2];
        ls[2]= x_lsset[i][2];
      } 
            else 
            {
        g[2]=0.0;
        ls[2]=0;
      }
    }
  }
  
  // If chosen set produced a negative value, make it zero and
  // calculate direction that corresponds  to these new
  // gain values. This happens when the virtual source is outside of
  // all loudspeaker sets. 
  
  //
    gains_modified=0;
    for(i=0;i<dim;i++)
        if(g[i]<-0.01){
            g[i]=0.0001;
            gains_modified=1;
        }   
    if(gains_modified==1){
        new_cartdir[0] =  x_set_matx[winner_set][0] * g[0] 
                        + x_set_matx[winner_set][1] * g[1]
                        + x_set_matx[winner_set][2] * g[2];
        new_cartdir[1] =  x_set_matx[winner_set][3] * g[0] 
                        + x_set_matx[winner_set][4] * g[1] 
                        + x_set_matx[winner_set][5] * g[2];
        if(dim==3){
            new_cartdir[2] =  x_set_matx[winner_set][6] * g[0] 
                            + x_set_matx[winner_set][7] * g[1]
                            + x_set_matx[winner_set][8] * g[2];
        } else new_cartdir[2] = 0;
        cart_to_angle(new_cartdir,new_angle_dir);
        x_azi = (long) (new_angle_dir[0] + 0.5);
        x_ele = (long) (new_angle_dir[1] + 0.5);
     }
  //}
  
  power=sqrt(g[0]*g[0] + g[1]*g[1] + g[2]*g[2]);
  g[0] /= power;
  g[1] /= power;
  g[2] /= power;
}



/*--------------------------------------------------------------------------*/
// Initialize the speaker positions from ls-directions list
// void initContent_ls_directions(t_def_ls *x,int ac,t_atom*av) // original
void initContent_ls_directions(int ac, float *av)
{
    x_ls_read = 0;
    
    long d = 0;
    // if (av[0].a_type == A_LONG) d = av[0].a_w.w_long;
    // else if(av[0].a_type == A_FLOAT) d = (long)av[0].a_w.w_float;
    // else { object_error((t_object*)x, "define-loudspeakers: dimension NaN"); return; }
    d = (long)av[0];

    if (d==2 || d==3)
    {
         x_def_ls_dimension= d;
         x_ls_read = 1;
    } 
    else
    {
        x_def_ls_dimension= 0;
        // object_error((t_object*)x, "define-loudspeakers: Dimension has to be 2 or 3!");
        printf("ERROR: define-loudspeakers: Dimension has to be 2 or 3!\n");
        return;
    }
        
    int pointer = 1;
    x_def_ls_amount= (ac-1) / (x_def_ls_dimension - 1);

    // read loudspeaker direction angles  
    for(int i=0; i < x_def_ls_amount;i++)
    {
        float azi = 0;
        // if(av[pointer].a_type == A_LONG) azi = (float) av[pointer].a_w.w_long;
        // else if(av[pointer].a_type == A_FLOAT) azi = av[pointer].a_w.w_float;
        // else { object_error((t_object*)x, "define-loudspeakers: direction angle #%d NaN",i+1); x->x_ls_read = 0; return; }
        azi = av[pointer];

        x_ls[i].azi = azi;
        
        pointer++;

        float ele = 0; // in 2d elevation is zero
        if(x_def_ls_dimension == 3)
        {  // 3-D 
            // if(av[pointer].a_type == A_LONG) ele = (float) av[pointer].a_w.w_long;
            // else if(av[pointer].a_type == A_FLOAT) ele = av[pointer].a_w.w_float;
            // else { object_error((t_object*)x, "define-loudspeakers: elevation #%d NaN",i+1);  x->x_ls_read = 0; return; }
            ele = av[pointer];

            pointer++;
        } 
        x_ls[i].ele = ele;
    }
    
    if(x_ls_read == 1)
    {
        for(int i=0;i<x_def_ls_amount;i++)
        {
            ls_angles_to_cart(&x_ls[i]); 
        }
    }
    x_triplets_specified=0;
    x_ls_set = NULL;
}



/*--------------------------------------------------------------------------*/
// selects the loudspeaker pairs, calculates the inversion matrices and stores the data to a global array
// void choose_ls_tuplets(t_def_ls *x) // original
void choose_ls_tuplets()
{
  //float atorad = (2 * 3.1415927 / 360) ;
  int i,j;
  //float w1,w2;
  //float p1,p2;
  int sorted_lss[MAX_LS_AMOUNT];
  int exist[MAX_LS_AMOUNT];   
  int amount=0;
  float inv_mat[MAX_LS_AMOUNT][4];  // In 2-D ls amount == max amount of LS pairs
  float mat[MAX_LS_AMOUNT][4];
  //float *ptr;   
  //float *ls_table;
  t_ls *lss = x_ls;
  long ls_amount= x_def_ls_amount;
  long list_length;
  // t_atom *at;
  float *at;
  long pointer;
  
  for(i=0;i<MAX_LS_AMOUNT;i++){
    exist[i]=0;
  }

  /* sort loudspeakers according their aximuth angle */
  sort_2D_lss(x_ls,sorted_lss,ls_amount);

  /* adjacent loudspeakers are the loudspeaker pairs to be used.*/
  for(i=0;i<(ls_amount-1);i++){
    if((lss[sorted_lss[i+1]].azi - 
        lss[sorted_lss[i]].azi) <= (180 - 10)){
      if (calc_2D_inv_tmatrix( lss[sorted_lss[i]].azi, 
                               lss[sorted_lss[i+1]].azi, 
                               inv_mat[i],mat[i]) != 0){
        exist[i]=1;
        amount++;
      }
    }
  }

  if(((360 - lss[sorted_lss[ls_amount-1]].azi) 
      +lss[sorted_lss[0]].azi) <= (180 -  10)) {
    if(calc_2D_inv_tmatrix(lss[sorted_lss[ls_amount-1]].azi, 
                           lss[sorted_lss[0]].azi, 
                           inv_mat[ls_amount-1],mat[ls_amount-1]) != 0) { 
        exist[ls_amount-1]=1;
        amount++;
    } 
  }
  
  
  // Output
  list_length= amount * 10  + 2;
  // at= (t_atom *) sysmem_newptr (list_length*sizeof(t_atom));
  at = (float *) malloc(list_length * sizeof(float));
  
  // atom_setlong(&at[0], x->x_def_ls_dimension);
  at[0] = (float) x_def_ls_dimension;
  // atom_setlong(&at[1], x->x_def_ls_amount);
  at[1] = (float) x_def_ls_amount;
  pointer=2;
  
  for (i=0;i<ls_amount - 1;i++){
    if(exist[i] == 1) {
        // atom_setlong(&at[pointer], sorted_lss[i]+1);
        at[pointer] = (float) sorted_lss[i]+1;
        pointer++;
        // atom_setlong(&at[pointer], sorted_lss[i+1]+1);
        at[pointer] = (float) sorted_lss[i+1]+1;
        pointer++;
        for(j=0;j<4;j++) {
            // atom_setfloat(&at[pointer], inv_mat[i][j]);
            at[pointer] = inv_mat[i][j];
            pointer++;
        }
       for(j=0;j<4;j++) {
            // atom_setfloat(&at[pointer], mat[i][j]);
            at[pointer] = mat[i][j];
            pointer++;
        }
    }
  }
  if(exist[ls_amount-1] == 1) {
    // atom_setlong(&at[pointer], sorted_lss[ls_amount-1]+1);
    at[pointer] = (float) sorted_lss[ls_amount-1]+1;
    pointer++;
    // atom_setlong(&at[pointer], sorted_lss[0]+1);
    at[pointer] = (float) sorted_lss[0]+1;
    pointer++;
    for(j=0;j<4;j++) {
        // atom_setfloat(&at[pointer], inv_mat[ls_amount-1][j]);
        at[pointer] = inv_mat[ls_amount-1][j];
        pointer++;
    }
    for(j=0;j<4;j++) {
        // atom_setfloat(&at[pointer], mat[ls_amount-1][j]);
        at[pointer] = mat[ls_amount-1][j];
        pointer++;
    }
  }
    // sendLoudspeakerMatrices(x,list_length, at);
  vbap_matrix(list_length, at); // we jump to vbap_matrix directly
  //outlet_anything(x->x_outlet0, gensym("loudspeaker-matrices"), list_length, at);
  // sysmem_freeptr(at);
  free(at);
}

/*--------------------------------------------------------------------------*/
// sort loudspeakers according to azimuth angle
void sort_2D_lss(t_ls lss[MAX_LS_AMOUNT], int sorted_lss[MAX_LS_AMOUNT], int ls_amount)
{
  float tmp, tmp_azi;
//  float rad2ang = 360.0f / ( 2.0f * M_PI );

  //float x,y;
  /* Transforming angles between -180 and 180 */
  for (int i=0;i<ls_amount;i++) 
    {
    ls_angles_to_cart(&lss[i]);
    lss[i].azi = acos( lss[i].x) * rad2ang;
    if (fabs(lss[i].y) <= 0.001)
        tmp = 1.0;
    else
        tmp = lss[i].y / fabs(lss[i].y);
    lss[i].azi *= tmp;
  }
  for (int i=0;i<ls_amount;i++)
    {
    tmp = 2000;
        int index = 0;
    for (int j=0 ; j<ls_amount;j++)
        {
      if (lss[j].azi <= tmp)
            {
        tmp=lss[j].azi;
        index = j;
      }
    }
    sorted_lss[i]=index;
    tmp_azi = (lss[index].azi);
    lss[index].azi = (tmp_azi + (float) 4000.0);
  }
  for (int i=0;i<ls_amount;i++) 
    {
    tmp_azi = (lss[i].azi);
    lss[i].azi = (tmp_azi - (float) 4000.0);
  }
}

/*--------------------------------------------------------------------------*/
// calculate inverse 2x2 matrix
int calc_2D_inv_tmatrix(float azi1,float azi2, float inv_mat[4],float mat[4])
{
  float x1,x2,x3,x4; /* x1 x3 */
  //float y1,y2,y3,y4; /* x2 x4 */
  float det;
  
  mat[0]=x1 = cos(azi1 / rad2ang);
  mat[1]=x2 = sin(azi1 / rad2ang);
  mat[2]=x3 = cos(azi2 / rad2ang);
  mat[3]=x4 = sin(azi2 / rad2ang);
  det = (x1 * x4) - ( x3 * x2 );
   if(fabsf(det) <= 0.001) {

    inv_mat[0] = 0.0;
    inv_mat[1] = 0.0;
    inv_mat[2] = 0.0;
    inv_mat[3] = 0.0;
    return 0;
  } else {
    inv_mat[0] =   (x4 / det);
    inv_mat[1] =  (-x3 / det);
    inv_mat[2] =   (-x2 / det);
    inv_mat[3] =    (x1 / det);
    return 1;
  }
}



/*--------------------------------------------------------------------------*/
/* Selects the loudspeaker triplets, and calculates the inversion matrices
    for each selected triplet. A line (connection) is drawn between each
    loudspeaker. The lines denote the sides of the triangles. The triangles 
    should not be intersecting. All crossing connections are searched and the
    longer connection is erased. This yields non-intesecting triangles,
    which can be used in panning. 
   See theory in paper Pulkki, V. Lokki, T. "Creating Auditory Displays with
    Multiple Loudspeakers Using VBAP: A Case Study with DIVA Project" in
    International Conference on Auditory Displays -98. */
void choose_ls_triplets()
{
  int i,j,k,l,/*m,li,*/ table_size;
  //int *i_ptr;
  //t_ls vb1,vb2,tmp_vec;
  int connections[MAX_LS_AMOUNT][MAX_LS_AMOUNT];
  //float angles[MAX_LS_AMOUNT];
  //int sorted_angles[MAX_LS_AMOUNT];
  float distance_table[((MAX_LS_AMOUNT * (MAX_LS_AMOUNT - 1)) / 2)];
  int distance_table_i[((MAX_LS_AMOUNT * (MAX_LS_AMOUNT - 1)) / 2)];
  int distance_table_j[((MAX_LS_AMOUNT * (MAX_LS_AMOUNT - 1)) / 2)];
  float distance;
  t_ls_set *trip_ptr, *prev, *tmp_ptr;
  int ls_amount = x_def_ls_amount;
  t_ls *lss = x_ls;
  if (ls_amount == 0) { 
    // object_post((t_object*)x, "define-loudspeakers: Number of loudspeakers is zero");
    printf("ERROR: define-loudspeakers: Number of loudspeakers is zero\n");
    return; 
  }
 
  for(i=0;i<ls_amount;i++)
    for(j=i+1;j<ls_amount;j++)
      for(k=j+1;k<ls_amount;k++)
            {
        if(vol_p_side_lgth(i,j,k, x_ls) > MIN_VOL_P_SIDE_LGTH)
                {
          connections[i][j]=1;
          connections[j][i]=1;
          connections[i][k]=1;
          connections[k][i]=1;
          connections[j][k]=1;
          connections[k][j]=1;
          add_ldsp_triplet(i,j,k);
        }
      }
   

  /*calculate distancies between all lss and sorting them*/
  table_size =(((ls_amount - 1) * (ls_amount)) / 2); 
  for(i=0;i<table_size; i++)
    distance_table[i] = 100000.0;
  for(i=0;i<ls_amount;i++)
    { 
    for(j=(i+1);j<ls_amount; j++)
        { 
          if(connections[i][j] == 1) 
          {
        distance = fabs(vec_angle(lss[i],lss[j]));
        k=0;
        while(distance_table[k] < distance)
          k++;
        for(l=(table_size - 1);l > k ;l--)
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
            {
        table_size--;
            }
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
      for(j=0; j<ls_amount ; j++)
            {
        for(k=j+1; k<ls_amount; k++)
                {
          if( (j!=fst_ls) && (k != sec_ls) && (k!=fst_ls) && (j != sec_ls))
                    {
            if(lines_intersect(fst_ls, sec_ls, j,k,x_ls) == 1)
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
  trip_ptr = x_ls_set;
  prev = NULL;
  while (trip_ptr != NULL)
    {
    i = trip_ptr->ls_nos[0];
    j = trip_ptr->ls_nos[1];
    k = trip_ptr->ls_nos[2];
    if(connections[i][j] == 0 || 
       connections[i][k] == 0 || 
       connections[j][k] == 0 ||
             any_ls_inside_triplet(i,j,k,x_ls,ls_amount) == 1 )
        {
      if(prev != NULL) 
            {
        prev->next = trip_ptr->next;
        tmp_ptr = trip_ptr;
        trip_ptr = trip_ptr->next;
        free(tmp_ptr);
      } 
            else 
            {
        x_ls_set = trip_ptr->next;
        tmp_ptr = trip_ptr;
        trip_ptr = trip_ptr->next;
        free(tmp_ptr);
      }
    } 
        else 
        {
      prev = trip_ptr;
      trip_ptr = trip_ptr->next;
    }
  }
  x_triplets_specified=1;
}

/*--------------------------------------------------------------------------*/
/* calculate volume of the parallelepiped defined by the loudspeaker
    direction vectors and divide it with total length of the triangle sides. 
    This is used when removing too narrow triangles. */
float vol_p_side_lgth(int i, int j,int k, t_ls lss[MAX_LS_AMOUNT])
{

  float volper, lgth;
  t_ls xprod;
  ls_cross_prod(lss[i], lss[j], &xprod);
  volper = fabsf(vec_prod(xprod, lss[k]));
  lgth = (fabsf(vec_angle(lss[i],lss[j])) 
          + fabsf(vec_angle(lss[i],lss[k])) 
          + fabsf(vec_angle(lss[j],lss[k])));
  if(lgth>0.00001)
    return volper / lgth;
  else
    return 0.0;
}

/*--------------------------------------------------------------------------*/
// adds i,j,k triplet to structure
void add_ldsp_triplet(int i, int j, int k)
{
  struct t_ls_set *trip_ptr, *prev;
  trip_ptr = x_ls_set;
  prev = NULL;
  while (trip_ptr != NULL)
    {
    prev = trip_ptr;
    trip_ptr = trip_ptr->next;
  }
  trip_ptr = (struct t_ls_set*) malloc (sizeof (struct t_ls_set));
  if(prev == NULL)
    x_ls_set = trip_ptr;
  else 
    prev->next = trip_ptr;
  trip_ptr->next = NULL;
  trip_ptr->ls_nos[0] = i;
  trip_ptr->ls_nos[1] = j;
  trip_ptr->ls_nos[2] = k;
}

/*--------------------------------------------------------------------------*/
// returns 1 if there is loudspeaker(s) inside given ls triplet
int any_ls_inside_triplet(int a, int b, int c,t_ls lss[MAX_LS_AMOUNT],int ls_amount)
{
  float invdet;
  t_ls *lp1, *lp2, *lp3;
  float invmx[9];
  int i,j;
  float tmp;
  int any_ls_inside, this_inside;

  lp1 =  &(lss[a]);
  lp2 =  &(lss[b]);
  lp3 =  &(lss[c]);

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
  for(i=0; i< ls_amount; i++) 
    {
    if (i != a && i!=b && i != c)
        {
      this_inside = 1;
      for(j=0; j< 3; j++)
            {
        tmp = lss[i].x * invmx[0 + j*3];
        tmp += lss[i].y * invmx[1 + j*3];
        tmp += lss[i].z * invmx[2 + j*3];
        if(tmp < -0.001)
          this_inside = 0;
      }
      if(this_inside == 1)
        any_ls_inside=1;
    }
  }
  return any_ls_inside;
}

/*--------------------------------------------------------------------------*/
// Calculates the inverse matrices for 3D
void  calculate_3x3_matrixes()
{  
  float invdet;
  t_ls *lp1, *lp2, *lp3;
  float *invmx;
  //float *ptr;
  struct t_ls_set *tr_ptr = x_ls_set;
  int triplet_amount = 0, /*ftable_size,*/i,pointer,list_length=0;
  float *at;
  t_ls *lss = x_ls;
  
  if (tr_ptr == NULL)
    {
    // object_error((t_object*)x, "define-loudspeakers: Not valid 3-D configuration\n");
    printf("ERROR: define-loudspeakers: Not valid 3-D configuration\n");
    return;
  }
    
  /* counting triplet amount */
  while(tr_ptr != NULL)
    {
    triplet_amount++;
    tr_ptr = tr_ptr->next;
  }
  tr_ptr = x_ls_set;
  list_length= triplet_amount * 21 + 3;
  at= (float *) malloc(list_length*sizeof(float));
  
  // atom_setlong(&at[0], x->x_def_ls_dimension);
  at[0] = (float) x_def_ls_dimension;
  // atom_setlong(&at[1], x->x_def_ls_amount);
  at[1] = (float) x_def_ls_amount;
  pointer=2;
  
  while(tr_ptr != NULL){
    lp1 =  &(lss[tr_ptr->ls_nos[0]]);
    lp2 =  &(lss[tr_ptr->ls_nos[1]]);
    lp3 =  &(lss[tr_ptr->ls_nos[2]]);

    /* matrix inversion */
    invmx = tr_ptr->inv_mx;
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

    for(i=0;i<3;i++){
        // atom_setlong(&at[pointer], tr_ptr->ls_nos[i]+1);
        at[pointer] = (float) (tr_ptr->ls_nos[i]+1);
        pointer++;
    }
    for(i=0;i<9;i++){
        // atom_setfloat(&at[pointer], invmx[i]);
        at[pointer] = invmx[i];
        pointer++;
    }
    // atom_setfloat(&at[pointer], lp1->x);
    at[pointer] = lp1->x; pointer++;
    // atom_setfloat(&at[pointer], lp2->x);
    at[pointer] = lp2->x; pointer++;
    // atom_setfloat(&at[pointer], lp3->x);
    at[pointer] = lp3->x; pointer++;
    // atom_setfloat(&at[pointer], lp1->y);
    at[pointer] = lp1->y; pointer++;
    // atom_setfloat(&at[pointer], lp2->y);
    at[pointer] = lp2->y; pointer++;
    // atom_setfloat(&at[pointer], lp3->y);
    at[pointer] = lp3->y; pointer++;
    // atom_setfloat(&at[pointer], lp1->z);
    at[pointer] = lp1->z; pointer++;
    // atom_setfloat(&at[pointer], lp2->z);
    at[pointer] = lp2->z; pointer++;
    // atom_setfloat(&at[pointer], lp3->z);
    at[pointer] = lp3->z; pointer++;
 
    tr_ptr = tr_ptr->next;
  }
    // sendLoudspeakerMatrices(x,list_length, at);
  vbap_matrix(list_length, at); // we jump to vbap_matrix directly
//  outlet_anything(x->x_outlet0, gensym("loudspeaker-matrices"), list_length, at);
  free(at);
}



/*--------------------------------------------------------------------------*/
// read in loudspeaker matrices
// void vbap_matrix(t_vbap *x, t_symbol *s, int ac, t_atom *av) // original
void vbap_matrix(int ac, float *av)
{
    int datapointer = 0; 
    if(ac>0) 
    {
        int d = 0;
        // if(av[datapointer].a_type == A_LONG) d = av[datapointer++].a_w.w_long;
        // else if(av[datapointer].a_type == A_FLOAT) d = (long)av[datapointer++].a_w.w_float;
        // else { object_error((t_object*)x, "vbap: Dimension NaN"); x->x_lsset_available=0; return; }
        d = (long)av[datapointer++];

        if (d!=2 && d!=3) { 
            // object_error((t_object*)x, "vbap %s: Dimension can be only 2 or 3",s->s_name);
            printf("ERROR: vbap: Dimension can be only 2 or 3\n");
            x_lsset_available=0;
            return;
        }

        x_dimension = d;
        x_lsset_available=1;
    }
    else {
        // object_error((t_object*)x, "vbap %s: bad empty parameter list",s->s_name);
        printf("ERROR: vbap: bad empty parameter list\n");
        x_lsset_available=0;
        return;
    }

    if(ac>1) 
    {
        long a = 0;
        // if(av[datapointer].a_type == A_LONG) a = av[datapointer++].a_w.w_long;
        // else if(av[datapointer].a_type == A_FLOAT) a = (long) av[datapointer++].a_w.w_float;
        // else { object_error((t_object*)x, "vbap: ls_amount NaN");  x->x_lsset_available=0; return; }
        a = (long) av[datapointer++];

        x_ls_amount = a;
    }
    
    long counter = (ac - 2) / ((x_dimension * x_dimension*2) + x_dimension);
    x_lsset_amount=counter;

    if(counter==0) {
        // object_error((t_object*)x, "vbap %s: not enough parameters",s->s_name);
        printf("ERROR: vbap: not enough parameters\n");
        x_lsset_available=0;
        return;
    }
    
    long setpointer=0;
    long i;
 
    while(counter-- > 0)
    {
        for(i=0; i < x_dimension; i++)
        {
            // if(av[datapointer].a_type == A_LONG)
            // {
            //      x->x_lsset[setpointer][i]=av[datapointer++].a_w.w_long;
            // }
            // else { object_error((t_object*)x, "vbap %s: param %d is not an int",s->s_name,datapointer); x->x_lsset_available=0; return; }
            x_lsset[setpointer][i]=(long) av[datapointer++];
        }

        for(i=0; i < x_dimension*x_dimension; i++)
        {
            // if(av[datapointer].a_type == A_FLOAT)
            // {
            //     x->x_set_inv_matx[setpointer][i]=av[datapointer++].a_w.w_float;
            // }
            // else { object_error((t_object*)x, "vbap %s: param %d is not a float",s->s_name,datapointer); x->x_lsset_available=0; return; }
            x_set_inv_matx[setpointer][i]=av[datapointer++];
        }
        
        for(i=0; i < x_dimension*x_dimension; i++)
        {
            // if(av[datapointer].a_type == A_FLOAT)
            // {
            //     x->x_set_matx[setpointer][i]=av[datapointer++].a_w.w_float;
            // }
            // else { object_error((t_object*)x, "vbap %s: param %d is not a float",s->s_name,datapointer); x->x_lsset_available=0; return; }
            x_set_matx[setpointer][i]=av[datapointer++];
        }
    
        setpointer++;
    }
    // if (_enable_trace) object_post((t_object*)x, "vbap: Loudspeaker setup configured!");
}
