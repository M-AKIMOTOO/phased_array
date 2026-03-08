#ifndef __TOOL__
#define __TOOL__

#include <fftw3.h>
#include <xmmintrin.h>

#define PHS_DIV 256

float cos_table[PHS_DIV],sin_table[PHS_DIV];

void gico3_init();
void gico3_tone (__v4sf *vis_x,              __v4sf *sum,int num);
void gico3_cross(__v4sf *vis_x,__v4sf *vis_y,__v4sf *sum,int num);
void gico3_track(fftwf_complex *in,fftwf_complex *out,unsigned int phs0,unsigned int phs,int num);

void float_01bit_01ch(unsigned int *raw_data,float *sample,int num,int ch,float level[4]);
void float_01bit_02ch(unsigned int *raw_data,float *sample,int num,int ch,float level[4]);
void float_01bit_04ch(unsigned int *raw_data,float *sample,int num,int ch,float level[4]);
void float_01bit_08ch(unsigned int *raw_data,float *sample,int num,int ch,float level[4]);
void float_01bit_16ch(unsigned int *raw_data,float *sample,int num,int ch,float level[4]);
// Added by KT 19st Jan, 2015
void float_01bit_32ch(unsigned int *raw_data,float *sample,int num,int ch,float level[4]);

void float_02bit_01ch(unsigned int *raw_data,float *sample,int num,int ch,float level[4]);
void float_02bit_02ch(unsigned int *raw_data,float *sample,int num,int ch,float level[4]);
void float_02bit_04ch(unsigned int *raw_data,float *sample,int num,int ch,float level[4]);
void float_02bit_08ch(unsigned int *raw_data,float *sample,int num,int ch,float level[4]);
void float_02bit_16ch(unsigned int *raw_data,float *sample,int num,int ch,float level[4]);

void vssp32_01bit_01ch(unsigned int *raw_data,float *sample,int num,int ch,float level[4]);
void vssp32_01bit_04ch(unsigned int *raw_data,float *sample,int num,int ch,float level[4]);
void vssp32_02bit_01ch(unsigned int *raw_data,float *sample,int num,int ch,float level[4]);
void vssp32_02bit_04ch(unsigned int *raw_data,float *sample,int num,int ch,float level[4]);

#endif
