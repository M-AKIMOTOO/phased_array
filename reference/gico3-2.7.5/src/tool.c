#include <math.h>
#include <fftw3.h>
#include <xmmintrin.h>
#include "tool.h"


static __v4sf sign_table  ={+1.0,-1.0,+1.0,-1.0};//複素演算でADDSUB用の符号設定ベクトル

void gico3_init(){int phs;for(phs=0;phs<PHS_DIV;phs++){cos_table[phs]=cos(2*M_PI*phs/PHS_DIV);sin_table[phs]=sin(2*M_PI*phs/PHS_DIV);}}

void gico3_cross(__v4sf *vis_x,__v4sf *vis_y,__v4sf *sum,int num)
{
  int s;
  for(s=0;s<num/2;s++){
    sum[s]=__builtin_ia32_addps(sum[s],__builtin_ia32_mulps(__builtin_ia32_mulps( vis_y[s],sign_table  ),__builtin_ia32_shufps(vis_x[s],vis_x[s],160)));
    sum[s]=__builtin_ia32_addps(sum[s],__builtin_ia32_mulps(__builtin_ia32_shufps(vis_y[s],vis_y[s],177),__builtin_ia32_shufps(vis_x[s],vis_x[s],245)));
  }
}
void gico3_tone (__v4sf *vis_x,              __v4sf *sum,int num)
{
  int s;
  for(s=0;s<num/2;s++){
    sum[s]=__builtin_ia32_addps(sum[s],vis_x[s]); // sum = sum + vis

  }
  //printf("%d,test\n",num);
}

void gico3_track0(__v4sf *in,__v4sf *out,unsigned long long phs0,unsigned long long phs,int num)
{
  int s;__v4sf vector_cos,vector_sin;
  for(s=0;s<num/2;s++,phs0+=2*phs){
    vector_sin=__builtin_ia32_mulps(in[s],_mm_set_ps1(sin_table[phs0>>56]));vector_cos=__builtin_ia32_mulps(in[s],_mm_set_ps1(cos_table[phs0>>56]));
    vector_sin=__builtin_ia32_mulps(vector_sin,sign_table);vector_sin=__builtin_ia32_shufps(vector_sin,vector_sin,177);out[s]=__builtin_ia32_addps(vector_cos,vector_sin);
  }
}

void gico3_track(fftwf_complex *in,fftwf_complex *out,unsigned int phs0,unsigned int phs,int num)
{
  int s;
  for(s=0;s<num;s++){
    out[s][0]=+in[s][0]*cos_table[phs0>>24]-in[s][1]*sin_table[phs0>>24];
    out[s][1]=+in[s][0]*sin_table[phs0>>24]+in[s][1]*cos_table[phs0>>24];
    phs0+=phs;
  }
}

void float_01bit_01ch(unsigned int *raw_data,float *sample,int num,int ch,float level[4]){int s,t;for(s=0;s<num/32;s+=1) for(t=0;t<32;t++) sample[32*s+t]=level[(raw_data[s]>>( 1*t+32*ch)) & 0x01];} 
void float_01bit_02ch(unsigned int *raw_data,float *sample,int num,int ch,float level[4]){int s,t;for(s=0;s<num/16;s+=1) for(t=0;t<16;t++) sample[16*s+t]=level[(raw_data[s]>>( 1*t+16*ch)) & 0x01];} 
void float_01bit_04ch(unsigned int *raw_data,float *sample,int num,int ch,float level[4]){int s,t;for(s=0;s<num/ 8;s+=1) for(t=0;t< 8;t++) sample[ 8*s+t]=level[(raw_data[s]>>( 1*t+ 8*ch)) & 0x01];} 
void float_01bit_08ch(unsigned int *raw_data,float *sample,int num,int ch,float level[4]){int s,t;for(s=0;s<num/ 4;s+=1) for(t=0;t< 4;t++) sample[ 4*s+t]=level[(raw_data[s]>>( 1*t+ 4*ch)) & 0x01];} 
void float_01bit_16ch(unsigned int *raw_data,float *sample,int num,int ch,float level[4]){int s,t;for(s=0;s<num/ 2;s+=1) for(t=0;t< 2;t++) sample[ 2*s+t]=level[(raw_data[s]>>( 1*t+ 2*ch)) & 0x01];} 
// Added by KT 19st Jan, 2015
void float_01bit_32ch(unsigned int *raw_data,float *sample,int num,int ch,float level[4]){int s,t;for(s=0;s<num/ 1;s+=1) for(t=0;t< 1;t++) sample[ 1*s+t]=level[(raw_data[s]>>( 1*t+ 1*ch)) & 0x01];} 

void float_02bit_01ch(unsigned int *raw_data,float *sample,int num,int ch,float level[4]){int s,t;for(s=0;s<num/16;s+=1) for(t=0;t<16;t++) sample[16*s+t]=level[(raw_data[s]>>( 2*t+32*ch)) & 0x03];} 
void float_02bit_02ch(unsigned int *raw_data,float *sample,int num,int ch,float level[4]){int s,t;for(s=0;s<num/ 8;s+=1) for(t=0;t< 8;t++) sample[ 8*s+t]=level[(raw_data[s]>>( 2*t+16*ch)) & 0x03];}
void float_02bit_04ch(unsigned int *raw_data,float *sample,int num,int ch,float level[4]){int s,t;for(s=0;s<num/ 4;s+=1) for(t=0;t< 4;t++) sample[ 4*s+t]=level[(raw_data[s]>>( 2*t+ 8*ch)) & 0x03];} 
void float_02bit_08ch(unsigned int *raw_data,float *sample,int num,int ch,float level[4]){int s,t;for(s=0;s<num/ 2;s+=1) for(t=0;t< 2;t++) sample[ 2*s+t]=level[(raw_data[s]>>( 2*t+ 4*ch)) & 0x03];} 
void float_02bit_16ch(unsigned int *raw_data,float *sample,int num,int ch,float level[4]){int s,t;for(s=0;s<num/ 1;s+=1) for(t=0;t< 1;t++) sample[ 1*s+t]=level[(raw_data[s]>>( 2*t+ 2*ch)) & 0x03];} 

void vssp32_01bit_01ch(unsigned int *raw_data,float *sample,int num,int ch,float level[4]){int s,t;for(s=0;s<num/32;s+=1) for(t=0;t<32;t++) sample[32*s+t]=level[(raw_data[s]>>( 1*t+32*ch)) & 0x01];}
void vssp32_01bit_04ch(unsigned int *raw_data,float *sample,int num,int ch,float level[4]){int s,t;for(s=0;s<num/ 8;s+=1) for(t=0;t< 8;t++) sample[ 8*s+t]=level[(raw_data[s]>>( 4*t+ 1*ch)) & 0x01];}
void vssp32_02bit_01ch(unsigned int *raw_data,float *sample,int num,int ch,float level[4]){int s,t;for(s=0;s<num/16;s+=1) for(t=0;t<16;t++) sample[16*s+t]=level[(raw_data[s]>>( 2*t+32*ch)) & 0x03];}
void vssp32_02bit_04ch(unsigned int *raw_data,float *sample,int num,int ch,float level[4]){int s,t;for(s=0;s<num/ 4;s+=1) for(t=0;t< 4;t++) sample[ 4*s+t]=level[(raw_data[s]>>( 8*t+ 2*ch)) & 0x03];}


