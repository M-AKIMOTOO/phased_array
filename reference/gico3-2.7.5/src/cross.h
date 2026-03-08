#include <cutil.h>
__global__ void   cross1x1(float2 *A,          float2 *vis,int fft,int batch);
__global__ void   cross2x2(float2 *A,          float2 *vis,int fft,int batch);
__global__ void   cross3x3(float2 *A,          float2 *vis,int fft,int batch);
__global__ void   cross4x4(float2 *A,          float2 *vis,int fft,int batch);
__global__ void   cross5x5(float2 *A,          float2 *vis,int fft,int batch);
__global__ void   cross6x6(float2 *A,          float2 *vis,int fft,int batch);
__global__ void   cross7x7(float2 *A,          float2 *vis,int fft,int batch);
__global__ void   cross8x8(float2 *A,          float2 *vis,int fft,int batch);
__global__ void  matrix1x1(float2 *A,float2 *B,float2 *vis,int fft,int batch);
__global__ void  matrix2x2(float2 *A,float2 *B,float2 *vis,int fft,int batch);
__global__ void  matrix3x3(float2 *A,float2 *B,float2 *vis,int fft,int batch);
__global__ void  matrix4x4(float2 *A,float2 *B,float2 *vis,int fft,int batch);
__global__ void  matrix5x5(float2 *A,float2 *B,float2 *vis,int fft,int batch);
__global__ void  matrix6x6(float2 *A,float2 *B,float2 *vis,int fft,int batch);
__global__ void  matrix7x7(float2 *A,float2 *B,float2 *vis,int fft,int batch);
__global__ void  matrix8x8(float2 *A,float2 *B,float2 *vis,int fft,int batch);
