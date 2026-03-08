#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <unistd.h>
#include <pthread.h>
#include <time.h>
#include <sys/time.h>
#include <xmmintrin.h>
#include "header.h"
#include "schedule.h"
#include "GICO3.h"

#include "tool.h"

void *func_float(Terminal *terminal)
{
  if(strstr(terminal->name,"VSSP32")!=NULL){
    if((terminal->bit==1) && (terminal->channel== 1)) return vssp32_01bit_01ch;if((terminal->bit==2) && (terminal->channel== 1)) return vssp32_02bit_01ch;
    if((terminal->bit==1) && (terminal->channel== 4)) return vssp32_01bit_04ch;if((terminal->bit==2) && (terminal->channel== 4)) return vssp32_02bit_04ch;
  }
  if((terminal->bit==1) && (terminal->channel== 1)) return float_01bit_01ch;if((terminal->bit==2) && (terminal->channel== 1)) return float_02bit_01ch;
  if((terminal->bit==1) && (terminal->channel== 2)) return float_01bit_02ch;if((terminal->bit==2) && (terminal->channel== 2)) return float_02bit_02ch;
  if((terminal->bit==1) && (terminal->channel== 4)) return float_01bit_04ch;if((terminal->bit==2) && (terminal->channel== 4)) return float_02bit_04ch;
  if((terminal->bit==1) && (terminal->channel== 8)) return float_01bit_08ch;if((terminal->bit==2) && (terminal->channel== 8)) return float_02bit_08ch;
  if((terminal->bit==1) && (terminal->channel==16)) return float_01bit_16ch;if((terminal->bit==2) && (terminal->channel==16)) return float_02bit_16ch;
// Added by KT 19st Jan, 2015
  if((terminal->bit==1) && (terminal->channel==32)) return float_01bit_32ch;// if((terminal->bit==2) && (terminal->channel==32)) return float_02bit_32ch;

  fprintf(stdout,"GICO3 :  unsupported terminal mode such as bit=%d channel=%d <%s:%d>\n",terminal->bit,terminal->channel,__FILE__,__LINE__);exit(0);
}

void GICO3(int dev)
{
  static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
  fftwf_plan fft[MAX_STATION];void (*func[MAX_STATION])(unsigned int *raw_data,float *sample,int num,int ch,float level[4]);
  int station_index,index,s,ant1,ant2,block,points;int  offset;unsigned int phase0,phase1,phase2,phs;double w;
  float *sample,*temp;fftwf_complex *fr,*freq[MAX_STATION],*vis[MAX_STATION][MAX_STATION];
  int page;

  points=stream[STREAM_INDEX].points;page=CURRENT_EPOCH%MAX_BUFFER;
  //  fprintf (stderr,"points : %d", points);

  temp  =fftwf_malloc(8*points+32);                                  if(                                           temp==NULL) {fprintf(stderr,"GICO3 : fftwf_malloc('%d') <%s:%d>",8*points,__FILE__,__LINE__);exit(0);}
  sample=fftwf_malloc(8*points+32);                                  if(                                         sample==NULL) {fprintf(stderr,"GICO3 : fftwf_malloc('%d') <%s:%d>",8*points,__FILE__,__LINE__);exit(0);}
  fr    =fftwf_malloc(8*points+32);                                  if(                                             fr==NULL) {fprintf(stderr,"GICO3 : fftwf_malloc('%d') <%s:%d>",8*points,__FILE__,__LINE__);exit(0);}
  for(station_index=0;station_index<STATION;station_index++)         if((freq[station_index]=fftwf_malloc(8*points+32))==NULL) {fprintf(stderr,"GICO3 : fftwf_malloc('%d') <%s:%d>",8*points,__FILE__,__LINE__);exit(0);}
  for(ant1=0;ant1<STATION;ant1++) for(ant2=ant1;ant2<STATION;ant2++) if((    vis[ant1][ant2]=fftwf_malloc(8*points+32))==NULL) {fprintf(stderr,"GICO3 : fftwf_malloc('%d') <%s:%d>",8*points,__FILE__,__LINE__);exit(0);}
  for(ant1=0;ant1<STATION;ant1++) for(ant2=ant1;ant2<STATION;ant2++) memset(vis[ant1][ant2],0,8*points);
  if(complex==0) {pthread_mutex_lock(&mutex);for(station_index=0;station_index<STATION;station_index++) fft[station_index]=fftwf_plan_dft_r2c_1d(points,(        float *)sample,fr,             FFTW_ESTIMATE);pthread_mutex_unlock(&mutex);}
  if(complex==1) {pthread_mutex_lock(&mutex);for(station_index=0;station_index<STATION;station_index++) fft[station_index]=fftwf_plan_dft_1d    (points,(fftwf_complex *)sample,fr,FFTW_FORWARD,FFTW_ESTIMATE);pthread_mutex_unlock(&mutex);}
  if(iq     ==1) {pthread_mutex_lock(&mutex);for(station_index=0;station_index<STATION;station_index++) fft[station_index]=fftwf_plan_dft_1d    (points,(fftwf_complex *)sample,fr,FFTW_FORWARD,FFTW_ESTIMATE);pthread_mutex_unlock(&mutex);}
  for(station_index=0;station_index<STATION;station_index++) func[station_index]=func_float(&terminal[station_index]);
  for(block=0;prm3[dev][0][block].offset!=0x7FFFFFFF;block++){
    for(station_index=0;station_index<STATION;station_index++){
      offset=prm3[dev][station_index][block].offset;phase0=prm3[dev][station_index][block].phase0;phase1=prm3[dev][station_index][block].phase1;phase2=prm3[dev][station_index][block].phase2;      //      printf("offset=%08x raw_offset=%08x\n",offset,raw_offset[page][station_index]);
      //  printf("sizeof(fftwf) %d,%d\n",sizeof(fftwf_complex),sizeof(float)); = 8,4
      if(offset==0x80000000) memset(sample,0,sizeof(fftwf_complex)*points);
      else{	func[station_index](&raw_data[page][station_index][offset-raw_offset[page][station_index]],sample,points,special[station_index][STREAM_INDEX].channel-1,terminal[station_index].level);      }
      if(special[station_index][STREAM_INDEX].sideband==LSB ){
	//      if(special[station_index][STREAM_INDEX].sideband==LSB || special[station_index][STREAM_INDEX].sideband==IQ){
	if(terminal[station_index].bit*terminal[station_index].channel!=32) for(index=       0;index<points;index+=2) sample[index]=-sample[index];
	else                                                                for(index=offset&1;index<points;index+=2) sample[index]=-sample[index];
      }
      if(iq==1 ){ // for complex I/Q
	if( special[station_index][STREAM_INDEX].sideband==IQ){
	  for(index=0;index<points;index+=2 ){temp[2*index+0]=sample[index];temp[2*index+1]=sample[index+1];temp[2*index+2]=0.0;temp[2*index+3]=0.0;} memcpy(sample,temp,sizeof(fftwf_complex)*points);//printf("haystack\n");
	}else{ // for real 
	  for(index=0;index<points;index++ ){temp[2*index+0]=sample[index];  temp[2*index+1]=0.0;                                  } 
	  for(index=0 ; index<points; index++){ /*/printf("index: %d,%f,%f\n",index,sample[index],temp[index]); */}
	  memcpy(sample,temp,sizeof(fftwf_complex)*points);//printf("Kashimax\n");
	}//ifelse
      }//if IQ
      //printf("phs: %f\n",phs) ;
      if(complex==1){	phs=phase0;for(index=0;index<points;index++){ temp[2*index+0]=sample[index]*cos_table[phs>>24];temp[2*index+1]=sample[index]*sin_table[phs>>24];phs+=phase1;} 	memcpy(sample,temp,sizeof(fftwf_complex)*points);      }
      fftwf_execute(fft[station_index]);
      if(complex==0) gico3_track(fr,freq[station_index],phase0,phase2,points/2);
      else gico3_track(fr,freq[station_index],     0,phase2,points/2);

    }
    if(pcal==0) for(ant1=0;ant1<STATION;ant1++) for(ant2=ant1  ;ant2<STATION;ant2++) {gico3_cross((__v4sf *)freq[ant1],(__v4sf *)freq[ant2],(__v4sf *)vis[ant1][ant2],points/2   );}
    if(pcal==1) for(ant1=0;ant1<STATION;ant1++) for(ant2=ant1+1;ant2<STATION;ant2++) {gico3_cross((__v4sf *)freq[ant1],(__v4sf *)freq[ant2],(__v4sf *)vis[ant1][ant2],points/2);}
    //modified by KT  
  if(pcal==1) for(ant1=0;ant1<STATION;ant1++) for(ant2=ant1  ;ant2==ant1  ;ant2++) {gico3_tone ((__v4sf *)freq[ant1],                     (__v4sf *)vis[ant1][ant2],points/2);}
  
 }
  pthread_mutex_lock(&mutex);

  w=1/(0.5*1.25*pow(points,2.0)*BLOCK*cpu);
  //   printf("%f,%f,%d]\n", w,BLOCK,cpu);
    for(ant1=0;ant1<STATION;ant1++) for(ant2=ant1;ant2<STATION;ant2++) for(s=0;s<points/2  ;s++){total[ant1][ant2][s][0]+=vis[ant1][ant2][s][0]*w;total[ant1][ant2][s][1]+=vis[ant1][ant2][s][1]*w;}
  //  for(ant1=0;ant1<STATION;ant1++) for(ant2=ant1;ant2<STATION;ant2++) for(s=0;s<points  ;s++){printf("s:%d,%d, %d,%d,%f,%f\n",ant1,ant2,points,s,vis[ant1][ant2][s][0],vis[ant1][ant2][s][1]);}
  //  for(ant1=0;ant1<STATION;ant1++) for(ant2=ant1;ant2<STATION;ant2++) for(s=0;s<points/2  ;s++){total[ant1][ant2][s][0]+=(vis[ant1][ant2][2*s][0]+vis[ant1][ant2][2*s+1][0])*w;total[ant1][ant2][s][1]+=(vis[ant1][ant2][2*s][1]+vis[ant1][ant2][2*s+1][1])*w;}
  fftwf_free(temp);fftwf_free(sample);fftwf_free(fr);for(ant1=0;ant1<STATION;ant1++){fftwf_destroy_plan(fft[ant1]);fftwf_free(freq[ant1]);for(ant2=ant1;ant2<STATION;ant2++) fftwf_free(vis[ant1][ant2]);}
  pthread_mutex_unlock(&mutex);
  return;
}
