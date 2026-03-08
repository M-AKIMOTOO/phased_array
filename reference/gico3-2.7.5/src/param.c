#define _XOPEN_SOURCE 600
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "param.h"
#include "config.h"
#include "common.h"

inline long double clock_delay(Clock *geo_data,time_t sec,long double tim)
{
  long double dt1,dt2,dt3,dt4;
  dt1=((sec)-(geo_data->sec))+((tim)-(0.000000001*geo_data->nsec));dt2=dt1*dt1/2;dt3=dt2*dt1/3;dt4=dt3*dt1/4;
  return (geo_data->delay)+dt1*(geo_data->rate)+dt2*(geo_data->acel)+dt3*(geo_data->jerk)+dt4*(geo_data->snap);
}
inline long double clock_rate(Clock *geo_data,time_t sec,long double tim)
{
  long double dt1,dt2,dt3,dt4;
  dt1=((sec)-(geo_data->sec))+((tim)-(0.000000001*geo_data->nsec));dt2=dt1*dt1/2;dt3=dt2*dt1/3;dt4=dt3*dt1/4;
  return (geo_data->rate)+dt1*(geo_data->acel)+dt2*(geo_data->jerk)+dt3*(geo_data->snap);
}

Prm3 param3(Clock *gdata,Clock *cdata,Terminal *terminal,Stream *stream,double rotation,time_t epoch,double tim)
{
  Prm3 prm3;long double delay,rate,sample,sample0,sample1,phs0,phs1,phase0,phase1,phase2,phase3;int spw;
  spw=32/((terminal->bit)*(terminal->channel));
  if(fixdelay==0){ 
    delay =clock_delay(gdata,epoch,tim)+clock_delay(cdata,epoch,tim);
    rate  =clock_rate (gdata,epoch,tim)+clock_rate (cdata,epoch,tim);
  }else{// added by KT for zenith observation between kas-ishioka 2015aug27
    delay =clock_delay(cdata,epoch,tim);//clock_delay(gdata,epoch,tim)
    rate  =clock_rate (cdata,epoch,tim);//clock_rate (gdata,epoch,tim)+clock_rate (cdata,epoch,tim);
  }
  
  sample=terminal->sps*(tim+delay);sample0=spw*rintl(sample/spw);sample1=sample-sample0;prm3.offset=sample0/spw;prm3.phase3=0;
  phs0=rotation*epoch;phs0-=floorl(phs0);phs0+=rotation*sample0/terminal->sps;phs1=rotation/terminal->sps;
  //printf("Test %lf %Lf %lf sample=%.16Lf sample0=%Lf sample1=%Lf\n",tim,phs0,rotation,sample,sample0,sample1);
  phase0=stream->frequency*delay              +phs0;prm3.phase0=pow(2,32)*(phase0-floorl(phase0));
  //NG lost fringe : attempt by KT  phase1=stream->frequency*rate +phs1;prm3.phase1=pow(2,32)*(phase1-floorl(phase1));
  if(message>=4)  printf("Test1 phase0[cycle] %Lf freq %lf gdata[us] %Lf cdata[us] %LF rate[Hz] %.16Lf %f \n",phase0, stream->frequency, clock_delay(gdata,epoch,tim)*1e6,clock_delay(cdata,epoch,tim)*1e6, rate,rotation);

  //original 
  phase1=stream->frequency*rate /terminal->sps+phs1;prm3.phase1=pow(2,32)*(phase1-floorl(phase1)); // phase velocity per 1 sample unit
  phase2=sample1/stream->points                    ;prm3.phase2=pow(2,32)*(phase2-floorl(phase2));

  return prm3;
}

void PARAM3(int dev)
{
  unsigned int random;int index,page,block,unit;Prm3 v;struct tm tm;char time_code[32];long double tim2;double tim;
  gmtime_r(&CURRENT_EPOCH,&tm);strftime(time_code,sizeof(time_code),"%Y/%j %H:%M:%S",&tm);
  page=CURRENT_EPOCH%MAX_BUFFER;
  for(block=0;block<BLOCK;block++){
    tim=(PP+(double)dev/cpu+(double)block/(BLOCK*cpu))/stream[STREAM_INDEX].Hz;if(pcal==0) random=genrand_int32();else random=0;
    for(index=0;index<STATION;index++){
      prm3[dev][index][block]=param3(&gdata[index],&cdata[index],&terminal[index],&stream[STREAM_INDEX],special[index][STREAM_INDEX].rotation,CURRENT_EPOCH,tim);
      //prm3[dev][index][block].phase2+=random;
      unit=stream[STREAM_INDEX].points/(32/(terminal[index].channel*terminal[index].bit));
      if(prm3[dev][index][block].offset+0*unit<raw_offset[page][index]+0*raw_length[page][index]) prm3[dev][index][block].offset=0x80000000;
      if(prm3[dev][index][block].offset+1*unit>raw_offset[page][index]+1*raw_length[page][index]) prm3[dev][index][block].offset=0x80000000;
      if(message>=3){v=prm3[dev][index][block];fprintf(stdout,"GICO3 : [level-2] param3[%08s][%s.%09d]={%08x,%08x,%08x,%08x} %08x %08x\n",
						       station[index].name,time_code,(int)rint(1000000000*tim),v.offset,v.phase0,v.phase1,v.phase2,raw_offset[page][index],raw_length[page][index]);}
    }
  }
  for(index=0;index<STATION;index++) prm3[dev][index][block].offset=0x7FFFFFFF;
}
  
