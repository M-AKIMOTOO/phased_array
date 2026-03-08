#define _XOPEN_SOURCE 600
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/time.h>
#include <errno.h>

#include "GICO.h"
#include "GICO3.h"
#include "param.h"
#include "apri.h"
#include "read.h"
#include "common.h"

int WRITE()
{
  FILE *fp;Sector sector;int ant1,ant2,ret;char time_code[32],file_name[256];
  time_t  timetmp = process.epoch+(time_t)process.skip; // add by KT 2015 11 jun
  strftime(time_code,sizeof(time_code),"%Y%j%H%M%S",gmtime(&timetmp));memset(&sector,0,sizeof(sector));
  for(ant1=0;ant1<STATION;ant1++) for(ant2=ant1;ant2<STATION;ant2++){
    if((strlen(unit[ant1])==0)&&(strlen(unit[ant2])==0)) sprintf(file_name,"%s/%s_%s_%s_%s.cor"      ,COR_FILE,station[ant1].name,station[ant2].name,time_code,stream[STREAM_INDEX].label                       );
    if((strlen(unit[ant1])!=0)&&(strlen(unit[ant2])!=0)) sprintf(file_name,"%s/%s_%s_%s_%s_%s-%s.cor",COR_FILE,station[ant1].name,station[ant2].name,time_code,stream[STREAM_INDEX].label,unit[ant1], unit[ant2]);
    if((strlen(unit[ant1])!=0)&&(strlen(unit[ant2])==0)) sprintf(file_name,"%s/%s_%s_%s_%s_%s-%s.cor",COR_FILE,station[ant1].name,station[ant2].name,time_code,stream[STREAM_INDEX].label,"#"       , unit[ant2]);
    if((strlen(unit[ant1])==0)&&(strlen(unit[ant2])!=0)) sprintf(file_name,"%s/%s_%s_%s_%s_%s-%s.cor",COR_FILE,station[ant1].name,station[ant2].name,time_code,stream[STREAM_INDEX].label,unit[ant1],"#"        );
    fp =fopen(file_name, "a");if(fp==NULL){fprintf(stdout,"GICO3 : [ERROR]  fopen('%s') -> %s <%s:%d>\n",file_name,strerror(errno),__FILE__,__LINE__);exit(0);}
    sector.st_time=CURRENT_EPOCH;sector.st_nsec=(1000000000LL*(PP+0))/stream[STREAM_INDEX].Hz;
    sector.et_time=CURRENT_EPOCH;sector.et_nsec=(1000000000LL*(PP+1))/stream[STREAM_INDEX].Hz;    
    //sector.amp[0]=           hypot(pcal[ant1][special[ant1][STREAM_INDEX].channel][1],pcal[ant1][special[ant1][STREAM_INDEX].channel][0]);
    //sector.phs[0]=18000/M_PI*atan2(pcal[ant1][special[ant1][STREAM_INDEX].channel][1],pcal[ant1][special[ant1][STREAM_INDEX].channel][0]);
    //sector.amp[1]=           hypot(pcal[ant2][special[ant2][STREAM_INDEX].channel][1],pcal[ant2][special[ant2][STREAM_INDEX].channel][0]);
    //sector.phs[1]=18000/M_PI*atan2(pcal[ant2][special[ant2][STREAM_INDEX].channel][1],pcal[ant2][special[ant2][STREAM_INDEX].channel][0]);
    if(fixdelay==0){    sector.geo_data[0]=gdata[ant1];    sector.geo_data[1]=gdata[ant2]; }
    sector.effective=(double)(Max(cpu,gpu)*BLOCK*stream[STREAM_INDEX].points)/terminal[0].sps;
    ret=fwrite(&sector,sizeof(sector),1,fp);if(ret!=1) {fprintf(stdout,"GICO3 : [ERROR]  write('%s') -> %s <%s:%d>\n",file_name,strerror(errno),__FILE__,__LINE__);exit(0);}
    ret=fwrite(total[ant1][ant2],sizeof(fftwf_complex)*stream[STREAM_INDEX].points/2,1,fp);if(ret!=1) {fprintf(stdout,"GICO3 : [ERROR]  write('%s') -> %s <%s:%d>\n",file_name,strerror(errno),__FILE__,__LINE__);exit(0);}

    fclose(fp);
  }
  return 0;
}

void OPEN()
{
  FILE *fp;Header header;int station_index,ant1,ant2,sps,ret;char time_code[32],file_name[256];
  time_t  timetmp = process.epoch+(time_t)process.skip; // add by KT 2015 11 jun
memset(&header,0,sizeof(header));
  header.header_magic    =__GICO_HEADER_MAGIC__  ;
  header.header_version  =__GICO_HEADER_VERSION__;
  //  header.scan_index      =PROCESS_INDEX;
  header.scan_index      =stream[0].Hz;
  
  strftime(time_code,sizeof(time_code),"%Y%j%H%M%S",gmtime(&timetmp ));
  sps=0;for(station_index=0;station_index<STATION;station_index++) sps=Max(sps,terminal[station_index].sps);
  for(ant1=0;ant1<STATION;ant1++) for(ant2=ant1;ant2<STATION;ant2++){
    for(STREAM_INDEX=0;STREAM_INDEX<STREAM;STREAM_INDEX++){
      header.sample_speed=sps;
      header.points=stream[STREAM_INDEX].points   ;
      header.frequency=stream[STREAM_INDEX].frequency;
      header.sector_number=stream[STREAM_INDEX].Hz*(process.length-process.skip); // add process.skip by KT 2015 11 jun
      header.station[0]=station[ant1];
      header.station[1]=station[ant2];
      header.source=source[STREAM_INDEX];
      header.clk_data[0]=cdata[ant1];
      header.clk_data[1]=cdata[ant2];
      if((strlen(unit[ant1])==0)&&(strlen(unit[ant2])==0)) sprintf(file_name,"%s/%s_%s_%s_%s.cor"      ,COR_FILE,station[ant1].name,station[ant2].name,time_code,stream[STREAM_INDEX].label                       );
      if((strlen(unit[ant1])!=0)&&(strlen(unit[ant2])!=0)) sprintf(file_name,"%s/%s_%s_%s_%s_%s-%s.cor",COR_FILE,station[ant1].name,station[ant2].name,time_code,stream[STREAM_INDEX].label,unit[ant1], unit[ant2]);
      if((strlen(unit[ant1])!=0)&&(strlen(unit[ant2])==0)) sprintf(file_name,"%s/%s_%s_%s_%s_%s-%s.cor",COR_FILE,station[ant1].name,station[ant2].name,time_code,stream[STREAM_INDEX].label,"#"       , unit[ant2]);
      if((strlen(unit[ant1])==0)&&(strlen(unit[ant2])!=0)) sprintf(file_name,"%s/%s_%s_%s_%s_%s-%s.cor",COR_FILE,station[ant1].name,station[ant2].name,time_code,stream[STREAM_INDEX].label,unit[ant1],"#"        );
      fp =fopen(file_name, "w");
      if(fp==NULL){fprintf(stdout,"GICO3 : [ERROR]  fopen('%s') -> %s <%s:%d>\n",file_name,strerror(errno),__FILE__,__LINE__);exit(0);}
      ret=fwrite(&header,sizeof(header),1,fp);if(ret!=1)  {fprintf(stdout,"GICO3 : [ERROR] fwrite('%s') -> %s <%s:%d>\n",file_name,strerror(errno),__FILE__,__LINE__);exit(0);}
      fclose(fp);
    }
  }
  return;
}


void GeoInit()
{
  int station_index,stream_index;FILE *fp;struct stat file_stat;char file_name[256],time_code[32];
  strftime(time_code,sizeof(time_code),"%Y/%j %H:%M:%S",gmtime(&process.epoch));
  for(station_index=0;station_index<STATION;station_index++) for( stream_index=0;stream_index<STREAM;stream_index++){
    sprintf(file_name,"%s/%s_%s_%s.geo",GEO_FILE,station[station_index].name,      time_code,stream[stream_index].label);if(stat(file_name,&file_stat)==0) goto skip;
    sprintf(file_name,"%s/%s_%s_%s.geo",GEO_FILE,station[station_index].name,      time_code,source[stream_index].name );if(stat(file_name,&file_stat)==0) goto skip;
    sprintf(file_name,"%s/%s_%s_%s.geo",GEO_FILE,station[station_index].name,"YYYYDDDHHMMSS",stream[stream_index].label);if(stat(file_name,&file_stat)==0) goto skip;
    sprintf(file_name,"%s/%s_%s_%s.geo",GEO_FILE,station[station_index].name,"YYYYDDDHHMMSS",source[stream_index].name );if(stat(file_name,&file_stat)==0) goto skip;
    sprintf(file_name,"%s/%s_%s.geo"   ,GEO_FILE,station[station_index].name,      time_code                           );if(stat(file_name,&file_stat)==0) goto skip;
    sprintf(file_name,"%s/%s_%s.geo"   ,GEO_FILE,station[station_index].name,"YYYYDDDHHMMSS"                           );if(stat(file_name,&file_stat)==0) goto skip;
    fprintf(stdout,"GICO3 : [ERROR] fopen('%s') -> %s <%s:%d>\n",file_name,strerror(errno),__FILE__,__LINE__);goto next;
  skip:
    if((fp=fopen(file_name,"r"))==NULL) {fprintf(stdout,"GICO3 : [ERROR]  fopen('%s') -> %s <%s:%d>\n",file_name ,strerror(errno),__FILE__,__LINE__);goto next;}
    geo_data[station_index][stream_index]=malloc(file_stat.st_size+64);fread(geo_data[station_index][stream_index],file_stat.st_size,1,fp);geo_length[station_index][stream_index]=file_stat.st_size/sizeof(Clock);fclose(fp);
  }  
  return;
 next:
  exit(0);
}

void GeoFree()
{
  int station_index,stream_index;
  for(station_index=0;station_index<STATION;station_index++) for( stream_index=0;stream_index<STREAM;stream_index++) free(geo_data[station_index][stream_index]);
}

Clock find(int station_index,int stream_index,time_t epoch)
{
  int s,s0,s1;Clock gdata;char time_code[32];s0=0;s1=geo_length[station_index][stream_index]-1;
  if(geo_data[station_index][stream_index][s0].sec+0.000000001*geo_data[station_index][stream_index][s0].nsec>epoch) return geo_data[station_index][stream_index][s0];
  if(geo_data[station_index][stream_index][s1].sec+0.000000001*geo_data[station_index][stream_index][s1].nsec<epoch) return geo_data[station_index][stream_index][s1];
  while(s1-s0>1){
    s=(s1+s0)/2;
    if(geo_data[station_index][stream_index][s].sec+0.000000001*geo_data[station_index][stream_index][s].nsec>=epoch) s1=s;
    if(geo_data[station_index][stream_index][s].sec+0.000000001*geo_data[station_index][stream_index][s].nsec<=epoch) s0=s;
  }
  if(fabs(geo_data[station_index][stream_index][s0].sec+0.000000001*geo_data[station_index][stream_index][s0].nsec-epoch)<=fabs(geo_data[station_index][stream_index][s1].sec+0.000001*geo_data[station_index][stream_index][s1].nsec-epoch)) s=s0;
  if(fabs(geo_data[station_index][stream_index][s0].sec+0.000000001*geo_data[station_index][stream_index][s0].nsec-epoch)>=fabs(geo_data[station_index][stream_index][s1].sec+0.000001*geo_data[station_index][stream_index][s1].nsec-epoch)) s=s1;
  gdata=geo_data[station_index][stream_index][s];epoch=gdata.sec;strftime(time_code,sizeof(time_code),"%Y/%j %H:%M:%S",gmtime(&epoch));
  if(message>=2) fprintf(stdout,"GICO3 : [level-2] apri[%8s@%s.%09d %+.7e,%+.7e,%+.7e,%+.7e,%+.7e]\n",station[station_index].name,time_code,gdata.nsec,gdata.delay,gdata.rate,gdata.acel,gdata.jerk,gdata.snap);
  return geo_data[station_index][stream_index][s];
}

void GICO(Process process)
{
  struct timeval s0,s1,s2;double t1,t2,busy=0,idle=0;
  int index,point,ant1,ant2;char time_code[32];
  int dev;time_t epoch;
  pthread_t read_thread,pthread[MAX_DEVICE];
  point=0        ;for(STREAM_INDEX=0;STREAM_INDEX<STREAM;STREAM_INDEX++) point=Max(stream[STREAM_INDEX].points,point);
  for(dev=0;dev<MAX_DEVICE;dev++)                                    for(index=0;index<STATION;index++) posix_memalign((void *)&prm3[dev][index]             ,512,sizeof(Prm3)*1000000+64);
  for(dev=0;dev<MAX_DEVICE;dev++)                                    for(index=0;index<STATION;index++) posix_memalign((void *)&prm4[dev][index]             ,512,sizeof(Prm4)*1000000+64);
  for(ant1=0;ant1<STATION;ant1++) for(ant2=ant1;ant2<STATION;ant2++)                                    posix_memalign((void *)&total[ant1][ant2]            ,512,sizeof(float)*point +64);
  if(strlen(GEO_FILE)!=0) GeoInit();
  
  //CURRENT_EPOCH=process.epoch;READ_EPOCH=process.epoch;pthread_create(&read_thread,NULL,READ1SEC, process.length);
  if(elecs){
    CURRENT_EPOCH=process.epoch+process.skip;READ_EPOCH=process.epoch;pthread_create(&read_thread,NULL,Read1SecElecs, NULL);
  }else if(vdif){
    CURRENT_EPOCH=process.epoch+process.skip;READ_EPOCH=process.epoch;pthread_create(&read_thread,NULL,Read1SecVDIF, NULL);
  }else{
 CURRENT_EPOCH=process.epoch+process.skip;READ_EPOCH=process.epoch;pthread_create(&read_thread,NULL,Read1Sec, NULL);
  }
  for(CURRENT_EPOCH=process.epoch+process.skip;CURRENT_EPOCH<process.epoch+process.length;CURRENT_EPOCH++){
    gettimeofday(&s0,NULL);while(READ_EPOCH<=CURRENT_EPOCH) usleep(100000);gettimeofday(&s1,NULL);
    for(STREAM_INDEX=0;STREAM_INDEX<STREAM;STREAM_INDEX++) for(PP=0;PP<stream[STREAM_INDEX].Hz;PP++){
	BLOCK=terminal[0].sps/(stream[STREAM_INDEX].points*stream[STREAM_INDEX].Hz*Max(cpu,gpu));
	double pp = PP ; double hz = stream[STREAM_INDEX].Hz;
	double msec = pp/hz; //printf("msec : %6.1lf",msec);
	
	// fprintf(stdout,"cpu%f, %d,%d,%d,%d,%d]\n", BLOCK,cpu,terminal[0].sps,stream[STREAM_INDEX].points,stream[STREAM_INDEX].Hz,PP);
	if(BLOCK < 1 ){        printf("Error Block is too small! reduce the cpu point \n"); exit(1); }
	for(ant1=0;ant1<STATION;ant1++) for(ant2=ant1;ant2<STATION;ant2++)  memset(total[ant1][ant2],0,sizeof(float)*point);
	for(index=0;index<STATION;index++) if(strlen(GEO_FILE)==0) gdata[index]=apri(station[index],source[STREAM_INDEX],CURRENT_EPOCH, msec);
	for(index=0;index<STATION;index++) if(strlen(GEO_FILE)!=0) gdata[index]=find(         index,       STREAM_INDEX ,CURRENT_EPOCH);
	for(dev=0;dev<cpu;dev++) pthread_create(&pthread[dev],NULL,(void *)PARAM3,(void *)dev);for(dev=0;dev<cpu;dev++) pthread_join(pthread[dev],NULL);
	for(dev=0;dev<cpu;dev++) pthread_create(&pthread[dev],NULL,(void *)GICO3 ,(void *)dev);for(dev=0;dev<cpu;dev++) pthread_join(pthread[dev],NULL);
	WRITE();
      } 
    gettimeofday(&s2,NULL);epoch=CURRENT_EPOCH;strftime(time_code,sizeof(time_code),"%Y/%j %H:%M:%S",gmtime(&epoch));
    t1=(s1.tv_sec-s0.tv_sec)+0.000001*(s1.tv_usec-s0.tv_usec);t2=(s2.tv_sec-s1.tv_sec)+0.000001*(s2.tv_usec-s1.tv_usec);busy+=t2;idle+=t1;
    if(message>=1) fprintf(stdout,"GICO3 : [level-1] 観測データ[%s,%04d<sec>]の処理に要した時間 %6.3lf秒 <計算=%6.3lf秒 待機=%6.3lf秒>\n",time_code,1,t1+t2,t2,t1);fflush(stdout);
  }
  pthread_join(read_thread,NULL);
  if(strlen(GEO_FILE)!=0) GeoFree();
  time_t  timetmp = process.epoch+(time_t)process.skip; // add by KT 2015 11 jun
  strftime(time_code,sizeof(time_code),"%Y/%j %H:%M:%S",gmtime(&timetmp));
  fprintf(stdout,"GICO3 : [LEVEL-0] 観測データ[%s,%04d<sec>]の処理に要した時間 %6.1lf秒 <計算=%6.1lf秒 待機=%6.1lf秒>\n",time_code,(process.length-process.skip),busy+idle,busy,idle);fflush(stdout);
  for(dev=0;dev<MAX_DEVICE;dev++) for(index=0;index<STATION;index++) free( prm3[dev][index]);
  for(dev=0;dev<MAX_DEVICE;dev++) for(index=0;index<STATION;index++) free( prm4[dev][index]);
  for(ant1=0;ant1<STATION;ant1++) for(ant2=ant1;ant2<STATION;ant2++) free(total[ant1][ant2]);
}
