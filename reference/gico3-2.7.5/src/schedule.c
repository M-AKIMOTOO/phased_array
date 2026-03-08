#define _XOPEN_SOURCE 
#include "common.h"
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>
#include <errno.h>
#include <math.h>
#include <time.h>

#include <libxml/tree.h>
#include <libxml/parser.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>

#include "schedule.h"
#include "config.h"
#include "header.h"

typedef enum {EPOCH,LENGTH,SKIP,OBJECT,STATIONS,NAME,TERMINAL,UNIT,POS_X,POS_Y,POS_Z,SHUFFLE,RAW_NAME,SPEED,CHANNEL,BIT,LEVEL,LABEL,SOURCE,FREQUENCY,FFT_POINT,OUTPUT,RA,DEC,DELAY,RATE,ACEL,JERK,SNAP,OFFSET,SIDEBAND,ROTATION};

#define RA_to_Radian(      hh,mm,ss)  (((ra_hh)+(ra_mm)/60.0+(ra_ss)/3600.0)/12.0*M_PI)
#define DEC_to_Radian(sign,dd,mm,ss)  ((sign=='+') ? (+((dec_dd)+(dec_mm)/60.0+(dec_ss)/3600.0)/180.0*M_PI) : (-((dec_dd)+(dec_mm)/60.0+(dec_ss)/3600.0)/180.0*M_PI))

int lookup_process(xmlXPathContext *context,const char *exp,int pp,char *ret)
{
  char key[256];int num=0;xmlXPathObject *res=NULL;
  if(ret!=NULL) strcpy(ret,"");
  sprintf(key,"/schedule/process[%d]/%s",pp,exp);res=xmlXPathEvalExpression((const xmlChar *)key,context);
  if(xmlXPathNodeSetIsEmpty(res->nodesetval)==0){num=res->nodesetval->nodeNr;if(ret) sscanf((char *)res->nodesetval->nodeTab[0]->children->content, " %[^\n]",ret);goto next;}
  sprintf(key,"/schedule/%s"               ,exp);res=xmlXPathEvalExpression((const xmlChar *)key,context);
  if(xmlXPathNodeSetIsEmpty(res->nodesetval)==0){num=res->nodesetval->nodeNr;if(ret) sscanf((char *)res->nodesetval->nodeTab[0]->children->content, " %[^\n]",ret);goto next;}
 next:
  if(res) xmlXPathFreeObject(res);
  if(ret!=NULL) while(ret[strlen(ret)-1]==' ') ret[strlen(ret)-1]=0;
  return num;
}

int lookup_context(xmlXPathContext *context,const char *exp,char *ret)
{
  int num=0;xmlXPathObject *res=NULL;
  if(ret!=NULL) strcpy(ret,"");
  res=xmlXPathEvalExpression((const xmlChar *)exp,context);if(res==NULL) goto skip;
  if(xmlXPathNodeSetIsEmpty(res->nodesetval)==0) num=res->nodesetval->nodeNr;if(num==0) goto skip;
  if((res->nodesetval->nodeTab[0]->children->content)==NULL) goto skip;
  if(ret!=NULL) sscanf((char *)res->nodesetval->nodeTab[0]->children->content," %[^\n] ",ret);
 skip:
  if(ret!=NULL) while(ret[strlen(ret)-1]==' ') ret[strlen(ret)-1]=0;
  if(res) xmlXPathFreeObject(res);
  return num;
}

void put_process(Process *process)
{
  time_t epoch;char time_code[32];int station_index,stream_index;
  strftime(time_code,sizeof(time_code),"%Y/%j %H:%M:%S",gmtime(&process->epoch));
  printf("-------------------------------------------------------------------------------------------------------------------------\n");
  printf("      Process  :  epoch='%s' length=%04d skip=%04d object=%s stations=%s <ut1utc=%+.6lf cpu=%02d pcal=%d>\n",time_code,process->length,process->skip,process->object,process->stations,ut1utc,cpu,pcal);  
  for(station_index=0;(station_index<MAX_STATION)&&(station_index<strlen(process->stations));station_index++){
    Station s=station[station_index];
    printf("   Station[%c]  :  name=%-8s position=<%+15.6lf %+15.6lf %+15.6lf> unit=%s\n",process->stations[station_index],s.name,s.pos_x,s.pos_y,s.pos_z,unit[station_index]);
  }
  for(station_index=0;(station_index<MAX_STATION)&&(station_index<strlen(process->stations));station_index++){
    Clock   s=cdata[station_index];epoch=s.sec;strftime(time_code,sizeof(time_code),"%Y/%j %H:%M:%S",gmtime(&epoch));
    printf("     Clock[%c]  :  epoch='%s' delay=%+.2e rate=%+.2e acel=%+.2e jerk=%+.2e snap=%+.2e\n",process->stations[station_index],time_code,s.delay,s.rate,s.acel,s.jerk,s.snap);
  }
  for(station_index=0;(station_index<MAX_STATION)&&(station_index<strlen(process->stations));station_index++){
    Terminal s=terminal[station_index];
    if(s.bit==1) printf("  Terminal[%c]  :  name=%-8s speed=%04d[Msps] bit=%02d channel=%02d level=<%+4.1lf,%+4.1lf>                \n",process->stations[station_index],s.name,s.sps/1000000,s.bit,s.channel,s.level[0],s.level[1]                      );
    if(s.bit==2) printf("  Terminal[%c]  :  name=%-8s speed=%04d[Msps] bit=%02d channel=%02d level=<%+4.1lf,%+4.1lf,%+4.1lf,%+4.1lf>\n",process->stations[station_index],s.name,s.sps/1000000,s.bit,s.channel,s.level[0],s.level[1],s.level[2],s.level[3]);
  }
  for(station_index=0;(station_index<MAX_STATION)&&(station_index<strlen(process->stations));station_index++){
    int *b=shuffle_list[station_index];
    printf("  Shuffle[%c]   :  %02d,%02d,%02d,%02d,%02d,%02d,%02d,%02d,%02d,%02d,%02d,%02d,%02d,%02d,%02d,%02d,%02d,%02d,%02d,%02d,%02d,%02d,%02d,%02d,%02d,%02d,%02d,%02d,%02d,%02d,%02d,%02d\n",process->stations[station_index],
	   b[31],b[30],b[29],b[28],b[27],b[26],b[25],b[24],b[23],b[22],b[21],b[20],b[19],b[18],b[17],b[16],b[15],b[14],b[13],b[12],b[11],b[10],b[ 9],b[ 8],b[ 7],b[ 6],b[ 5],b[ 4],b[ 3],b[ 2],b[ 1],b[ 0]);
  }
  for(station_index=0;(station_index<MAX_STATION)&&(station_index<strlen(process->stations));station_index++) printf("  Raw-File[%c]  :  directory=%s name=%26s\n",process->stations[station_index],RAW_FILE,raw_name[station_index]);
  for(stream_index=0;(stream_index<MAX_STREAM)&&(strlen(stream[stream_index].label)!=0);stream_index++){
    Stream s=stream[stream_index];
    printf("  Stream[%03d]  :  label=%-8s source=%-8s frequency=%+13.9lf[GHz] channel=%02d fft=%06d output=%04d[Hz]\n",stream_index+1,s.label,source[stream_index].name,s.frequency/1000000000.0,s.channel,s.points,s.Hz);
    for(station_index=0;(station_index<MAX_STATION)&&(station_index<strlen(process->stations));station_index++){
      Special s=special[station_index][stream_index];if((s.offset==0)&&(s.channel==stream[stream_index].channel)&&(s.sideband==USB)&&(s.rotation==0.0)) continue;
      if(s.sideband==USB) printf("  ->Special[%c] :  offset=%+06d[point] channel=%02d sideband=USB,rotation=%+.6lf[MHz]\n",process->stations[station_index],s.offset,s.channel,s.rotation/1000000);
      if(s.sideband==LSB) printf("  ->Special[%c] :  offset=%+06d[point] channel=%02d sideband=LSB,rotation=%+.6lf[MHz]\n",process->stations[station_index],s.offset,s.channel,s.rotation/1000000);
      if(s.sideband==IQ)  printf("  ->Special[%c] :  offset=%+06d[point] channel=%02d sideband=IQ ,rotation=%+.6lf[MHz]\n",process->stations[station_index],s.offset,s.channel,s.rotation/1000000);
    }
  }
 
}

int read_process(xmlXPathContext *context,int index,Process *process,int *STATION,int *STREAM,char *remove,double add_delay) // modified KT 2013 DEC 6
{
  int *b;
  char key,exp[256],label[32][256],ans[256],time_code[32];int s,station_index,stream_index;struct tm tm;float *v;int ra_hh,ra_mm;double ra_ss;char sign;int dec_dd,dec_mm;double dec_ss;
  memset(process,0,sizeof(Process));memset(station,0,sizeof(station));memset(cdata,0,sizeof(cdata));memset(terminal,0,sizeof(terminal));memset(source,0,sizeof(source));memset(stream,0,sizeof(stream));memset(unit,0,sizeof(unit));
  sprintf(exp,"/schedule/process[%d]/epoch"   ,index);if(lookup_context(context,exp,label[EPOCH   ])!=1) {sprintf(ans,"GICO3 : [ERROR] process[%04d]/epoch    -> not assigned  <%s:%d>\n",index,__FILE__,__LINE__);goto exit;}
  sprintf(exp,"/schedule/process[%d]/length"  ,index);if(lookup_context(context,exp,label[LENGTH  ])!=1) {sprintf(ans,"GICO3 : [ERROR] process[%04d]/length   -> not assigned  <%s:%d>\n",index,__FILE__,__LINE__);goto exit;}
  sprintf(exp,"/schedule/process[%d]/skip"    ,index);if(lookup_context(context,exp,label[SKIP    ])==1) {process->skip=atoi(label[SKIP]);}// add by KT 2015 11 jun
  if(strptime(label[EPOCH],"%Y/%j %H:%M:%S",&tm)==0) {sprintf(ans,"GICO3 : process[%04d]/epoch    -> illigal format<%s:%d>\n",index,__FILE__,__LINE__);goto exit;}
  process->epoch=mktime(&tm);process->length=atoi(label[LENGTH]);
  sprintf(exp,"/schedule/process[%d]/object"  ,index);if(lookup_context(context,exp,label[OBJECT  ])!=1) {sprintf(ans,"GICO3 : [ERROR] process[%04d]/object   -> not assigned  <%s:%d>\n",index,__FILE__,__LINE__);goto exit;}
  sprintf(exp,"/schedule/process[%d]/stations",index);if(lookup_context(context,exp,label[STATIONS])!=1) {sprintf(ans,"GICO3 : [ERROR] process[%04d]/stations -> not assigned  <%s:%d>\n",index,__FILE__,__LINE__);goto exit;}
  strcpy(process->object,label[OBJECT]);strcpy(process->stations,label[STATIONS]);

  strftime(time_code,sizeof(time_code),"%Y%j%H%M%S",gmtime(&process->epoch));
  for(s=0;s<strlen(remove);s++) if(strchr(process->stations,remove[s])!=NULL) strcpy(strchr(process->stations,remove[s]),strchr(process->stations,remove[s])+1);
  for(station_index=0;(station_index<MAX_STATION)&&(key=process->stations[station_index]);station_index++){
    sprintf(exp,"station[@key='%c']"         ,key);if(lookup_process(context,exp,index,           NULL)!=1){sprintf(ans,"GICO3 : [ERROR] station[@key='%c']           -> not assigned <%s:%d>\n",key,__FILE__,__LINE__);goto exit;}
    sprintf(exp,"station[@key='%c']/name"    ,key);if(lookup_process(context,exp,index,label[    NAME])!=1){sprintf(ans,"GICO3 : [ERROR] station[@key='%c']/name      -> not assigned <%s:%d>\n",key,__FILE__,__LINE__);goto exit;}
    sprintf(exp,"station[@key='%c']/terminal",key);if(lookup_process(context,exp,index,label[TERMINAL])!=1){sprintf(ans,"GICO3 : [ERROR] station[@key='%c']/terminal  -> not assigned <%s:%d>\n",key,__FILE__,__LINE__);goto exit;}
    sprintf(exp,"station[@key='%c']/unit"    ,key);if(lookup_process(context,exp,index,label[    UNIT])==1) strcpy(unit[station_index],label[UNIT]);else unit[station_index][0]=0;
    sprintf(exp,"station[@key='%c']/pos-x"   ,key);if(lookup_process(context,exp,index,label[   POS_X])==1) station[station_index].pos_x=atof(label[POS_X]);
    sprintf(exp,"station[@key='%c']/pos-y"   ,key);if(lookup_process(context,exp,index,label[   POS_Y])==1) station[station_index].pos_y=atof(label[POS_Y]);
    sprintf(exp,"station[@key='%c']/pos-z"   ,key);if(lookup_process(context,exp,index,label[   POS_Z])==1) station[station_index].pos_z=atof(label[POS_Z]);
    sprintf(exp,"station[@key='%c']/file"    ,key);if(lookup_process(context,exp,index,label[RAW_NAME])==1)                                          sprintf(raw_name[station_index],"%s"         ,label[RAW_NAME]                          );
    if(elecs){// elecs file
      sprintf(exp,"station[@key='%c']/file"    ,key);if(lookup_process(context,exp,index,label[RAW_NAME])==0) if(strstr(label[TERMINAL],"VSSP")==NULL) if(label[UNIT][0]==0) sprintf(raw_name[station_index],"%s_%s.bin"   ,label[    NAME],&time_code[0]);
      sprintf(exp,"station[@key='%c']/file"    ,key);if(lookup_process(context,exp,index,label[RAW_NAME])==0) if(strstr(label[TERMINAL],"VSSP")==NULL) if(label[UNIT][0]!=0) sprintf(raw_name[station_index],"%s_%s_%s.bin",label[    NAME],&time_code[0],label[UNIT]);
    }else if(vdif){// vdif file
      sprintf(exp,"station[@key='%c']/file"    ,key);if(lookup_process(context,exp,index,label[RAW_NAME])==0) if(strstr(label[TERMINAL],"VSSP")==NULL) if(label[UNIT][0]==0) sprintf(raw_name[station_index],"%s_%s.vdif"   ,label[    NAME],&time_code[0]);
      sprintf(exp,"station[@key='%c']/file"    ,key);if(lookup_process(context,exp,index,label[RAW_NAME])==0) if(strstr(label[TERMINAL],"VSSP")==NULL) if(label[UNIT][0]!=0) sprintf(raw_name[station_index],"%s_%s_%s.vdif",label[    NAME],&time_code[0],label[UNIT]);
    }else{ // rawfile
      sprintf(exp,"station[@key='%c']/file"    ,key);if(lookup_process(context,exp,index,label[RAW_NAME])==0) if(strstr(label[TERMINAL],"VSSP")==NULL) if(label[UNIT][0]==0) sprintf(raw_name[station_index],"%s_%s.raw"   ,label[    NAME],&time_code[0]);
      sprintf(exp,"station[@key='%c']/file"    ,key);if(lookup_process(context,exp,index,label[RAW_NAME])==0) if(strstr(label[TERMINAL],"VSSP")==NULL) if(label[UNIT][0]!=0) sprintf(raw_name[station_index],"%s_%s_%s.raw",label[    NAME],&time_code[0],label[UNIT]);
    }

    sprintf(exp,"station[@key='%c']/file"    ,key);if(lookup_process(context,exp,index,label[RAW_NAME])==0) if(strstr(label[TERMINAL],"VSSP")!=NULL) sprintf(raw_name[station_index],"%c%s%s.dat" ,key,            &time_code[4],label[UNIT]);
    {int ret;struct stat buf;char tmp[256];sprintf(tmp,"%s/%s",RAW_FILE,raw_name[station_index]);ret=stat(tmp,&buf);if(ret==-1) {sprintf(ans,"File : stat('%s') -> %s <%s:%d>\n",tmp,strerror(errno),__FILE__,__LINE__);goto exit;}}
    sprintf(exp,"terminal[@name='%s']"        ,label[TERMINAL]);if(lookup_process(context,exp,index,          NULL)!=1) {sprintf(ans,"GICO3 : [ERROR] terminal[%s]         -> not assigned <%s:%d>\n",label[TERMINAL],__FILE__,__LINE__);goto exit;}
    sprintf(exp,"terminal[@name='%s']/speed"  ,label[TERMINAL]);if(lookup_process(context,exp,index,label[SPEED  ])!=1) {sprintf(ans,"GICO3 : [ERROR] terminal[%s]/speed   -> not assigned <%s:%d>\n",label[TERMINAL],__FILE__,__LINE__);goto exit;}
    sprintf(exp,"terminal[@name='%s']/channel",label[TERMINAL]);if(lookup_process(context,exp,index,label[CHANNEL])!=1) {sprintf(ans,"GICO3 : [ERROR] terminal[%s]/channel -> not assigned <%s:%d>\n",label[TERMINAL],__FILE__,__LINE__);goto exit;}
    sprintf(exp,"terminal[@name='%s']/bit"    ,label[TERMINAL]);if(lookup_process(context,exp,index,label[BIT    ])!=1) {sprintf(ans,"GICO3 : [ERROR] terminal[%s]/bit     -> not assigned <%s:%d>\n",label[TERMINAL],__FILE__,__LINE__);goto exit;}
    sprintf(exp,"terminal[@name='%s']/level"  ,label[TERMINAL]);if(lookup_process(context,exp,index,label[LEVEL  ])!=1) {sprintf(ans,"GICO3 : [ERROR] terminal[%s]/level   -> not assigned <%s:%d>\n",label[TERMINAL],__FILE__,__LINE__);goto exit;}
    strcpy(station[station_index].name,label[NAME]);station[station_index].key=key;strcpy(terminal[station_index].name,label[TERMINAL]);
    terminal[station_index].sps=atoi(label[SPEED]);terminal[station_index].channel=atoi(label[CHANNEL]);terminal[station_index].bit=atoi(label[BIT]);v=terminal[station_index].level;
    if(terminal[station_index].bit==1) if(sscanf(label[LEVEL],"      %f , %f     ",&v[0],&v[1]            )!=2) {sprintf(ans,"GICO3 : [ERROR] terminal[%s]/bit  -> not assigned <%s:%d>\n",label[TERMINAL],__FILE__,__LINE__);goto exit;}
    if(terminal[station_index].bit==2) if(sscanf(label[LEVEL]," %f , %f , %f , %f",&v[0],&v[1],&v[2],&v[3])!=4) {sprintf(ans,"GICO3 : [ERROR] terminal[%s]/bit  -> not assigned <%s:%d>\n",label[TERMINAL],__FILE__,__LINE__);goto exit;}
    sprintf(exp,"clock[@key='%c']/epoch",key);if(lookup_process(context,exp,index,label[EPOCH])==1){strptime(label[EPOCH],"%Y/%j %H:%M:%S",&tm);cdata[station_index].sec=mktime(&tm);} else cdata[station_index].sec=process->epoch;
    sprintf(exp,"clock[@key='%c']/delay",key);if(lookup_process(context,exp,index,label[DELAY])==1) if(station_index == 0 ) {cdata[station_index].delay=atof(label[DELAY])+add_delay;} else{ cdata[station_index].delay=atof(label[DELAY])+add_delay;} 
    sprintf(exp,"clock[@key='%c']/rate" ,key);if(lookup_process(context,exp,index,label[RATE ])==1) cdata[station_index].rate =atof(label[ RATE]);
    sprintf(exp,"clock[@key='%c']/acel" ,key);if(lookup_process(context,exp,index,label[ACEL ])==1) cdata[station_index].acel =atof(label[ ACEL]);
    sprintf(exp,"clock[@key='%c']/jerk" ,key);if(lookup_process(context,exp,index,label[JERK ])==1) cdata[station_index].jerk =atof(label[ JERK]);
    sprintf(exp,"clock[@key='%c']/snap" ,key);if(lookup_process(context,exp,index,label[SNAP ])==1) cdata[station_index].snap =atof(label[ SNAP]);
    b=shuffle_list[station_index];
    sprintf(exp,"shuffle[@key='%c']"    ,key);if(lookup_process(context,exp,index,label[SHUFFLE])!=1) {int ss;for(ss=0;ss<32;ss++) b[ss]=ss;}
    sprintf(exp,"shuffle[@key='%c']"    ,key);if(lookup_process(context,exp,index,label[SHUFFLE])==1) {
      if(sscanf(label[SHUFFLE],"%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d",
		&b[31],&b[30],&b[29],&b[28],&b[27],&b[26],&b[25],&b[24],&b[23],&b[22],&b[21],&b[20],&b[19],&b[18],&b[17],&b[16],
		&b[15],&b[14],&b[13],&b[12],&b[11],&b[10],&b[ 9],&b[ 8],&b[ 7],&b[ 6],&b[ 5],&b[ 4],&b[ 3],&b[ 2],&b[ 1],&b[ 0])!=32) {sprintf(ans,"GICO3 : [ERROR] shuffle[%c]  -> not assigned <%s:%d>\n",key,__FILE__,__LINE__);goto exit;}
    }
    {unsigned int ss,byte,bit;memset(shuffle_table[station_index],0,sizeof(shuffle_table));for(byte=0;byte<4;byte++) for(ss=0;ss<256;ss++) for(bit=00;bit<32;bit++) if(((ss<<(8*byte))>>b[bit])&0x01) shuffle_table[station_index][byte][ss]|=(1<<bit);}
  }
  for(stream_index=0;(stream_index<MAX_STREAM)&&(stream_index<lookup_process(context,"stream",index,NULL));stream_index++){
    sprintf(exp,"stream[%d]/label"    ,stream_index+1);if(lookup_process(context,exp,index,label[LABEL    ])!=1) {sprintf(ans,"GICO3 : [ERROR] stream[%02d]/label     -> not assigned <%s:%d>\n",stream_index+1,__FILE__,__LINE__);goto exit;}
    sprintf(exp,"stream[%d]/source"   ,stream_index+1);if(lookup_process(context,exp,index,label[SOURCE   ])!=1) strcpy(label[SOURCE   ],"None");
    sprintf(exp,"stream[%d]/frequency",stream_index+1);if(lookup_process(context,exp,index,label[FREQUENCY])!=1) strcpy(label[FREQUENCY],"0.0");
    sprintf(exp,"stream[%d]/channel"  ,stream_index+1);if(lookup_process(context,exp,index,label[CHANNEL  ])!=1) {sprintf(ans,"GICO3 : [ERROR] stream[%02d]/channel   -> not assigned <%s:%d>\n",stream_index+1,__FILE__,__LINE__);goto exit;}
    sprintf(exp,"stream[%d]/fft"      ,stream_index+1);if(lookup_process(context,exp,index,label[FFT_POINT])!=1) {sprintf(ans,"GICO3 : [ERROR] stream[%02d]/fft       -> not assigned <%s:%d>\n",stream_index+1,__FILE__,__LINE__);goto exit;}
    sprintf(exp,"stream[%d]/output"   ,stream_index+1);if(lookup_process(context,exp,index,label[OUTPUT   ])!=1) {sprintf(ans,"GICO3 : [ERROR] stream[%02d]/output    -> not assigned <%s:%d>\n",stream_index+1,__FILE__,__LINE__);goto exit;}
    strcpy(stream[stream_index].label,label[LABEL]);if(strcmp(label[OBJECT],"Multi")!=0) strcpy(source[stream_index].name,label[OBJECT]);else strcpy(source[stream_index].name,label[SOURCE]);
    stream[stream_index].frequency=atof(label[FREQUENCY]);stream[stream_index].channel=atoi(label[CHANNEL]);stream[stream_index].points=atoi(label[FFT_POINT]);stream[stream_index].Hz=atoi(label[OUTPUT]);
    sprintf(exp,"source[@name='%s']"    ,source[stream_index].name);if(lookup_process(context,exp,index,      NULL)!=1){sprintf(ans,"GICO3 : [ERROR] source[%8s]     -> not assigned <%s:%d>\n",source[stream_index].name,__FILE__,__LINE__);goto exit;}
    sprintf(exp,"source[@name='%s']/ra" ,source[stream_index].name);if(lookup_process(context,exp,index,label[ RA])!=1){sprintf(ans,"GICO3 : [ERROR] source[%8s]/ra  -> not assigned <%s:%d>\n",source[stream_index].name,__FILE__,__LINE__);goto exit;}
    sprintf(exp,"source[@name='%s']/dec",source[stream_index].name);if(lookup_process(context,exp,index,label[DEC])!=1){sprintf(ans,"GICO3 : [ERROR] source[%8s]/dec -> not assigned <%s:%d>\n",source[stream_index].name,__FILE__,__LINE__);goto exit;}
    sscanf(label[RA],"%02dh%02dm%lf",&ra_hh,&ra_mm,&ra_ss);sscanf(label[DEC],"%c%02dd%02d'%lf",&sign,&dec_dd,&dec_mm,&dec_ss);
    source[stream_index].right_ascension=RA_to_Radian(ra_hh,ra_mm,ra_ss);source[stream_index].declination=DEC_to_Radian(sign,dec_dd,dec_mm,dec_ss);
    for(s=0;(s<station_index)&&(key=process->stations[s]);s++){
      special[s][stream_index].channel=stream[stream_index].channel;special[s][stream_index].sideband=USB;if(stream[stream_index].frequency<0) special[s][stream_index].sideband=LSB;
      sprintf(exp,"stream[%d]/special[@key='%c']/offset"  ,stream_index+1,key);if(lookup_process(context,exp,index,label[  OFFSET])==1) special[s][stream_index].offset =atoi(label[OFFSET ]);
      sprintf(exp,"stream[%d]/special[@key='%c']/channel" ,stream_index+1,key);if(lookup_process(context,exp,index,label[ CHANNEL])==1) special[s][stream_index].channel=atoi(label[CHANNEL]);
      sprintf(exp,"stream[%d]/special[@key='%c']/sideband",stream_index+1,key);if(lookup_process(context,exp,index,label[SIDEBAND])==1) if(strcmp(label[SIDEBAND],"LSB")==0) special[s][stream_index].sideband=LSB; if(strcmp(label[SIDEBAND],"IQ")==0) special[s][stream_index].sideband=IQ;
      sprintf(exp,"stream[%d]/special[@key='%c']/rotation",stream_index+1,key);if(lookup_process(context,exp,index,label[ROTATION])==1) special[s][stream_index].rotation=atof(label[ROTATION]);
    }
    stream[stream_index].frequency=fabs(atof(label[FREQUENCY]));
  }
  *STATION=station_index;
  *STREAM = stream_index;
  return SUCCEED;
 exit:
  fputs(ans,stderr);
  return FAILURE;
}
