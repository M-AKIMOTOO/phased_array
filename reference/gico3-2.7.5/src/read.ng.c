

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <time.h>
#include <pthread.h>
#include <errno.h>
#include "common.h"

/*

void READ1SEC(int index)
{ 
  int ss;unsigned int table[4][256];unsigned char *adr;static int select=0;
  static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;time_t epoch;FILE *fp[MAX_STATION];struct stat buf;struct tm tm;char time_code[32],name[256];int ret,page;off_t len[MAX_STATION];long long s0,s1;struct timeval st,et;double sec,speed;
  sprintf(name,"%s/%s",RAW_FILE,raw_name[index]);file_epoch[index]=process.epoch;if(strchr(raw_name[index],'_')!=NULL) if(strptime(strchr(raw_name[index],'_'),"_%Y%j%H%M%S",&tm)!=NULL) file_epoch[index]=mktime(&tm);
  ret      = stat(name,&buf);if(      ret==  -1) {fprintf(stdout,"GICO3 : [ERROR]   stat('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);} else len[index]=buf.st_size/4;
  fp[index]=fopen(name, "r");if(fp[index]==NULL) {fprintf(stdout,"GICO3 : [ERROR]  fopen('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);}
  for(page=0;page<MAX_BUFFER;page++) if((raw_data[page][index]=malloc(4*Speed(terminal[index])*(1+2*MAX_DELAY)+65536))==NULL) {fprintf(stdout,"GICO3 : [ERROR] malloc('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);}
  for(epoch=process.epoch;epoch<process.epoch+process.length;epoch++){
    gmtime_r(&epoch,&tm);strftime(time_code,sizeof(time_code),"%Y/%j %H:%M:%S",&tm);page=epoch%MAX_BUFFER;while(raw_epoch[page][index]!=0) usleep(10000);
    s0=Min(Max(Speed(terminal[index])*((epoch-file_epoch[index]-MAX_DELAY+0)+cdata[index].delay),0),len[index]);//読み出し開始アドレス
    s1=Min(Max(Speed(terminal[index])*((epoch-file_epoch[index]+MAX_DELAY+1)+cdata[index].delay),0),len[index]);//読み出し終了アドレス
    if(parallel==0) pthread_mutex_lock(&mutex);
    gettimeofday(&st,NULL);
    if(fseek(            fp[index],4*    s0  ,SEEK_SET)!=    0) {fprintf(stdout,"GICO3 : [ERROR]  fseek('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);}
    if(fread(raw_data[page][index],4,(s1-s0),fp[index])!=s1-s0) {fprintf(stdout,"GICO3 : [ERROR]  fread('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);}
    adr=raw_data[page][index];memcpy(table,shuffle_table[index],sizeof(shuffle_table[index]));for(ss=0;ss<s1-s0;ss++) raw_data[page][index][ss]=(table[0][adr[4*ss+0]] | table[1][adr[4*ss+1]] | table[2][adr[4*ss+2]] | table[3][adr[4*ss+3]]);
    gettimeofday(&et,NULL);
    if(parallel==0) pthread_mutex_unlock(&mutex);
    raw_offset[page][index]=s0-Speed(terminal[index])*(epoch-file_epoch[index]);raw_length[page][index]=s1-s0;raw_epoch[page][index]=epoch;sec=(et.tv_sec-st.tv_sec)+0.000001*(et.tv_usec-st.tv_usec);speed=32*(s1-s0)/1000000.0/sec;
    if(message>=2) fprintf(stdout,"GICO3 : [level-2] fread[File='%26s',epoch='%s',offset=%12Ld,length=%8Ld BitRate=%4.0lfMbps]\n",raw_name[index],time_code,s0,s1-s0,speed);
  }
  fclose(fp[index]);for(page=0;page<MAX_BUFFER;page++) {while(raw_epoch[page][index]!=0) usleep(10000);free(raw_data[page][index]);}
}

*/


void Read1Sec(int index)
{   
  time_t raw_epoch[MAX_STATION];size_t size;int ss;unsigned int table[4][256];unsigned char *adr;
  FILE *fp[MAX_STATION];struct stat buf;struct tm tm;char time_code[32],name[256];int ret,page;off_t len[MAX_STATION];long long s0,s1;struct timeval st,et;double sec,speed;
  for(index=0;(index<MAX_STATION)&&(strlen(station[index].name));index++){
    sprintf(name,"%s/%s",RAW_FILE,raw_name[index]);if(strchr(raw_name[index],'_')!=NULL) if(strptime(strchr(raw_name[index],'_'),"_%Y%j%H%M%S",&tm)!=NULL) raw_epoch[index]=mktime(&tm);raw_epoch[index]=process.epoch;
    ret      = stat(name,&buf);if(      ret==  -1) {fprintf(stdout,"GICO3 : [ERROR]   stat('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);} else len[index]=buf.st_size/4;
    fp[index]=fopen(name, "r");if(fp[index]==NULL) {fprintf(stdout,"GICO3 : [ERROR]  fopen('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);}
    size=4*Speed(terminal[index])*(1+2*MAX_DELAY)+4096;
    for(page=0;page<MAX_BUFFER;page++) if((raw_data[page][index]=malloc(size))==NULL) {fprintf(stdout,"GICO3 : [ERROR] malloc('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);}
  }
  // add skip by KT 2015 11 jun
  for(READ_EPOCH=process.epoch+process.skip;READ_EPOCH<process.epoch+process.length;READ_EPOCH++){
    for(index=0;(index<MAX_STATION)&&(strlen(station[index].name));index++){
      s0=Min(Max(Speed(terminal[index])*((READ_EPOCH-raw_epoch[index]-MAX_DELAY+0)+cdata[index].delay),0),len[index]);//読み出し開始アドレス
      s1=Min(Max(Speed(terminal[index])*((READ_EPOCH-raw_epoch[index]+MAX_DELAY+1)+cdata[index].delay),0),len[index]);//読み出し終了アドレス
      gmtime_r(&READ_EPOCH,&tm);strftime(time_code,sizeof(time_code),"%Y/%j %H:%M:%S",&tm);page=READ_EPOCH%MAX_BUFFER;while(CURRENT_EPOCH+MAX_BUFFER<=READ_EPOCH) usleep(10000);gettimeofday(&st,NULL);
      if(fseek(            fp[index],4*    s0  ,SEEK_SET)!=    0) {fprintf(stdout,"GICO3 : [ERROR]  fseek('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);}
      if(fread(raw_data[page][index],4,(s1-s0),fp[index])!=s1-s0) {fprintf(stdout,"GICO3 : [ERROR]  fread('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);}

      /* if(message>=2 ){ */
      /* 	long se=0;for( se = 0 ; se< s1-s0 ; se+=320){ */
      /* 	  fprintf(stdout,"GICO3:[level-2]VSI[index= %2d epoch= %s s0= %15ld s1= %15ld pos= %15ld  data%d= %x [byte]]\n", index,time_code,s0*4,s1*4,se,index,raw_data[page][index][se]  ); */
      /* 	} */
      /* } */
      adr=raw_data[page][index];memcpy(table,shuffle_table[index],sizeof(table));for(ss=0;ss<s1-s0;ss++) raw_data[page][index][ss]=(table[0][adr[4*ss+0]] | table[1][adr[4*ss+1]] | table[2][adr[4*ss+2]] | table[3][adr[4*ss+3]]);
      gettimeofday(&et,NULL);raw_offset[page][index]=s0-Speed(terminal[index])*(READ_EPOCH-raw_epoch[index]);raw_length[page][index]=s1-s0;sec=(et.tv_sec-st.tv_sec)+0.000001*(et.tv_usec-st.tv_usec);speed=32*(s1-s0)/1000000.0/sec;
      if(message>=2) fprintf(stdout,"GICO3 : [level-2] fread[File='%26s',epoch='%s',offset=%12Ld,length=%8Ld BitRate=%4.0lfMbps]\n",raw_name[index],time_code,s0,s1-s0,speed);
    }
  }
  while(CURRENT_EPOCH!=READ_EPOCH) usleep(10000);for(index=0;(index<MAX_STATION)&&(strlen(station[index].name));index++) {fclose(fp[index]);for(page=0;page<MAX_BUFFER;page++) free(raw_data[page][index]);}
}

void Read1SecElecs(int index) // currently only for elecs VdIF format 1ch
{   
  time_t raw_epoch[MAX_STATION];size_t size;long long ss;unsigned int table[4][256];unsigned char *adr;
  FILE *fp[MAX_STATION];struct stat buf;struct tm tm;char time_code[32],name[256];int ret,page;off_t len[MAX_STATION];long long s0,s1;struct timeval st,et;double sec,speed;

  long long framenum,frame, cluster;long posinframe;size_t vdifsize;// 2015 jul 10 add by KT
  long cluster_size = 102300; //
  VDIF vdif_header; long total_jump[MAX_STATION]; int t ; for(t = 0 ; t<MAX_STATION ; t++) total_jump[t]=0;
  unsigned int		initial_sec   = 0;  

 for(index=0;(index<MAX_STATION)&&(strlen(station[index].name));index++){
    sprintf(name,"%s/%s",RAW_FILE,raw_name[index]);if(strchr(raw_name[index],'_')!=NULL) if(strptime(strchr(raw_name[index],'_'),"_%Y%j%H%M%S",&tm)!=NULL) raw_epoch[index]=mktime(&tm);raw_epoch[index]=process.epoch;
    ret      = stat(name,&buf);if(      ret==  -1) {fprintf(stdout,"GICO3 : [ERROR]   stat('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);} else len[index]=buf.st_size/4;
    fp[index]=fopen(name, "r");if(fp[index]==NULL) {fprintf(stdout,"GICO3 : [ERROR]  fopen('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);}
    size     = 4*Speed(terminal[index])*(1+2*MAX_DELAY)*2+4096;
    vdifsize = 4*Speed(terminal[index])*(1+2*MAX_DELAY)*(1312.0/1280.0)+128*2+4096; // vdif frame size / original vsi size = 1312 / 1280
    //size=268804096 vdifsize=275524352
    for(page=0;page<MAX_BUFFER;page++) {
      if((raw_data[page][index]=malloc(size))==NULL ) {	fprintf(stdout,"GICO3 : [ERROR] malloc('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);      } 
      if((raw_buffer[page][index]=malloc(vdifsize))==NULL) {	fprintf(stdout,"GICO3 : [ERROR] malloc('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);      } 
      // 2015 jul 10 add by KT
    }// for page
  }// for index
  // add skip by KT 2015 11 jun

  for(READ_EPOCH=process.epoch+process.skip;READ_EPOCH<process.epoch+process.length;READ_EPOCH++){
    unsigned long		jump	      =	0 ;       
    for(index=0;(index<MAX_STATION)&&(strlen(station[index].name));index++){

      s0=Min(Max(Speed(terminal[index])*((READ_EPOCH-raw_epoch[index]-MAX_DELAY+0)+cdata[index].delay),0),len[index]);//読み出し開始アドレス
      s1=Min(Max(Speed(terminal[index])*((READ_EPOCH-raw_epoch[index]+MAX_DELAY+1)+cdata[index].delay),0),len[index]);//読み出し終了アドレス
      gmtime_r(&READ_EPOCH,&tm);strftime(time_code,sizeof(time_code),"%Y/%j %H:%M:%S",&tm);page=READ_EPOCH%MAX_BUFFER;while(CURRENT_EPOCH+MAX_BUFFER<=READ_EPOCH) usleep(10000);gettimeofday(&st,NULL);
      
      // 1  decide cluster and frame number  
      framenum = (320.0+s0)/320.0;   posinframe = (s0%320)*4; cluster = framenum/cluster_size; // cluster = (4.0*s0)/(8.0*pow(16,6)) ;
      if(message>=2 ) fprintf(stdout,"GICO3 : [level-2] VDIF[index=%2d ,cluster= %6d,frameNUM=%6d, posinframe=%6d s0by=%15ld s1by=%15ld seek=%15ld unit]\n",index,cluster,framenum,posinframe,4*s0,4*s1,((framenum)*1312+cluster*128));
      // 2  read frame unit f0 and f1 
      if(fseek(            fp[index],framenum*1312+cluster*128  ,SEEK_SET)!=    0) {fprintf(stdout,"GICO3 : [ERROR]  fseek('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);}
      
      frame = (320.0+(s1-s0))/320.0;       cluster = frame/cluster_size;  long long readsize = (frame*1312+cluster*128)/4 ;
      if(message>=2) fprintf(stdout,"GICO3 : [level-2] VDIF[index=%2d, sub cluster= %6d,frame=%6d, read=%12d [unit]]\n",index, cluster,frame,(frame*1312+cluster*128));
      
      if(fread(raw_buffer[page][index],4,readsize,fp[index]) != readsize ) {fprintf(stdout,"GICO3 : [ERROR]  fread('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);}
      
      // 3  extract data without header ( offset byte of f0 and f1 

      //      time_t		initial_sec   = process.epoch+process.skip;
      if(initial_sec == 0 ){ initial_sec =  ((char*)raw_buffer[page][index][0]); initial_sec = initial_sec & 0x3FFFFFFF; }
      //      unsigned int current_sec	      =	READ_EPOCH;  
      unsigned int	current_sec   = ((char*)raw_buffer[page][index][0]); current_sec = current_sec & 0x3FFFFFFF; 
      const  int	header_size   = 8 ;	// 4byte unit 32 byte
      const int		frame_size    = 328;	// 4byte unit 1312byte
      const  int	vsi_size      = frame_size-header_size;	// 320 unit
      const  int	dummy	      = 32;	// 4byte unit 128byte blanc data in elecs vdif per 102300 frames
      const int		frame_pos     = posinframe/4.0;
      const long	diff_vsi_size = vsi_size - frame_pos;
      unsigned long pos		      =	0, p1=0;	//frame_pos;  
      unsigned int	next_frame    =	0 ;
      unsigned int	current_frame;	// current_frame is 200000 unit
      current_frame		      = ((char*)raw_buffer[page][index][1]); current_frame = current_frame & 0x00FFFFFF;
      unsigned long	current_jump  =	0;
      jump			      =	0 ;       

      if(process.skip == 0 && (current_sec-initial_sec)==0){       
	current_jump =  (current_frame-framenum) +1 ; // +1 means blanc the first frame of elecs data
      }else if (process.skip != 0 ) {
	current_jump = (current_sec-initial_sec+process.skip-1)*200000 - (framenum-current_frame) +1 ; // +1 means blanc the first frame of elecs data
      }else{
	current_jump = (current_sec-initial_sec)*200000 - (framenum-current_frame) +1 ; // +1 means blanc the first frame of elecs data
      }

      if(current_jump != total_jump[index]){ // in case jump occurs before 1sec data, gico3 reads a bit before 1sec
	total_jump[index] =  current_jump;
      }else{
	pos += (current_jump*vsi_size); // shift due to packet loss data for next second data   
      }
      if(message>=2){ fprintf(stdout,"GICO3:[level-2]VDIF[ vdifsec = %d cur_frame = %d total_frame= %10d curr_jump= %10d total_jump= %10d pos=%10d, sec= %10d ]\n" ,(current_sec-initial_sec), current_frame, framenum, current_jump,total_jump[index],pos,current_sec) ;}

      if (current_jump<0 || current_jump > 200000 ){       fprintf(stdout,"ERROR out of framenum (0<x<200000)! elecs data has some errors! ]\n" );	     exit(1);	    }      
      //////////////////////////////////////////////////////////
      for(ss=header_size+frame_pos ; ss<readsize ; ss++){
	if( p1 == diff_vsi_size ){ 	  ss +=	header_size;   framenum++;
	  next_frame = ((char*)raw_buffer[page][index][ss-7]); next_frame = next_frame & 0x00FFFFFF;	  
	  if( (framenum % cluster_size)	     ==	0 ){  ss+=dummy;  // elecs data 102300 frame blank 32*4 byte data
	    next_frame = ((char*)raw_buffer[page][index][ss-7]); next_frame = next_frame & 0x00FFFFFF;
	  }// if framenum %
	  
	  jump=(next_frame-1)-current_frame; // calc jump
	  
	  if(current_frame+1 == next_frame || current_frame+1 == next_frame+200000 ){
	  }else{ // data loss  due to recording error
	    if( next_frame > current_frame ){
	      if(message>=2 ){ fprintf(stdout,"GICO3:[level-2]VDIF[framejump+ frame= %10d current= %10d next=%10d jump=%10d pos=%10d pos=%12d]\n" ,framenum,current_frame,next_frame,jump,pos,next_frame-current_frame);}
	      pos += (jump*vsi_size); // shift due to jump
	      total_jump[index] += jump; // accumlate jump packet loss
	    }
	    if( next_frame < current_frame ){	     
	      if(message>=2 ){ fprintf(stdout,"GICO3:[level-2]VDIF[overflow frame= %10d cur= %10d next=%10d jump=%10d dif(f-n)=%10d dif(n-c)=%12d]\n" ,framenum,current_frame,next_frame,jump,framenum-current_frame,next_frame-current_frame);}
	      fprintf(stdout,"ERROR elecs data has some errors! ]\n" );	     exit(1);	    }
	  }// if else
	  current_frame=next_frame;  
	}// p1=diff_vsi_size
	if( p1==vsi_size ){	  p1 = 0; 	}
      raw_data[page][index][pos++] = raw_buffer[page][index][ss]; p1++;
      }// for ss
      adr=raw_data[page][index];memcpy(table,shuffle_table[index],sizeof(table));for(ss=0;ss<s1-s0;ss++) raw_data[page][index][ss]=(table[0][adr[4*ss+0]] | table[1][adr[4*ss+1]] | table[2][adr[4*ss+2]] | table[3][adr[4*ss+3]]);
      gettimeofday(&et,NULL);raw_offset[page][index]=s0-Speed(terminal[index])*(READ_EPOCH-raw_epoch[index]);raw_length[page][index]=s1-s0;sec=(et.tv_sec-st.tv_sec)+0.000001*(et.tv_usec-st.tv_usec);speed=32*(s1-s0)/1000000.0/sec;
      if(message>=2) fprintf(stdout,"GICO3 : [level-2] fread[File='%26s',epoch='%s',offset=%12Ld,length=%8Ld BitRate=%4.0lfMbps]\n",raw_name[index],time_code,s0,s1-s0,speed);
      //      fprintf(stdout, "      Process  :  epoch='%s' jump=%6d file=%s \n",time_code, total_jump[index], raw_name[index]);  
      if(total_jump[index] !=0){      fprintf(stdout, "      Process  :  epoch='%s' framejump=%6d file=%s \n",time_code, total_jump[index], raw_name[index]);  }
    } // for index
  }// for READ_EPOCH
  
  while(CURRENT_EPOCH!=READ_EPOCH) usleep(10000);for(index=0;(index<MAX_STATION)&&(strlen(station[index].name));index++) {fclose(fp[index]);for(page=0;page<MAX_BUFFER;page++){ free(raw_data[page][index]);free(raw_buffer[page][index]);}}
}

void Read1SecElecs_NG(int index) // currently only for elecs VdIF format 1ch
{   
  time_t raw_epoch[MAX_STATION];size_t size;long long ss;unsigned int table[4][256];unsigned char *adr;
  FILE *fp[MAX_STATION];struct stat buf;struct tm tm;char time_code[32],name[256];int ret,page;off_t len[MAX_STATION];long long s0,s1;struct timeval st,et;double sec,speed;

  long long framenum,frame, cluster;long posinframe;size_t vdifsize;// 2015 jul 10 add by KT
  long cluster_size = 102300; //
  VDIF vdif_header; long total_jump[MAX_STATION]; int t ; for(t = 0 ; t<MAX_STATION ; t++) total_jump[t]=0;
  unsigned int		initial_sec   = 0;  

 for(index=0;(index<MAX_STATION)&&(strlen(station[index].name));index++){
    sprintf(name,"%s/%s",RAW_FILE,raw_name[index]);if(strchr(raw_name[index],'_')!=NULL) if(strptime(strchr(raw_name[index],'_'),"_%Y%j%H%M%S",&tm)!=NULL) raw_epoch[index]=mktime(&tm);raw_epoch[index]=process.epoch;
    ret      = stat(name,&buf);if(      ret==  -1) {fprintf(stdout,"GICO3 : [ERROR]   stat('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);} else len[index]=buf.st_size/4;
    fp[index]=fopen(name, "r");if(fp[index]==NULL) {fprintf(stdout,"GICO3 : [ERROR]  fopen('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);}
    size     = 4*Speed(terminal[index])*(1+2*MAX_DELAY)*2+4096;
    vdifsize = 4*Speed(terminal[index])*(1+2*MAX_DELAY)*(1312.0/1280.0)+128*2+4096; // vdif frame size / original vsi size = 1312 / 1280
    //size=268804096 vdifsize=275524352
    for(page=0;page<MAX_BUFFER;page++) {
      if((raw_data[page][index]=malloc(size))==NULL ) {	fprintf(stdout,"GICO3 : [ERROR] malloc('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);      } 
      if((raw_buffer[page][index]=malloc(vdifsize))==NULL) {	fprintf(stdout,"GICO3 : [ERROR] malloc('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);      } 
      // 2015 jul 10 add by KT
    }// for page
  }// for index
  // add skip by KT 2015 11 jun
 

  for(READ_EPOCH=process.epoch+process.skip;READ_EPOCH<process.epoch+process.length;READ_EPOCH++){

  skip: ;  // due to packet overflow 

    long		jump	      =	0 ;       
    for(index=0;(index<MAX_STATION)&&(strlen(station[index].name));index++){

      s0=Min(Max(Speed(terminal[index])*((READ_EPOCH-raw_epoch[index]-MAX_DELAY+0)+cdata[index].delay),0),len[index]);//読み出し開始アドレス
      s1=Min(Max(Speed(terminal[index])*((READ_EPOCH-raw_epoch[index]+MAX_DELAY+1)+cdata[index].delay),0),len[index]);//読み出し終了アドレス
      gmtime_r(&READ_EPOCH,&tm);strftime(time_code,sizeof(time_code),"%Y/%j %H:%M:%S",&tm);page=READ_EPOCH%MAX_BUFFER;while(CURRENT_EPOCH+MAX_BUFFER<=READ_EPOCH) usleep(10000);gettimeofday(&st,NULL);
      
      // 1  decide cluster and frame number  
      framenum = (320.0+s0)/320.0;   posinframe = (s0%320)*4; cluster = framenum/cluster_size; // cluster = (4.0*s0)/(8.0*pow(16,6)) ;
      if(message>=2 ) fprintf(stdout,"GICO3 : [level-2] VDIF[index=%2d ,cluster= %6d,frameNUM=%6d, posinframe=%6d s0by=%15ld s1by=%15ld seek=%15ld unit]\n",index,cluster,framenum,posinframe,4*s0,4*s1,((framenum)*1312+cluster*128));
      // 2  read frame unit f0 and f1 
      if(fseek(            fp[index],framenum*1312+cluster*128  ,SEEK_SET)!=    0) {fprintf(stdout,"GICO3 : [ERROR]  fseek('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);}
      
      frame = (320.0+(s1-s0))/320.0;       cluster = frame/cluster_size;  long long readsize = (frame*1312+cluster*128)/4 ;
      if(message>=2) fprintf(stdout,"GICO3 : [level-2] VDIF[index=%2d, sub cluster= %6d,frame=%6d, read=%12d [unit]]\n",index, cluster,frame,(frame*1312+cluster*128));
      
      if(fread(raw_buffer[page][index],4,readsize,fp[index]) != readsize ) {fprintf(stdout,"GICO3 : [ERROR]  fread('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);}
      
      // 3  extract data without header ( offset byte of f0 and f1 
      //      time_t		initial_sec   = process.epoch+process.skip;
      if(initial_sec == 0 ){ initial_sec =  ((char*)raw_buffer[page][index][0]); initial_sec = initial_sec & 0x3FFFFFFF; }
      //      unsigned int current_sec	      =	READ_EPOCH;  
      unsigned int	current_sec   = ((char*)raw_buffer[page][index][0]); current_sec = current_sec & 0x3FFFFFFF; 
      const  int	header_size   = 8 ;	// 4byte unit 32 byte
      const int		frame_size    = 328;	// 4byte unit 1312byte
      const  int	vsi_size      = frame_size-header_size;	// 320 unit
      const  int	dummy	      = 32;	// 4byte unit 128byte blanc data in elecs vdif per 102300 frames
      const int		frame_pos     = posinframe/4.0;
      const long	diff_vsi_size = vsi_size - frame_pos;
      unsigned long pos		      =	0, p1=0;	//frame_pos;  
      unsigned int	next_frame    =	0 ;
      unsigned int	current_frame;	// current_frame is 200000 unit
      current_frame		      = ((char*)raw_buffer[page][index][1]); current_frame = current_frame & 0x00FFFFFF;
      unsigned long	current_jump  =	0;
      jump			      =	0 ;       
      int flag=0;
      if(process.skip == 0 && (current_sec-initial_sec)==0){       
	current_jump =  (current_frame-framenum) +1 ; // +1 means blanc the first frame of elecs data
      }else if (process.skip != 0 ) {
	current_jump = (current_sec-initial_sec+process.skip-1)*200000 - (framenum-current_frame) +1 ; // +1 means blanc the first frame of elecs data
      }else{
	current_jump = (current_sec-initial_sec)*200000 - (framenum-current_frame) +1 ; // +1 means blanc the first frame of elecs data
      }

      if(current_jump != total_jump[index]){ // in case jump occurs before 1sec data, gico3 reads a bit before 1sec
	total_jump[index] =  current_jump;
      }else{
	pos += (current_jump*vsi_size); // shift due to packet loss data for next second data   
      }
      if(message>=2){ fprintf(stdout,"GICO3:[level-2]VDIF[ vdifsec = %d cur_frame = %d total_frame= %10d curr_jump= %10d total_jump= %10d pos=%10d, sec= %10d ]\n" ,(current_sec-initial_sec), current_frame, framenum, current_jump,total_jump[index],pos,current_sec) ;}
      
      if (current_jump<0 || current_jump > 200000 ){       fprintf(stdout,"ERROR out of framenum (0<x<200000)! elecs data has some errors! ]\n" );	     exit(1);	    }      
      //////////////////////////////////////////////////////////
      for(ss=header_size+frame_pos ; ss<readsize ; ss++){
	if( p1 == diff_vsi_size ){ 	  ss +=	header_size;   framenum++;
	  next_frame = ((char*)raw_buffer[page][index][ss-7]); next_frame = next_frame & 0x00FFFFFF;	  
	  if( (framenum % cluster_size)	     ==	0 ){  ss+=dummy;  // elecs data 102300 frame blank 32*4 byte data
	    next_frame = ((char*)raw_buffer[page][index][ss-7]); next_frame = next_frame & 0x00FFFFFF;
	  }// if framenum %
	  
	  jump=(next_frame-1)-current_frame; // calc jump
	  if( next_frame < current_frame){// || next_sec-current_sec == 2){	      
	      jump += 200000; 
	      //	      current_sec += 1; 
	    }

	  //  if(jump < 0 || next_sec == current_sec+1 ){ jump += 200000; }
	  unsigned int next_sec   = ((char*)raw_buffer[page][index][ss-8]); next_sec = next_sec & 0x3FFFFFFF; 
	  //	  if(message>=2 && flag==1){ fprintf(stdout,":[level-2]VDIF[framejump+ frame= %10d current= %10d next=%10d jump=%10d pos=%10d sec=%3d]\n" ,framenum,current_frame,next_frame,jump,pos,next_sec-current_sec);}

	  
	  if(current_frame+1 == next_frame || current_frame+1 == next_frame+200000 ){
	    current_frame=next_frame;  
	  }else{ // data loss  due to recording error
	    if( next_frame < current_frame ){	     
	      if(framenum*vsi_size > pos+(jump*vsi_size)){
		pos += (jump*vsi_size); // shift due to jump
		total_jump[index] += jump; // accumlate jump packet loss
		
	      }else{ // over read size
		ss = readsize;
		current_frame=next_frame;  
	      }
	      if(message>=2 ){ fprintf(stdout,":[level-2]VDIF[jump+ frame= %10d current= %10d next=%10d jump=%10d pos=%10d sec=%3d]\n" ,framenum,current_frame,next_frame,jump,pos,next_sec-current_sec);} 
	      //	    }
	    }
	    if( next_frame < current_frame ){	     

	      fprintf(stdout,"ERROR elecs data has some errors! ]\n" );	     exit(1);	     
	      /* fprintf(stdout,"GICO3:[level-2]VDIF[overflow frame= %10d cur= %10d next=%10d jump=%10d dif(f-n)=%10d dif(n-c)=%12d index=%d]\n" ,framenum,current_frame,next_frame,jump,framenum-current_frame,next_frame-current_frame,index); */
	      /* for(t=-16 ; t<0 ; t++){ */
	      /* 	//		fprintf(stdout,"GICO3:[level-2]VDIF[overflow framenum= %10d next= %10d raw=%10x t= %4d]\n" ,framenum,  ((char*)raw_buffer[page][index][ss+t]) ,  ((char*)raw_buffer[page][index][ss+t]), t); */
	      /* }//for */
	      /* ss+=((current_frame - next_frame) * frame_size ); */
	      /* readsize -= ((current_frame - next_frame) * frame_size ); */
	      /* framenum +=(current_frame - next_frame); */
	      /* current_frame++; */
	      /* //next_sec   = ((char*)raw_buffer[page][index][ss-8]); next_sec = next_sec & 0x3FFFFFFF;  */
	      /* if(next_sec - current_sec >= 1 ){ 	      */
	      /* 	fprintf(stdout,"Second jump! in elecs data %d! ]\n",next_sec-current_sec );	 */
	      /* 	READ_EPOCH += (next_sec  - current_sec); */
	      /* 	//		goto skip ;  */
	      /* }   //  exit(1);	     */
	      /* fprintf(stdout,"GICO3:[level-2]VDIF[overflow frame= %10d cur= %10d next=%10d jump=%10d dif(f-n)=%10d dif(n-c)=%12d index=%d]\n" ,framenum,current_frame,next_frame,jump,framenum-current_frame,next_frame-current_frame,index); */
	      /* for(t=-32 ; t<16 ; t++){ */
	      /* 	fprintf(stdout,"GICO3:[level-2]VDIF[overflow framenum= %10d next= %10d raw=%10x t= %4d]\n" ,framenum,  ((char*)raw_buffer[page][index][ss+t]) ,  ((char*)raw_buffer[page][index][ss+t]), t); */
	      /* }//for */
	      
	    }
	  }// if else
	  
	}// p1=diff_vsi_size
	if( p1==vsi_size ){	  p1 = 0; 

	  //	  if(message>=2 && flag==1){ fprintf(stdout,":[level-2]VDIF[framejump+ frame= %10d current= %10d next=%10d jump=%10d pos=%10d sec=%3d]\n" ,framenum,current_frame,p1);}
	}
	raw_data[page][index][pos++] = raw_buffer[page][index][ss]; p1++;
      }// for ss
      adr=raw_data[page][index];memcpy(table,shuffle_table[index],sizeof(table));for(ss=0;ss<s1-s0;ss++) raw_data[page][index][ss]=(table[0][adr[4*ss+0]] | table[1][adr[4*ss+1]] | table[2][adr[4*ss+2]] | table[3][adr[4*ss+3]]);
      gettimeofday(&et,NULL);raw_offset[page][index]=s0-Speed(terminal[index])*(READ_EPOCH-raw_epoch[index]);raw_length[page][index]=s1-s0;sec=(et.tv_sec-st.tv_sec)+0.000001*(et.tv_usec-st.tv_usec);speed=32*(s1-s0)/1000000.0/sec;
      if(message>=2) fprintf(stdout,"GICO3 : [level-2] fread[File='%26s',epoch='%s',offset=%12Ld,length=%8Ld BitRate=%4.0lfMbps]\n",raw_name[index],time_code,s0,s1-s0,speed);
      //      fprintf(stdout, "      Process  :  epoch='%s' jump=%6d file=%s \n",time_code, total_jump[index], raw_name[index]);  
      if(total_jump[index] !=0){      fprintf(stdout, "      Process  :  epoch='%s' framejump=%6d file=%s \n",time_code, total_jump[index], raw_name[index]);  }
    } // for index
  }// for READ_EPOCH
  
  while(CURRENT_EPOCH!=READ_EPOCH) usleep(10000);for(index=0;(index<MAX_STATION)&&(strlen(station[index].name));index++) {fclose(fp[index]);for(page=0;page<MAX_BUFFER;page++){ free(raw_data[page][index]);free(raw_buffer[page][index]);}}
}

void Read1SecVDIF(int index){ // currently only for vdif software made by Miyauchi-san and Sekido-san

  time_t raw_epoch[MAX_STATION];size_t size;long long ss;unsigned int table[4][256];unsigned char *adr;
  FILE *fp[MAX_STATION];struct stat buf;struct tm tm;char time_code[32],name[256];int ret,page;off_t len[MAX_STATION];long long s0,s1;struct timeval st,et;double sec,speed;

  long long framenum=0,frame=0, cluster=0;unsigned long posinframe=0;size_t vdifsize=0;// 2015 jul 10 add by KT
  long cluster_size = 10; //
  VDIF vdif_header; //long total_jump[MAX_STATION]; int t ; for(t = 0 ; t<index ; t++) total_jump[index]=0;

  unsigned  int		header_size,frame_size,vsi_size,current_frame, current_sec, initial_sec;

  for(index=0;(index<MAX_STATION)&&(strlen(station[index].name));index++){
    
    sprintf(name,"%s/%s",RAW_FILE,raw_name[index]);if(strchr(raw_name[index],'_')!=NULL) if(strptime(strchr(raw_name[index],'_'),"_%Y%j%H%M%S",&tm)!=NULL) raw_epoch[index]=mktime(&tm);raw_epoch[index]=process.epoch;
    ret      = stat(name,&buf);if(      ret==  -1) {fprintf(stdout,"GICO3 : [ERROR]   stat('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);} else len[index]=buf.st_size/4;
    fp[index]=fopen(name, "r");if(fp[index]==NULL) {fprintf(stdout,"GICO3 : [ERROR]  fopen('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);}
    size     = 4*Speed(terminal[index])*(1+2*MAX_DELAY)*2+4096;
    vdifsize = 4*Speed(terminal[index])*(1+2*MAX_DELAY)*2+4096; // 
    //size=268804096 vdifsize=275524352
    for(page=0;page<MAX_BUFFER;page++) {
      if((raw_data[page][index]=malloc(size))==NULL ) {	fprintf(stdout,"GICO3 : [ERROR] malloc('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);      } 
      if((raw_buffer[page][index]=malloc(vdifsize))==NULL) {	fprintf(stdout,"GICO3 : [ERROR] malloc('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);      } 

      // 2015 jul 10 add by KT
    }// for page
  }// for index
  // add skip by KT 2015 11 jun
  /////////////////////////////////////////////////
  // decide date and frame value
  unsigned int *buf_data[MAX_STATION];
  for(index=0;(index<MAX_STATION)&&(strlen(station[index].name));index++){
    if((buf_data[index]=malloc(320))==NULL ) {	fprintf(stdout,"GICO3 : [ERROR] malloc('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);      } 
    if(fseek(   fp[index],0  ,SEEK_SET)!=    0) {fprintf(stdout,"GICO3 : [ERROR]  fseek('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);}
    if(fread(buf_data[index],4,64,fp[index]) !=64 ) {fprintf(stdout,"GICO3 : [ERROR]  fread('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);}
    
    initial_sec      = ((char*)buf_data[index][0]); initial_sec = initial_sec & 0x3FFFFFFF;
    header_size      = ((char*)buf_data[index][0]); header_size = (2 - (header_size & 0x40000000) )*4; // 4byte unit
    current_frame    = ((char*)buf_data[index][1]); current_frame = current_frame & 0x00FFFFFF;
    frame_size       = ((char*)buf_data[index][2]); frame_size = (frame_size & 0x00FFFFFF)*2; // 4byte unit
    vsi_size      = frame_size-header_size;// 4byte unit
    if(message>=2){
      fprintf(stdout,"GICO3:[level-2]VDIF[ vdifsec= %10d  cur_frame= %10d fram_size= %10d vsi_size= %10d %x %x header_size=%d (byte) ]\n" ,initial_sec,current_frame,frame_size*4,vsi_size*4,current_sec,current_frame,header_size*4) ;}
  }// for check vdif header 

///////////////////////////////////////////////
  for(READ_EPOCH=process.epoch+process.skip;READ_EPOCH<process.epoch+process.length;READ_EPOCH++){
    for(index=0;(index<MAX_STATION)&&(strlen(station[index].name));index++){
      s0=Min(Max(Speed(terminal[index])*((READ_EPOCH-raw_epoch[index]-MAX_DELAY+0)+cdata[index].delay),0),len[index]);//読み出し開始アドレス
      s1=Min(Max(Speed(terminal[index])*((READ_EPOCH-raw_epoch[index]+MAX_DELAY+1)+cdata[index].delay),0),len[index]);//読み出し終了アドレス
      gmtime_r(&READ_EPOCH,&tm);strftime(time_code,sizeof(time_code),"%Y/%j %H:%M:%S",&tm);page=READ_EPOCH%MAX_BUFFER;while(CURRENT_EPOCH+MAX_BUFFER<=READ_EPOCH) usleep(10000);gettimeofday(&st,NULL);
      
      // 1  decide cluster and frame number  
      framenum = s0/vsi_size;    posinframe = (s0%vsi_size)*4;
      if(message>=2 ) fprintf(stdout,"GICO3 : [level-2] VDIF[index=%2d ,frameNUM=%6d, posinframe=%6d s0by=%15ld s1by=%15ld ]\n",index,framenum,posinframe,4*s0,4*s1 );
      // 2  read frame unit f0 and f1 
      if(fseek(            fp[index],framenum*frame_size ,SEEK_SET)!=    0) {fprintf(stdout,"GICO3 : [ERROR]  fseek('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);}
      
      frame = (s1-s0)/vsi_size;
      long long readsize = frame*frame_size;
      if(message>=2) fprintf(stdout,"GICO3 : [level-2] VDIF[index=%2d, sub cluster= %6d,frame=%6d, read=%12d [unit]]\n",index, cluster,frame,readsize );
      
      if(fread(raw_buffer[page][index],4,readsize,fp[index]) != readsize ) {fprintf(stdout,"GICO3 : [ERROR]  fread('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);}

      unsigned int		frame_pos     = posinframe/4.0;
      unsigned long		diff_vsi_size = vsi_size - frame_pos;
      unsigned long pos			      =	0, p1=0;	//frame_pos;  
      unsigned long		jump	      =	0 ;
      unsigned long current_jump=0;
      unsigned int	next_frame    =	0 ;
      current_sec      = ((char*)buf_data[index][0]); current_sec = current_sec & 0x3FFFFFFF;
      current_frame		      = ((char*)raw_buffer[page][index][1]); current_frame = current_frame & 0x00FFFFFF;
          
      //      pos += (current_jump*vsi_size); // shift due to packet loss data for next second data 
      //////////////////////////////////////////////////////////
      for(ss=header_size+frame_pos ; ss<readsize ; ss++){
	if( p1 == diff_vsi_size ){ 	  ss +=	header_size;   framenum++;
	  next_frame = ((char*)raw_buffer[page][index][ss-7]); next_frame = next_frame & 0x00FFFFFF;	  
	  
	  jump=(next_frame-1)-current_frame; // calc jump
	  
	  if(current_frame+1 == next_frame || current_frame+1 == next_frame+200000 ){
	  }else{ // data loss  due to recording error
	    if( next_frame > current_frame ){
	      pos += (jump*vsi_size); // shift due to jump
	    }
	    if( next_frame < current_frame ){	      fprintf(stdout,"ERROR vdif data has some errors! ]\n" );	     exit(1);	    }
	  }// if else
	  current_frame=next_frame;  
	}// p1=diff_vsi_size
	if( p1==vsi_size ){	  p1 = 0; 	}
	raw_data[page][index][pos++] = raw_buffer[page][index][ss]; p1++;
      }
      adr=raw_data[page][index];memcpy(table,shuffle_table[index],sizeof(table));for(ss=0;ss<s1-s0;ss++) raw_data[page][index][ss]=(table[0][adr[4*ss+0]] | table[1][adr[4*ss+1]] | table[2][adr[4*ss+2]] | table[3][adr[4*ss+3]]);
      
      gettimeofday(&et,NULL);raw_offset[page][index]=s0-Speed(terminal[index])*(READ_EPOCH-raw_epoch[index]);raw_length[page][index]=s1-s0;sec=(et.tv_sec-st.tv_sec)+0.000001*(et.tv_usec-st.tv_usec);speed=32*(s1-s0)/1000000.0/sec;
      if(message>=2) fprintf(stdout,"GICO3 : [level-2] fread[File='%26s',epoch='%s',offset=%12Ld,length=%8Ld BitRate=%4.0lfMbps]\n",raw_name[index],time_code,s0,s1-s0,speed);
    } // for index
  }// for READ_EPOCH
  
  while(CURRENT_EPOCH!=READ_EPOCH) usleep(10000);for(index=0;(index<MAX_STATION)&&(strlen(station[index].name));index++) {fclose(fp[index]);for(page=0;page<MAX_BUFFER;page++){ free(raw_data[page][index]);free(raw_buffer[page][index]);}}
  
}













