#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <assert.h>
#include <time.h>
#include <pthread.h>
#include <errno.h>
#include "common.h"

xmlKeepBlanksDefault(0); 
time_t  ConvertEpoch(  unsigned int counter, unsigned long sec){
  
  uint year = (counter/2)+2000-1900;
  uint month;  if((counter % 2)){ month=6;}else{month=0;}

  // 2000 Jan 1st
  struct tm epoch;
  epoch.tm_sec = 0;       /* 秒 */
  epoch.tm_min = 0;       /* 分 */
  epoch.tm_hour = 0;      /* 時 */
  epoch.tm_mday = 1;      /* 日 */
  epoch.tm_mon = month;       /* 月 ( 1月＝0 ) */
  epoch.tm_year = year;    /* 西暦年 - 1900 */
  epoch.tm_isdst = -1;    /* サマータイムフラグ */
  time_t time = mktime(&epoch); 
  time+=sec;
  //    struct tm *out;
  //    out = localtime(&time);
  //    epoch=*out;
  //    printf("%s,%d\n",asctime(&epoch),sec);
  return time;
}

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
      if(message>=2) fprintf(stdout,"GICO3 : [level-2] fread[File='%26s',epoch='%s',offset=%12Ld,length=%8Ld BitRate=%4.0lfMbps delay=%f]\n",raw_name[index],time_code,s0,s1-s0,speed,cdata[index].delay);
    }
  }
  while(CURRENT_EPOCH!=READ_EPOCH) usleep(10000);for(index=0;(index<MAX_STATION)&&(strlen(station[index].name));index++) {fclose(fp[index]);for(page=0;page<MAX_BUFFER;page++) free(raw_data[page][index]);}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////

void Read1SecElecs(int index) // currently only for elecs VdIF format 1ch
{   
  time_t raw_epoch[MAX_STATION];size_t size;long long ss;unsigned int table[4][256];unsigned char *adr;
  FILE *fp[MAX_STATION];struct stat buf;struct tm tm;char time_code[32],name[256];int ret,page;off_t len[MAX_STATION];long long s0,s1;struct timeval st,et;double sec,speed;

  long long framenum,frame, cluster;long posinframe;size_t vdifsize;// 2015 jul 10 add by KT
  long cluster_size = 102300; //
  long total_jump[MAX_STATION]; int t ; for(t = 0 ; t<MAX_STATION ; t++) total_jump[t]=0;
  unsigned int		initial_sec   = 0;  

 for(index=0;(index<MAX_STATION)&&(strlen(station[index].name));index++){
    sprintf(name,"%s/%s",RAW_FILE,raw_name[index]);if(strchr(raw_name[index],'_')!=NULL) if(strptime(strchr(raw_name[index],'_'),"_%Y%j%H%M%S",&tm)!=NULL) raw_epoch[index]=mktime(&tm);raw_epoch[index]=process.epoch;
    ret      = stat(name,&buf);if(      ret==  -1) {fprintf(stdout,"GICO3 : [ERROR]   stat('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);} else len[index]=buf.st_size/4;
    fp[index]=fopen(name, "r");if(fp[index]==NULL) {fprintf(stdout,"GICO3 : [ERROR]  fopen('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);}
    size     = 4*Speed(terminal[index])*(1+1.2*MAX_DELAY)*2+4096;
    vdifsize = 4*Speed(terminal[index])*(1+2*MAX_DELAY)*(1312.0/1280.0)*2+128*2+4096; // vdif frame size / original vsi size = 1312 / 1280, add *2 for huge packet loss case gv15226
    //size=268804096 vdifsize=275524352
    for(page=0;page<MAX_BUFFER;page++) {
      if((raw_data[page][index]=malloc(size))==NULL ) {	fprintf(stdout,"GICO3 : [ERROR] malloc('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);      } 
      if((raw_buffer[page][index]=malloc(vdifsize))==NULL) {	fprintf(stdout,"GICO3 : [ERROR] malloc('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);      } 
      // 2015 jul 10 add by KT
      //  printf("page : %x %x\n", raw_data[page][index],raw_buffer[page][index]);
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
      framenum = (320.0+s0)/320.0 -total_jump[index];  assert(framenum>0); posinframe = (s0%320)*4; cluster = framenum/cluster_size; // cluster = (4.0*s0)/(8.0*pow(16,6)) ;
      if(message>=2 ) fprintf(stdout,"GICO3 : [level-2] VDIF[index=%2d ,cluster= %6d,frameNUM=%6d (jump%6d), posinframe=%6d s0by=%15ld s1by=%15ld seek=%15ld unit]\n",index,cluster,framenum,total_jump[index],posinframe,4*s0,4*s1,((framenum)*1312+cluster*128));
      // 2  read frame unit f0 and f1 
      if(fseek(            fp[index],framenum*1312+cluster*128  ,SEEK_SET)!=    0) {fprintf(stdout,"GICO3 : [ERROR]  fseek('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);}
      
      frame = (320.0+(s1-s0))/320.0;       cluster = frame/cluster_size;  long long readsize = (frame*1312+cluster*128)/4 ;
      if(message>=2) fprintf(stdout,"GICO3 : [level-2] VDIF[index=%2d, sub cluster= %6d,frame=%6d, read=%12d [unit]]\n",index, cluster,frame,(frame*1312+cluster*128));
      
      if(fread(raw_buffer[page][index],4,readsize*2,fp[index]) != readsize*2 ) {fprintf(stdout,"GICO3 : [ERROR]  fread('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);}
      
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
      long long pos		      =	0, p1=0;	//frame_pos;  
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
	pos+=( current_jump-total_jump[index] )*vsi_size; // ATTENTION pos must be minus !!
	total_jump[index] =  current_jump;
      }else{
	//	pos += (current_jump*vsi_size); // shift due to packet loss data for next second data   
      }
      if(message>=2){ fprintf(stdout,"GICO3:[level-2]VDIF[ vdifsec = %d cur_frame = %d total_frame= %10d curr_jump= %10d total_jump= %10d pos=%10d, sec= %10d ]\n" ,(current_sec-initial_sec), current_frame, framenum, current_jump,total_jump[index],pos,current_sec) ;}

      if (current_jump<0 ){       fprintf(stdout,"ERROR out of framenum < 0! elecs data has some errors! ]\n" );	     exit(1);	    }      
      //////////////////////////////////////////////////////////
      for(ss=header_size+frame_pos ; ss<readsize ; ss++){
	if( p1 == diff_vsi_size ){ 	  ss +=	header_size;   framenum++;
	  next_frame = ((char*)raw_buffer[page][index][ss-7]); next_frame = next_frame & 0x00FFFFFF;	  
	  if( (framenum % cluster_size)	     ==	0 ){  ss+=dummy;  // elecs data 102300 frame blank 32*4 byte data
	    next_frame = ((char*)raw_buffer[page][index][ss-7]); next_frame = next_frame & 0x00FFFFFF;
	  }// if framenum %
	  while(next_frame < current_frame ) { next_frame+=200000;}	  
	  jump=(next_frame-1)-current_frame; // calc jump
	  
	  if(current_frame+1 == next_frame || current_frame+1 == next_frame+200000 ){
	  }else{ // data loss  due to recording error
	    if( next_frame > current_frame ){
	      //      if(message>=2 &&index==0 ){ fprintf(stdout,"GICO3:[level-2]VDIF[framejump+ frame= %10d current= %10d next=%10d jump=%10d pos=%10d total_jump=%12d]\n" ,framenum,current_frame,next_frame,jump,pos,total_jump[index]);}
	      pos += (jump*vsi_size); // shift due to jump
	      total_jump[index] += jump; // accumlate jump packet loss
	    }
	    if( next_frame < current_frame ){	      fprintf(stdout,"ERROR elecs data has some errors! ]\n" );	     exit(1);	    }
	  }// if else
	  current_frame=next_frame;  
	}// p1=diff_vsi_size
	if( p1==vsi_size ){	  p1 = 0; 	}
	if(pos<0){pos++;p1++;}
	else{      raw_data[page][index][pos++] = raw_buffer[page][index][ss]; p1++;}
      }// for ss
      adr=raw_data[page][index];memcpy(table,shuffle_table[index],sizeof(table));for(ss=0;ss<s1-s0;ss++) raw_data[page][index][ss]=(table[0][adr[4*ss+0]] | table[1][adr[4*ss+1]] | table[2][adr[4*ss+2]] | table[3][adr[4*ss+3]]);
      gettimeofday(&et,NULL);raw_offset[page][index]=s0-Speed(terminal[index])*(READ_EPOCH-raw_epoch[index]);raw_length[page][index]=s1-s0;sec=(et.tv_sec-st.tv_sec)+0.000001*(et.tv_usec-st.tv_usec);speed=32*(s1-s0)/1000000.0/sec;
      if(message>=2) fprintf(stdout,"GICO3 : [level-2] fread[File='%26s',epoch='%s',offset=%12Ld,length=%8Ld BitRate=%4.0lfMbps]\n",raw_name[index],time_code,s0,s1-s0,speed);
      //      fprintf(stdout, "      Process  :  epoch='%s' jump=%6d file=%s \n",time_code, total_jump[index], raw_name[index]);  
      if(total_jump[index] !=0){      fprintf(stdout, "      Process  :  epoch='%s' framejump=%6d file=%s \n",time_code, total_jump[index], raw_name[index]);  }
    } // for index
  }// for READ_EPOCH

  while(CURRENT_EPOCH!=READ_EPOCH) usleep(1000);
  for(index=0;(index<MAX_STATION)&&(strlen(station[index].name));index++) {

    for(page=0;page<MAX_BUFFER;page++){
      //   printf("page %d %d %x %x \n",page,(index<MAX_STATION)&&(strlen(station[index].name)),raw_data[page][index],raw_buffer[page][index]); 
	   free(raw_buffer[page][index]);
      //      assert(raw_data[page][index] != NULL );
      free(raw_data[page][index]); 
      
    }
    fclose(fp[index]);
  }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
void Read1SecVDIF(int index){ // currently only for vdif same frame size 2018 Apr 20

  time_t raw_epoch[MAX_STATION];size_t size;long long ss;unsigned int table[4][256];unsigned char *adr;
  FILE *fp[MAX_STATION];struct stat buf;struct tm tm;char time_code[32],name[256];int ret,page;off_t len[MAX_STATION];long long s0,s1;struct timeval st,et;double sec,speed;
  double fractional_sec=0;  double diff_sec=0;
  long long framenum=0, cluster=0;unsigned long offsetframe=0;size_t vdifsize=0;// 2015 jul 10 add by KT
  unsigned long Frames1sec; // 1sec frame size
  VDIF header[MAX_STATION]; 
  long total_jump[MAX_STATION]; int t ; for(t = 0 ; t<MAX_STATION ; t++) total_jump[t]=0;
  long		header_size,frame_size,vsi_size,current_frame, current_sec, initial_sec; 
  long jump_frame=0,current_jump=0;

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
  // decide date and frame value header read
  unsigned int *buf_data[MAX_STATION];
  for(index=0;(index<MAX_STATION)&&(strlen(station[index].name));index++){
    if((buf_data[index]=malloc(32))==NULL ) {	fprintf(stdout,"GICO3 : [ERROR] malloc('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);      }
    if(fseek(   fp[index],0  ,SEEK_SET)!=    0) {fprintf(stdout,"GICO3 : [ERROR]  fseek('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);}
    if(fread(buf_data[index],4,8,fp[index]) !=8 ) {fprintf(stdout,"GICO3 : [ERROR]  fread('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);}

    memcpy(&header[index],(void*)buf_data[index],32); // convert
    initial_sec      =  ConvertEpoch(  header[index].Ref_Epoch_hy_counter,header[index].Sec_from_Ref);
	//	printf("init %d", initial_sec);
    current_frame    = header[index].Frame_Number; // start from 0
    // fix values 
    header_size      =  8-header[index].Legacy_flag*4;// 4byte unit
    frame_size       = header[index].Frame_Length*2; // 4byte unit  OK
    vsi_size         = frame_size-header_size;// 4byte unit
    Frames1sec       = terminal[index].sps*terminal[index].bit*terminal[index].channel/8/(vsi_size*4);
    fractional_sec = (double)current_frame / Frames1sec ;
    diff_sec = (long)((double)(initial_sec +fractional_sec- process.epoch - process.skip  ) );
    //        total_jump[index] =  (long)((double)((process.epoch + process.skip)-(initial_sec +fractional_sec)  )  *Frames1sec);

    if(message>=2){
      fprintf(stdout,"GICO3:[level-2]VDIF[%s vdifsec= %10d.%f diffsec=%10f jumpframe=%d cur_frame= %10d fram_size= %10d vsi_size= %10d %x %x header_size=%d Frames1sec=%d ]\n" ,station[index].name,initial_sec,fractional_sec,diff_sec, jump_frame,current_frame,frame_size*4,vsi_size*4,current_sec,current_frame,header_size*4,Frames1sec) ;}
  }// for check vdif header

///////////////////////////////////////////////
  for(READ_EPOCH=process.epoch+process.skip;READ_EPOCH<process.epoch+process.length;READ_EPOCH++){
    for(index=0;(index<MAX_STATION)&&(strlen(station[index].name));index++){
      s0=Min(Max(Speed(terminal[index])*((READ_EPOCH-raw_epoch[index]-MAX_DELAY+0)+cdata[index].delay),0),len[index]);//読み出し開始アドレス
      s1=Min(Max(Speed(terminal[index])*((READ_EPOCH-raw_epoch[index]+MAX_DELAY+1)+cdata[index].delay),0),len[index]);//読み出し終了アドレス
      if(message>=2 ) fprintf(stdout,"GICO3 : [level-2] MAX[ %s delay=%f len[index]=%ld,s0by=%15ld s1by=%15ld ]\n",station[index].name,(cdata[index].delay),len[index],4*s0,4*s1 );

      gmtime_r(&READ_EPOCH,&tm);strftime(time_code,sizeof(time_code),"%Y/%j %H:%M:%S",&tm);page=READ_EPOCH%MAX_BUFFER;while(CURRENT_EPOCH+MAX_BUFFER<=READ_EPOCH) usleep(10000);gettimeofday(&st,NULL);
      
      // 1  decide cluster and frame number  
      framenum = s0/vsi_size;    offsetframe = (s0%vsi_size);
      // 2  read frame unit f0 and f1 with jump, jump is enabled next second
      
      if(fseek(            fp[index],((framenum-total_jump[index])*frame_size)*4 /* 4byte*/ ,SEEK_SET)!=    0) {fprintf(stdout,"GICO3 : [ERROR]  fseek('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);}
      unsigned long readsize = ((s1-s0)*frame_size)/vsi_size; // how many frames in  read data
      if(message>=2) fprintf(stdout,"GICO3 : [level-2] VDIF[%s, framenum= %6d, offset4byte=%d, read=%12d [unit]]\n",station[index].name, framenum,offsetframe,readsize );
      // read 1.5times volumn data
      fread(raw_buffer[page][index],4,readsize,fp[index]);
	    //	    if(fread(raw_buffer[page][index],4,readsize*1.5,fp[index])) != readsize*1.5 ) {fprintf(stdout,"GICO3 : [ERROR]  fread('%s/%s') -> %s <%s:%d>\n",RAW_FILE,raw_name[index],strerror(errno),__FILE__,__LINE__);exit(0);}

      unsigned long	diff_vsi_size = vsi_size - offsetframe;
      unsigned long     pos	      =	0, p1=0;	//frame_pos;  
      
//////////////////////////////////////////////////////////

      for(ss=offsetframe+header_size ; ss<readsize ; ss++){ // 4byte unit
	if( p1 == diff_vsi_size ){ 	 
	  memcpy(&header[index],(void*)(raw_buffer[page][index]+ss),32); // convert 
	  
	  ss +=	header_size;   framenum++; 
	  current_sec      =  ConvertEpoch(  header[index].Ref_Epoch_hy_counter,header[index].Sec_from_Ref);
	  current_frame    = header[index].Frame_Number;
	  current_jump = (current_sec-initial_sec+process.skip)*Frames1sec - (framenum-current_frame)  ; 

	  if( framenum != (current_sec-initial_sec+process.skip)*Frames1sec +current_frame){
	    if(current_jump > 0 ) {
	    pos=( current_jump-total_jump[index] )*vsi_size; // ATTENTION pos must be minus !!
	    }else{// if recording start from schedule epoch
	      ss+=( (-1)*current_jump )*vsi_size; 
	    }
	    total_jump[index] +=  current_jump; 
	    framenum += current_jump;
	    printf("%s : Frame jumped %d [frames] total %d [frames] \n",station[index].name, current_jump,total_jump[index]);
	  } 

	  if( message>=2){ fprintf(stdout,"GICO3 : VDIF2[%s,cur_sec=%ld, frame#= %6d diff=%d diff_frame=%d %s\n",station[index].name,current_sec,current_frame,p1,  (current_sec-initial_sec+process.skip)*Frames1sec - (framenum-current_frame)  , header[index].sync_data); 
	  }
	}// p1=diff_vsi_size
	if( p1==vsi_size ){	  
	  if( message>=2){ fprintf(stdout,"GICO3 : [level-2] VDIFIVS[%s,cur_sec=%ld, frame#= %6d diff=%d %s\n",station[index].name,current_sec,current_frame,p1,header[index].sync_data);}
	  p1 = 0; 	
	}
	raw_data[page][index][pos++] = raw_buffer[page][index][ss]; p1++;
      }
      adr=raw_data[page][index];memcpy(table,shuffle_table[index],sizeof(table));for(ss=0;ss<s1-s0;ss++) raw_data[page][index][ss]=(table[0][adr[4*ss+0]] | table[1][adr[4*ss+1]] | table[2][adr[4*ss+2]] | table[3][adr[4*ss+3]]);
      
      gettimeofday(&et,NULL);raw_offset[page][index]=s0-Speed(terminal[index])*(READ_EPOCH-raw_epoch[index]);raw_length[page][index]=s1-s0;sec=(et.tv_sec-st.tv_sec)+0.000001*(et.tv_usec-st.tv_usec);speed=32*(s1-s0)/1000000.0/sec;
      if(message>=2) fprintf(stdout,"GICO3 : [level-2] fread[File='%26s',epoch='%s',offset=%12Ld,length=%8Ld BitRate=%4.0lfMbps]\n",raw_name[index],time_code,s0,s1-s0,speed);
    } // for index
  }// for READ_EPOCH
  
  while(CURRENT_EPOCH!=READ_EPOCH) usleep(10000);for(index=0;(index<MAX_STATION)&&(strlen(station[index].name));index++) {fclose(fp[index]);for(page=0;page<MAX_BUFFER;page++){ free(raw_data[page][index]);free(raw_buffer[page][index]);}}
  
}













