

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

      if(message>=2 ){
      	long se=0;for( se = 0 ; se< s1-s0 ; se+=320){
      	  fprintf(stdout,"GICO3:[level-2]VSI[index= %2d epoch= %s s0= %15ld s1= %15ld pos= %15ld  data%d= %x [byte]]\n", index,time_code,s0*4,s1*4,se,index,raw_data[page][index][se]  );
      	}
      }
      adr=raw_data[page][index];memcpy(table,shuffle_table[index],sizeof(table));for(ss=0;ss<s1-s0;ss++) raw_data[page][index][ss]=(table[0][adr[4*ss+0]] | table[1][adr[4*ss+1]] | table[2][adr[4*ss+2]] | table[3][adr[4*ss+3]]);
      gettimeofday(&et,NULL);raw_offset[page][index]=s0-Speed(terminal[index])*(READ_EPOCH-raw_epoch[index]);raw_length[page][index]=s1-s0;sec=(et.tv_sec-st.tv_sec)+0.000001*(et.tv_usec-st.tv_usec);speed=32*(s1-s0)/1000000.0/sec;
      if(message>=2) fprintf(stdout,"GICO3 : [level-2] fread[File='%26s',epoch='%s',offset=%12Ld,length=%8Ld BitRate=%4.0lfMbps]\n",raw_name[index],time_code,s0,s1-s0,speed);
    }
  }
  while(CURRENT_EPOCH!=READ_EPOCH) usleep(10000);for(index=0;(index<MAX_STATION)&&(strlen(station[index].name));index++) {fclose(fp[index]);for(page=0;page<MAX_BUFFER;page++) free(raw_data[page][index]);}
}


void Read1SecVDIF(int index) // currently only for elecs VdIF format 1ch
{   
  time_t raw_epoch[MAX_STATION];size_t size;long long ss;unsigned int table[4][256];unsigned char *adr;
  FILE *fp[MAX_STATION];struct stat buf;struct tm tm;char time_code[32],name[256];int ret,page;off_t len[MAX_STATION];long long s0,s1;struct timeval st,et;double sec,speed;

  long long framenum,frame, cluster;long posinframe;size_t vdifsize;// 2015 jul 10 add by KT
  long cluster_size = 102300; //
  VDIF vdif_header; long total_jump[MAX_STATION]; int t ; for(t = 0 ; t<index ; t++) total_jump[index]=0;
  
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
      unsigned int current_sec=READ_EPOCH;// current_sec =  ((char*)raw_buffer[page][index][0]); current_sec = current_sec & 0x3FFFFFFF;  */
      time_t initial_sec = process.epoch+process.skip;
      int header_size = 8 ; // 4byte unit 32 byte
      int frame_size = 328; // 4byte unit 1312byte
      int vsi_size = frame_size-header_size; // 320 unit
      int dummy = 32; // 4byte unit 128byte
      int frame_pos = posinframe/4.0;
      long diff_vsi_size  = vsi_size - frame_pos;
      long pos=0, p1=0;//frame_pos;  
      unsigned int flag=0;
      long jump=0 ;
      unsigned int next_frame=0 ;
      unsigned int current_frame;// current_frame is 200000 unit
      current_frame= ((char*)raw_buffer[page][index][1]); current_frame = current_frame & 0x00FFFFFF;

      //      total_jump[index] = (current_sec-initial_sec)*200000+(framenum-current_frame); // +1 means blanc the first frame of elecs data
      if(message>=2){ fprintf(stdout,"GICO3:[level-2]VDIF[ sec = %d current_frame = %d total_frame= %10d jump= %10d]\n" ,(current_sec-initial_sec), current_frame, framenum, total_jump[index]) ;}
         pos += (total_jump[index]*vsi_size);
      ///////////////////////////////////////////////////////////
      //      if(current_frame != framenum) seekframe() or 0 padding; // seek before reading
      //////////////////////////////////////////////////////////
      for(ss=header_size+frame_pos ; ss<readsize ; ss++){
	
	if( p1 == diff_vsi_size ){ 	  ss +=	header_size;   framenum++;
	  next_frame = ((char*)raw_buffer[page][index][ss-7]); next_frame = next_frame & 0x00FFFFFF;	  
	  //  if(message>=2 ){ fprintf(stdout,"GICO3:[level-2]VDIF[ ss= %12d pos= %12d diff= %8d framenum= %12d p1= %8d dat= %x [byte]]\n", ss  ,pos, (ss-pos)/8,framenum,p1,raw_buffer[page][index][ss-diff_vsi_size]);  }
	  if( (framenum % cluster_size)	     ==	0 ){  ss+=dummy;  
	    /* for(pos=-16 ; pos<0 ; pos++){   */
	    /*   if(message>=2 ){ fprintf(stdout,"GICO3:[level-2]VDIF[framenum= %10d next= %10d raw=%10x pos= %4d]\n" ,framenum,  ((char*)raw_buffer[page][index][ss+pos]) ,  ((char*)raw_buffer[page][index][ss+pos]), pos);} */
	    /* }//for */
	    next_frame = ((char*)raw_buffer[page][index][ss-7]); next_frame = next_frame & 0x00FFFFFF;
	    
	  }// if framenum %
	  
	  jump=(next_frame-1)-current_frame;
	  if(message>=2 && framenum%10000==0){ fprintf(stdout,"GICO3:[level-2]VDIF[ frame= %10d cur= %10d next=%10d jump=%10d dif(f-n)=%10d dif(n-c)=%12d]\n" ,framenum,current_frame,next_frame,jump,framenum-current_frame,next_frame-current_frame);}
	  if(current_frame+1 == next_frame || current_frame+1 == next_frame+200000 ){//||current_frame+2 == next_frame || current_frame+2 == next_frame+200000 ){
	  }else{ // data loss  due to recording error
	    if( next_frame > current_frame ){
	      if(message>=2 ){	    fprintf(stdout,"GICO3:[level-2]VDIF[ framejump= %10d cur= %10d next=%10d jump=%10d total_jump=%10d dif(n-c)=%12d]\n" ,framenum,current_frame,next_frame,jump,total_jump[index],next_frame-current_frame);}
	      pos += (jump*vsi_size);
	      total_jump[index] += jump;
	      //framenum+=jump;
	      flag=1;
	    }
	    if( next_frame < current_frame ){
	      fprintf(stdout,"ERROR elecs data has some errors! ]\n" );	   //   exit(1);
	    }
	  }// if else
	  current_frame=next_frame;  
	}// p1=diff_vsi_size
	
	
	if( (p1==vsi_size) ){	    // %==0 is nessesarry in case of p1==0
	  //  if(message>=2 && index==0 ){ fprintf(stdout,"GICO3:[level-2]VDIF[ ss= %12d pos= %12d diff= %8d framenum= %12d p1= %8d data= %x epoch=%s [byte]]\n", ss  ,pos, (ss-pos)/8,framenum,p1,raw_buffer[page][index][ss-320],time_code);}
	  p1 = 0; 
	}

	//	if(message>=2 && index==0 && framenum >= 0) fprintf(stdout,"GICO3 : [level-2] VDIF[ss= %12d pos= %12d diff = %8d framenum=%12d p1=%8d vsi= %x dat=%x [byte]]\n", ss  ,pos, (ss-pos)/8,framenum,p1,raw_data[page][index][pos-1],raw_buffer[page][index][ss]);
	raw_data[page][index][pos++]  = raw_buffer[page][index][ss]; p1++;
      }
      //if(message>=2) fprintf(stdout,"GICO3 : [level-2] VDIF[index=%2d pos= %12d s1-s0 = %12d diff=%12d index=%2d data=%x %x[byte]]\n", index, pos*4  ,(s1-s0)*4, (pos-s1+s0)*4, index, raw_data[page][index][0],raw_data[page][index][60000000] );
      // 4 copy read data to raw_data structure
      
      adr=raw_data[page][index];memcpy(table,shuffle_table[index],sizeof(table));for(ss=0;ss<s1-s0;ss++) raw_data[page][index][ss]=(table[0][adr[4*ss+0]] | table[1][adr[4*ss+1]] | table[2][adr[4*ss+2]] | table[3][adr[4*ss+3]]);
      
      gettimeofday(&et,NULL);raw_offset[page][index]=s0-Speed(terminal[index])*(READ_EPOCH-raw_epoch[index]);raw_length[page][index]=s1-s0;sec=(et.tv_sec-st.tv_sec)+0.000001*(et.tv_usec-st.tv_usec);speed=32*(s1-s0)/1000000.0/sec;
      
    } // for index
  }// for READ_EPOCH
  
  while(CURRENT_EPOCH!=READ_EPOCH) usleep(10000);for(index=0;(index<MAX_STATION)&&(strlen(station[index].name));index++) {fclose(fp[index]);for(page=0;page<MAX_BUFFER;page++){ free(raw_data[page][index]);free(raw_buffer[page][index]);}}
  
}













