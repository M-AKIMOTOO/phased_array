#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <getopt.h>
#include <errno.h>

#include <libxml/tree.h>
#include <libxml/parser.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>
#include "GICO.h"
#include "tool.h"
#include "schedule.h"

char *strptime(const char *s, const char *format, struct tm *tm);

static struct option options[]={
  {name:"schedule",has_arg:required_argument,flag:NULL,val:   1},
  {name:"UT1UTC"  ,has_arg:required_argument,flag:NULL,val:   1},
  {name:"raw-file",has_arg:required_argument,flag:NULL,val:   1},
  {name:"geo-file",has_arg:required_argument,flag:NULL,val:   1},
  {name:"cor-file",has_arg:required_argument,flag:NULL,val:   1},
  {name:"epoch"   ,has_arg:required_argument,flag:NULL,val:   1},
  {name:"length"  ,has_arg:required_argument,flag:NULL,val:   1},
  {name:"remove"  ,has_arg:required_argument,flag:NULL,val:   1},
  {name:"cpu"     ,has_arg:required_argument,flag:NULL,val:   1},
  {name:"pcal"    ,has_arg:      no_argument,flag:NULL,val:   1},
  //  {name:"complex" ,has_arg:      no_argument,flag:NULL,val:   1}, removal by KT on 24 oct 2016
  {name:"short" ,has_arg:      no_argument,flag:NULL,val:   1}, //added by KT on 24 oct 2016
  {name:"iq"      ,has_arg:      no_argument,flag:NULL,val:   1}, // added by KT for Haystack IQ 2015 3,Feb
  {name:"elecs"      ,has_arg:      no_argument,flag:NULL,val:   1}, // added by KT for vdif data input
  {name:"vdif"      ,has_arg:      no_argument,flag:NULL,val:   1}, // added by KT for vdif data input
  {name:"message" ,has_arg:required_argument,flag:NULL,val:   1},
  {name:"delay"   ,has_arg:required_argument,flag:NULL,val:   1},
  {name:"fixdelay",has_arg:      no_argument,flag:NULL,val:   1},// added by KT for zenith observation between kas-ishioka 2015aug27
  {name:"check"   ,has_arg:      no_argument,flag:NULL,val:   1},
  {name:"realtime",has_arg:      no_argument,flag:NULL,val:   1},
  {name:"history" ,has_arg:      no_argument,flag:NULL,val:   1},
  {name:"help"    ,has_arg:      no_argument,flag:NULL,val:   1},
  {name:"version" ,has_arg:      no_argument,flag:NULL,val:   1},
  {name:      NULL,has_arg:                0,flag:NULL,val:   0},
};
void put_help()
{
  fprintf(stdout,"名前\n\t gico3      - スケジュールファイルに従い相関処理を実行します\n");
  fprintf(stdout,"書式\n\t gico3       --schedule=string --raw-file-string            \n");
  fprintf(stdout,"      \t             --geo-file=string --cor-file=string --core=int \n");
  fprintf(stdout,"      \t             --epoch='YYYY/DDD HH:MM:SS' --length=int       \n");
  fprintf(stdout,"      \t             --pcal --short --realtime --delay=double                  \n");
  fprintf(stdout,"      \t             --remove=string [--history] [--help] [--version]\n");
  fprintf(stdout,"オプション                                                          \n");
  fprintf(stdout,"\t--schedule \tGICO用の相関スケジュールファイルを指定します         \n");
  fprintf(stdout,"\t--raw-file \t観測ファイルの格納ディレクトリ名を指定します         \n");
  fprintf(stdout,"\t--geo-file \t遅延ファイルの格納ディレクトリ名を指定します         \n");
  fprintf(stdout,"\t--cor-file \t相関ファイルの格納ディレクトリ名を指定します         \n");
  fprintf(stdout,"\t--cpu      \t相関処理で使用するCPUコア数を指定します              \n");
  fprintf(stdout,"\t--epoch    \t[epoch:epoch+length]間の相関処理のみを実行します     \n");
  fprintf(stdout,"\t--length   \t[epoch:epoch+length]間の相関処理のみを実行します     \n");
  fprintf(stdout,"\t--pcal     \t特殊相関処理(周波数分解能毎のPCAL信号を検出します)   \n");
  //  fprintf(stdout,"\t--complex  \tFFT前にフリンジ回転処理をおこないます                \n");
  fprintf(stdout,"\t--short  \tFFT前にフリンジ回転処理をおこなわない for short baseline   \n");
  fprintf(stdout,"\t--iq  \t へいすたっくのIQデータ処理をおこないます            \n");
  fprintf(stdout,"\t--delay    \t Add delay only for fringe test              \n");
  fprintf(stdout,"\t--fixdelay    \t Fix delay for zenith observation              \n");
  fprintf(stdout,"\t--elecs    \t apply vdif file of elecs sampler processing              \n");
  fprintf(stdout,"\t--vdif    \t apply vdif file for sampler processing              \n");
  fprintf(stdout,"\t--remove   \t指定した観測局を除いて相関処理を実行します           \n");
  fprintf(stdout,"\t--realtime \t各スキャンの終了時刻まで相関処理を待機します         \n");
  fprintf(stdout,"\t--history  \t更新履歴を表示します　                               \n");
  fprintf(stdout,"\t--help     \tこの画面を表示します                                 \n");
  fprintf(stdout,"\t--version  \tソフトのバージョンとコンパイルされた日時を表示します \n");
  exit(0);
}
void put_version()
{
  fprintf(stdout,"gico3 2.7.5 for small memory<%s %s>\n",__DATE__,__TIME__);
  fprintf(stdout,"Copyright (C) 2007- Moritaka Kimura\n");
  fprintf(stdout,"Copyright (C) 2014- Modified by Kazuhiro Takefuji\n");
  exit(0);
}

void put_history()
{
  fprintf(stdout,"gico3 2.6.0 : 公開版のプロトタイプ\n");
  fprintf(stdout,"gico3 2.6.8 : マルチユニットに対応\n");
  fprintf(stdout,"gico3 2.6.9 : あぷりおりを高精度化。光行差＋大気導入\n");
  fprintf(stdout,"gico3 2.7.0 : Add skip[sec] mode. length[sec] should be same as schedule epoch\n");
  fprintf(stdout,"gico3 2.7.1 : Apply vdif file of elecs sampler input \n");
  fprintf(stdout,"gico3 2.7.2.1 : Apply vdif file  input \n");
  fprintf(stdout,"gico3 2.7.3 : Added option of fix delay for zenith observation between kashima - ishioka \n");
  fprintf(stdout,"gico3 2.7.5 on 24 Oct 2017: Added option of short and remove complex description \n");
  fprintf(stdout,"gico3 2.7.5.1 on 20 Apr 2018: Bug fixed for vdif file input \n");

  exit(0);
}

int main(int argc,char **argv)
{
  xmlDoc *doc;xmlXPathContext *context;
  int index,length;struct tm tm;double add_delay;
  char remove[32]="";int check;
  gico3_init();message=0;shortbase=0;cpu=1;gpu=0;ut1utc=0;pcal=0;check=0;realtime=0;iq=0;vdif=0;elecs=0;fixdelay=0;
  setenv("TZ","GMT",1);tzset();strcpy(RAW_FILE,"./raw-file");strcpy(COR_FILE,"./cor-file");strcpy(GEO_FILE,"");
  strptime("1970/001 12:00:00","%Y/%j %H:%M:%S",&tm);length=0x7fffffff;//相関処理開始時刻の初期値
  {FILE *fp;fp=popen("cat /proc/cpuinfo |grep processor|wc","r");fscanf(fp,"%d",&cpu);pclose(fp);}
  if(argc==1) put_help();
  while(getopt_long(argc,argv,"",options,&index)!=-1){
    if(strcmp(options[index].name,"schedule")==0) strcpy(SCHEDULE,optarg);
    if(strcmp(options[index].name,"ut1utc"  )==0) ut1utc=atof(optarg);
    if(strcmp(options[index].name,"raw-file")==0) strcpy(RAW_FILE,optarg);
    if(strcmp(options[index].name,"cor-file")==0) strcpy(COR_FILE,optarg);
    if(strcmp(options[index].name,"geo-file")==0) strcpy(GEO_FILE,optarg);
    if(strcmp(options[index].name,"remove"  )==0) strcpy(remove  ,optarg);
    if(strcmp(options[index].name,"epoch"   )==0) strptime(optarg,"%Y/%j %H:%M:%S",&tm);
    if(strcmp(options[index].name,"length"  )==0) length    =atoi(optarg);
    if(strcmp(options[index].name,"cpu"     )==0) cpu       =atoi(optarg);
    //    if(strcmp(options[index].name,"complex" )==0) complex   =           1;
    if(strcmp(options[index].name,"short" )==0) shortbase   =           1;

    if(strcmp(options[index].name,"iq" )==0) iq   =           1;
    if(strcmp(options[index].name,"fixdelay" )==0) fixdelay   =           1;
    if(strcmp(options[index].name,"message" )==0) message   =atoi(optarg);
    if(strcmp(options[index].name,"delay"   )==0) add_delay   =atof(optarg);
    if(strcmp(options[index].name,"pcal"    )==0) pcal=1;
    if(strcmp(options[index].name,"elecs"    )==0) elecs=1;  
    if(strcmp(options[index].name,"vdif"    )==0) vdif=1;
    if(strcmp(options[index].name,"check"   )==0) check=1;
    if(strcmp(options[index].name,"realtime")==0) realtime=1;
    if(strcmp(options[index].name,"history" )==0) put_history();
    if(strcmp(options[index].name,"help"    )==0) put_help();
    if(strcmp(options[index].name,"version" )==0) put_version();
  }
  xmlInitParser();doc=xmlReadFile(SCHEDULE,NULL,0);context=xmlXPathNewContext(doc);
 loop:
  for(PROCESS_INDEX=1;PROCESS_INDEX<=lookup_context(context,"/schedule/process",NULL);PROCESS_INDEX++){
    if(read_process(context,PROCESS_INDEX,&process,&STATION,&STREAM,remove,add_delay)==FAILURE) continue;
    if(process.epoch+process.length<mktime(&tm)) continue;
    while((process.epoch+0*process.length<mktime(&tm)+0*length)&&(process.length>0)){process.epoch++;process.length--;if(process.length==0) break;}
    while((process.epoch+1*process.length>mktime(&tm)+1*length)&&(process.length>0)){               ;process.length--;if(process.length==0) break;}
    if((process.epoch<mktime(&tm))||(process.length<=0)) continue;
    if(check==1) {process.epoch=process.epoch+process.length/2;process.length=1;}
    if(realtime==1) while(process.epoch+process.length+4>time(NULL)) sleep(1);
    put_process(&process);OPEN();GICO(process);
  }
  return 0;
}
