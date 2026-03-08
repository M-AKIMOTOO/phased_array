#ifndef __COMMON__
#define __COMMON__

#include <fftw3.h>
#include "header.h"
#include "config.h"
#include "param.h"

_Bool       pcal; //  0:PCAL検出をしない 　　　　1:PCAL検出をする
//_Bool    complex; //  0:実数型のFFTをする 　     1:複素型のFFTをする
_Bool    shortbase; //  0:複素型のFFTをする        1:実数型のFFTをする 　     

_Bool    iq   ; //  0:実数型のFFTをする 　     1:複素型のFFTをする
_Bool      check; //  0:全時間の処理をする 　    1:中央時間のみ処理をする
_Bool   realtime; //  
_Bool   fixdelay; //  
_Bool   vdif; // input file is vdif // add by kt 2015 jul 11  
_Bool   elecs; // input file is vdif // add by kt 2015 jul 11  a

int          cpu; //  相関処理に使用するCPUコア数
int          gpu; //  相関処理に使用するGPUコア数                                                                                                                                                                                                                          
int      message; //  メッセージ表示のレベル     

unsigned int              *raw_data[MAX_BUFFER][MAX_STATION];//VLBI観測データの読み込みバッファ
unsigned int              *raw_buffer[MAX_BUFFER][MAX_STATION];//VLBI観測データの読み込みバッファ input to raw_data 2015 jul 10 add by KT
signed   int                   raw_offset[MAX_BUFFER][MAX_STATION];//VLBI観測データの観測時刻
signed   int                   raw_length[MAX_BUFFER][MAX_STATION];//VLBI観測データの観測時刻
unsigned int                       shuffle_list[MAX_STATION][32];
unsigned int                      shuffle_table[MAX_STATION][4][256];

Prm3              *prm3[MAX_DEVICE][MAX_STATION];
Prm4              *prm4[MAX_DEVICE][MAX_STATION];

volatile time_t CURRENT_EPOCH;//    相関処理時間
volatile time_t    READ_EPOCH;
int    PROCESS_INDEX;
int     STREAM_INDEX;//　ストリーム番号
int               PP;//          PP番号
int          STATION;//      総観測局数
int           STREAM;//  総ストリーム数
float            BLOCK;//    各ブロック数 (change float for pulsar gating) 2014 aug by KT
double        ut1utc;

char   RAW_FILE[256];
char   GEO_FILE[256];
char   COR_FILE[256];
char   SCHEDULE[256];

Process                process;
Station   station[MAX_STATION];
Terminal terminal[MAX_STATION];
Clock       cdata[MAX_STATION];
Clock       gdata[MAX_STATION];
Stream     stream[MAX_STREAM ];
Source     source[MAX_STREAM ];
Clock   *geo_data[MAX_STATION][MAX_STREAM ];
int    geo_length[MAX_STATION][MAX_STREAM ];
char     raw_name[MAX_STATION][256];
Special   special[MAX_STATION][MAX_STREAM ];
char         unit[MAX_STATION][16];
fftwf_complex  *total[MAX_STATION][MAX_STATION];

#endif


