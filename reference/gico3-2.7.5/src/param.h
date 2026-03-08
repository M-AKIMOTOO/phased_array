#ifndef __PARAM__
#define __PARAM__

#include "header.h"

typedef struct {
  signed   int offset;//  時間領域でのオフセット数 
  unsigned int phase0;//  時間領域での初期位相        
  unsigned int phase1;//  時間領域での位相速度        
  unsigned int phase2;//  時間領域での初期位相
  unsigned int phase3;//  時間領域での位相速度
} Prm3;

typedef struct {
  signed   int offset;//  時間領域でのオフセット数
  float        phase0;//  時間領域での初期位相
  float        phase1;//  時間領域での位相速度
  float        phase2;//  周波数領域での初期位相
} Prm4;

typedef struct {
  time_t        epoch; //パラメータ作成開始時刻のUNIX秒
  int            nsec; //パラメータ作成開始時刻のナノ秒
  double           dt; //1ブロック長に相当する時間[秒]
  int          blocks; //総ブロック数
} Param;

void PARAM0(         );
void PARAM3(int index);
void PARAM2(         );


#endif
