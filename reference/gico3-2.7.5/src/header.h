#ifndef __HEADER__
#define __HEADER__

#include <time.h>
#include <sys/types.h>
 
// histrory 
// add VDIF header 2015 july 16 by KT

#define __GICO_HEADER_MAGIC__     0x3ea2f983
#define __GICO_HEADER_VERSION__   0x01030000
#define SUCCEED                   0x00000000
#define FAILURE                   0xFFFFFFFF


typedef struct {
  // unsigned = 4byte = 32bits
  unsigned Sec_from_Ref:30;   /** Secondsn from reference Epoch **/ 
  unsigned  Legacy_flag:1;   /** Flag for legacy compatibility, 0:32 Bytes header, 1:16 Bytes header*/
  unsigned  Invalid_flag:1;   /** Invalid flag 1: Invalid, 0: Valid */
  
  unsigned  Frame_Number:24;   /** Data Frame Number */
  unsigned  Ref_Epoch_hy_counter:6;   /** Reference Epoch half year counter since 2000.0 */
  unsigned  unused:2;   /** unused  */

  unsigned  Frame_Length:24;   /** Data Frame length -1, unit: 8 bytes  */
  unsigned  Lg_Nch:5;   /** Log2(nch)*/
  unsigned  Version_No:3;   /** VDIF Version Number */

  unsigned char Station_ID[2];   /** Station ID*/
  unsigned  Thread_Id:10;   /** Thread ID*/
  unsigned  Bit:5;   /** Quantization Bit/sample -1 */
  unsigned  Complex_data_flag:1;   /** Comlex data flag. 1:Complex, 0: Real data */
  // NICT format
  unsigned srate:23;   /** Extended user data 1 */
  unsigned uflag:1;
  unsigned Extended_Data_Version:8;   /** Extended Data Version */

  unsigned char sync_data[4];   /** Extended user data 2 */

  unsigned char user_data2[8];   /** Extended user data 2 */
} VDIF;

typedef struct {
  unsigned int                     sec;
  unsigned int                    nsec;
  double                         delay;
  double                          rate;
  double                          acel;
  double                          jerk;
  double                          snap;
} Clock;

typedef struct {
  char                        name[16];
  unsigned int                     sps;
  unsigned int                     bit;
  unsigned int                 channel;
  float                       level[4];
} Terminal;

typedef struct {
  char                        name[16];
  double                         pos_x;
  double                         pos_y;
  double                         pos_z;
  char                             key;
  char                        dummy[7];
} Station;

typedef struct {
  char                        name[16];
  double               right_ascension;
  double                   declination;
} Source;

typedef struct {
  unsigned int            header_magic;
  unsigned int          header_version;
  unsigned int              scan_index;
  unsigned int            sample_speed; 
  double                     frequency;
  unsigned int                  points;
  int                    sector_number;
  Station                   station[2];
  Source                        source;
  Clock                    clk_data[2];
} Header;

typedef struct {
  unsigned int                 st_time;
  unsigned int                 st_nsec;
  unsigned int                 et_time;
  unsigned int                 et_nsec;
  Clock                    geo_data[2];
  float                      effective;
  float                         amp[2];
  short int                     phs[2];//互換性維持のため2バイト型で格納　1unit=0.01度
} Sector;

typedef struct {
  char                       label[16];
  double                     frequency;
  int                          channel;
  int                           points;
  int                               Hz;
} Stream;

typedef struct {
  int                          offset;
  int                         channel;
  int                        sideband;
  double                     rotation;
} Special;

typedef struct {
  int   offset;//  時間領域でのオフセット数
  float phase0;//  時間領域での初期位相
  float phase1;//  時間領域での位相速度 
  float phase2;//  周波数領域での位相速度
} Parameter4;

typedef struct {
  time_t              epoch;
  off_t              offset;
  size_t             length;
  unsigned int        *data;
} Raw_Buffer;

typedef struct {
  unsigned int                  number;
  time_t                         epoch;
  unsigned int                  length;
  unsigned int                  skip; // add by KT 2015 11 jun
  char                      object[16];
  char                    stations[64];
} Process;

#define Speed(terminal) ((double)(((terminal.sps)*(terminal.bit)*(terminal.channel))/32))

#endif
