#ifndef __READ__
#define __READ__

time_t ConvertEpoch(  unsigned int counter, unsigned long sec);
void READ1SEC(int index);
void Read1Sec(     void);
void Read1SecElecs(     void);
void Read1SecVDIF(     void);
void CheckVDIFHeader (void);
#endif
