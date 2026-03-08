#ifndef __CONFIG__
#define __CONFIG__

#define MAX_STATION   64
#define MAX_STREAM    64
#define MAX_DEVICE    64
#define MAX_DELAY  0.025
#define MAX_BUFFER     4

#ifndef Max
#define Max(a,b) ( (a) > (b) ? (a) : (b) )
#endif

#ifndef Min
#define Min(a,b) ( (a) < (b) ? (a) : (b) )
#endif

#endif
