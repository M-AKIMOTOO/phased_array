#ifndef __SCHEDULE__
#define __SCHEDULE__
#define USB 0x00000000
#define LSB 0x00000001
#define  IQ 0x00000002

#include <libxml/tree.h>
#include <libxml/parser.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>
#include "header.h"
int lookup_context(xmlXPathContext *context,const char *exp,char *ret);
void put_process(Process *process);
int read_process(xmlXPathContext *context,int index,Process *process,int *STATION,int *STREAM,char *remove, double add_delay); // modified by KT 2013 DEC 6


#endif
