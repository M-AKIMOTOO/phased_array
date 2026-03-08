bin_PROGRAMS  = gico3 

gico3_CFLAGS  = `xml2-config --cflags` -g -O3  -msse -Wall -I/usr/local/sofa/include -I/usr/local/include/ 
gico3_LDFLAGS = `xml2-config --libs` -L/usr/local/sofa/lib -L/usr/local/lib 
gico3_LDADD   =  -lsofa_c  -lpthread -lfftw3f -lm 
gico3_SOURCES = GICO.c GICO3.c apri.c main.c mt19937ar.c param.c read.c schedule.c tool.c
