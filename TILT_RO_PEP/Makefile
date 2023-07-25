CFLAGS= -g -lm
CFILES= 
HFILES= get_hel_ulr.h ParseLine.h
IFLAGS= -I. -I/usr/include
CPPFILES= get_hel_ulr.cpp

IFLAGS= -I. -I/usr/include
READ_DCD_CPPFILES=  ReadDcd.cpp ReadPdb.cpp ReadCM.cpp
READ_DCD_HFILES= ReadPdb.h ReadDcd.h dcdStruct.h
TILT_CPPFILES= get_tilt_ulr.cpp
NTILT_CPPFILES= get_tilt_ulr_neg.cpp




get_hel_ulr: $(CPPFILES) $(HFILES) Makefile
	g++ -g3 -o get_hel_ulr $(CPPFILES) $(IFLAGS) -lpthread -lm

get_tilt_ulr: $(TILT_CPPFILES) $(READ_DCD_CPPFILES) $(READ_DCD_HFILES) $(HFILES) Makefile
	g++ -g3 -o get_tilt_ulr $(TILT_CPPFILES) $(READ_DCD_CPPFILES) $(IFLAGS) -lpthread -lm

get_tilt_ulr_neg: $(NTILT_CPPFILES) $(READ_DCD_CPPFILES) $(READ_DCD_HFILES) $(HFILES) Makefile
	g++ -g3 -o get_tilt_ulr_neg $(NTILT_CPPFILES) $(READ_DCD_CPPFILES) $(IFLAGS) -lpthread -lm


