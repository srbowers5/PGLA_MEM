CFLAGS= -g -lm
CFILES= 
HFILES= dcdStruct.h  ReadDcd.h  ReadPdb.h get_com_xsc.h
IFLAGS= -I. -I/usr/include
READ_DCD_CPPFILES=  ReadDcd.cpp ReadPdb.cpp ReadCM.cpp
READ_DCD_HFILES= ReadPdb.h ReadDcd.h dcdStruct.h
CPPFILES= get_com_xsc.cpp
get_com_xsc: $(CPPFILES) $(READ_DCD_CPPFILES) $(READ_DCD_HFILES) $(HFILES) Makefile
	g++ -g3 -o get_com_xsc $(CPPFILES) $(READ_DCD_CPPFILES) $(IFLAGS) -lpthread -lm
