CFLAGS= -g -lm
CFILES= 
HFILES= dcdStruct.h  ReadDcd.h  ReadPdb.h ReadCM.h get_scd.h
IFLAGS= -I. -I/usr/include
READ_DCD_CPPFILES=  ReadDcd.cpp ReadPdb.cpp ReadCM.cpp
READ_DCD_HFILES= ReadPdb.h ReadDcd.h dcdStruct.h
CPPFILES= get_scd.cpp
MEMBCPPFILES= get_scd_memb.cpp
get_scd: $(CPPFILES) $(READ_DCD_CPPFILES) $(READ_DCD_HFILES) $(HFILES) Makefile
	g++ -g3 -o get_scd $(CPPFILES) $(READ_DCD_CPPFILES) $(IFLAGS) -lpthread -lm


get_scd_memb: $(MEMBCPPFILES) $(READ_DCD_CPPFILES) $(READ_DCD_HFILES) $(HFILES) Makefile
	g++ -g3 -o get_scd_memb $(MEMBCPPFILES) $(READ_DCD_CPPFILES) $(IFLAGS) -lpthread -lm


