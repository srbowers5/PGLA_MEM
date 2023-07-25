CFLAGS= -g -lm
CFILES=
HFILES= dcdStruct.h  ReadDcd.h  ReadPdb.h system_conf.h
IFLAGS= -I. -I/usr/include
READ_DCD_CPPFILES=  ReadDcd.cpp ReadPdb.cpp
READ_DCD_HFILES= ReadPdb.h ReadDcd.h dcdStruct.h get_water_region.h
CPPFILES= get_water_region.cpp
get_water_region: $(CPPFILES) $(READ_DCD_CPPFILES) $(READ_DCD_HFILES) $(HFILES) Makefile
	g++ -g3 -o get_water_region $(CPPFILES) $(READ_DCD_CPPFILES) $(IFLAGS) -lpthread -lm

