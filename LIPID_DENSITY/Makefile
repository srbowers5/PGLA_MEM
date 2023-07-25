CFLAGS= -g -lm
CFILES= 
HFILES= dcdStruct.h  ReadDcd.h  ReadPdb.h ReadCM.h get_lipid_den.h system_conf.h
IFLAGS= -I. -I/usr/include
READ_DCD_CPPFILES=  ReadDcd.cpp ReadPdb.cpp ReadCM.cpp
READ_DCD_HFILES= ReadPdb.h ReadDcd.h dcdStruct.h
CPPFILES= get_lipid_den.cpp
get_lipid_den: $(CPPFILES) $(READ_DCD_CPPFILES) $(READ_DCD_HFILES) $(HFILES) Makefile
	g++ -g3 -o get_lipid_den $(CPPFILES) $(READ_DCD_CPPFILES) $(IFLAGS) -lpthread -lm

test: test.cpp
	g++ -g3 -o test test.cpp
