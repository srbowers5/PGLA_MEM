CFLAGS= -g -lm
CFILES=
HFILES= dcdStruct.h  ReadDcd.h  ReadPdb.h system_conf.h get_water_wire.h
IFLAGS= -I. -I/usr/include
READ_DCD_CPPFILES=  ReadDcd.cpp ReadPdb.cpp
READ_DCD_HFILES= ReadPdb.h ReadDcd.h dcdStruct.h
CPPFILES= get_water_wire.cpp
get_water_wire: $(CPPFILES) $(READ_DCD_CPPFILES) $(READ_DCD_HFILES) $(HFILES) Makefile
	g++ -g3 -o get_water_wire $(CPPFILES) $(READ_DCD_CPPFILES) $(IFLAGS) -lpthread -lm

test: test.cpp
	g++ -g3 -o test test.cpp $(IFLAGS) -lpthread -lm
test2: test2.cpp
	g++ -g3 -o test2 test2.cpp $(IFLAGS) -lpthread -lm
