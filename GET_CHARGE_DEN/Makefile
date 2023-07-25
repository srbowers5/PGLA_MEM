CFLAGS= -g -lm
CFILES= 
HFILES= dcdStruct.h  ReadDcd.h  ReadPdb.h ReadCM.h get_charge_den.h
IFLAGS= -I. -I/usr/include
READ_DCD_CPPFILES=  ReadDcd.cpp ReadPdb.cpp ReadPsf.cpp
READ_DCD_HFILES= ReadPdb.h ReadDcd.h dcdStruct.h ReadPsf.h
CPPFILES= get_charge_den.cpp get_area_xy.cpp
CPPFILES_Z= get_charge_den_z.cpp get_area_xy.cpp
get_charge_den: $(CPPFILES) $(READ_DCD_CPPFILES) $(READ_DCD_HFILES) $(HFILES) Makefile
	g++ -g3 -o get_charge_den $(CPPFILES) $(READ_DCD_CPPFILES) $(IFLAGS) -lpthread -lm

get_charge_den_z: $(CPPFILES_Z) $(READ_DCD_CPPFILES) $(READ_DCD_HFILES) $(HFILES) Makefile
	g++ -g3 -o get_charge_den_z $(CPPFILES_Z) $(READ_DCD_CPPFILES) $(IFLAGS) -lpthread -lm

test: test.cpp
	g++ -g3 -o test test.cpp
