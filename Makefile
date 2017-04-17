

# Requires libsspt, use the following to point to its exported header files
# and library
PKG_INC = -I$(BUILD_DIR)/include
PKG_LIB = -L$(BUILD_DIR)/lib -lsspt

#NC_INC = -I$(HOME)/include
#NC_LIB = -Wl,-Bstatic -lnetcdf -Wl,-Bdynamic
#NC_LIB =  -L/usr/local/lib -lnetcdf 
NC_INC = `nc-config --cflags`
NC_LIB = `nc-config --libs` # -lnetcdf 




INCS =  $(PKG_INC) $(NC_INC) 
LIBS =  $(PKG_LIB) $(NC_LIB) -lz -lm -lpthread


# -pg
CFLAGS = -g -Wall -D_REENTRANT 


DEPRECATED_OBJS = vcf33.o vcf33-translator.o   vcf40-translator.o   
OBJS =    vcf40field-translator.o vcfvariable.o datasetdescription.o vcf40.o utilsnetcdf.o utilstext.o

PROGS = vcf2nc


all: $(PROGS)

C++ = g++

%.o : %.cpp %.h
	$(C++) -c $(CFLAGS) $(INCS) $<


vcf2nc: vcf-translator.cpp $(OBJS)
	g++ -o $@ $(CFLAGS) $(INCS) $^ $(LIBS)

#vcfinfo: vcfinfo.cpp $(OBJS)
#	g++ -o $@ $(CFLAGS) $(INCS) $^ $(LIBS)




install: $(PROGS)
	cp $(PROGS) $(BUILD_DIR)/bin

clean:
	rm -f *.o *.a *~ $(PROGS)  core
