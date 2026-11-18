# choose static/dynamic library
LIBEXTENSION=a	
#LIBEXTENSION=so

# main target
LIBNAME=libresample
LIB=$(LIBNAME).$(LIBEXTENSION)

CXXFLAGS+=-Wall -ggdb3 -pipe -O3 -pg
CXXFLAGS+=-I. -fPIC

ROOTCXXFLAGS+=$(shell root-config --cflags)
ROOTLDFLAGS+=$(shell root-config --libs)

OBJECTS=CommonBlockParticleCONEX.o  CommonBlockWrapperCONEX.o \
	CommonBlockParticleCRMC.o  CommonBlockWrapperCRMC.o \
        CommonBlockParticleSIBYLL.o  CommonBlockWrapperSIBYLL.o \
        CommonBlockWrapperSIBYLL_LAB.o \
        ParticleBlock.o  resample.o

.PHONY: all

all: $(LIB)

# dependencies
$(LIB): CommonBlockParticleCONEX.h  CommonBlockWrapperCONEX.h  CommonBlockCONEX.h \
	CommonBlockParticleCRMC.h  CommonBlockWrapperCRMC.h  CommonBlockCRMC.h \
        CommonBlockParticleSIBYLL.h  CommonBlockWrapperSIBYLL.h  CommonBlockSIBYLL.h \
	CommonBlockWrapperSIBYLL_LAB.h  CommonBlockSIBYLL_LAB.h \
        ParticleBlockEntry.h  ParticleBlock.h  resample.h Verbosity.h ResamplingMode.h \
	CommonBlockParticleCONEX.cc  CommonBlockWrapperCONEX.cc \
	CommonBlockParticleCRMC.cc  CommonBlockWrapperCRMC.cc \
        CommonBlockParticleSIBYLL.cc  CommonBlockWrapperSIBYLL.cc \
	CommonBlockWrapperSIBYLL_LAB.cc  \
        ParticleBlock.cc  resample.cc

resample.o: resample.cc
	$(CXX) -c -o $@ $(CXXFLAGS) $(ROOTCXXFLAGS) $^

$(LIBNAME).a: $(OBJECTS)
	$(AR) r $@ $(OBJECTS)

$(LIBNAME).so: $(OBJECTS) 
	$(CXX) -o $@ --shared $(LDFLAGS) $^

clean:
	rm -f *.o
	rm -f $(LIBNAME).*
