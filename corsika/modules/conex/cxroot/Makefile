# conex version 
VER := 7.801
# conex extension
# activate this, if you want to manipulate cross-sections and 
# multiplicities etc., or you want to write/read particle lists
# during the simulation
#CONEX_EXTENSIONS+=-DCONEX_EXTENSIONS  # 
# system type and directory layout
SYSTEM := $(shell uname -s)

ifdef CONEXBINROOT
  CONEX_PREFIX := $(CONEXBINROOT)
else
  CONEX_PREFIX := $(CURDIR)
endif
TABDIR_PRE := $(CONEX_PREFIX)/tabs
ifdef CORSIKA_DATA
  TABDIR_PRE := $(CORSIKA_DATA)
else
  ifdef CONEXTABROOT
    TABDIR_PRE := $(CONEXTABROOT)
  endif
endif
TABDIR := $(TABDIR_PRE)
ifneq ($(shell echo $(TABDIR_PRE) | egrep .+[^\/]$), )
  TABDIR := $(TABDIR_PRE)/
endif

SRCDIR = $(CURDIR)/src
OBJDIR = $(CONEX_PREFIX)/obj/$(SYSTEM)
LIBDIR = $(CONEX_PREFIX)/lib/$(SYSTEM)
BINDIR = $(CONEX_PREFIX)/bin/
CFGDIR = $(CONEX_PREFIX)/cfg
# general compiling options
MC3D   = 
#MC3D   = -D__MC3D__ #-D__COAST__
#LOWMEMORY = -D__SORT_FOR_ENERGY__
LOWMEMORY= -D__SAVEMEMO__  -D__SORT_FOR_ENERGY__ -D__CXDEBUG__ -D__GHEISHA__
CPPFLAGS = -D__URQMD__ -D__CXSUB__  $(LOWMEMORY) $(CONEX_EXTENSIONS) $(MC3D) -fPIC
# compiler selection
COMP := `which gfortran 2>&1 | grep -v which`
COMPIL := "$(shell echo $(COMP))"
#COMPIL = ""
ifeq ($(COMPIL), "")
 FC  = g77  # for g77
 G2C = g2c  # for g77
else
 FC = gfortran  # for gfortran
 G2C = gfortran  # for gfortran
 FFLAGS += -std=legacy -fPIC -Wall
endif
ifndef CXX
  CXX = g++
endif
ifndef CC
  CC = gcc
endif
LD = $(CXX)

ROOTCXXFLAGS := 
ROOTLDFLAGS :=
CCFLAGS := 
ifndef CX_NO_ROOT
  ROOTCXXFLAGS := -I$(ROOTSYS)/include $(shell root-config --cflags)
  ROOTLDFLAGS := -Wl,--no-as-needed $(shell root-config --libs)
  CCFLAGS+= $(ROOTCFLAGS)
endif
PRESHOW_OBJS := 
ifdef CORSIKA_8
  # use HE/LE model code from C8, not from CONEX
  CPPFLAGS += -D__CORSIKA8__
else
  # do not include preshower in C8
  CPPFLAGS += -D__PRESHOW__
  PRESHOW_OBJS := $(OBJDIR)/preshw.o $(OBJDIR)/utils.o $(OBJDIR)/veto.o
endif

ifeq ("$(SYSTEM)", "Darwin")
   LDFLAGS += -bind_at_load -dynamic
endif
ifeq ("$(SYSTEM)", "Linux")
   TRAPFPE =  src/trapfpe.c 
   SOFLAGS = -fPIC -ggdb3 -Wall -shared -l$(G2C)
   SOEXT = so
else
  ifeq ("$(SYSTEM)", "Darwin")
    TRAPFPE =  src/trapfpe_darwin.c 
    SOFLAGS = -fPIC -ggdb3 -Wall -dynamiclib -l$(G2C) 
    SOEXT = dylib
    ifndef CX_NO_ROOT
      SOFLAGS += $(ROOTLDFLAGS)
    endif
  endif
endif

LDFLAGS += -lbz2
SOFLAGS += -lbz2
ifneq ("x$(BZIP2_INCLUDE_DIRS)", "x")
 CXXFLAGS += -I${BZIP2_INCLUDE_DIRS}
endif
CXXFLAGS += -std=c++17 -I$(TABDIR)/readLib/externals/compress_stream/include


CXXFLAGS += $(CONEX_EXTENSIONS)
FFLAGS += -w -fno-second-underscore -fno-automatic -Isrc/ -Isrc/urqmd13/
CXXFLAGS += -g -fPIC $(MC3D) -Isrc 
CCFLAGS += -g -fPIC $(MC3D) 

DEBUG = -ggdb3
OPTIMZE = -O2
FFLAGS += $(OPTIMZE) $(DEBUG)
CFLAGS += $(OPTIMZE) $(DEBUG)
CXXFLAGS += $(OPTIMZE) $(DEBUG)

OBJS := $(OBJDIR)/conex_f.o #$(OBJDIR)/leadingInteractionsData.o

ifdef CONEX_EXTENSIONS
 ifndef CX_NO_ROOT
   OBJS += $(OBJDIR)/particle.o $(OBJDIR)/particleDict.o
   OBJS += $(OBJDIR)/conexExtensionsF.o $(OBJDIR)/conexExtensions.o 
 endif
 OBJS += $(OBJDIR)/conexResampling.o 
 OBJS += $(OBJDIR)/Modifier.o
 OBJS += $(LIBDIR)/libresample.a
 CXXFLAGS += -Iresampler
 LIBRESAMPLE_SOURCES := $(wildcard resampler/*.cc) $(wildcard resampler/*.h)
endif

ifeq ("$(SYSTEM)", "Darwin")
 # RU So 31. Mai 12:27:23 CEST 2020 I am not sure any of this is really helpful and needed:
 OBJS += $(OBJDIR)/ConexDynamicInterface.o
 ifndef CX_NO_ROOT
   OBJS += $(OBJDIR)/CxRoot.o
 endif
endif

CPPFILES=al.cpp

OBJS += $(OBJDIR)/Interface.o $(PRESHOW_OBJS)
OBJSS := $(OBJS) $(OBJDIR)/leadingInteractionsData.o
BZ2=".bz2"
INTERFACE_SOURCE := $(TABDIR)/readLib/source/Interface.bz2.cc

PARAMS := $(patsubst src/%.paramin, $(CFGDIR)/%.param, $(wildcard src/*.paramin))


.PHONY: std epos qII qIII dpm dpmjetIII qgsjetIII qgsjetII qgsjet sibyll \
        param_files config_auto \
        commons rebuild dirs help first_target all \
        binary dynamic note

first_target: help

all: epos sibyll qgsjet qgsjetII dpmjetIII qgsjetIII

ifdef CX_NO_ROOT
binary:
else
binary: $(BINDIR)/conex2r
endif

ifdef CORSIKA_8
epos: commons binary dynamic $(LIBDIR)/libCONEXepos.a
qgsjet: commons binary dynamic $(LIBDIR)/libCONEXqgsjet.a
qgsjetII: commons binary dynamic $(LIBDIR)/libCONEXqgsjetII.a
qgsjetIII: commons binary dynamic $(LIBDIR)/libCONEXqgsjetIII.a
dpmjetIII: commons binary dynamic $(LIBDIR)/libCONEXdpmjetIII.a
dpm: commons binary dynamic $(LIBDIR)/libCONEXdpmjetIII.a
qII: commons binary dynamic $(LIBDIR)/libCONEXqgsjetII.a
qIII: commons binary dynamic $(LIBDIR)/libCONEXqgsjetIII.a
sibyll: commons binary dynamic $(LIBDIR)/libCONEXsibyll.a
.NOTPARALLEL: $(LIBDIR)/libCONEXsibyll.a $(LIBDIR)/libCONEXqgsjetII.a $(LIBDIR)/libCONEXdpmjetIII.a
else
epos: commons binary dynamic $(LIBDIR)/libCONEXepos.$(SOEXT) note
qgsjet: commons binary dynamic $(LIBDIR)/libCONEXqgsjet.$(SOEXT) note
qgsjetII: commons binary dynamic $(LIBDIR)/libCONEXqgsjetII.$(SOEXT) note
qgsjetIII: commons binary dynamic $(LIBDIR)/libCONEXqgsjetIII.$(SOEXT) note
dpmjetIII: commons binary dynamic $(LIBDIR)/libCONEXdpmjetIII.$(SOEXT) note
dpm: commons binary dynamic $(LIBDIR)/libCONEXdpmjetIII.$(SOEXT) note
qII: commons binary dynamic $(LIBDIR)/libCONEXqgsjetII.$(SOEXT) note
qIII: commons binary dynamic $(LIBDIR)/libCONEXqgsjetIII.$(SOEXT) note
sibyll: commons binary dynamic $(LIBDIR)/libCONEXsibyll.$(SOEXT) note
.NOTPARALLEL: $(LIBDIR)/libCONEXsibyll.$(SOEXT) $(LIBDIR)/libCONEXqgsjetII.$(SOEXT) $(LIBDIR)/libCONEXdpmjetIII.$(SOEXT)
endif

dynamic: $(LIBDIR)/libCONEXdynamic.a

crossSection:
	$(MAKE) -f Makefile.crossSection

commons: dirs param_files config_auto
param_files: dirs $(PARAMS)

.NOTPARALLEL: dirs param_files config_auto binary dynamic note

# because of the automatic creation of src/conexConfigAuto.h, the dependent target all
# build every time...
config_auto: dirs
	@echo "#ifndef _include_conexconfigauto_h__" > $(OBJDIR)/.conexConfigAuto.h
	@echo "#define _include_conexconfigauto_h__" >> $(OBJDIR)/.conexConfigAuto.h
ifdef CONEX_EXTENSIONS
	@echo "#define CONEX_EXTENSIONS 1" >> $(OBJDIR)/.conexConfigAuto.h
endif
ifdef CONEX_ROOT
	@echo "#define CONEX_CONFIG_PATH \"$(CONEX_ROOT)\"" >> $(OBJDIR)/.conexConfigAuto.h
else
	@echo "#define CONEX_CONFIG_PATH \"$(CURDIR)/cfg\"" >> $(OBJDIR)/.conexConfigAuto.h
endif
	@echo "#define CONEX_DATA_PATH \"$(TABDIR)\"" >> $(OBJDIR)/.conexConfigAuto.h
	@echo "#endif" >> $(OBJDIR)/.conexConfigAuto.h
	@cmp -s $(OBJDIR)/.conexConfigAuto.h src/conexConfigAuto.h; \
	RETVAL=$$?; \
	if [ $$RETVAL -ne 0 ]; then \
            cp -vf $(OBJDIR)/.conexConfigAuto.h src/conexConfigAuto.h; \
        fi

rebuild: clean all

$(CFGDIR)/%.param: src/%.paramin
	@(echo "==[make]==> generating $@ from : $^")
	@cat $^ | sed s%@BZ2@%$(BZ2)%g | sed s%@CONEXTABROOT@%$(TABDIR)%g > $@ ;

$(OBJDIR)/ConexDynamicInterface.o: src/ConexDynamicInterface.cc src/ConexDynamicInterface.h src/conexConfig.h
	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS)  -D_CONEX_PREFIX=\"$(CONEX_PREFIX)\" -D_CONEX_SYSTEM=\"$(SYSTEM)\"  src/ConexDynamicInterface.cc -c -o $@

# the optional compressed data reading 
$(OBJDIR)/Interface.o: $(INTERFACE_SOURCE) 
	$(CXX) -I$(TABDIR)/readLib $(CXXFLAGS) $^ -c -o $@

$(OBJDIR)/CxRoot.o: src/CxRoot.cc src/CxRoot.h src/conexConfig.h src/leadingInteractionsData.cc src/leadingInteractionsData.h
	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS)  -D_CONEX_PREFIX=\"$(CONEX_PREFIX)\" -D_CONEX_SYSTEM=\"$(SYSTEM)\" src/CxRoot.cc -c -o $@

$(BINDIR)/conex2r: src/conex.cc $(OBJDIR)/ConexDynamicInterface.o $(OBJDIR)/CxRoot.o $(OBJDIR)/leadingInteractionsData.o $(OBJDIR)/trapfpe.o src/conexConfig.h
	@(echo "==[make]==> compiling CONEX-ROOT interface")
	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS) $(ROOTLDFLAGS) $(LDFLAGS) -D_CONEX_PREFIX=\"$(CONEX_PREFIX)\" -D_CONEX_SYSTEM=\"$(SYSTEM)\" $(OBJDIR)/ConexDynamicInterface.o $(OBJDIR)/CxRoot.o $(OBJDIR)/leadingInteractionsData.o -l$(G2C) src/conex.cc -o $@

$(OBJDIR)/conex_f.o: src/conex2r.F src/conexConfig.h
	$(FC) $(CPPFLAGS) $(FFLAGS) -c src/conex2r.F -o $(OBJDIR)/conex_f.o

$(OBJDIR)/Modifier.o: src/Modifier.h src/Modifier.cc src/conexConfig.h
	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS) -c src/Modifier.cc -o $(OBJDIR)/Modifier.o

$(OBJDIR)/leadingInteractionsData.o: src/leadingInteractionsData.h src/leadingInteractionsData.cc src/conexConfig.h
	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS) -c src/leadingInteractionsData.cc -o $(OBJDIR)/leadingInteractionsData.o

$(OBJDIR)/conexExtensions.o: src/conexExtensions.h src/conexExtensions.cc src/conexConfig.h
	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS) -c src/conexExtensions.cc -o $(OBJDIR)/conexExtensions.o

$(OBJDIR)/conexExtensionsF.o: src/conexExtensionsF.h src/conexExtensionsF.cc src/conexConfig.h
	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS) -c src/conexExtensionsF.cc -o $(OBJDIR)/conexExtensionsF.o

$(OBJDIR)/conexResampling.o: src/conexResampling.h src/conexResampling.cc src/conexConfig.h
	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS) -c src/conexResampling.cc -o $(OBJDIR)/conexResampling.o

$(OBJDIR)/particleDict.cc: src/particle.h src/particleLinkDef.h
	rootcint -f $@ -c $(CXXFLAGS) $(ROOTCXXFLAGS) $^ 

$(OBJDIR)/particleDict.o: $(OBJDIR)/particleDict.cc
	$(CXX) $(CXXFLAGS) -I. $(ROOTCXXFLAGS) -c $^ -o $@

$(OBJDIR)/particle.o: src/particle.cc
	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS) -c $^ -o $@

$(OBJDIR)/preshw.o: src/preshw.c
	$(CC) $(CCFLAGS)  -D_GNU_SOURCE -c $^ -o $@

$(OBJDIR)/utils.o: src/utils.c
	$(CC) $(CCFLAGS)  -D_GNU_SOURCE -c $^ -o $@

$(OBJDIR)/veto.o: src/veto.c
	$(CC) $(CCFLAGS)  -D_GNU_SOURCE -c $^ -o $@

$(OBJDIR)/eposmemory.o: src/eposmemory.c
	$(CC) $(CCFLAGS)  -D_GNU_SOURCE -c $^ -o $@

$(OBJDIR)/trapfpe.o: $(TRAPFPE)
	$(CC) $(CCFLAGS)  -D_GNU_SOURCE -c $^ -o $@

$(CPPFILES:%.cpp=$(OBJDIR)%.o) : $(OBJDIR)%.o : src/%.cpp  
	$(CXX) -O2 $(CXXFLAGS) -o $@ -c $<

$(OBJDIR)/conex_sub_sibyll.o: src/conex_sub.F  src/conexConfig.h
	@(echo "==[make]==> compiling main CONEX routines for SIBYLL (this may take a while ...)")
	$(FC) $(FFLAGS) $(CPPFLAGS) -D__SIBYLL21__ -c src/conex_sub.F -o $@
$(OBJDIR)/conex_sub_epos.o: src/conex_sub.F src/conexConfig.h
	@(echo "==[make]==> compiling main CONEX routines for EPOS (this may take a while ...)")
	$(FC) $(FFLAGS)  $(CPPFLAGS) -D__EPOS__ -c src/conex_sub.F -o $@
$(OBJDIR)/conex_sub_qgsjetII.o: src/conex_sub.F src/conexConfig.h
	@(echo "==[make]==> compiling main CONEX routines for QGSJETII (this may take a while ...)")
	$(FC) $(FFLAGS) $(CPPFLAGS) -D__QGSJETII__ -c src/conex_sub.F -o $@
$(OBJDIR)/conex_sub_qgsjetIII.o: src/conex_sub.F src/conexConfig.h
	@(echo "==[make]==> compiling main CONEX routines for QGSJETIII (this may take a while ...)")
	$(FC) $(FFLAGS) $(CPPFLAGS) -D__QGSJETIII__ -c src/conex_sub.F -o $@
$(OBJDIR)/conex_sub_dpmjetIII.o: src/conex_sub.F src/conexConfig.h
	@(echo "==[make]==> compiling main CONEX routines for DPMJETIII (this may take a while ...)")
	$(FC) $(FFLAGS) $(CPPFLAGS) -D__DPMJET__ -c src/conex_sub.F -o $@
$(OBJDIR)/conex_sub_qgsjet.o: src/conex_sub.F src/conexConfig.h
	@(echo "==[make]==> compiling main CONEX routines for QGSJET01 (this may take a while ...)")
	$(FC) $(FFLAGS) $(CPPFLAGS) -D__QGSJET__ -c src/conex_sub.F -o $@

$(OBJDIR)/conex_mod_sibyll.o: src/conex_mod.F  src/conexConfig.h
	@(echo "==[make]==> compiling main CONEX routines for SIBYLL (this may take a while ...)")
	$(FC) $(FFLAGS) $(CPPFLAGS) -D__SIBYLL21__ -c src/conex_mod.F -o $@
$(OBJDIR)/conex_mod_epos.o: src/conex_mod.F src/conexConfig.h
	@(echo "==[make]==> compiling main CONEX routines for EPOS (this may take a while ...)")
	$(FC) $(FFLAGS)  $(CPPFLAGS) -D__EPOS__ -c src/conex_mod.F -o $@
$(OBJDIR)/conex_mod_qgsjetII.o: src/conex_mod.F src/conexConfig.h
	@(echo "==[make]==> compiling main CONEX routines for QGSJETII (this may take a while ...)")
	$(FC) $(FFLAGS) $(CPPFLAGS) -D__QGSJETII__ -c src/conex_mod.F -o $@
$(OBJDIR)/conex_mod_qgsjetIII.o: src/conex_mod.F src/conexConfig.h
	@(echo "==[make]==> compiling main CONEX routines for QGSJETIII (this may take a while ...)")
	$(FC) $(FFLAGS) $(CPPFLAGS) -D__QGSJETIII__ -c src/conex_mod.F -o $@
$(OBJDIR)/conex_mod_dpmjetIII.o: src/conex_mod.F src/conexConfig.h
	@(echo "==[make]==> compiling main CONEX routines for DPMJETIII (this may take a while ...)")
	$(FC) $(FFLAGS) $(CPPFLAGS) -D__DPMJET__ -c src/conex_mod.F -o $@
$(OBJDIR)/conex_mod_qgsjet.o: src/conex_mod.F src/conexConfig.h
	@(echo "==[make]==> compiling main CONEX routines for QGSJET01 (this may take a while ...)")
	$(FC) $(FFLAGS) $(CPPFLAGS) -D__QGSJET__ -c src/conex_mod.F -o $@

$(LIBDIR)/libCONEXdynamic.a: $(OBJDIR)/ConexDynamicInterface.o
	@(echo "==[make]==> compiling particle resampling routines ...")
	$(AR) r $@ $^

$(LIBDIR)/libresample.a: $(LIBRESAMPLE_SOURCES)
	@(echo "==[make]==> compiling particle resampling routines ...")
	$(MAKE) -C resampler
	ln -fs ../../resampler/libresample.a $@

clean:
	@rm -rf $(CONEX_PREFIX)/obj
	@rm -rf $(CONEX_PREFIX)/lib
	@rm -rf $(BINDIR)
	@rm -rf $(CFGDIR)
	$(MAKE) clean -C resampler

dirs:
	@mkdir -p $(OBJDIR)
	@mkdir -p $(BINDIR)
	@mkdir -p $(LIBDIR)
	@mkdir -p $(CFGDIR)


# the dynamic libaries
$(LIBDIR)/libCONEXqgsjetII.$(SOEXT): $(OBJS) $(OBJDIR)/conex_sub_qgsjetII.o $(OBJDIR)/conex_mod_qgsjetII.o
	@(echo "")
	@(echo "==[make]==> creating QGSJETII shared library")
	$(LD) $^ -o $@ $(SOFLAGS)
$(LIBDIR)/libCONEXqgsjetIII.$(SOEXT): $(OBJS) $(OBJDIR)/conex_sub_qgsjetIII.o $(OBJDIR)/conex_mod_qgsjetIII.o
	@(echo "")
	@(echo "==[make]==> creating QGSJETIII shared library")
	$(LD) $^ -o $@ $(SOFLAGS)
$(LIBDIR)/libCONEXdpmjetIII.$(SOEXT): $(OBJS) $(OBJDIR)/conex_sub_dpmjetIII.o $(OBJDIR)/conex_mod_dpmjetIII.o
	@(echo "")
	@(echo "==[make]==> creating DPMJETIII shared library")
	$(LD) $^ -o $@ $(SOFLAGS)
$(LIBDIR)/libCONEXqgsjet.$(SOEXT): $(OBJS) $(OBJDIR)/conex_sub_qgsjet.o $(OBJDIR)/conex_mod_qgsjet.o
	@(echo "")
	@(echo "==[make]==> creating QGSJET01 shared library")
	$(LD)  $^ -o $@ $(SOFLAGS)
$(LIBDIR)/libCONEXsibyll.$(SOEXT): $(OBJS) $(OBJDIR)/conex_sub_sibyll.o $(OBJDIR)/conex_mod_sibyll.o
	@(echo "")
	@(echo "==[make]==> creating SIBYLL shared library")
	$(LD) $^ -o $@ $(SOFLAGS)
$(LIBDIR)/libCONEXepos.$(SOEXT): $(OBJS) $(OBJDIR)/conex_sub_epos.o $(OBJDIR)/conex_mod_epos.o $(CPPFILES:%.cpp=$(OBJDIR)%.o) $(OBJDIR)/eposmemory.o
	@(echo "")
	@(echo "==[make]==> creating EPOS shared library")
	$(LD) $^ -o $@ $(SOFLAGS)

# the static libraries for CORSIKA8
$(LIBDIR)/libCONEXqgsjetII.a: $(OBJSS) $(OBJDIR)/conex_sub_qgsjetII.o $(OBJDIR)/conex_mod_qgsjetII.o
	@(echo "")
	@(echo "==[make]==> creating QGSJETII static library")
	$(AR) r $@ $^
$(LIBDIR)/libCONEXqgsjetIII.a: $(OBJSS) $(OBJDIR)/conex_sub_qgsjetIII.o $(OBJDIR)/conex_mod_qgsjetIII.o
	@(echo "")
	@(echo "==[make]==> creating QGSJETIII static library")
	$(AR) r $@ $^
$(LIBDIR)/libCONEXdpmjetIII.a: $(OBJSS) $(OBJDIR)/conex_sub_dpmjetIII.o $(OBJDIR)/conex_mod_dpmjetIII.o
	@(echo "")
	@(echo "==[make]==> creating DPMJETIII static library")
	$(AR) r $@ $^
$(LIBDIR)/libCONEXqgsjet.a: $(OBJSS) $(OBJDIR)/conex_sub_qgsjet.o $(OBJDIR)/conex_mod_qgsjet.o
	@(echo "")
	@(echo "==[make]==> creating QGSJET01 static library")
	$(AR) r $@ $^
$(LIBDIR)/libCONEXsibyll.a: $(OBJSS) $(OBJDIR)/conex_sub_sibyll.o $(OBJDIR)/conex_mod_sibyll.o
	@(echo "")
	@(echo "==[make]==> creating SIBYLL static library")
	$(AR) r $@ $^
$(LIBDIR)/libCONEXepos.a: $(OBJSS) $(OBJDIR)/conex_sub_epos.o $(OBJDIR)/conex_mod_epos.o $(CPPFILES:%.cpp=$(OBJDIR)%.o) $(OBJDIR)/eposmemory.o
	@(echo "")
	@(echo "==[make]==> creating EPOS static library")
	$(AR) r $@ $^


help:
	@echo " "
	@echo "  CONEX installation Makefile "
	@echo " "
	@echo "  Select interaction models you want to compile "
	@echo "  Options are: 'epos', 'sibyll', 'dpmjetIII', 'qgsjet', 'qgsjetII', 'qgsjetIII', "
	@echo "                'all' selects all of them"
	@echo " "
	@echo "  Edit 'src/conexConfig.h' for special options"
	@echo " "
	@echo "  Set CORSIKA_DATA to the location of the data tables."
	@echo "  See README for more information "
	@echo " "

note:
	@echo " "
	@echo "  CONEX installed successfully "
	@echo " "
	@echo "  Run with: bin/conex2r -h "
	@echo " "
	@echo "  Data tables assumed in: $(TABDIR)"
	@echo "  If you want to change that, use CORSIKA_DATA variable and re-compile. "
	@echo " "
	@echo "  See README for more information "
	@echo " "

dist:
	@echo "==[preparing files]==> this will take a while ...";
	@rm -fr conex2r$(VER)
	@mkdir conex2r$(VER)
	@cp -r CMakeLists.txt Makefile README macros src tabs resampler conex2r$(VER)/
	@rm -rf conex2r$(VER)/src/*.h.in
	@find conex2r$(VER) -name "tmp" -type d -print | xargs rm -rf
	@find conex2r$(VER) -name "CVS" -type d -print | xargs rm -rf
	@find conex2r$(VER) -name ".svn" -type d -print | xargs rm -rf
	@find conex2r$(VER) -name ".git*" -type d -print | xargs rm -rf
	@find conex2r$(VER) -name "\#*" -type f -print | xargs rm -f
	@find conex2r$(VER) -name "*~" -type f -print | xargs rm -f
	@find conex2r$(VER) -name "*.o" -type f -print | xargs rm -f 
	@find conex2r$(VER) -name "*.so" -type f -print | xargs rm -f 
	@find conex2r$(VER) -name "*.a" -type f -print | xargs rm -f 
	@rm -fr conex2r$(VER)/tabs/QGSJetII # too large
	@rm -fr conex2r$(VER)/tabs/QGSJetIII # too large
	@rm -fr conex2r$(VER)/tabs/CONEX/*.*flu* # not needed publicly
	@rm -fr conex2r$(VER)/tabs/CONEX/*.*21* # not needed publicly
	@rm -fr conex2r$(VER)/tabs/CONEX/*.*22* # not needed publicly
	@rm -fr conex2r$(VER)/tabs/CONEX/*.*23* # not needed publicly
	@rm -fr conex2r$(VER)/tabs/CONEX/*.*24* # not needed publicly
	@rm -fr conex2r$(VER)/tabs/CONEX/*.*25* # not needed publicly
	@rm -fr conex2r$(VER)/tabs/CONEX/*.*26* # not needed publicly
	@rm -fr conex2r$(VER)/tabs/CONEX/*.*27* # not needed publicly
	@rm -fr conex2r$(VER)/tabs/CONEX/*.*28* # not needed publicly
	@rm -fr conex2r$(VER)/tabs/CONEX/*.*29* # not needed publicly
	@rm -fr conex2r$(VER)/tabs/CONEX/*.p2* # not needed publicly
	@echo "==[packing files]==> this will take another while ...";
	@tar -zcf conex2r$(VER).tgz conex2r$(VER)
	@echo "==[finished]==>";
