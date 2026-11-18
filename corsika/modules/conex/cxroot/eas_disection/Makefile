LDFLAGS+=-Wl,--no-as-needed 

# coast
LDFLAGS += -fPIC -L$(COAST_DIR)/lib -L$(COAST_DIR)/lib/unknown -lCorsikaFileIO -lCorsikaIntern 
CXXFLAGS += -I$(COAST_DIR)/include

CXXFLAGS+=$(shell root-config --cflags) -ggdb3 
LDFLAGS+=$(shell root-config --libs) -ggdb3

LDFLAGS+=-lmpfr -lgfortran


FC=gfortran

all: overlay forSTACKIN stackinPlots

#%.o: %.cc
#	$(CXX) -c $^ $(CXXFLAGS) $(LDFLAGS) -o $@

forSTACKIN: idtrafocx.o
#forSTACKIN: forSTACKIN.o

stackinPlots: stackinPlots.cc
	$(CXX) $^ $(CXXFLAGS) $(LDFLAGS) -o $@

clean:
	rm -f overlay forSTACKIN stackinPlots *~ *.o
