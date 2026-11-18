CXXFLAGS+=$(shell root-config --cflags)
LDFLAGS+=-Wl,--no-as-needed $(shell root-config --libs)

all: merge

merge: merge.C
	$(CXX) -o $@ $(CXXFLAGS) $(LDFLAGS) $^

clean:
	rm -f merge *~ *.so *.d
