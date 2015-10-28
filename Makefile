# Include configure output
include config.mk
CF = $(CPPFLAGS) $(CXXFLAGS)

progs =
tests =

# Default rule
%.o: %.cc
	$(CXX) $(CF) -MMD -c $< -o $@

# Base source files
base_src=Geom2d.cc QuadEdge.cc Delaunay.cc Delaunay_IO.cc
delaunay_src=$(base_src) DoDelaunay.cc
-include $(delaunay_src:.cc=.d)

# Delaunay prog
progs += DoDelaunay
DoDelaunay: $(delaunay_src:.cc=.o) Serial.o
	$(CXX) $(LDFLAGS) $^ -o $@

# Tests
tests += QuadEdge_test
QuadEdge_test: QuadEdge_test.o QuadEdge.o Geom2d.o
	$(CXX) $(LDFLAGS) $^ -o $@

ifdef MPICXX
# Parallel source
Parallel.o: Parallel.cc
	$(MPICXX) $(CF) -DUSE_PARALLEL -MMD -c $< -o $@
-include Parallel.d

# Parallel Delaunay prog
progs += pDoDelaunay
pDoDelaunay: $(delaunay_src:.cc=.o) Parallel.o
	$(MPICXX) $(LDFLAGS) $^ -o $@
endif

all: $(progs)
.DEFAULT_GOAL=all

# check: Build and execute tests
SHELL=/bin/bash
.PHONY: check
check: $(tests)
	status=0; \
	for t in $(tests); do ./$$t || let status += 1; done; exit $$status;


.PHONY: clean
clean:
	rm -f *.o *.d $(progs)
