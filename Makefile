# Include configure output
include config.mk
CF = $(CPPFLAGS) $(CXXFLAGS)

progs =
tests =

# Default rule
%.o: %.C
	$(CXX) $(CF) -MMD -c $< -o $@

# Base source files
base_src=Geom2d.C QuadEdge.C Delaunay.C Delaunay_IO.C
delaunay_src=$(base_src) DoDelaunay.C
-include $(delaunay_src:.C=.d)

# Delaunay prog
progs += DoDelaunay
DoDelaunay: $(delaunay_src:.C=.o) Serial.o
	$(CXX) $(LDFLAGS) $^ -o $@

# Tests
tests += QuadEdge_test
QuadEdge_test: QuadEdge_test.o QuadEdge.o Geom2d.o
	$(CXX) $(LDFLAGS) $^ -o $@
-include QuadEdge_test.d

tests += Delaunay_test
Delaunay_test: Delaunay_test.o Delaunay.o Delaunay_IO.o Serial.o QuadEdge.o Geom2d.o
	$(CXX) $(LDFLAGS) $^ -o $@
-include Delaunay_test.d

ifdef MPICXX
# Parallel source
Parallel.o: Parallel.C
	$(MPICXX) $(CF) -DUSE_PARALLEL -MMD -c $< -o $@
-include Parallel.d

# Parallel Delaunay prog
progs += pDoDelaunay
pDoDelaunay: $(delaunay_src:.C=.o) Parallel.o
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

.PHONY: lint
lint:
	uncrustify -c uncrustify.cfg --no-backup $(shell find . -name "*.C" -or -name "*.H")


.PHONY: clean
clean:
	rm -f *.o *.d $(progs) $(tests)
