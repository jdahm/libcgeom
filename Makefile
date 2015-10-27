include config.mk

all: DoDelaunay pDoDelaunay

%.o: %.cc
	$(CPPC) $(CFLAGS) -MMD -c $< -o $@

Parallel.o: Parallel.cc
	$(CPPC) $(CFLAGS) -DUSE_PARALLEL -MMD -c $< -o $@
-include Parallel.d

BASE_SRC=Geom2d.cc QuadEdge.cc Delaunay.cc Delaunay_IO.cc
-include $(BASE_SRC:.cc=.d)

DoDelaunay: DoDelaunay.cc $(BASE_SRC:.cc=.o) Serial.o
	$(CPPC) $(CFLAGS) $(LFLAGS) $^ -o $@

pDoDelaunay: DoDelaunay.cc $(BASE_SRC:.cc=.o) Parallel.o
	$(MPICPPC) $(CFLAGS) $(LFLAGS) $^ -o $@

.PHONY: clean
clean:
	rm -f *.o *.d
