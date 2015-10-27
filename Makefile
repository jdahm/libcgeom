include config.mk

all: test_serial test_parallel

%.o: %.cc
	$(CPPC) $(CFLAGS) -MMD -c $< -o $@

Parallel.o: Parallel.cc
	$(CPPC) $(CFLAGS) -DUSE_PARALLEL -MMD -c $< -o $@

BASE_SRC=Geom2d.cc QuadEdge.cc Delaunay.cc Delaunay_IO.cc

-include $(SRC:.cc=.d)

test_serial: test.cc $(BASE_SRC:.cc=.o) Serial.o
	$(CPPC) $(CFLAGS) $(LFLAGS) $^ -o $@

test_parallel: test.cc $(BASE_SRC:.cc=.o) Parallel.o
	$(MPICPPC) $(CFLAGS) $(LFLAGS) $^ -o $@

.PHONY: clean
clean:
	rm -f *.o *.d
