include config.mk
CF = $(CPPFLAGS) $(CXXFLAGS) -I.

progs =
tests =

# Default rule
%.o: %.cpp
	$(MPICXX) $(CF) -MMD -c $< -o $@

# Base source
cgl_base=cgl/geom2d.o cgl/quad_edge.o cgl/point_set.o cgl/subdivision.o cgl/triangulation.o cgl/load_balance.o
ifdef USE_ZOLTAN
cgl_base += cgl/load_balance_zoltan.o
else
cgl_base += cgl/load_balance_nozoltan.o
endif
-include $(cgl_base:.o=.d)

par_base=par/communicator.o par/environment.o par/status.o par/request.o
-include $(par_base:.o=.d)

progs += do_delaunay
do_delaunay: do_delaunay.o $(cgl_base) $(par_base)
	$(MPICXX) $^ -o $@ $(LDFLAGS)
-include do_delaunay.d

progs += generate_test
generate_test: generate_test.o
	$(MPICXX) $^ -o $@ $(LDFLAGS)
-include generate_test.d

# Tests
tests += geom2d_test
geom2d_test: cgl/geom2d_test.o cgl/geom2d.o
	$(MPICXX) $^ -o $@ $(LDFLAGS)
-include cgl/geom2d_test.d

tests += quad_edge_test
quad_edge_test: cgl/quad_edge_test.o cgl/quad_edge.o cgl/geom2d.o
	$(MPICXX) $^ -o $@ $(LDFLAGS)
-include cgl/quad_edge_test.d

tests += point_set_test
point_set_test: cgl/point_set_test.o cgl/point_set.o cgl/geom2d.o cgl/load_balance.o $(par_base)
	$(MPICXX) $^ -o $@ $(LDFLAGS)
-include cgl/point_set_test.d

tests += subdivision_test
subdivision_test: cgl/subdivision_test.o cgl/subdivision.o cgl/geom2d.o cgl/quad_edge.o $(par_base)
	$(MPICXX) $^ -o $@ $(LDFLAGS)
-include cgl/subdivision_test.d

tests += triangulation_test
triangulation_test: cgl/triangulation_test.o cgl/triangulation.o cgl/subdivision.o cgl/point_set.o cgl/geom2d.o cgl/quad_edge.o cgl/load_balance.o $(par_base)
	$(MPICXX) $^ -o $@ $(LDFLAGS)
-include cgl/triangulation_test.d


test: $(tests)

# check: Build and execute tests
SHELL=/bin/bash
.PHONY: check
check: test
	status=0; \
	for t in $(tests); do ./$$t || let status += 1; done; exit $$status;

.PHONY: clean
clean:
	rm -f cgl/*.o cgl/*.d test/* $(progs) $(tests) *.pyc

all: $(progs)
.DEFAULT_GOAL=all

