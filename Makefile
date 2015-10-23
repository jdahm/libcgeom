CPPC=g++
CPPFLAGS=-Wall -Wextra -Werror -g -O2

%.o: %.cc
	$(CPPC) $(CPPFLAGS) -MMD -c $< -o $@

SRC=QuadEdge.cc Delaunay.cc Subdivision.cc

-include $(SRC:.cc=.d)

test: test.cc $(SRC:.cc=.o)
	$(CPPC) $(CPPFLAGS) $^ -o $@

.PHONY: clean
clean:
	rm -f *.o
