CPPC=g++
CPPFLAGS=-std=c++14 -Wall -Wextra -Werror -g -O2

all: test

%.o: %.cc
	$(CPPC) $(CPPFLAGS) -MMD -c $< -o $@

SRC=Geom2d.cc QuadEdge.cc Delaunay.cc

-include $(SRC:.cc=.d)

test: test.cc $(SRC:.cc=.o)
	$(CPPC) $(CPPFLAGS) $^ -o $@

.PHONY: clean
clean:
	rm -f *.o *.d
