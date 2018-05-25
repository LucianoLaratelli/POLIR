CC = g++

SRCS = $(wildcard *.cpp)

CFLAGS = -Wall -Wextra -Werror -ggdb3 -O2 -funroll-loops -std=c++11

PROGS = $(patsubst %.cpp,%,$(SRCS))

all: $(PROGS)

%: %.cpp
	$(CC) $(CFLAGS) -o $@ $< 

polir: ls4equil.f90 r4.f ps_intra_polir_b.f90
	pgf90 -fast ls4equil.f90 r4.f ps_intra_polir_b.f90 -o water.equil

debug_polir: ls4equil.f90 r4.f ps_intra_polir_b.f90
	pgf90 -g ls4equil.f90 r4.f ps_intra_polir_b.f90 -o water.equil

clean:
	rm -rf *.o water.equil
	rm -f $(PROGS)
