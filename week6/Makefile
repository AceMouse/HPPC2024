LIB = -L/opt/openmpi-4.1.6/lib -Wl,-rpath -Wl,/opt/openmpi-4.1.6/lib -Wl,--enable-new-dtags -lmpi
INC = -I/opt/openmpi-4.1.6/include

FLAGS := -O3 -DNDEBUG -g -Wall

#FLAGS := -Os -DNDEBUG -g -Wall -pg

FLAGS := $(FLAGS) -fopenmp

.PHONY: clean all

all: fwc_sequential fwc_parallel fwc_parallel_slab

fwc_sequential: fwc_sequential.cpp
	h5c++ $(FLAGS) $(LIB) $(INC) -o fwc_sequential fwc_sequential.cpp

fwc_parallel: fwc_parallel.cpp
	h5c++ $(FLAGS) $(LIB) $(INC) -o fwc_parallel fwc_parallel.cpp

fwc_parallel_slab: fwc_parallel_slab.cpp
	h5c++ $(FLAGS) $(LIB) $(INC) -o fwc_parallel_slab fwc_parallel_slab.cpp

clean:
	rm -f fwc_sequential fwc_parallel fwc_parallel_slab *.o
