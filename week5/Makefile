CXX := nvc++


#######################################################
# Optimization flags are chosen as the last definition.
# Comment out using "#" at the begining of the line or rearrange according to your needs.
#
# Fastest executable (-ffast-math removes checking for NaNs and other things)
OPT=-O3 -fast

# Add profiling to code
#OPT=-O1 -pg

# Faster compilation time
#OPT=-O1

INCLUDE := -I../include

CXXFLAGS := $(INCLUDE) $(OPT) -Wall -march=native -g -std=c++14

# use on CPU for testing and development
#ACC := -acc=multicore -Minfo=acc

# use on GPU for benchmarking and real tests of code
ACC := -acc -gpu=cuda11.6 -Minfo=acc

.PHONY: clean all

all: sw_sequential sw_parallel sw_par sw_par_async

sw_sequential: sw_sequential.cpp
	$(CXX) $(CXXFLAGS) sw_sequential.cpp -o sw_sequential

sw_parallel: sw_parallel.cpp
	$(CXX) $(CXXFLAGS) $(ACC) sw_parallel.cpp -o sw_parallel

sw_par: sw_par.cpp
	$(CXX) $(CXXFLAGS) $(ACC) sw_par.cpp -o sw_par
sw_par_async: sw_par_async.cpp
	$(CXX) $(CXXFLAGS) $(ACC) sw_par_async.cpp -o sw_par_async
clean:
	rm -f sw_sequential sw_parallel sw_par sw_par_async
