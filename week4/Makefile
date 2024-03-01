CXX := clang++-12
CXX := g++ 

#######################################################
# Optimization flags are chosen as the last definition.
# Comment out using "#" at the begining of the line or rearrange according to your needs.
#
# Fastest executable (-ffast-math removes checking for NaNs and other things)
OPT=-O3 -ffast-math

# Add thread checking to code
#OPT=-O1 -fsanitize=thread -fPIE -pie

# Faster compilation time
#OPT=-O1

CXXFLAGS := $(OPT) -Wall -march=native -g -pg -std=c++14

default: seq mp mpnc mpps

seq: seismogram_seq.cpp
	$(CXX) seismogram_seq.cpp $(CXXFLAGS) -o seq

mp: seismogram_omp.cpp
	$(CXX) seismogram_omp.cpp $(CXXFLAGS) -fopenmp -o mp

mpnc: seismogram_no_complex_omp_vec.cpp
	$(CXX) seismogram_no_complex_omp_vec.cpp $(CXXFLAGS) -fopenmp -o mpnc

mpps: seismogram_par_sec.cpp
	$(CXX) seismogram_par_sec.cpp $(CXXFLAGS) -fopenmp -o mpps 

clean:
	rm -fr seq mp mpnc mpps 

