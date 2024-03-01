/*
  Assignment: Make an OpenMP parallelised wave propagation
  model for computing the seismic repsonse for a wave
  propagating through a horizontally stratified medium
*/
#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <chrono>
#include <thread>
#include <complex>
#include <cmath>
#include <omp.h>


// ======================================================
// The number of frequencies sets the cost of the problem
const long NTHREADS=1;            // number of threads
const long NFREQ=1024*1024;         // number of frequencies per core
const long nfreq=NFREQ*NTHREADS;  // frequencies in spectrum

// ======================================================
template <class T> class NUMA_Allocator {
public:
  typedef T* pointer;
  typedef const T* const_pointer;
  typedef T& reference;
  typedef const T& const_reference;
  typedef size_t size_type;
  typedef T value_type;

  NUMA_Allocator() { }
  NUMA_Allocator(const NUMA_Allocator& _r) { }
  ~NUMA_Allocator() { }

  // allocate raw memory including page placement
  pointer allocate(size_t numObjects, const void *localityHint=0) {
    size_t len = numObjects * sizeof(T);
    void *m = std::malloc(len);
    char *p = static_cast<char*>(m);
    if(!omp_in_parallel()) {
#pragma omp parallel for schedule(static)
      for(size_t i=0; i<len; i+=sizeof(T)) {
        for(size_t j=0; j<sizeof(T); ++j) {
	  p[i+j]=0;
	}
      }
    }
    
    return static_cast<T*>(m);
  }

  // free raw memory
  void deallocate(pointer ptrToMemory, size_t numObjects) {
    std::free(ptrToMemory);
  }

  // construct object at given address
  void construct(pointer p, const value_type& x) {
    new(p) value_type(x);
  }

  // destroy object at given address
  void destroy(pointer p) {
    p->~value_type();
  }

private:
  void operator=(const NUMA_Allocator&) {}
};

// shorthand name for complex number type definition below
typedef std::complex<double> Complex;

// shorthand names for vector types

// Use NUMA-aware first touch allocator
//typedef std::vector<Complex, NUMA_Allocator<Complex>> ComplexVector;
//typedef std::vector<double, NUMA_Allocator<double>> DoubleVector;

// Use standard allocator
//typedef std::vector<Complex> ComplexVector;
typedef std::vector<double> DoubleVector;

// Initialize Basic Constants
const double dT=0.001;     // sampling distance
const long nsamp=2*nfreq;  // samples in seismogram

// Frequency resolution (frequency sampling distance)
double dF = 1/(nsamp*dT);

// read data file with one number per line
std::vector<double> read_txt_file(std::string fname) {
    std::vector<double> data;  // vector of data points
    std::string line;          // string to read in each line    
    std::ifstream file(fname); // open file
    while (std::getline(file, line))     // loop over lines until end of file
        data.push_back(std::stod(line)); // convert string to double and store in array
    return data;
}

// Cooley–Tukey FFT (in-place computation)
void fft(std::vector<double>& x_re_vec, std::vector<double>& x_im_vec)
{
	const long N = x_re_vec.size();
	if (N <= 1) return;

    double* x_re = &x_re_vec[0];
    double* x_im = &x_im_vec[0];
	// divide
	std::vector<double> even_re_vec(N/2),even_im_vec(N/2), odd_re_vec(N/2), odd_im_vec(N/2);
    double* even_re = &even_re_vec[0];
    double* even_im = &even_im_vec[0];
    double* odd_re = &odd_re_vec[0];
    double* odd_im = &odd_im_vec[0];

//    #pragma omp parallel for if(N>1000)
    for (long i=0; i<N/2; i++) {//vectorised
        even_im[i] = x_im[2*i];
        even_re[i] = x_re[2*i];
        odd_im[i]  = x_im[2*i+1];
        odd_re[i]  = x_re[2*i+1];
    }
	// conquer
//    #pragma omp task if(N>1000) shared(even_re_vec,even_im_vec)
    {
        fft(even_re_vec, even_im_vec);
    }
//    #pragma omp task if(N>1000) shared(odd_re_vec,odd_im_vec)
    {
        fft(odd_re_vec, odd_im_vec);
    }
    // combine
//    #pragma omp taskwait
//    #pragma omp parallel for if(N>1000)
    for (long k = 0; k < N/2; k++) { //not vectorised: reason = cos and sin
        double theta = -2 * M_PI * k / N;
        double p_re = cos(theta);
        double p_im = sin(theta);
        double t_re = p_re * odd_re[k]-p_im*odd_im[k];
        double t_im = p_re * odd_im[k]+p_im*odd_re[k];
        x_re[k    ] = even_re[k] + t_re;
        x_im[k    ] = even_im[k] + t_im;
        x_re[k+N/2] = even_re[k] - t_re;
        x_im[k+N/2] = even_im[k] - t_im;
    }
}

// inverse fft (in-place)
void ifft(std::vector<double>& x_re_vec, std::vector<double>& x_im_vec)
{
	const long N = x_re_vec.size();

    double* x_re = &x_re_vec[0];
    double* x_im = &x_im_vec[0];
    double inv_size = 1.0 / N;
//#pragma omp parallel
    {
//#pragma omp for
        for (int i = 0; i<N; i++){ //vectorised
            x_im[i] = -x_im[i]; 
        }
        fft(x_re_vec, x_im_vec);  	   // forward fft
//        #pragma GCC ivdep
//#pragma omp for
        for (int i = 0; i<N; i++) { //vectorised
            x_im[i] = -x_im[i] * inv_size;
            x_re[i] =  x_re[i] * inv_size; 
        }
    }
}

// Main routine: propgate wave through layers and compute seismogram
DoubleVector propagator(std::vector<double> wave,
                        std::vector<double> density,
                        std::vector<double> velocity) {
    const long nlayers = density.size();
    std::vector<double> imp(nlayers);      // impedance
    std::vector<double> ref_vec(nlayers-1);    // reflection coefficient
    double* ref = &ref_vec[0];
    std::vector<double> half_filter_re(nfreq/2+1,1); // half filter
    std::vector<double> filter_re(nfreq+1);  // full filter
    DoubleVector half_wave(nfreq+1,0); // half wave
    std::vector<double> wave_spectral_re(nsamp); // FFT(wave)
    std::vector<double> wave_spectral_im(nsamp); // FFT(wave)
    std::vector<double> U_re(nfreq+1,0);     // Upgoing waves
    std::vector<double> U_im(nfreq+1,0);     // Upgoing waves
    std::vector<double> Upad_re(nsamp,0);    // FFT(seismogram)
    std::vector<double> Upad_im(nsamp,0);    // FFT(seismogram)
    DoubleVector seismogram(nsamp); // final seismogram

    long n_wave = wave.size();             // size of wave array
    long lc = std::lround(std::floor(nfreq*0.01)); // low-cut indices
    double mean_wave = 0.;                 // wave zero point

    std::chrono::time_point<std::chrono::system_clock> tstart1,tstart2,tend1,tend2;

    auto tstart = std::chrono::high_resolution_clock::now(); // start time (nano-seconds)
    
    // Compute seismic impedance
    #pragma GCC ivdep
    for (long i=0; i < nlayers; i++){ //vectorised
        imp[i] = density[i] * velocity[i];
    }
    
    // Reflection coefficients at the base of the layers :
    for (long i=0; i < nlayers-1; i++){ //vectorised
        ref[i] = (imp[i+1] - imp[i])/(imp[i+1] + imp[i]);
    }
            // spectrum U of upgoing waves just below the surface.
            // See eq. (43) and (44) in Ganley (1981).
    #pragma omp parallel 
    {
        
        #pragma omp single nowait
        {
            // Spectral window (both low- and high cut)
            for (long i=0; i < lc+1; i++) {//not vectorised: reason = sin 
                half_filter_re[i]= (sin(M_PI*(2*i-lc)/(2*lc)))/2+0.5;
            }

            for (long i=0; i < nfreq/2+1; i++) { //vectorised
                filter_re[i] = half_filter_re[i];
            }

            filter_re[nfreq/2+1] = 1;

            for (long i=nfreq/2+2; i < nfreq+1; i++) { //vectorised
                filter_re[i] = half_filter_re[nfreq+1-i];
            }
        }
        #pragma omp single nowait
        {
            for (long i=0; i < n_wave/2; i++) { //vectorised
                half_wave[i] = wave[n_wave/2-1+i];
            }

            for (long i=0; i < nfreq; i++) { //vectorised
                wave_spectral_re[i] = half_wave[i];
                mean_wave += wave_spectral_re[i];
            }

            for (long i=nfreq; i < 2*nfreq; i++) { //vectorised
                wave_spectral_re[i] = half_wave[2*nfreq-i];
                mean_wave += wave_spectral_re[i];
            }

            mean_wave = mean_wave / nsamp;

            for (long i=0; i < 2*nfreq; i++){ //vectorised
                wave_spectral_re[i] -= mean_wave;
            }
            tstart1 = std::chrono::high_resolution_clock::now(); // start time (nano-seconds)
            // Fourier transform waveform to frequency domain
            fft(wave_spectral_re, wave_spectral_im);
            tend1 = std::chrono::high_resolution_clock::now(); // end time (nano-seconds)
        }
        #pragma omp for
        for (long i=0; i < nfreq+1; i++) { // not vectorised: reason = inner loop
            Complex omega{0, 2*M_PI*i*dF};
            Complex exp_omega = exp( - dT * omega);
            Complex Y = 0;
            for (long n=nlayers-2; n > -1; n--){ // not vectorised: reason = Y_n+1 depends on Y_n
                Y = exp_omega * (ref[n] + Y) / (1.0 + ref[n]*Y);
            }
            U_re[i] = Y.real();
            U_im[i] = Y.imag();
        }
        #pragma omp barrier
        // Compute seismogram
        #pragma omp for 
        for (long i=0; i < nfreq+1; i++) { //vectorised
        //(ac - bd) + i(ad + bc)
            U_re[i] = U_re[i]*filter_re[i]; 
            U_im[i] = U_im[i]*filter_re[i]; 
        }
        #pragma omp for
        for (long i=0; i < nfreq+1; i++) { //vectorised
            Upad_re[i] = U_re[i];
            Upad_im[i] = U_im[i];
        }

        #pragma omp for
        for (long i=nfreq+1; i < nsamp; i++){ //vectorised
            Upad_re[i] = Upad_re[nsamp - i];
            Upad_im[i] = -Upad_im[nsamp - i];
        }

        #pragma omp for
        for (long i=0; i < nsamp; i++){ //vectorised
            Upad_re[i] =Upad_re[i]*wave_spectral_re[i]-Upad_im[i]*wave_spectral_im[i];
            Upad_im[i] =Upad_re[i]*wave_spectral_im[i]+Upad_im[i]*wave_spectral_re[i];
        }
        #pragma omp barrier
        #pragma omp single
        {
            // Fourier transform back again
            tstart2 = std::chrono::high_resolution_clock::now(); // start time (nano-seconds)
            ifft(Upad_re, Upad_im);
            tend2 = std::chrono::high_resolution_clock::now(); // end time (nano-seconds)
        } 
    //#pragma GCC ivdep
        #pragma omp for
        for (long i=0; i < nsamp; i++){ //vectorised
            seismogram[i] = Upad_re[i];
        }
        
    }
    


    auto tend = std::chrono::high_resolution_clock::now(); // end time (nano-seconds)

    std::cout <<  "Wave zero-point        : "  << std::setw(9) << std::setprecision(5) 
              << mean_wave<< "\n";    
    std::cout <<  "Seismogram first coeff : "  << std::setw(9) << std::setprecision(5) 
              << seismogram[0] << ", " << seismogram[1] << ", " << seismogram[2] << ", " << seismogram[3] <<"\n";    
    std::cout <<  "Elapsed time for FFTs  :" << std::setw(9) << std::setprecision(4)
              << (tend1 - tstart1 + tend2 - tstart2).count()*1e-9 << "\n";
    std::cout <<  "Elapsed time without FFTs:" << std::setw(9) << std::setprecision(4)
              << (tend - tstart - (tend1 - tstart1 + tend2 - tstart2)).count()*1e-9 << "\n";
    std::cout <<  "Elapsed time:" << std::setw(9) << std::setprecision(4)
              << (tend - tstart).count()*1e-9 << "\n";
    
    return seismogram;
}

//======================================================================================================
//======================== Main function ===============================================================
//======================================================================================================
int main(int argc, char* argv[]){    
    // Load the wave profile and the density and velocity structure of the rock from text files
    std::vector<double> wave = read_txt_file("wave_data.txt");         // input impulse wave in medium
    std::vector<double> density = read_txt_file("density_data.txt");   // density as a function of depth
    std::vector<double> velocity = read_txt_file("velocity_data.txt"); // seismic wave velocity as a function of depth

    // Propagate wave
    auto seismogram = propagator(wave,density,velocity);
    
    // write output and make checksum
    double checksum=0;
    std::ofstream file("seismogram.txt"); // open file
    for (long i=0; i < nsamp; i++) {
        file << seismogram[i] << '\n';
        checksum += abs(seismogram[i]);
    }
    std::cout <<  "Checksum    :" << std::setw(20) << std::setprecision(15)
              << checksum << "\n";
}

