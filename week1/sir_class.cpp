#include <iostream>
#include <fstream>
#include <cstring>


class Epidemic {
    public:
        // population
        int N = 1000;
        double I = 1;
        double S = N-I;
        double R = N - S - I;

        // step/step size
        double const step_size = 0.5;

        // params
        double const beta = 0.2;
        double const gamma = 0.1;

        void update() {
            S += dS() * step_size;
            I += dI() * step_size;
            R += dR() * step_size;
        };

    private:
        // diff equations
        double dS() {
            return -beta*I*(S/N);
        };
        double dI() {
            return beta*I*(S/N) - gamma*I;
        };
        double dR() {
            return gamma*I;
        };
};


int main() {

    // init epidemic
    Epidemic epi;

    // init file and header
    std::ofstream myfile;
    myfile.open("sir_output.txt");    
    myfile << "S\t I\t R\n";

    // write content 
    do {
        myfile << epi.S << "\t" << epi.I << "\t" << epi.R << "\n";
        epi.update();
    }
    while (epi.I >= 1);
    myfile.close();

    return 0;
};