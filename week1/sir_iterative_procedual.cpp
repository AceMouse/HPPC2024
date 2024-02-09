#include <iostream>
#include <fstream>
#include <stdio.h>

double I = 1;          // Infected
double N = 1000;       // Population
double S = N - I;      // Suspectible but not yet infected
double R = N - S - I;  // Recovered and now immune

double infection_rate = 1./5.;  // Beta
double recovery_rate = 1./10.;  // Gamma

int d = 1000; //steps pr day. 
double step_size = 1./(double)d;
int steps = 10000000;

int main (int argc, char *argv[]) {
    FILE *outFile;
    outFile = fopen("sir_output.txt", "w");

    fprintf(outFile, "S I R\n");
    double is; 
    double ir; 
    double frac = 1/N; //save on divisions
    while(steps-- && I>=1){
        if ((steps % d) == 0) { // only write to file onece pr day. 
            fprintf(outFile, "%.3f %.3f %.3f\n", S, I, R);
        }
        is = infection_rate * I * S * frac;
        ir = recovery_rate * I;
        S -= is * step_size;
        I += (is - ir) * step_size;
        R += ir * step_size;
    }
    fclose(outFile);
    return 0;
}
