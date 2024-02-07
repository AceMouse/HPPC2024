#include <iostream>
#include <fstream>
#include <stdio.h>

double I = 1;          // Infected
double N = 1000;       // Population
double S = N - I;      // Suspectible but not yet infected
double R = N - S - I;  // Recovered and now immune

double infection_rate = 1./5.;  // Beta
double recovery_rate = 1./10.;  // Gamma

double step_size = .5;
double steps = 10000;

// Forward Euler Method
void sir(FILE *outFile, double S, double I, double R, int32_t steps) {
  if (steps == 0 || I < 1) return;

  fprintf(outFile, "%.3f %.3f %.3f\n", S, I, R);

  double reduction_rate = S/N;
  double dSdt = -infection_rate * I * reduction_rate;
  double S_new = S + dSdt * step_size;
  double dIdt = infection_rate * I * reduction_rate - recovery_rate * I;
  double I_new = I + dIdt * step_size;
  double dRdt = recovery_rate * I;
  double R_new = R + dRdt * step_size;
  sir(outFile, S_new, I_new, R_new, steps - 1);
}

int main (int argc, char *argv[]) {
  FILE *outFile;
  outFile = fopen("sir_output.txt", "w");

  // Write Header
  fprintf(outFile, "S I R\n");

  sir(outFile, S, I, R, steps);
  fclose(outFile);

  return 0;
}
