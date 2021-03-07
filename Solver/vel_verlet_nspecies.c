/* TODO:
1. the function takes x, t, parameters
2. Two species problem - x and V
3. z = x, V at the same time
4. generalize to n species

1. add interface as ode solver, don't add order yet though

*/

#include <stdio.h>
#include <stdlib.h>

#define MAX (1000)

// how the function works = it takes arrays of x, t, parameters, and array to store returning values
// NOTE that number is the number of pairs x+V


void velocity_verlet(int number, double* z, double* t, void (*dzdt)(double*,double*,double*,double*), 
       double* parameters, double dt, int steps);

/* Euler step to compute x and V
void step(double* results_z, double* z, double* farray, double dt){
  for(int i; i<2; i++){
    *(results_z+k) = *(z+k) + *(farray+k) * dt;
  }
} */

void f(double* z, double* t, double* parameters, double* fcomp){
  fcomp[0] = -z[0];
  fcomp[1] = -z[2];
}

int main(){
  int i,k;
  int number = 2;         // number of species - meaning coordinates x and y
  double n[2*number*MAX], t[MAX];
  double dt=0.5;
  double fparams[] = {1.0, 1.0, 1.0}; // tc

  n[0]=0.0, n[1]=10.0;
  n[2]=10.0, n[3]=0.0;
  t[0]=0.0;

  velocity_verlet(number,n,t,&f,fparams,dt,MAX);

  for(i=0; i<MAX; i++ ) {
    for(k=0; k<2*number; k++){
      printf("%f\t", n[2*number*i+k]);
    }
    printf("%f\n", t[i]);
  }

  return 0;
}

void velocity_verlet(int number, double* z, double* t, void (*dzdt)(double*,double*,double*,double*), 
       double* parameters, double dt, int steps){
  int i,k;                  // counter for loop
  double f[2*number];       // to store two f at a time for each coordinate
  double halfV[number];     // half V

  // compute right-hand function at x0, t0
  // pass pointers to arrays

  // compute the velocity at t = 1/2
  // z[0] is initial x, z[1] is initial V

  for(i=0; i<steps-1; i++){
    // compute all f at initial x
    (*dzdt)(z + 2*i*number, t + i, parameters, f);

    for(k=0; k<number; k++){
      halfV[k] = *(z+ 2*i*number + 2*k + 1) + 0.5 * f[k] * dt;
    }
    
    // x_n+1 = x_n + h * V_(n+1/2)
    for(k=0; k<number; k++){
      *(z + 2*(i+1)*number + 2*k) = *(z + 2*i*number + 2*k) + halfV[k] * dt;
    }

    // F(x_(n+1))
    (*dzdt)(z + 2*(i+1)*number, t + 2*(i+1), parameters,f+number);

    for(k=0; k<number; k++){
      // V_n+1 = V_(n+1/2) + 0.5 * h * F(x_(n+1))
      *(z + 2*(i+1)*number + 2*k + 1) = halfV[k] + 0.5 * dt * f[number+k];
    }

    *(t+i+1) = *(t+i) + dt;
  }

}