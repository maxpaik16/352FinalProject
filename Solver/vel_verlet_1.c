/* TODO:
1. the function takes x, t, parameters
2. Two species problem - x and V
3. z = x, V at the same time
4. 
IT WORKS
*/

#include <stdio.h>
#include <stdlib.h>

#define MAX (10000)


// how the function works = it takes arrays of x, t, parameters, and array to store returning values

void velocity_verlet(double* z, double* t, void (*dzdt)(double*,double*,double*,double*), 
       double* parameters, double dt, int steps);

/* Euler step to compute x and V
void step(double* results_z, double* z, double* farray, double dt){
  for(int i; i<2; i++){
    *(results_z+k) = *(z+k) + *(farray+k) * dt;
  }
} */

void f(double* z, double* t, double* parameters, double* fcomp){
  fcomp[0] = -z[0];
}

int main(){
  int i,k;
  double n[2*MAX], t[MAX];
  double dt=0.5;
  double fparams[] = {1.0, 1.0, 1.0}; // tc

  n[0]=0.0, n[1]=10.0;
  t[0]=0.0;

  velocity_verlet(n,t,&f,fparams,dt,MAX);

  for(i=0; i<MAX; i++ ) {
    for(k=0; k<2; k++){
      printf("%f\t", n[2*i+k]);
    }
    printf("%f\n", t[i]);
  }

  return 0;
}

void velocity_verlet(double* z, double* t, void (*dzdt)(double*,double*,double*,double*), 
       double* parameters, double dt, int steps){
  int i;            // counter for loop
  double f[2];      // to store two f at a time
  double halfV;        // half V

  // compute right-hand function at x0, t0
  // pass pointers to arrays

  // compute the velocity at t = 1/2
  // z[0] is initial x, z[1] is initial V

  for(i=0; i<steps-1; i++){
    (*dzdt)(z + 2*i, t + 2*i, parameters, f);
    halfV = *(z+ 2*i + 1) + 0.5 * f[0] * dt;
    
    // x_n+1 = x_n + h * V_(n+1/2)
    *(z + 2*(i+1)) = *(z + 2*i) + halfV * dt;

    // F(x_(n+1))
    (*dzdt)(z + 2*(i+1), t + 2*(i+1), parameters,f+1);

    // V_n+1 = V_(n+1/2) + 0.5 * h * F(x_(n+1))
    *(z + 2*(i+1) + 1) = halfV + 0.5 * dt * f[1];
    *(t+i+1) = *(t+i) + dt;
  }

}