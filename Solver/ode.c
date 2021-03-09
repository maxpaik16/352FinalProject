/*******************************************************************
 *
 *                    generic integrator
 *                    to compile:
 *                    gcc -shared -O2 ode.c -o libode.so
 *       
 *******************************************************************
*/
#include <stdio.h>
#include <stdlib.h>

/****** function prototypes begin ******/
int solve_ode(double *pt, double *px, double dt, int nsteps, int nvar, int order, double *params,
              void (*dxdt)(double t, double *px, double *params, double *derivs));

void stepEuler(double dt, double t, double a[], double anew[], int nvar,
               double params[], double *buf,
               void (*dxdt)(double t, double a[], double params[], double derivs[]));

void stepEulercromer(double dt, double t, double a[], double anew[], int nvar,
                     double params[], double *buf,
                     void (*dxdt)(double t, double a[], double params[], double derivs[]));

void stepRK2(double dt, double t, double a[], double anew[], int nvar,
             double params[], double *buf,
             void (*dxdt)(double t, double a[], double params[], double derivs[])); 

void stepRK4(double dt, double t, double a[], double anew[], int nvar,
             double params[], double *buf,
             void (*dxdt)(double t, double a[], double params[], double derivs[]));

void stepVelocityVerlet(double dt, double t, double a[], double anew[], int nvar,
                     double params[], double *buf,
                     void (*dxdt)(double t, double a[], double params[], double derivs[]));

/****** function prototypes end ******/



/*****************************
*
* ODE library functions start 
*
******************************/

// Generic ODE solver:
// pt and px should hold nsteps+1 elements 
// nvar  = number of dependent variables 
// order = order of integration: 1 = Euler, -1 = Euler-Cromer, 2 = RK2, 4 = RK4, 13 = Velocity Verlet
int solve_ode(double *pt, double *px, double dt, int nsteps, int nvar, int order, double *params,
               void (*dxdt)(double t, double *px, double *params, double *derivs))
{
  double *buf;
  //pointer to time stepper function
  void (*step_ode)(double dt, double t, double a[], double anew[], int nvar,
                     double params[], double *buf,
                     void (*dxdt)(double t, double a[], double params[], double derivs[]));
  int step;
  
  //allocate memory for storing temporary variables
  buf = calloc( 5*nvar, sizeof(double) );
  
  if(NULL == buf) {
    fprintf(stderr, "solve_ode(): could not allocate memory\n");
    return 1;
  }
  
  //choose solver
  switch(order){
    case 1:
      step_ode = &stepEuler;
      break;
    case -1: 
      if(2!=nvar) {
        fprintf(stderr, "Euler-Cromer only works for 2 variables\n");
        return -1;
      }
      step_ode = &stepEulercromer;
      break;
    case 2: 
      step_ode = &stepRK2;
      break;
    case 4:
      step_ode = &stepRK4;
      break;
    case 13:
    default:
      //use Velocity Verlet stepping by default
      step_ode = &stepVelocityVerlet;      
  }

  for(step = 0; step < nsteps; step++) {
    pt[step+1] = pt[0] + (step+1)* dt;    
    /* take step */
    (*step_ode)(dt, pt[step], &px[step*nvar], &px[(step+1)*nvar], nvar, params, buf, dxdt); 
  }
  //free memory
  free(buf);
  buf = NULL;
  return 0;
}

/******** step  *********************/
/*** Take a step with the Euler method, output anew ***/
/*** Calls function derivs to get derivative ***/
void stepEuler(double dt, double t, double a[], double anew[], int nvar,
               double params[], double *buf,
               void (*dxdt)(double t, double a[], double params[], double derivs[]))
{
    int i ;
    double *dadt = buf ; 
    (*dxdt)(t, a, params, dadt) ; 
    for (i=0 ; i<nvar; i++) anew[i] = a[i] + dt*dadt[i] ; 
}

void stepEulercromer(double dt, double t, double a[], double anew[], int nvar,
                     double params[], double *buf,
                     void (*dxdt)(double t, double a[], double params[], double derivs[]))
{
    int i ;
    double *dadt = buf;
    (*dxdt)(t, a, params, dadt) ; 
    anew[1] = a[1] + dt*dadt[1] ; 
    anew[0] = a[0] + anew[1]*dt ; 
    
}

void stepRK2(double dt, double t, double a[], double anew[], int nvar,
             double params[], double *buf,
             void (*dxdt)(double t, double a[], double params[], double derivs[]))
{
    int i ;
    double *f1 = buf, *f2 = f1+nvar; 
    double *a2 = f2 + nvar; 
    (*dxdt)(t, a, params, f1); 
    for (i=0 ; i<nvar ; i++) a2[i] = a[i] + (0.5*dt)*f1[i]; 
    (*dxdt)(t+0.5*dt, a2, params, f2) ; 
    for (i=0 ; i<nvar ; i++) anew[i] = a[i] + dt*f2[i]; 
}

void stepRK4(double dt, double t, double a[], double anew[], int nvar,
             double params[], double *buf,
             void (*dxdt)(double t, double a[], double params[], double derivs[]))
{
    int i ;
    double *f1 = buf, *f2 = f1+nvar, *f3 = f2+nvar, *f4 = f3+nvar, *ax = f4+nvar; 
    const double sixth = 0.16666666666666666667; //save 1./6. to avoid recomputing division at every step
    
    (*dxdt)(t, a, params, f1); 
    for (i=0 ; i<nvar ; i++) ax[i] = a[i] + (0.5*dt)*f1[i] ; 
    (*dxdt)(t+0.5*dt, ax, params, f2); 
    for (i=0 ; i<nvar ; i++) ax[i] = a[i] + (0.5*dt)*f2[i] ; 
    (*dxdt)(t+0.5*dt, ax, params, f3); 
    for (i=0 ; i<nvar ; i++) ax[i] = a[i] + (dt)*f3[i] ; 
    (*dxdt)(t+dt, ax, params, f4); 
    
    for (i=0 ; i<nvar ; i++) 
      anew[i] = a[i] + sixth*(f1[i]+2.0*f2[i]+2.0*f3[i]+f4[i])*dt; 
}

void stepVelocityVerlet(double dt, double t, double a[], double anew[], int nvar,
                     double params[], double *buf,
                     void (*dxdt)(double t, double a[], double params[], double derivs[])){
  int i;                      // counter for loop

  /* allocate memory using buf */
  double *f = buf, *halfV = f+nvar;
  // to store two f at a time for each coordinate
  // half V

  // compute the velocity at t = 1/2
  // z[0] is initial x, z[1] is initial V
  // compute right-hand function at x0, t0
  (*dxdt)(t, a, params, f);

  for(i=0; i<(nvar/2); i++){
    halfV[i] = *(a + 2*i + 1) + 0.5 * f[i] * dt;
  }
    
  // x_n+1 = x_n + h * V_(n+1/2)
  for(i=0; i<(nvar/2); i++){
    *(anew + 2*i) = *(a + 2*i) + halfV[i] * dt;
  }

  // F(x_(n+1))
  (*dxdt)(t + dt, anew, params, f+nvar/2);

  for(i=0; i<(nvar/2); i++){
    // V_n+1 = V_(n+1/2) + 0.5 * h * F(x_(n+1))
    *(anew + 2*i + 1) = halfV[i] + 0.5 * dt * f[nvar/2+i];
  }

  /* for(k=0; k<number; k++){
    f[k] = f[number + k];
  } */
}

/**************************
*
* ODE library functions end 
* 
**************************/