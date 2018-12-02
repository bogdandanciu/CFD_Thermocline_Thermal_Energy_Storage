#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

#define PI 3.14159265358979323846

using namespace std; 

//Structure that contains the inputs 
typedef struct variables
{
    double H;           // height
    double D;           // diameter
    int    N;           // number of cells 
    double Ti;          // initial temperature 
    double Tf;          // final temperature 

    double t_charge;    // duration of charging state
    double t_discharge; // duration of discharginh state
    double t_idle;      // duration of idle state between charge and discharge
    int    n_cycles;    // number of cycles 
    int    n_timeSteps; // number of tine steps per cycle 
    
    double u_f;         //fluid speed
    double u_s;         //solid speed
    double u_d;         //idle speed

    double T_bcl;
    double T_bcr;

    double k_f;
    double k_s;
    double epsilon;
    double delta_t;
    double Cp_f;
    double C_s;
    double rho_f;
    double rho_s; 
    double h_v;

}variables;


//Read and user inputs 
void read_inputs(variables* inputs);

//Write state data
void write_data(const int N, double Ti);

//Write state data
void write_state(double t_charge, double t_discharge, double t_idle, const int n_cycles, const int n_timeSteps);

//Charging Equations
void charging_equation(variables* inputs, double alpha_f, double alpha_s, double delta_t, double h, double T_old[][3], double T_new[][3]);

//Idle Equation 
void idle_equation(variables* inputs, double alpha_f, double alpha_s, double delta_t, double h, double T_old[][3], double T_new[][3]);

//Discharge Equation 
void discharge_equation(variables* inputs, double alpha_f, double alpha_s, double delta_t, double h, double T_old[][3], double T_new[][3]);

//LU decomposition 
void luDecomposition(double A[][2], double b[][1], double x[2]);

//Method of Manufacture solutions
double MMS(int n, double x, double L, int state);

//Solver    
void solver(variables* inputs);








