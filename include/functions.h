#include <iostream>
#include <fstream>
#include <cstring>
#include <math.h>
#include <cassert>
#include <iomanip>

#define PI 3.14159265358979323846

using namespace std; 

//Structure that contains the variables required to run the simulation 
typedef struct 
{
    double H;           // height of the storage unit 
    double D;           // diameter
    int    N;           // number of cells 
    double Ti;          // initial temperature 
    double Tf;          // final temperature 

    double t_charge;    // duration of charging state
    double t_discharge; // duration of discharginh state
    double t_idle;      // duration of idle state between charge and discharge
    int    n_cycles;    // number of cycles 
    int    n_timeSteps; // number of tine steps per cycle 
    
    double u_f;         // fluid speed
    double u_s;         // solid speed
    double u_d;         // idle speed

    double T_bcl;       // temperature at the left BC
    double T_bcr;       // temperature at the right BC

    double k_f;         // thermal conductivity of the fluid phase 
    double k_s;         // thermal conductivity of the solid phase 
    double epsilon;     // porosity 
    float  delta_t;     // time step 
    double Cp_f;        // specific heat of the fluid phase at constant pressure
    double C_s;         // specific heat of solid phase  
    double rho_f;       // density of the fluid
    double rho_s;       // density of the solid 
    double h_v;         // volumetric heat transfer coefficient 
    
}variables;


//Get user inputs 
void read_inputs(variables* inputs, int choice);

//Write data
void write_temperature(double **T, int n, int time_step, float delta_t, ofstream &tempFile);

//Write state data
void write_state(int state, int time_step, float delta_t, ofstream &stateFile);

//Write error data
void write_error(float h, double err, double err_avg, float Pe, int n_wave, ofstream &errorFile);

//Charging Equations
void charging_equation(variables* inputs, double alpha_f, double alpha_s, double delta_t, double h, double **T_old, double **T_new);

//Idle Equation 
void idle_equation(variables* inputs, double alpha_f, double alpha_s, double delta_t, double h, double **T_old, double **T_new);

//Discharge Equation 
void discharge_equation(variables* inputs, double alpha_f, double alpha_s, double delta_t, double h, double **T_old, double **T_new);

//LU decomposition 
void luDecomposition(double A[][2], double b[][1], double x[2]);

//Method of Manufactured solutions function 
double MMS(int n, double x, double L, int state);

//Solver    
int solver(variables* inputs);








