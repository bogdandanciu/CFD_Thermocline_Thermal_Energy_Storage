#include <iostream>
#include <fstream>
#include <cstring>
#include <math.h>
#include <cassert>
#include <iomanip>

#define PI 3.14159265358979323846

using namespace std; 

//Structure that contains the parameters required to run the simulation 
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
    double u_f;         // fluid speed

    double T_c;         // temperature at charge
    double T_d;         // temperature at the discharge

    double k_f;         // thermal conductivity of the fluid phase 
    double k_s;         // thermal conductivity of the solid phase 
    double epsilon;     // porosity 
    float  delta_t;     // time step 
    double Cp_f;        // specific heat of the fluid phase at constant pressure
    double C_s;         // specific heat of solid phase  
    double rho_f;       // density of the fluid
    double rho_s;       // density of the solid 
	double m_f;         // mass flow 
    double h_v;         // volumetric heat transfer coefficient 
    
	int save_file;       //this number dictates how often output files will be saved 
}parameters;


//Get user inputs 
void read_inputs(parameters* inputs, int choice);

//Write state 
void write_state(int state, int cycle, int time_step, float delta_t, ofstream &stateFile);

//Write data
void write_temperature(double **T, int n, int cycle, int time_step, float delta_t, ofstream &tempFile);

//Write data
void write_data(double exergy_flux, double thermal_energy, double themral_energy_max, int cycle, int time_step, float, ofstream &dataFile);

//Write error 
void write_error(float h, double err, double E1, double E2, float Pe, int n_wave, ofstream &errorFile);

//Charging Equations
void charging_equation(parameters* inputs, double alpha_f, double alpha_s, double delta_t, double h, double **T_old, double **T_new);

//Idle Equation 
void idle_equation(parameters* inputs, double alpha_f, double alpha_s, double delta_t, double h, double **T_old, double **T_new);

//Discharge Equation 
void discharge_equation(parameters* inputs, double alpha_f, double alpha_s, double delta_t, double h, double **T_old, double **T_new);

//LU decomposition 
void luDecomposition(double A[][2], double b[][1], double x[2]);

//Method of Manufactured solutions function 
double MMS_func(int n, double x, double L, int state);

//Source term of MMS
double MMS_source(double x,double alpha_f,double alpha_s,double u, float H, int n, double h_v_f, double h_v_s, int state);

//Main Solver    
int solver(parameters* inputs);








