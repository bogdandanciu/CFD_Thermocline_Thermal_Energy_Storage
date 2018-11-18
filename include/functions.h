#include <iostream>
#include <fstream>
#include <string>

using namespace std; 

//Structure that contains the inputs 
typedef struct variables
{
    double H;  // height
    double D;  // diameter
    int    N;  // number of cells 
    double Ti; // initial temperature 
    double Tf; // final temperature 

    double t_charge;    // duration of charging state
    double t_discharge; // duration of discharginh state
    double t_idle;      // duration of idle state between charge and discharge
    int    n_cycles;    // number of cycles 
    int    n_timeSteps; // number of tine steps per cycle 




}variables;


//Read and user inputs 
void read_inputs(variables* inputs);

//Write state data
void write_data(const int N, double Ti);

//Write state data
void write_state(double t_charge, double t_discharge, double t_idle, const int n_cycles, const int n_timeSteps);

//Discrete Equations
void charging_equation

void discharging_equation

void idle_equation
