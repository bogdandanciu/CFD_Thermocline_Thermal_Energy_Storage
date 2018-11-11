#include <iostream>
#include <fstream>
#include <string>

using namespace std; 

//Structure that contains the inputs 
typedef struct variables
{
    double H;  //height
    double D;  //diameter
    int    N;  //number of cells 
    double Ti; //initial temperature 
    double Tf; //final temperature 

}variables;


//Read and user inputs 
void read_inputs(variables* inputs);

//Write simulation data
void write2file(const int N, double Ti);
