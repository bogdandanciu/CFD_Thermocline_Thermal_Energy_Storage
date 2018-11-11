#include "functions.h" 

using namespace std; 

//%%%FUNC%%%////////////////////////////////////////////////////////////////////
// Function Name: read_inputs 
//
// Description:   This functions reads the inputs necessary to run the 
//                simulation.
//
////////////////////////////////////////////////////////////////////////////////
//
// Parameters:
//   IN     : <None>
//   OUT    : <None>
//   RETURN : <None>
//
////////////////////////////////////////////////////////////////////////////////
//
// Global Variables used:
//
//%%%FUNC%%%////////////////////////////////////////////////////////////////////

void read_inputs(variables* inputs)
{
    cout << "Please input the diameter of the storage: ";
    cin  >> inputs->D;
    cout << "Please input the height of the storage: ";
    cin  >> inputs->H;
    cout << "Please input the number of cells: ";
    cin  >> inputs->N;
    cout << "Please input the initial temperature: "; 
    cin  >> inputs->Ti; 

    cout << "Please enter the duration of the charge state: ";
    cin  >> inputs->t_charge;
    cout << "Please enter the duration of the discharge state: ";
    cin  >> inputs->t_discharge;
    cout << "Please enter the duration of the idle state: ";
    cin  >> inputs->t_idle;
    cout << "Please enter the number of cycles";
    cin  >> inputs->n_cycles;
    cout << "Please enter the number of time steps";
    cin  >> inputs->n_timeSteps;
}

//%%%FUNC%%%////////////////////////////////////////////////////////////////////
// Function Name: write_data 
//
// Description:   This functions reads the inputs necessary to run the 
//                simulation.
//
////////////////////////////////////////////////////////////////////////////////
//
// Parameters:
//   IN     : <None>
//   OUT    : <None>
//   RETURN : <None>
//
////////////////////////////////////////////////////////////////////////////////
//
// Global Variables used:
//
//%%%FUNC%%%////////////////////////////////////////////////////////////////////

void write_data(const int N, double Ti)
{
    ofstream simFile("simData.dat");
    if (simFile.is_open())
    {
        for (int i = 0; i < N; i++)
        {
            simFile << i << " " << Ti << endl;
        }
        simFile.close();
    }
    else cout << "Unable to open file" << endl; 

}

//%%%FUNC%%%////////////////////////////////////////////////////////////////////
// Function Name: write_state 
//
// Description:   This functions reads the inputs necessary to run the 
//                simulation.
//
////////////////////////////////////////////////////////////////////////////////
//
// Parameters:
//   IN     : <None>
//   OUT    : <None>
//   RETURN : <None>
//
////////////////////////////////////////////////////////////////////////////////
//
// Global Variables used:
//
//%%%FUNC%%%////////////////////////////////////////////////////////////////////

void write_state(double t_charge, double t_discharge, double t_idle, const int n_cycles, const int n_timeSteps)
{
//    ofstrean stateFile("stateData.dat");
//    if (stateFile.is_open())
//    {
//        stateFile 
//    }
}
