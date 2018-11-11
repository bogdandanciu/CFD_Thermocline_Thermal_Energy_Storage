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

void write2file(const int N, double Ti)
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

