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




//*******Charging equation 
void charging_equation(variables *inputs, 
                       double alpha_f, double alpha_s, 
                       double delta_t, 
                       double h, 
                       double T_old[][3], double T_new[][3])
{
    //Boundary condition on the left 
    T_new[0][0] = h;
    T_new[0][1] = T_old[0][1] + alpha_s*(delta_t/(h*h)) * (T_old[1][1] - T_old[0][1]);
    T_new[0][2] = T_old[0][2] - inputs->u_f*(delta_T/h) * (T_old[0][2] - inputs->T_bcl) 
                              + alpha_f*(delta_t/(h*h)) * (T_old[1][2] - T_old[0][2]); 


    //Main body of computation 
    for (int i = 1; i < inputs->N; i++)
    {
        T_new[i][0] = (i+1)*h;
        T_new[i][1] = T_old[i][1] + alpha_s*(delta_t/(h*h)) * (T_old[i+1][1] - 2*T_old[i][1] + T_old[i-1][1]);
        T_new[i][2] = T_old[i][2] - inputs->u_f*(delta_T/h) * (T_old[i+1][2] - T_old[i-1][2]) 
                                  + alpha_f*(delta_t/(h*h)) * (T_old[i+1][2] - 2*T_old[i][2] + T_old[i-1][2]); 
    }

    //Boundary conditions on the right 
    T_new[inputs->N - 1][0] = inputs->N*h;
    T_new[inputs->N - 1][1] = T_old[inputs->N - 1][1] + alpha_s*(delta_t/(h*h)) * (T_old[inputs->N - 2][1] - T_old[inputs->N - 1][1]);
    T_new[inputs->N - 1][2] = T_old[inputs->N - 1][2] - inputs->u_f*(delta_T/h) * (T_old[inputs->N - 2][2] - T_old[inputs->N - 2][2]) 
                                                      + alpha_f*(delta_t/(h*h)) * (T_old[inputs->N - 2][2] - T_old[inputs->N - 1][2]); 
}


//*******Idle equation 
void idle_equation(variables *inputs, 
                   double alpha_f, double alpha_s, 
                   double delta_t, 
                   double h, 
                   double T_old[][3], double T_new[][3])
{
    //Boundary condition on the left 
    T_new[0][0] = h;
    T_new[0][1] = T_old[0][1] + alpha_s*(delta_t/(h*h)) * (T_old[1][1] - T_old[0][1]);
    T_new[0][2] = T_old[0][2] - alpha_f*(delta_T/(h*h)) * (T_old[1][2] - T_old[0][2]) 


    //Main body of computation 
    for (int i = 1; i < inputs->N; i++)
    {
        T_new[i][0] = (i+1)*h;
        T_new[i][1] = T_old[i][1] + alpha_s*(delta_t/(h*h)) * (T_old[i+1][1] - 2*T_old[i][1] + T_old[i-1][1]);
        T_new[i][2] = T_old[i][2] + alpha_s*(delta_t/(h*h)) * (T_old[i+1][2] - 2*T_old[i][1] + T_old[i-1][2]);
    }

    //Boundary conditions on the right 
    T_new[inputs->N - 1][0] = inputs->N*h;
    T_new[inputs->N - 1][1] = T_old[inputs->N - 1][1] + alpha_s*(delta_t/(h*h)) * (T_old[inputs->N - 2][1] - T_old[inputs->N - 1][1]);
    T_new[inputs->N - 1][2] = T_old[inputs->N - 1][2] + alpha_f*(delta_T/(h*h)) * (T_old[inputs->N - 2][2] - T_old[inputs->N - 2][2]) 
}


//*******Discharge equation 
void discharge_equation(variables *inputs, 
                        double alpha_f, double alpha_s, 
                        double delta_t, 
                        double h, 
                        double T_old[][3], double T_new[][3])
{
    //Fluid velocity at discharge is the negative of the velocity at charge 
    double u_d = -1 * inputs->u_f;

    //Boundary condition on the left 
    T_new[0][0] = h;
    T_new[0][1] = T_old[0][1] + alpha_s*(delta_t/(h*h)) * (T_old[1][1] - T_old[0][1]);
    T_new[0][2] = T_old[0][2] - u_d*(delta_T/h) * (T_old[1][2] - T_old[0][2]) 
                              + alpha_f*(delta_t/(h*h)) * (T_old[1][2] - T_old[0][2]); 


    //Main body of computation 
    for (int i = 1; i < inputs->N; i++)
    {
        T_new[i][0] = (i+1)*h;
        T_new[i][1] = T_old[i][1] + alpha_s*(delta_t/(h*h)) * (T_old[i+1][1] - 2*T_old[i][1] + T_old[i-1][1]);
        T_new[i][2] = T_old[i][2] - u_d*(delta_T/h) * (T_old[i+1][2] - T_old[i][2]) 
                                  + alpha_f*(delta_t/(h*h)) * (T_old[i+1][2] - 2*T_old[i][2] + T_old[i-1][2]); 
    }

    //Boundary conditions on the right 
    T_new[inputs->N - 1][0] = inputs->N*h;
    T_new[inputs->N - 1][1] = T_old[inputs->N - 1][1] + alpha_s*(delta_t/(h*h)) * (T_old[inputs->N - 2][1] - T_old[inputs->N - 1][1]);
    T_new[inputs->N - 1][2] = T_old[inputs->N - 1][2] - inputs->u_d*(delta_T/h) * (T_old[inputs->N - 2][2] - T_old[inputs->N - 2][2]) 
                                                            + alpha_f*(delta_t/(h*h)) * (T_old[inputs->N - 2][2] - T_old[inputs->N - 1][2]); 
}










