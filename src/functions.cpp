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

void read_inputs(variables* inputs, int choice)
{
    
    //if choice == 0 the defaults values are read 
    if (choice == 0)
    {
        inputs->H           = 5.963810;
        inputs->D           = 8.0;
        inputs->t_charge    = 21600.0;
        inputs->t_discharge = 21600.0;
        inputs->t_idle      = 21600.0;
        inputs->Ti          = 293;
        inputs->T_bcl       = 873;
        inputs->T_bcr       = 293;
        inputs->epsilon     = 0.4;
        inputs->u_f         = 0.000271;
        inputs->rho_f       = 1835.6;
        inputs->Cp_f        = 1511.8;
        inputs->k_f         = 0.52;
        inputs->rho_s       = 2600;
        inputs->C_s         = 900;
        inputs->k_s         = 2;
        inputs->h_v         = 448.587788;
    }
    //if choice == 1 the values are manually inputed by the user  
    else if (choice == 1)
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
    
        cout << "Please enter the fluid velocity" << endl;
        cin  >> inputs->u_f;
    }

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
    T_new[0][0] = h; //grid spacing 
    T_new[0][1] = T_old[0][1] + alpha_s*(delta_t/(h*h)) * (T_old[1][1] - T_old[0][1]);
    T_new[0][2] = T_old[0][2] - inputs->u_f*(delta_t/h) * (T_old[0][2] - inputs->T_bcl) 
                              + alpha_f*(delta_t/(h*h)) * (T_old[1][2] - T_old[0][2]); 


    //Main body of computation 
    for (int i = 1; i < inputs->N; i++)
    {
        T_new[i][0] = (i+1)*h;
        T_new[i][1] = T_old[i][1] + alpha_s*(delta_t/(h*h)) * (T_old[i+1][1] - 2*T_old[i][1] + T_old[i-1][1]);
        T_new[i][2] = T_old[i][2] - inputs->u_f*(delta_t/h) * (T_old[i+1][2] - T_old[i-1][2]) 
                                  + alpha_f*(delta_t/(h*h)) * (T_old[i+1][2] - 2*T_old[i][2] + T_old[i-1][2]); 
    }

    //Boundary conditions on the right 
    T_new[inputs->N - 1][0] = inputs->N*h;
    T_new[inputs->N - 1][1] = T_old[inputs->N - 1][1] + alpha_s*(delta_t/(h*h)) * (T_old[inputs->N - 2][1] - T_old[inputs->N - 1][1]);
    T_new[inputs->N - 1][2] = T_old[inputs->N - 1][2] - inputs->u_f*(delta_t/h) * (T_old[inputs->N - 2][2] - T_old[inputs->N - 1][2]) 
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
    T_new[0][2] = T_old[0][2] - alpha_f*(delta_t/(h*h)) * (T_old[1][2] - T_old[0][2]);


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
    T_new[inputs->N - 1][2] = T_old[inputs->N - 1][2] + alpha_f*(delta_t/(h*h)) * (T_old[inputs->N - 2][2] - T_old[inputs->N - 1][2]);
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
    T_new[0][2] = T_old[0][2] - u_d*(delta_t/h) * (T_old[1][2] - T_old[0][2]) 
                              + alpha_f*(delta_t/(h*h)) * (T_old[1][2] - T_old[0][2]); 


    //Main body of computation 
    for (int i = 1; i < inputs->N; i++)
    {
        T_new[i][0] = (i+1)*h;
        T_new[i][1] = T_old[i][1] + alpha_s*(delta_t/(h*h)) * (T_old[i+1][1] - 2*T_old[i][1] + T_old[i-1][1]);
        T_new[i][2] = T_old[i][2] - u_d*(delta_t/h) * (T_old[i+1][2] - T_old[i][2]) 
                                  + alpha_f*(delta_t/(h*h)) * (T_old[i+1][2] - 2*T_old[i][2] + T_old[i-1][2]); 
    }

    //Boundary conditions on the right 
    T_new[inputs->N - 1][0] = inputs->N*h;
    T_new[inputs->N - 1][1] = T_old[inputs->N - 1][1] + alpha_s*(delta_t/(h*h)) * (T_old[inputs->N - 2][1] - T_old[inputs->N - 1][1]);
    T_new[inputs->N - 1][2] = T_old[inputs->N - 1][2] - inputs->u_d*(delta_t/h) * (T_old[inputs->N - 2][2] - T_old[inputs->N - 1][2]) 
                                                      + alpha_f*(delta_t/(h*h)) * (T_old[inputs->N - 2][2] - T_old[inputs->N - 1][2]); 
}



//*******LU decomposition 
void luDecomposition(double A[][2], double b[][1], double x[2])
{
    double lower[2][2], upper[2][2];
    double z1, z2, x1, x2;

//    memset(lower, 0, size(lower));
//    memset(upper, 0, size(upper));
    
    //Decomposition into Upper and Lower matrices 
    for (int i = 0; i < 2; i++)
    {
        //Upper Triangular Matrix 
        for (int k = i; k < 2; k++)
        {
            //Sum of L(i,j) * U(j,k)
            double sum = 0;
            for (int j = 0; j < i; j++) 
            {
                sum += (lower[i][j] * upper[j][k]);

            }
            //Evaluating upper matrix
            upper[i][k] = A[i][k] - sum;
        } 
        //Lower Triangular Matrix 
        for (int k = i; k < 2; k++)
        {
            if (i == k)
            //Diagonal is equal to 1 
                lower[i][i] = 1;   
            else
            {
                //Sum of L(k,j) * U(j,i)
                double sum = 0;
                for (int j = 1; j < i; j++)
                {
                    sum += (lower[k][j] * upper[j][i]);
                }       
                //Evaluating L(k,i)
                lower[k][i] = (A[k][i] - sum)/upper[i][i];
            }

        }
    }

    //Back substitution
    z1 = b[0][0];
    z2 = (b[1][0] - lower[1][0]*z1);
    x2 = z2/upper[1][1];
    x1 = (z1 - upper[0][1] * x2)/upper[0][0];
    x[0] = x1;
    x[1] = x2;

}


//******Order of verification study 
//
//*****Method of Manufactured solutions 
double MMS(int n, double x, double L, int state)
{
    double k = (2*PI*n)/L;

    if (state == 0) //fluid 
    {
        double T_sol_f = cos(k*x);

        return T_sol_f;
    }else if (state == 1) //solid
    {
        double T_sol_s = sin(k*x);

        return T_sol_s;
    }
}


//*********Solver
void solver(variables *inputs)
{
//    int         state;
//    double      h       = inputs->H/inputs->N;
//    const float delta_t = inputs->delta_t;
//    double      t_total = inputs->t_charge + inputs->t_discharge + 2*inputs->t_idle;
//    double      alpha_f = inputs->k_f / (inputs->epsilon * inputs->rho_f * inputs->Cp_f);
//    double      alpha_s = inputs->k_s / ((1-inputs->epsilon) * inputs->rho_s * inputs->C_s);
//    double      h_v_f   = inputs->h_v / (inputs->epsilon * inputs->rho_f * inputs->Cp_f);
//    double      h_v_s   = inputs->h_v / ((1-inputs->epsilon) * inputs->rho_s * inputs->C_s);



}

