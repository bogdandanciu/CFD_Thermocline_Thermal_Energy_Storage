#include "functions.h" 


//%%%FUNC%%%////////////////////////////////////////////////////////////////////
// Function Name: read_inputs 
//
// Description:   This function reads the inputs necessary to run the 
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
        inputs->t_charge    = 500.0;
        inputs->t_discharge = 500.0;
        inputs->t_idle      = 500.0;
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
// Description:   This function reads the inputs necessary to run the 
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
    else cout << "Unable to open file\n"; 

}


//%%%FUNC%%%////////////////////////////////////////////////////////////////////
// Function Name: write_error 
//
// Description:   This function reads the inputs necessary to run the 
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

void write_error(double err, double err_avg, float h, float Pe, ofstream &errorFile)
{
//    ofstream errorFile("errorData.dat");
    if (errorFile.is_open())
    {
        errorFile << h << " " << err << " " << err_avg << " " << Pe << endl;
    }
    else cout << "UNABLE TO OPEN FILE\n"; 

}

//%%%FUNC%%%////////////////////////////////////////////////////////////////////
// Function Name: write_state 
//
// Description:   This function reads the inputs necessary to run the 
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

void write_state(int state, int time_step, float delta_t, ofstream &stateFile)
{
//    ofstream stateFile("State_data.dat");
    if (stateFile.is_open())
    {
        stateFile << time_step*delta_t << " " << state << endl;
    }
    else cout << "UNABLE TO WRITE DATA IN FILE \n";
}




//*******Charging equation 
void charging_equation(variables *inputs, 
                       double alpha_f, double alpha_s, 
                       double delta_t, 
                       double h, 
                       double **T_old, double **T_new)
{
    //Boundary condition on the left 
    T_new[0][0] = h; //grid spacing 
    T_new[0][1] = T_old[0][1] + alpha_s*(delta_t/(h*h)) * (T_old[1][1] - T_old[0][1]);
    T_new[0][2] = T_old[0][2] - inputs->u_f*(delta_t/h) * (T_old[0][2] - inputs->T_bcl) 
                              + alpha_f*(delta_t/(h*h)) * (T_old[1][2] - T_old[0][2]); 


    //Main body of computation 
    for (int i = 1; i < inputs->N-1; i++)
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
                   double **T_old, double **T_new)
{
    //Boundary condition on the left 
    T_new[0][0] = h;
    T_new[0][1] = T_old[0][1] + alpha_s*(delta_t/(h*h)) * (T_old[1][1] - T_old[0][1]);
    T_new[0][2] = T_old[0][2] - alpha_f*(delta_t/(h*h)) * (T_old[1][2] - T_old[0][2]);


    //Main body of computation 
    for (int i = 1; i < inputs->N-1; i++)
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
                        double **T_old, double **T_new)
{
    //Fluid velocity at discharge is the negative of the velocity at charge 
    double u_d = -1 * inputs->u_f;

    //Boundary condition on the left 
    T_new[0][0] = h;
    T_new[0][1] = T_old[0][1] + alpha_s*(delta_t/(h*h)) * (T_old[1][1] - T_old[0][1]);
    T_new[0][2] = T_old[0][2] - u_d*(delta_t/h) * (T_old[1][2] - T_old[0][2]) 
                              + alpha_f*(delta_t/(h*h)) * (T_old[1][2] - T_old[0][2]); 


    //Main body of computation 
    for (int i = 1; i < inputs->N-1; i++)
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
    }
    else if (state == 1) //solid
    {
        double T_sol_s = sin(k*x);

        return T_sol_s;
    }

    return 0; 
}



//*************************************************


//*********Solver
int solver(variables *inputs)
{

    cout << "Plese input the number of cells: ";
    cin >> inputs->N;

    cout << "Please input the number of cycles: ";
    cin >> inputs->n_cycles;
    
    cout <<"Please input the time step delta_t: "; 
    cin >> inputs->delta_t;

    //Check for stability of inputs
    



    int         state;
    double      h       = inputs->H/inputs->N; //grid spacing 

    const float delta_t = inputs->delta_t;
    double      t_total = inputs->t_charge + inputs->t_discharge + 2*inputs->t_idle;
    double      alpha_f = inputs->k_f / (inputs->epsilon * inputs->rho_f * inputs->Cp_f);
    double      alpha_s = inputs->k_s / ((1-inputs->epsilon) * inputs->rho_s * inputs->C_s);
    double      h_v_f   = inputs->h_v / (inputs->epsilon * inputs->rho_f * inputs->Cp_f);
    double      h_v_s   = inputs->h_v / ((1-inputs->epsilon) * inputs->rho_s * inputs->C_s);
    double      Pe      = (inputs->u_f * inputs->H)/alpha_f;

    int         time_step;
    int         cycle = 0;


    int save_file = 1; 

    ofstream stateFile;
    stateFile.open("State_data.dat");

    ofstream errorFile; 
    errorFile.open("Error_data.dat");

    //Initialize the temperature domain 
    //Allocate memory 
//    double (*arr)[3];
//    arr = calloc(inputs->N, sizeof(*arr));
//    double (*T_old)[3] = arr; 
//    if (T_old == NULL)
//    {
//        cout << "Memory allocation failed! Aborting simulation!" << endl;
//        exit(EXIT_FAILURE);
//    }
    //Finish memory allocation 
//    T_old = new (nothrow) double[inputs->N][3];
//    assert(T_old != NULL);

    //Initialize the temperature domain 
//    double T_old[inputs->N][3];

    //Initialize the temperature domain 
    double** T_old = new double*[inputs->N];
    for (int i = 0; i < inputs->N; i++)
        T_old[i] = new double[3];

    for (int i = 0; i < inputs->N; i++)
    {
        T_old[i][0] = 1;
        T_old[i][1] = inputs->Ti;
        T_old[i][2] = inputs->Ti;
    }

    double max_error = 1.0;

    cout << "Start of the simulation\n"; 

    while (cycle < inputs->n_cycles)  //Main computation body  
    {
        time_step = 0;
        for (double simulation_time = 0; simulation_time <= t_total; simulation_time += delta_t)
        {
            //Intialize new temperature 
            double** T_new = new double*[inputs->N];
            for (int i = 0; i < inputs->N; i++)
                T_new[i]  = new double[3];

            //Charging phase 
            if (simulation_time <= inputs->t_charge) 
            {
                state = 1;
                if (time_step%save_file == 0)
                {
                    cout << "THERMOCLINE IS CHARGING!\n";
                    write_state(state, time_step, delta_t, stateFile);
                }
                charging_equation(inputs, alpha_f, alpha_s, delta_t, h, T_old, T_new);  
            }
            //Idle Phase
            else if (simulation_time > inputs->t_charge && simulation_time <= (inputs->t_charge + inputs->t_idle)) 
            {
                state = 0; 
                if (time_step%save_file == 0)
                {
                    cout << "THERMOCLINE IS IDLING! \n"; 
                    write_state(state, time_step, delta_t, stateFile);
                }   
                idle_equation(inputs, alpha_f, alpha_s, delta_t, h, T_old, T_new);
            }
            //Discharging phase
            else if (simulation_time > (inputs->t_charge + inputs->t_idle) && simulation_time <= (inputs->t_charge + inputs->t_idle + inputs->t_discharge))
            {
                state = 2;
                if (time_step%save_file == 0)
                {
                    cout << "THERMOCLINE IS DISCHARGING \n"; 
                    write_state(state, time_step, delta_t, stateFile);
                }
                discharge_equation(inputs, alpha_f, alpha_s, delta_t, h, T_old, T_new);
            }


            //Solve the linear system 
//            double A[2][2] = {{ 1 + (h_v_f*delta_t), (-1)*h_v_f*delta_t}, { (-1)*h_v_s*delta_t, 1 + h_v_s*delta_t}};
//            for (int i = 0; i < inputs->N; i++)
//            {
//                double x[2] = {};
//                double b[2][1] =  {{T_new[i][2]}, {T_new[i][1]}};
//                luDecomposition(A,b,x);
//                T_new[i][2] = x[0];
//                T_new[i][1] = x[1];
//            }


            //Error per iteration 
            double error_i;
            max_error = 0.0000;
            for (int i = 1; i < inputs->N; i++)
            {
                error_i = abs(T_new[i][2] - T_old[i][2]);
                if (error_i > max_error)
                    max_error = error_i;
            }


            //Convergence Check 
            double residual = 0;
            for (int i = 0; i < inputs->N; i++)
            {
                residual = residual + abs(T_new[i][2] - T_old[i][2]);
            }

            double error = residual/inputs->N;

            
            //OVS error
            if (error <= 1e-6)
            {
                double err = 0;
                double err_avg = 0;

                //Checking convergence for fluid 
                int fs_state = 0; //0 for fluid and 1 for solid
                
                for (int i = 1; i < inputs->N; i++)
                {
                    err = err + abs(T_new[i][2] - MMS(inputs->N, i*h, inputs->H, fs_state)) / (2 + MMS(inputs->N, i*h, inputs->H, fs_state));
                }

                err_avg = err/inputs->N;
                write_error(err, err_avg, h, Pe, errorFile);
            } 
            //End of Convergene check 

            
            //**Copy T_new to T_old
            //Free Memory for T_old 
            for (int i = 0; i < inputs->N; i++)
                delete [] T_old[i];
            delete [] T_old;
            //Initialize empty T_old
            double** T_old = new double*[inputs->N];
            for (int i = 0; i < inputs->N; i++)
                T_old[i] = new double[3];
            //Copy from T_new to T_old 
            memcpy(T_old,T_new, (inputs->N) * sizeof(**T_new));
            //Free T_new 
            for (int i = 0; i < inputs->N; i++)
                delete [] T_new[i];
            delete [] T_new;


            time_step++;

        }

        cycle++;
    }
    
//    cout << "Out of the main loop\n";

    //Free Memory for T_old 
//    for (int i = 0; i < inputs->N; i++)
//        delete [] T_old[i];
//    delete [] T_old;

    cout << "The solution convervged after " << cycle << " cycles"; 

    return 0;
}
