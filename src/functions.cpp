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
//   IN     : inputs, choice 
//   OUT    : inputs  
//   RETURN : <None>
//
////////////////////////////////////////////////////////////////////////////////
//
// Global Variables used: parameters
//
//%%%FUNC%%%////////////////////////////////////////////////////////////////////

void read_inputs(parameters* inputs, int choice)
{
    
    //if choice == 0 the defaults values are read 
    if (choice == 0)
    {
        inputs->H           = 23.8732;
        inputs->D           = 8;
        inputs->t_charge    = 21600.0;
        inputs->t_discharge = 21600.0;
        inputs->t_idle      = 21600.0;
        inputs->Ti          = 293;
        inputs->T_c         = 873;
        inputs->T_d         = 293;
        inputs->epsilon     = 0.4;
        inputs->u_f         = 0.0011;
        inputs->rho_f       = 1835.6;
		inputs->rho_s       = 2600;
        inputs->k_f         = 0.52;
        inputs->k_s         = 2;
		inputs->Cp_f        = 1511.8;
		inputs->C_s         = 900;
		inputs->m_f         = 10;
        inputs->h_v         = 1131.797;
		inputs->save_file   = 21600;
    }
    //if choice == 1 the values are manually inputed by the user  
    else if (choice == 1)
    {   
        cout << "Please input the height of the storage: ";
        cin  >> inputs->H;
        cout << "Please input the diameter of the storage: ";
        cin  >> inputs->D;
	    cout << "Please enter the duration of the charge state: ";
        cin  >> inputs->t_charge;
        cout << "Please enter the duration of the discharge state: ";
        cin  >> inputs->t_discharge;
        cout << "Please enter the duration of the idle state: ";
        cin  >> inputs->t_idle;
        cout << "Please input the initial temperature: "; 
        cin  >> inputs->Ti; 
	    cout << "Please input the charge temperature: "; 
        cin  >> inputs->T_c;
	    cout << "Please input the discharge temperature: "; 
        cin  >> inputs->T_d; 		
		cout << "Please enter the fluid velocity:" << endl;
        cin  >> inputs->epsilon;
        cout << "Please enter the fluid velocity:" << endl;
        cin  >> inputs->u_f;
	    cout << "Please enter the density of the fluid:" << endl;
        cin  >> inputs->rho_f;
	    cout << "Please enter the density of the solid:" << endl;
        cin  >> inputs->rho_s;
		cout << "Please enter the thermal conductivity of the fluid phase(k_f):" << endl;
        cin  >> inputs->k_f;
		cout << "Please enter the thermal conductivity of the solid phase(k_s):" << endl;
        cin  >> inputs->k_s;
		cout << "Please enter the specific heat of the fluid phase at constant pressure(Cp_f):" << endl;
        cin  >> inputs->Cp_f;
		cout << "Please enter the specific heat of solid phase(C_s):" << endl;
        cin  >> inputs->C_s;
		cout << "Please enter the mass flow (m_f):" << endl;
        cin  >> inputs->m_f;
		cout << "Please enter the volumetric heat transfer coefficient (h_v):" << endl;
        cin  >> inputs->h_v;
		cout << "Please enter the number of times you want to write data (the file will be written every xth time):" << endl;
        cin  >> inputs->save_file;
		
    }

}

//%%%FUNC%%%////////////////////////////////////////////////////////////////////
// Function Name: write_state 
//
// Description:   This function writes to a file the state at the current time.
//
////////////////////////////////////////////////////////////////////////////////
//
// Parameters:
//   IN     : state, cycle, time_step, delta_t, stateFile 
//   OUT    : stateFile
//   RETURN : <None>
//
////////////////////////////////////////////////////////////////////////////////
//
// Global Variables used: <None>
//
//%%%FUNC%%%////////////////////////////////////////////////////////////////////

void write_state(int state, int cycle, int time_step, float delta_t, ofstream &stateFile)
{
//    ofstream stateFile("State_data.dat");
    if (stateFile.is_open())
    {
        stateFile << (cycle-1)*86400 + time_step*delta_t << " " << state << endl;
    }
    else cout << "UNABLE TO WRITE DATA IN FILE \n";
}

//%%%FUNC%%%////////////////////////////////////////////////////////////////////
// Function Name: write_temperature 
//
// Description:   This function writes to a file the fluid and solid 
//                temepratures at the current time.
//
////////////////////////////////////////////////////////////////////////////////
//
// Parameters:
//   IN     : T, n, cycle, time_step, delta_t, tempFile
//   OUT    : tempFile
//   RETURN : <None>
//
////////////////////////////////////////////////////////////////////////////////
//
// Global Variables used: <None>
//
//%%%FUNC%%%////////////////////////////////////////////////////////////////////

void write_temperature(double **T, int n, int cycle, int time_step, float delta_t, ofstream &tempFile)
{
    if (tempFile.is_open())
    {
        tempFile << "Temperature profile at time: " << (cycle-1)*86400+time_step*delta_t << endl;
        for (int i = 0; i < n; i++)
        {    
            tempFile << fixed << setprecision(6) <<  T[i][0] << " " << T[i][1] 
                     << " " << T[i][2]  << endl;
        }
        tempFile << endl;
    }
    else cout << "Unable to open file\n"; 

}

//%%%FUNC%%%////////////////////////////////////////////////////////////////////
// Function Name: write_error 
//
// Description:   This function writes to a file the grid spacing the errors, 
//				  the Peclet number and the wave number.
//
////////////////////////////////////////////////////////////////////////////////
//
// Parameters:
//   IN     : h, err, E1, E2, Pe, n_wave
//   OUT    : errorFile
//   RETURN : <None>
//
////////////////////////////////////////////////////////////////////////////////
//
// Global Variables used: <None>
//
//%%%FUNC%%%////////////////////////////////////////////////////////////////////

void write_error(float h, double err, double E1, double E2, float Pe, int n_wave, ofstream &errorFile)
{
    if (errorFile.is_open())
    {
        errorFile << fixed << setprecision (4) << h << " " << err << " " 
                  << E1 << " " << E2 << " " << Pe << " " << n_wave << endl;
    }
    else cout << "UNABLE TO OPEN FILE\n"; 
}

//%%%FUNC%%%////////////////////////////////////////////////////////////////////
// Function Name: write_data 
//
// Description:   This function writes to a file the exergy flux, the
//                themral_energy and the themral_energy_max at the current time.
//
////////////////////////////////////////////////////////////////////////////////
//
// Parameters:
//   IN     : exergy, thermal_energy, thermal_energy_max, time_step, delta_t
//   OUT    : dataFile
//   RETURN : <None>
//
////////////////////////////////////////////////////////////////////////////////
//
// Global Variables used: <None>
//
//%%%FUNC%%%////////////////////////////////////////////////////////////////////
void write_data(double exergy, double thermal_energy, double themral_energy_max, int cycle,  int time_step, float delta_t, ofstream &dataFile)
{
//	dataFile << "Exergy flux " << "Thermal energy " << "Maximum Thermal Energy " << "Time" << endl;
    if (dataFile.is_open())
    {
        dataFile << exergy << " " << thermal_energy << " " << (cycle-1)*86400 + time_step*delta_t << endl;
    }
    else cout << "UNABLE TO WRITE DATA IN FILE \n";
}

//******Order of verification study 
//%%%FUNC%%%////////////////////////////////////////////////////////////////////
// Function Name: MMS_func
//
// Description:   This function returns the function value for the MMS function
//
////////////////////////////////////////////////////////////////////////////////
//
// Parameters:
//   IN     : n, x, L, state
//   OUT    : <None>
//   RETURN : T_sol_f
//
////////////////////////////////////////////////////////////////////////////////
//
// Global Variables used: <None>
//
//%%%FUNC%%%//////////////////////////////////////////////////////////////////// 
double MMS_func(int n, double x, double L, int state)
{
    double k = (2*PI*n)/L;
    double T_sol_f = cos(k*x);

    return T_sol_f;
}

//%%%FUNC%%%////////////////////////////////////////////////////////////////////
// Function Name: MMS_source
//
// Description:   This function provides the source term for the MMS function.
//
////////////////////////////////////////////////////////////////////////////////
//
// Parameters:
//   IN     : x, alpha_f, alpha_s, u, H, n, h_v_f, h_v_s, state
//   OUT    : <None>
//   RETURN : f_src or s_src
//
////////////////////////////////////////////////////////////////////////////////
//
// Global Variables used: <None>
//
//%%%FUNC%%%////////////////////////////////////////////////////////////////////
double MMS_source(double x,double alpha_f,double alpha_s,double u, float H, int n, double h_v_f, double h_v_s, int state)
{
    double f_src, s_src;
    if (state == 1) // fluid 
	{
        f_src = (( 4*alpha_f*n*n*PI*PI*cos((2*PI*n*x)/H))/(H*H) ) - ((2*n*u*PI*sin((2*PI*n*x)/H))/H) + h_v_f*(MMS_func(n,x,H,1) - MMS_func(n,x,H,0));
		
        return f_src;
    }
    else if (state == 0) // solid
	{
        s_src = (( 4*alpha_s*n*n*PI*PI*cos((2*PI*n*x)/H))/(H*H)) + h_v_s*(MMS_func(n,x,H,0) - MMS_func(n,x,H,1));
		
        return s_src;
    }
	
    return 0;
}

//%%%FUNC%%%////////////////////////////////////////////////////////////////////
// Function Name: charging_equation
//
// Description:   This function solves the discretized equation for the charging
//                phase. Also contains the source term for the OVS. 
//
////////////////////////////////////////////////////////////////////////////////
//
// Parameters:
//   IN     : inputs, alpha_f, alpha_s, delta_t, h, T_old, T_new
//   OUT    : <None>
//   RETURN : <None>
//
////////////////////////////////////////////////////////////////////////////////
//
// Global Variables used: parameters
//
//%%%FUNC%%%//////////////////////////////////////////////////////////////////// 
void charging_equation(parameters *inputs, 
                       double alpha_f, double alpha_s, 
                       double delta_t, 
                       double h,
                       double **T_old, double **T_new)
{
    //Boundary condition on the left 
    T_new[0][0] = h; //grid spacing 
    T_new[0][1] = T_old[0][1] + alpha_s*(delta_t/(h*h)) * (T_old[1][1] - T_old[0][1]);// + delta_t*MMS_source(h,alpha_f,alpha_s,inputs->u_f,inputs->H,n_wave,h_v_f,h_v_s,0) ;
    T_new[0][2] = T_old[0][2] - inputs->u_f*(delta_t/h) * (T_old[0][2] - inputs->T_c) 
                              + alpha_f*(delta_t/(h*h)) * (T_old[1][2] - T_old[0][2]);// + delta_t*MMS_source(h,alpha_f,alpha_s,inputs->u_f,inputs->H,n_wave,h_v_f,h_v_s,1);


    //Main body of computation 
    for (int i = 1; i < inputs->N-1; i++)
    {
        double l = (i+1)*h;
        T_new[i][0] = l;
        T_new[i][1] = T_old[i][1] + alpha_s*(delta_t/(h*h)) * (T_old[i+1][1] - 2*T_old[i][1] + T_old[i-1][1]);// + delta_t*MMS_source(h*(i+1),alpha_f,alpha_s,inputs->u_f,inputs->H,n_wave,h_v_f,h_v_s,0);
        T_new[i][2] = T_old[i][2] - inputs->u_f*(delta_t/h) * (T_old[i][2] - T_old[i-1][2]) 
                                  + alpha_f*(delta_t/(h*h)) * (T_old[i+1][2] - 2*T_old[i][2] + T_old[i-1][2]);// + delta_t*MMS_source(h*(i+1),alpha_f,alpha_s,inputs->u_f,inputs->H,n_wave,h_v_f,h_v_s,1);
    }

    //Boundary conditions on the right 
    T_new[inputs->N - 1][0] = inputs->N*h;
    T_new[inputs->N - 1][1] = T_old[inputs->N - 1][1] + alpha_s*(delta_t/(h*h)) * (T_old[inputs->N - 2][1] - T_old[inputs->N - 1][1]);// + delta_t*MMS_source((inputs->N)*h , alpha_f,alpha_s, inputs->u_f, inputs->H,n_wave,h_v_f,h_v_s,0);
    T_new[inputs->N - 1][2] = T_old[inputs->N - 1][2] - inputs->u_f*(delta_t/h) * (T_old[inputs->N - 1][2] - T_old[inputs->N - 2][2]) 
                                                      + alpha_f*(delta_t/(h*h)) * (T_old[inputs->N - 2][2] - T_old[inputs->N - 1][2]);// + delta_t*MMS_source((inputs->N)*h ,alpha_f,alpha_s,inputs->u_f,inputs->H,n_wave,h_v_f,h_v_s,1);
}

//%%%FUNC%%%////////////////////////////////////////////////////////////////////
// Function Name: idle_equation
//
// Description:   This function solves the discretized equation for the idle
//                phase. 
//
////////////////////////////////////////////////////////////////////////////////
//
// Parameters:
//   IN     : inputs, alpha_f, alpha_s, delta_t, h, T_old, T_new
//   OUT    : <None>
//   RETURN : <None>
//
////////////////////////////////////////////////////////////////////////////////
//
// Global Variables used: parameters
//
//%%%FUNC%%%//////////////////////////////////////////////////////////////////// 
void idle_equation(parameters *inputs, 
                   double alpha_f, double alpha_s, 
                   double delta_t, 
                   double h, 
                   double **T_old, double **T_new)
{
    //Boundary condition on the left 
    T_new[0][0] = h;
    T_new[0][1] = T_old[0][1] + alpha_s*(delta_t/(h*h)) * (T_old[1][1] - T_old[0][1]);
    T_new[0][2] = T_old[0][2] + alpha_f*(delta_t/(h*h)) * (T_old[1][2] - T_old[0][2]);


    //Main body of computation 
    for (int i = 1; i < inputs->N-1; i++)
    {
        double l = (i+1)*h;
        T_new[i][0] = l;
        T_new[i][1] = T_old[i][1] + alpha_s*(delta_t/(h*h)) * (T_old[i+1][1] - 2*T_old[i][1] + T_old[i-1][1]);
        T_new[i][2] = T_old[i][2] + alpha_s*(delta_t/(h*h)) * (T_old[i+1][2] - 2*T_old[i][1] + T_old[i-1][2]);
    }

    //Boundary conditions on the right 
    T_new[inputs->N - 1][0] = (inputs->N)*h;
    T_new[inputs->N - 1][1] = T_old[inputs->N - 1][1] + alpha_s*(delta_t/(h*h)) * (T_old[inputs->N - 2][1] - T_old[inputs->N - 1][1]);
    T_new[inputs->N - 1][2] = T_old[inputs->N - 1][2] + alpha_f*(delta_t/(h*h)) * (T_old[inputs->N - 2][2] - T_old[inputs->N - 1][2]);
}


//%%%FUNC%%%////////////////////////////////////////////////////////////////////
// Function Name: discharge_equation
//
// Description:   This function solves the discretized equation for the discharge
//                phase. 
//
////////////////////////////////////////////////////////////////////////////////
//
// Parameters:
//   IN     : inputs, alpha_f, alpha_s, delta_t, h, T_old, T_new
//   OUT    : <None>
//   RETURN : <None>
//
////////////////////////////////////////////////////////////////////////////////
//
// Global Variables used: parameters
//
//%%%FUNC%%%//////////////////////////////////////////////////////////////////// 
void discharge_equation(parameters *inputs, 
                        double alpha_f, double alpha_s, 
                        double delta_t, 
                        double h, 
                        double **T_old, double **T_new)
{
    //Fluid velocity at discharge is the negative of the velocity at charge 
    double u_d = (-1) * inputs->u_f;

    //Boundary condition on the left 
    T_new[0][0] = h;
    T_new[0][1] = T_old[0][1] + alpha_s*(delta_t/(h*h)) * (T_old[1][1] - T_old[0][1]);
    T_new[0][2] = T_old[0][2] - u_d*(delta_t/h) * (T_old[1][2] - T_old[0][2]) 
                              + alpha_f*(delta_t/(h*h)) * (T_old[1][2] - T_old[0][2]); 


    //Main body of computation 
    for (int i = 1; i < inputs->N-1; i++)
    {
        double l = (i+1)*h;
        T_new[i][0] = l;
        T_new[i][1] = T_old[i][1] + alpha_s*(delta_t/(h*h)) * (T_old[i+1][1] - 2*T_old[i][1] + T_old[i-1][1]);
        T_new[i][2] = T_old[i][2] - u_d*(delta_t/h) * (T_old[i+1][2] - T_old[i][2]) 
                                  + alpha_f*(delta_t/(h*h)) * (T_old[i+1][2] - 2*T_old[i][2] + T_old[i-1][2]); 
    }

    //Boundary conditions on the right 
    T_new[inputs->N - 1][0] = inputs->N*h;
    T_new[inputs->N - 1][1] = T_old[inputs->N - 1][1] + alpha_s*(delta_t/(h*h)) * (T_old[inputs->N - 2][1] - T_old[inputs->N - 1][1]);
    T_new[inputs->N - 1][2] = T_old[inputs->N - 1][2] - u_d*(delta_t/h) * (inputs->T_d - T_old[inputs->N - 1][2]) 
                                                      + alpha_f*(delta_t/(h*h)) * (T_old[inputs->N - 2][2] - T_old[inputs->N - 1][2]); 
}

//%%%FUNC%%%////////////////////////////////////////////////////////////////////
// Function Name: luDecomposition
//
// Description:   This provides a linear solver for a 2x2 system of equations 
//			      using the LU decompostion and back substitution. 
//
////////////////////////////////////////////////////////////////////////////////
//
// Parameters:
//   IN     : A[][2], b[][1], x[2] 
//   OUT    : <None>
//   RETURN : <None>
//
////////////////////////////////////////////////////////////////////////////////
//
// Global Variables used: <None>
//
//%%%FUNC%%%//////////////////////////////////////////////////////////////////// 
void luDecomposition(double A[][2], double b[][1], double x[2])
{
    double lower[2][2], upper[2][2];
    double z1, z2, x1, x2;
    
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

//%%%FUNC%%%////////////////////////////////////////////////////////////////////
// Function Name: solver
//
// Description:   This function represents the main solver used to compute the 
//                discretized two phase flow problem. It also outputs the results 
//				  (temperature, state, exergy_flux, themral_energy) to .dat 
// 				  files. 
//
////////////////////////////////////////////////////////////////////////////////
//
// Parameters:
//   IN     : inputs 
//   OUT    : <None>
//   RETURN : <None>
//
////////////////////////////////////////////////////////////////////////////////
//
// Global Variables used: parameters
//
//%%%FUNC%%%//////////////////////////////////////////////////////////////////// 
int solver(parameters *inputs)
{

    cout << "Plese input the number of cells: ";
    cin >> inputs->N;
    
    cout <<"Please input the time step delta_t: "; 
    cin >> inputs->delta_t;

    cout << "Please input the number of cycles: ";
    cin >> inputs->n_cycles;
	
    int         state; //1 for charging, 2 for discharging, 3 for idling,
    double      h       = inputs->H/inputs->N; //grid spacing 
    
//    double alpha_f = 2e-7;
//    double alpha_s = 9e-7;

    const float delta_t = inputs->delta_t;
    double      t_total = inputs->t_charge + inputs->t_discharge + 2*inputs->t_idle;
    double      alpha_f = inputs->k_f / (inputs->epsilon * inputs->rho_f * inputs->Cp_f);
    double      alpha_s = inputs->k_s / ((1-inputs->epsilon) * inputs->rho_s * inputs->C_s);
    double      h_v_f   = inputs->h_v / (inputs->epsilon * inputs->rho_f * inputs->Cp_f);
    double      h_v_s   = inputs->h_v / ((1-inputs->epsilon) * inputs->rho_s * inputs->C_s);
    double      Pe      = (inputs->u_f * inputs->H)/alpha_f;
    
    //Initialize output files
    ofstream stateFile;
    stateFile.open("State_data.dat");

//    ofstream errorFile; 
//    errorFile.open("Error_data.dat");

    ofstream tempFile;
    tempFile.open("Temperature_data.dat");

    ofstream dataFile; 
    dataFile.open("Exergy_ThermalEn.dat");
	dataFile << "Exergy flux; " << "Thermal energy; " << "Maximum Thermal Energy; " << "Time;" << endl;

//    ofstream exergyFile;
//    exergyFile.open("Exergy_data.dat");

    //Initialize the temperature domain 
    double** T_old = new double*[inputs->N];
    for (int i = 0; i < inputs->N; i++)
        T_old[i] = new double[3];

    for (int i = 0; i < inputs->N; i++)
    {
        double l = (i+1)*h;
        T_old[i][0] = l;
        T_old[i][1] = inputs->Ti;
        T_old[i][2] = inputs->Ti;
    }

    //Write first initial temperature
//    int time_step_init = 0; 
//	  int cycle_init     = 1;
//    write_temperature(T_old, inputs->N, cycle_init, time_step_init, delta_t, tempFile);
//    write_state(1, cycle_init,time_step_init, delta_t, stateFile);

    double max_error = 1.0;

	//Max capcity factor
	double thermal_energy_max = (inputs->epsilon*inputs->rho_f*inputs->Cp_f + (1-inputs->epsilon)*inputs->rho_s*inputs->C_s)*
								    PI/4*(inputs->D*inputs->D)*inputs->H*(inputs->T_c-inputs->T_d);
	

	int time_step;
    int cycle = 1;

	//Start of the main solver
    cout << "Start of the simulation\n"; 
    while (cycle <= inputs->n_cycles)  //Main computation body  
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
                if (time_step%inputs->save_file == 0)
                {
                    cout <<"Cycle " << cycle << ": THERMOCLINE IS CHARGING!\n";
                    write_state(state, cycle, time_step, delta_t, stateFile);
                }
                charging_equation(inputs, alpha_f, alpha_s, delta_t, h, T_old, T_new);  
            }
            //Idle Phase
            else if (simulation_time > inputs->t_charge && simulation_time <= (inputs->t_charge + inputs->t_idle)) 
            {
                state = 3; 
                if (time_step%inputs->save_file == 0)
                {
                    cout << "Cycle " << cycle << ": THERMOCLINE IS IDLING! \n"; 
                    write_state(state, cycle, time_step, delta_t, stateFile);
                }   
                idle_equation(inputs, alpha_f, alpha_s, delta_t, h, T_old, T_new);
            }
            //Discharging phase
            else if (simulation_time > (inputs->t_charge + inputs->t_idle) && simulation_time <= (inputs->t_charge + inputs->t_idle + inputs->t_discharge))
            {
                state = 2;
                if (time_step%inputs->save_file == 0)
                {
                    cout << "Cycle " << cycle  << ": THERMOCLINE IS DISCHARGING! \n"; 
                    write_state(state, cycle, time_step, delta_t, stateFile);
                }
                discharge_equation(inputs, alpha_f, alpha_s, delta_t, h, T_old, T_new);
            }
            //Idle Phase
            else if (simulation_time > (inputs->t_charge + inputs->t_idle + inputs->t_discharge)  && simulation_time <= t_total) 
            {
                state = 3; 
                if (time_step%inputs->save_file == 0)
                {
                    cout << "Cycle " << cycle << ": THERMOCLINE IS IDLING! \n"; 
                    write_state(state, cycle, time_step, delta_t, stateFile);
                }   
                idle_equation(inputs, alpha_f, alpha_s, delta_t, h, T_old, T_new);
            }

            //Solve the linear system 
            double A[2][2] = {{ 1 + (h_v_f*delta_t), (-1)*h_v_f*delta_t}, { (-1)*h_v_s*delta_t, 1 + h_v_s*delta_t}};
            for (int i = 0; i < inputs->N; i++)
            {
                double x[2] = {};
                double b[2][1] =  {{T_new[i][2]}, {T_new[i][1]}};
                luDecomposition(A,b,x);
                T_new[i][2] = x[0];
                T_new[i][1] = x[1];
            }

            //Error per iteration 
//            double error_i;
//            max_error = 0.0000;
//            for (int i = 1; i < inputs->N; i++)
//            {
//                error_i = fabs(T_new[i][2] - T_old[i][2]);
//                if (error_i > max_error)
//                    max_error = error_i;
//            }


            //Convergence Check 
//            double residual = 0;
//            for (int i = 0; i < inputs->N; i++)
//            {
//                residual = residual + fabs(T_new[i][2] - T_old[i][2]);
//            }
//
//            double error = residual/inputs->N;

            
            //OVS error
//           if (error <= 1e-5)
//            {
//                double E1 = 0;
//                double L1 = 0;
//
//                //Checking convergence for fluid 
//                int fs_state = 0; //0 for fluid and 1 for solid
//                
//                for (int i = 1; i < inputs->N; i++)
//                {
//                    //E1 = E1 + fabs(T_new[i][2] - MMS_func(n_wave, i*h, inputs->H, fs_state) / (2 + MMS_func(n_wave, i*h, inputs->H, fs_state)));
//                    E1 = E1 + fabs(T_new[i][2] - MMS_func(n_wave, i*h, inputs->H, fs_state));
//				      E2 = E2 + pow((T_new[i][2] - MMS_func(n_wave, i*h, inputs->H, fs_state),2);
//                }
//
//                E1 = E1/inputs->N;
//                E2 = sqrt(E2/inputs->N);
//                write_error(h, err, E1, E2, Pe, n_wave, errorFile);
//            } 
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
//            memcpy(T_old,T_new, (inputs->N) * sizeof(**T_new));
            for (int i = 0; i < inputs->N; i++)
            {
                T_old[i][0] = T_new[i][0];
                T_old[i][1] = T_new[i][1];
                T_old[i][2] = T_new[i][2];
            }
            //Free T_new 
            for (int i = 0; i < inputs->N; i++)
                delete [] T_new[i];
            delete [] T_new;


            //write temperature to file
            if (time_step%inputs->save_file==0)
            {
				cout << "Current time step: " << time_step << endl;
                write_temperature(T_old, inputs->N, cycle, time_step, delta_t, tempFile);
            }
      
			//Capacity factor and exergy at evey time step
			double T0 			   = 288.15;			
			double thermal_energy  = 0;
			double exergy_flux     = 0;
			for (int i = 0; i < inputs->N; i++)
			{
				thermal_energy = thermal_energy + PI/4*(inputs->D*inputs->D)*((inputs->epsilon*inputs->rho_f*inputs->Cp_f*(T_old[i][2] - inputs->T_d))*h 
							                                               + ((1-inputs->epsilon)*inputs->rho_s*inputs->C_s*(T_old[i][1] - inputs->T_d))*h);
																			 
			    exergy_flux = exergy_flux + inputs->m_f*inputs->Cp_f*(T_old[i][2]-T0-T0*log(T_old[i][2]/T0));
			}
			if (time_step%inputs->save_file==0)
			{
				write_data(exergy_flux,thermal_energy, thermal_energy_max, cycle, time_step, delta_t, dataFile);
			}
			
            time_step++;
		}

        cycle++;
    }

	cout << "End of simulation!";
	
    return 0;
}
