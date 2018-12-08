#include "functions.h"

int main()
{
    //Initialize the parameters used for the simulation
    variables inputs;

    //Read inputs
 //   read_inputs(&inputs);

    //Write simulation data
  //  write_data(inputs.N, inputs.Ti);
    
//    cout << "Test MMS function" << endl; 
//
//    double sol = MMS(10,4,5,0);
//
//    cout << "The solution is: " << sol << endl;
    
    //Linear solver test 
//    double A[2][2] = {{234,24},{2,353}};
//    double b[2][1] = {{25453},{3425}};
//    double x[2]    = {};
//
//    //Solve system
//    luDecomposition(A,b,x);
//    cout << "This is a linear sovler" << endl;
//    cout << "x1 " << x[0] << " and " << "x2 " << x[1];

    cout << "Start of the Simulation" << endl; 
    
    //Get required inputs to run the simulation 
    int choice;
    cout << "Choose 0 for default values or Choose 1 to input the manually\n";
    cin >> choice; 
    read_inputs(&inputs,choice); 

    //Run main solver
    solver(&inputs);



	return 0;
}
