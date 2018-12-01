#include "functions.h"

int main()
{
    //Initialize the parameters used for the simulation
    variables inputs;

    //Read inputs
 //   read_inputs(&inputs);

    //Write simulation data
  //  write_data(inputs.N, inputs.Ti);
    
    cout << "Test MMS function" << endl; 

    double sol = MMS(10,4,5,0);

    cout << "The solution is: " << sol << endl;

	return 0;
}
