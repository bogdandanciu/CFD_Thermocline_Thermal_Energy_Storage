#include "functions.h"

int main()
{
    //Initialize the parameters used for the simulation
    variables inputs;

    //Get required inputs to run the simulation 
    int choice;
    cout << "Choose 0 for default values or Choose 1 to input the values manually: ";
    cin >> choice; 
    read_inputs(&inputs,choice); 

    //Run main solver
    solver(&inputs);


	return 0;
}
