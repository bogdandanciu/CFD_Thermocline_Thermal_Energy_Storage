#include "functions.h"

int main()
{
    //Initialize the parameters used for the simulation
    variables inputs;

    //Read inputs
    read_inputs(&inputs);

    //Write simulation data
    write2file(inputs.N, inputs.Ti);

	return 0;
}
