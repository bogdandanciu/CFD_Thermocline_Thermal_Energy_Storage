#include <iostream>
#include <fstream>
#include <string> 

using namespace std;


int main()
{
	//Take user input routine 
	double diam, h, N, T_init;

//	cout << "Please input The diameter for the storage: ";
//        cin >> diam;	
//	cout << "Please input the height of the storage: ";
//        cin >> h;	
	cout << "Please input the number of cells: ";
	cin >> N;
	cout << "Please input the the initial temperature of the fluid and solid phases: ";
	cin >> T_init;

	//Write data to a text file
//        ofstream myfile("example.txt") ;		
//	if (myfile.is_open())
//	{
//		myfile << "This is a line.\n";
//		myfile << "This is another line.\n";
//		myfile.close();
//	}
//	else 
//	{
//		cout << "Unable to open file" << endl;	
//	}


	//Write simulations data
	ofstream simFile("simData.dat");
	if (simFile.is_open())
	{
		for (int i = 0; i < N; i++)
		{
			simFile << i << " " << T_init << endl;
		}
		simFile.close();
	}
	else cout << "Unable to open file" << endl; 
	

	return 0;

}
