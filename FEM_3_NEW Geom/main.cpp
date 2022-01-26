#include<iostream>
#include<Eigen/Dense>
#include<iomanip>
#include<fstream>
#include"material.h"
#include"element.h"

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/

using namespace std;
using namespace Eigen;

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/

void main() {

	cout << "Starting the Finite Element Program" << endl;
	//Creating csv files for storing stress, train and displacement data

	ofstream mydisplacefile("displacement.csv");
	ofstream mystressfile("stress.csv");
	ofstream mystrainfile("strain.csv");

	double ForceMax = 4700;			//Maximum Force
	double Area_1 = 12;					//csArea of left element
	double Area_2 = 4;				//csArea of right element
	double Length = 180;				//Length 
	
	double E = 120000;				//Youngs Modulus
	int steps_num = 10000;			//Number of Steps to solve the problem
	double steps = steps_num;
	double t_tot = 0.025;				//Total time for which force is applied
	double yield = 240;			//Yield Strength
	int elements = 25;				//Number of elements in which problem is to be divided
	double force = 0;					
	int gauss_weights = 2;
	double eta = 20;					//Viscosity
	MatrixXd Areas, Lengths;

	steps = t_tot / double(steps);

	//Calling Element Routine to get Areas & Lengths of all the elements

	tie(Areas, Lengths) = get_areas_lengths(elements, Area_1, Area_2, Length);	

	/*cout << "Areas: \n" << Areas << endl;
	cout << "Lengths: \n" << Lengths << endl;*/

	MatrixXd F_ext = MatrixXd::Zero(elements + 1,1);
	MatrixXd E_mat= E*MatrixXd::Ones(elements , 1);
	MatrixXd eps_pl = MatrixXd::Zero(elements , 1);
	MatrixXd epsilon = MatrixXd::Zero(elements , 1);
	MatrixXd u = MatrixXd::Zero(elements + 1, 1);
	MatrixXd Kt, F_int, strain, stress;
	
	int count = 0;

	//Iterating through all force steps

	while (force < ForceMax) {
		if (count != 0) { force = force + (steps / t_tot) * ForceMax; }		//Force Increment
		
		count++;
		
		//cout << "Force = " << force << endl;

		F_ext(int(elements / 2), 0) = force;								//Assigning the node to apply force in Global External Force Matrix

		//Entering the Newton-Raphson Method

		for (int NR = 0; NR < 6; NR++) {

			//Calling Material Routine to get the Global Tangent Stiffness Matrix, Global Internal Force Matrix, Plastic Strain, Total Strain & Stress Values of All Elements

			tie(Kt, F_int, eps_pl, strain, stress) = material_routine(elements, gauss_weights, Lengths, Areas, E_mat, u, steps, eta, eps_pl, yield, E);
			
			MatrixXd Kt_red = Kt.block(1, 1, Kt.cols() - 2, Kt.rows() - 2);		//Reducing the Tangent Stiffness
			MatrixXd inv_Kt_red = Kt_red.inverse();

			MatrixXd G = F_int - F_ext;
			//cout << "G=" << G << endl;
			MatrixXd G_red(elements - 1, 1),u_red(elements-1,1);					//Reducing the Residual 
			for (int c = 1; c < elements; c++) { G_red(c - 1, 0) = G(c, 0); }
			for (int c = 1; c < elements; c++) { u_red(c - 1, 0) = u(c, 0); }

			MatrixXd del_u_red = inv_Kt_red * G_red;								//Computing Delta U
			
			u_red = u_red - del_u_red;
			for (int c = 1; c < elements; c++) { u(c, 0) = u_red(c - 1, 0); }

		}
		
		//Appending values of Center Node/First Element to the csv files

		mydisplacefile << u((elements / 2), 0) << endl;
		mystressfile << stress((elements / 2) - 1, 0) << endl;
		mystrainfile << strain((elements / 2) - 1, 0) << endl;
		
	}

	mydisplacefile.close();
	mystressfile.close();
	mystrainfile.close();

	cout << "Computation Done ! Exiting the Program" << endl;
}

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/
