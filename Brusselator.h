#pragma once
#pragma once
#pragma once
# include<iostream>
# include <math.h>
#include <fstream>
#include<string.h>
using namespace std;
extern bool append;
float ran2(long& idum);
class Brusselator_reaction {
	// Class for different local reaction box
private:
	const static int N = 3; // Number of reactant species, we have Um,Vm,Uc, so N=3
	const static int M = 11;// Number of Reactions (including diffusioin between boxes)
	string filename; // filename for outputing files. (I set it to "Brusselator_System" in constructor, you can set it to the name you want.
	ofstream output; 
	int num; //index for reaction box. If we have 60 reaction box (each box have l=0.1, total length is 6), num should be in range [0,60) and it is integer
	Brusselator_reaction * next;
	Brusselator_reaction * previous; // linked list. Link to the box in front of it and box behind it. Use it to do diffusion simulation
	static double chemical_potential; // chemical potential, defined as ln(k_21+/k_21-)
	static int Um0, Uc0, Vm0; // initial species concentration.  You can decide it using homogeneous spatial concentration simulated by matlab code
	static double Nparticle; // Volume, we scale the total number of species by multiple all of them by V.
public:
	static double length; // length is the size of compartment (0.1) and L is the size of whole cell, L= Brusselator_reaction_system::numbox*length
	static double c[M + 1]; //  Reaction coefficient, we start from c[1], static member variable, initialize it in Brsselator.cpp (For detailed setting of Parameter, check Murray 2017 Nature Physics: Self-Organization and Positioning of bacterial Protein clusters
	static double DUm, DVm, DUc; // diffusion rate for species Um, Vm, Uc
	static double beta_12, k_21, k_12, k_13, k_32, k_31; // reactioin coefficient, same as defined in our note
	static double k_23, beta_21; //  beta_21 is reverse reaction coefficient of beta_12, k_23 is reverse reaction rate for k_32. Check our SI for detail
	friend class Brusselator_reaction_system; // Brusselator_reaction_System, whole reaction system containing many Brusselator_reaction
	int * S; //  List for different species
	double * a; // List for reaction rate
	int count; // count number of reactions

	// Below is constructor, destructor, copyfunction
	Brusselator_reaction() {
		S = new int[N];
		a = new double[M + 1];
	}; //default constructor, only assign memory for list
	Brusselator_reaction(int num); // constructor, getting  num from Brusselator_reaction_system
	~Brusselator_reaction(); //deconstructor
	void initialize_a(); // initialize reaction rate a
	Brusselator_reaction & operator =(const Brusselator_reaction & s); // & means reference, copy function
																	  
	void Gillespie_simulation(double &duration, double aa, double t, double &diffusion_W, double & reaction_W1, double & reaction_W2); // Function do simulation, when this box is chosen, it will start Gillespie simulation
};



class Brusselator_reaction_system
{
private:
	int sequence_number; // index indicating filenumber (exactly i defined in main.cpp used to initialize Brusselator_reaction_System)
	static int numbox;
	Brusselator_reaction * s; //List of Brusselator_reactions
	string filename;
	string path;
	ofstream system_output;
	ofstream diffusion_energy_file;  // file for output diffusion_energy
	ofstream reaction_energy_file1;  // file for output reaction_energy in small loop(main loop)
	ofstream reaction_energy_file2;  // file for output reaction_energy in other linkes linked with Uc
	ofstream information; // file store information about reaction coefficient and diffusion rate about this simulation.
	double t; // time of simulation
	double diffusion_W; // diffusion_energy variable
	double reaction_W1,reaction_W2; // reaction_energy variable
public:
	double timestep;  // timestep for outputting results to files
	Brusselator_reaction_system(string path, bool load, double timestep, int sequence_number1);
	// constructor, path is the path for storing files, load is bool variable to decide if we start simulation from scratch or continue, sequence_number is index of simulation ( i in main.cpp), sampling is the list used when we want to randomly sampling parameters (Set sampling=[0,0] if we want to turn this mode off and comment out lines as indicated in main.cpp and Brusselator_system.cpp)
	void System_evolve(double simulation_time); // function that do simulation.
	~Brusselator_reaction_system(); // destrctor

	void load_data(); // function that load data from save_data.txt
	void save_data(); // function that save data into save_data.txt

	int maxVm(); // return position for maximum concentration of Vm (The reason we use it is we check if we generate right pattern early by checking the maximum concentration's position at 200s)
	friend class Brusselator_reaction;
};