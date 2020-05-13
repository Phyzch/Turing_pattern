# include<iostream>
# include <math.h>
#include<fstream>
#include <sstream>
#include"Brusselator.h"
#include<time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <experimental/filesystem>
using namespace std;
namespace fs= std::experimental::filesystem;

// number of box in our simulation.  numbox * length = L (total length we are going to do simulation.)
int Brusselator_reaction_system::numbox=60;

// change size of box here.
double Brusselator_reaction::length = 0.1; // length stands for the size of box doing simulation.

//below are diffusion rates scaled with the length we set for each box(0.1)
double Brusselator_reaction::DVm = 0.012*1 / (Brusselator_reaction::length * Brusselator_reaction::length) ;
double Brusselator_reaction::DUm = 0.3*1 / (Brusselator_reaction::length * Brusselator_reaction::length);
double Brusselator_reaction::DUc = 0.3 / (Brusselator_reaction::length * Brusselator_reaction::length);

// change particle number N here
// Nparticle is particle number we can tune
double Brusselator_reaction::Nparticle =  10 *Brusselator_reaction::length *120/400;

//below are reaction coefficient beta_12, k_21, k_12, k_13, k_32, k_31, beta_21, k_23. remember we already tune beta12 to make it work with 10 times particle numbers, you can change it.
// Parameters only link with X_{1}, X_{2}
double Brusselator_reaction::beta_12 = 1.5*1e-4/(pow(10 * Brusselator_reaction::length *120/400,2)); // we stop rescale reaction parameters for particle number here.
double Brusselator_reaction::beta_21 = 0.1 *Brusselator_reaction::beta_12;
double Brusselator_reaction::k_21 = 3.6;
double Brusselator_reaction::k_12 =0.5;

// parameters link with X_{3}
double Brusselator_reaction::k_13 = log(2)/50;
double Brusselator_reaction::k_32 = log(2)/50;
double Brusselator_reaction::k_31 = 3 * Brusselator_reaction::k_13;
double Brusselator_reaction::k_23 = 0.001 *Brusselator_reaction::k_32;

// below are initial concentration of Um, Vm, Uc.
int Brusselator_reaction::Um0 = 123* Brusselator_reaction::Nparticle;
int Brusselator_reaction::Vm0 = 177* Brusselator_reaction::Nparticle;
int Brusselator_reaction::Uc0 = 100* Brusselator_reaction::Nparticle;

// total simulation time we are going to run in our program.
double simulation_time = 500;

//choice to append result to simulation result we already have or to do it from beginning.
bool append=false;

// early check if the pattern form 3 stripe pattern. (if not we will just stop simulation at t=200) (done by at t=200, we check the position of peak of pattern.)
bool early_check=false;

int main() {
    int i; // i means the file number, we cam generate many files simulated using different parameters or same parameters
	double timestep = 1;  //timestep for outputting result, set to 1 means we output concentration every 1 seconds
	bool load = true; // bool variable, if set to true, the programme will load information from save_data.txt and continue simulation. If load==false, it will delete all information and results and restart.
	//ifstream sample;
	char c;
	int Filenumber = 1;
	auto parent_path= "/home/phyzch/CLionProjects/fixed_position_Turing_pattern/save_data/test/";
	for (i = 1; i <= Filenumber; i++) {
		long timee = time(NULL);
		stringstream ss;
		string ch;
		ss << i;
		ss >> ch;
		string path = parent_path+ch;   // path_name, ch is the char form of number: "1","2","3",etc. Our result will be stored in these files.
        // check directory exists or not
		struct stat statbuf;
        int isDir = 0;
        if (stat(path.c_str(), &statbuf) != -1) {  // get permission to access directory
            if (S_ISDIR(statbuf.st_mode)) {
                isDir = 1;
            }
        }
        if (isDir==1){
            ;  // directory already exists
        }
        else {// create directory
            if (!mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)) {
                printf("File created!");
            } else {
                printf("Fail to create the file, directory may not exist.");
                exit(-1);
            }
        }

		Brusselator_reaction_system System(path, load, timestep, i);
		System.System_evolve(simulation_time);
	}
}