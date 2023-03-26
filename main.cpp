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


//choice to append result to simulation result we already have or to do it from beginning.
bool append=false;

// early check if the pattern form 3 stripe pattern. (if not we will just stop simulation at t=200) (done by at t=200, we check the position of peak of pattern.)
bool early_check=false;

bool load= false;// bool variable, if set to true, the programme will load information from save_data.txt and continue simulation. If load==false, it will delete all information and results and restart.

int main() {
    int i; // i means the file number, we cam generate many files simulated using different parameters or same parameters
	double timestep = 1;  //timestep for outputting result, set to 1 means we output concentration every 1 seconds
	//ifstream sample;
	char c;
	int Filenumber = 1;
	auto parent_path= "/home/phyzch/CLionProjects/fixed_position_Turing_pattern/save_data/test/";
	for (i = 1; i <= Filenumber; i++) {
		string ch;
        ch = std::to_string(i);
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
		System.System_evolve();
	}
}