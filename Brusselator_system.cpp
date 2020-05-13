#include"Brusselator.h"
#include<time.h>
#include<string.h>
#include<sstream>
#include<stdlib.h>
Brusselator_reaction_system::Brusselator_reaction_system(string path1, bool load, double timestep1, int sequence_number1, double * sampling) {
	// constructor
	srand((unsigned)time(NULL));
	t = 0;
	timestep = timestep1;
	diffusion_W = 0;
	reaction_W1 = 0;
	reaction_W2 = 0;
	sequence_number = sequence_number1;
	path = path1;
	// create memory space for Brusselator_reaction class entity, with number of numsystem(60 in our codes)
	void * rawMemory = operator new(numsystem * sizeof(Brusselator_reaction));
	s = reinterpret_cast<Brusselator_reaction *>(rawMemory);
	int i;
	for (i = 0; i < numsystem; i++) {
		new (&s[i]) Brusselator_reaction(i);  // constructor Brusselator_reaction_entity, with index i
	}
	long timee = time(NULL);
	/*Brusselator_reaction::beta1 = sampling[0] * Brusselator_reaction::beta;
	Brusselator_reaction::c[2] = Brusselator_reaction::beta1;
	Brusselator_reaction::alpha = sampling[1];
	Brusselator_reaction::c[4] = Brusselator_reaction::alpha;*/
	// uncomment above codes if you want to do simulation with random parameters stored in ./data/sampling/samplingparameter.txt

	if (load == true) {
		load_data();
	}
	// link different box together using link list
	if (numsystem > 1) {
		for (i = 0; i < numsystem; i++) {
			if (i == 0) {
				s[i].previous = NULL;
				s[i].next = &s[i + 1];
			}
			else if (i == numsystem - 1) {
				s[i].next = NULL;
				s[i].previous = &s[i - 1];
			}
			else {
				s[i].next = &s[i + 1];
				s[i].previous = &s[i - 1];
			}
		}
	}
	else {
		s[0].previous = s[0].next = NULL;
	}
    // initialize our total reaction rate a
	for (i = 0; i < numsystem; i++) {
		s[i].initialize_a();
	}
	stringstream ss;
	string ch;
	ss << sequence_number;
	ss >> ch;

	//���ղ�����ͬ�洢�ڲ�ͬsequence_number���ļ�������
	filename.assign("Brusselator_system.txt");
	if (load == true) {
		if (append == true) {
			system_output.open(path + "/" + filename,ios::app);
			diffusion_energy_file.open(path + "/diffusion_energy.txt",ios::app);
			reaction_energy_file1.open(path + "/reaction_energy1.txt",ios::app);
			reaction_energy_file2.open(path + "/reaction_energy2.txt",ios::app);
			information.open(path + "/information.txt");
		}
		else {
			system_output.open(path + "/" + filename);  // output our concentration to "Brusselator_system.txt"
			diffusion_energy_file.open(path + "/diffusion_energy.txt"); // output our energy dissipation to "diffusion_energy.txt"
			reaction_energy_file1.open(path + "/reaction_energy1.txt");
			reaction_energy_file2.open(path + "/reaction_energy2.txt");
			information.open(path + "/information.txt"); // output our reaction coefficient information to information.txt
		}
	}
	else {
		system_output.open(path + "/" + filename);  // output our concentration to "Brusselator_system.txt"
		diffusion_energy_file.open(path + "/diffusion_energy.txt"); // output our energy dissipation to "diffusion_energy.txt"
		reaction_energy_file1.open(path + "/reaction_energy1.txt");
		reaction_energy_file2.open(path + "/reaction_energy2.txt");
		information.open(path + "/information.txt"); // output our reaction coefficient information to information.txt
	}
};


Brusselator_reaction_system::~Brusselator_reaction_system() {
	// deconstructor function. This is required because we have to run System with different parameters everytime
	// so we have to make sure we have released memory when we deconstruct it, otherwise, we will see memory blow!
	int i;
	for (i = 0; i < numsystem; i++) {
		s[i].~Brusselator_reaction();
	}
	system_output.close();
	diffusion_energy_file.close();
	reaction_energy_file1.close();
	reaction_energy_file2.close();
}

int Brusselator_reaction_system::maxVm() { // output Vm's max concentration
	int maxVmparticlenumber = 0;
	int i;
	for (i = 0; i < numsystem; i++) {
		if (s[i].S[1] > maxVmparticlenumber) {
			maxVmparticlenumber = s[i].S[1];
		}
	}
	return maxVmparticlenumber;
}



void Brusselator_reaction_system::System_evolve() {
	int i;
	bool mark;
	double a;
	int maxVmparticlenumber;
	double random1, random2;
	long timee = long(time(NULL));
	double tau;
	double bb, aa;
	int miu;
	a = 0;
	bool mark1 = false;
	double starttime, endtime, duration, duration1;
	duration = 0;
	duration1 = 0;
	double deltaa = 0, deltanexta=0, deltaprevious=0;  // change of reaction rates a for box, previous box, next box
	information << Brusselator_reaction::beta << "," << Brusselator_reaction::beta1 << "," << Brusselator_reaction::gamma << "," <<
		Brusselator_reaction::alpha << "," << Brusselator_reaction::DUm << "," << Brusselator_reaction::DVm << "," << Brusselator_reaction::DUc << "," <<
		Brusselator_reaction::sigma << "," << Brusselator_reaction::sigma1 << "," << Brusselator_reaction::sigma2 << "," << Brusselator_reaction::ep << endl;

	information.close();

	for (i = 0; i < numsystem; i++) {// Randomly decide which box to simulate, using Gillespie algorithm
		a += s[i].a[0];
	}
	do random1 = ran2(timee); while (random1 == 0);
	tau = 1 / a * log(1 / random1);
	while (t < 5000) { // t: total time for simulation

		// Combining mark to output our results. every moment we approach integer times of timestep, we output our simulation results
		if ((t - int(t / timestep + 1)*timestep) > -5 * tau) {
			mark1 = false;
		}
		if ((t - int(t / timestep)*timestep) < 2 * tau && mark1 == false) {
			system_output << t;
			for (i = 0; i < numsystem; i++) {
				system_output << "," << s[i].S[0];
			}
			system_output << endl;
			system_output << t;
			for (i = 0; i < numsystem; i++) {
				system_output << "," << s[i].S[1];
			}
			system_output << endl;
			system_output << t;
			for (i = 0; i < numsystem; i++) {
				system_output << "," << s[i].S[2];
			}
			system_output << endl;
			diffusion_energy_file << t << "," << diffusion_W << endl;
			reaction_energy_file1 << t << "," << reaction_W1 << endl;
			reaction_energy_file2 << t << "," << reaction_W2 << endl;
			diffusion_W = 0; 
			reaction_W1 = 0;
			reaction_W2 = 0;
			mark1 = true;
		}

		// We make decision if we need to discard this simulation
		if (abs(t - 200) < 2 * tau  ) {
			maxVmparticlenumber = maxVm();
			if ((s[numsystem / 6].S[1] > 0.6*maxVmparticlenumber|| s[numsystem / 6-1].S[1] > 0.6*maxVmparticlenumber || s[numsystem / 6+1].S[1] > 0.6*maxVmparticlenumber)
				& (s[3 * numsystem / 6].S[1] > 0.6*maxVmparticlenumber || s[3 * numsystem / 6-1].S[1] > 0.6*maxVmparticlenumber || s[3 * numsystem / 6+1].S[1] > 0.6*maxVmparticlenumber) & (s[5 * numsystem / 6].S[1] > 0.6*maxVmparticlenumber || s[5 * numsystem / 6-1].S[1] > 0.6*maxVmparticlenumber || s[5 * numsystem / 6+1].S[1] > 0.6*maxVmparticlenumber)) {
				// we could continue our simulation
				cout << sequence_number << endl;
			}
			else {
				break; // pattern not right, discard this data.
			}
		}
		do random1 = ran2(timee); while (random1 == 0);
		do random2 = ran2(timee); while (random2 == 0);
		tau = 1 / a * log(1 / random1);
		bb = 0;
		aa = random2*a;
		//starttime = clock();
		// first decide which box to simulate ( miu )
		for (i = 0; i < numsystem; i++) {
			bb = bb + s[i].a[0];
			if (bb > aa) {
				miu = i;
				break;
			}
		}
		//endtime = clock();
		//duration += endtime - starttime;
		aa = aa - (bb - s[miu].a[0]);
		// now aa can be used in Brusselator_reaction::Gillespie_simulation() to decide which reaction the system take
 		if (aa > s[miu].a[0]) {
			system("pause"); // this means our code has something wrong here
		}
		deltaa = -s[miu].a[0]; // We record the original value of total reaction rate a in box miu
		if (miu != numsystem-1) {
			deltanexta = -s[miu + 1].a[0];  // original value of total reaction rate a in next box 
		}
		if (miu != 0) {
			deltaprevious = -s[miu - 1].a[0]; // original value of total reaction rate a in previous box
		}
		t = t + tau;
		s[miu].Gillespie_simulation(duration1, aa, t, diffusion_W, reaction_W1,reaction_W2);
		// update total reaction rate in whole system by only calculating change of reaction rate in our chosen box and box next to it
		deltaa += s[miu].a[0];
		a = a + deltaa;
		if (miu != numsystem - 1) {
			if (deltanexta != -s[miu + 1].a[0]) {
				deltanexta = deltanexta + s[miu + 1].a[0];
				a = a + deltanexta;
			}
		}
		if (miu != 0) {
			deltaprevious = deltaprevious + s[miu - 1].a[0];
			a = a + deltaprevious;
		}

	}


}

void Brusselator_reaction_system::load_data() {
	// load data from ./save_data.txt
	ifstream infile;
	infile.open("/home/phyzch/CLionProjects/fixed_position_Turing_pattern//save_data.txt");
	stringstream ss;
	int number;
	char c;
	int i = 0;
	if (infile.is_open()) {
		while (!infile.eof()) {
			infile >> number;
			s[i].S[0] = number;
			infile >> c;
			infile >> number;
			s[i].S[1] = number;
			infile >> c;
			infile >> number;
			s[i].S[2] = number;
			infile >> c;
			i++;
			if (i == numsystem) {
				break;
			}
		}
		infile >> t;
		infile >> c;
		if (append == false) {
			// we begin from beginning
			t = 0;
		}
		infile >> Brusselator_reaction::beta;
		infile >> c;
		infile >> Brusselator_reaction::beta1;
		infile >> c;
		infile >> Brusselator_reaction::gamma;
		infile >> c;
		infile >> Brusselator_reaction::alpha;
		infile >> c;
		infile >> Brusselator_reaction::DUm;
		infile >> c;
		infile >> Brusselator_reaction::DVm;
		infile >> c;
		infile >> Brusselator_reaction::DUc;
		infile >> c;
		infile >> Brusselator_reaction::sigma;
		infile >> c;
		infile >> Brusselator_reaction::sigma1;
		infile >> c;
		infile >> Brusselator_reaction::sigma2;
		infile >> c;
		infile >> Brusselator_reaction::ep;
		infile >> c;
		Brusselator_reaction::c[1] = Brusselator_reaction::beta;
		Brusselator_reaction::c[2] = Brusselator_reaction::beta1;
		Brusselator_reaction::c[3] = Brusselator_reaction::gamma;
		Brusselator_reaction::c[4] = Brusselator_reaction::alpha;
		Brusselator_reaction::c[5] = Brusselator_reaction::sigma;
		Brusselator_reaction::c[6] = Brusselator_reaction::ep;
		Brusselator_reaction::c[7] = Brusselator_reaction::sigma1;
		Brusselator_reaction::c[8] = Brusselator_reaction::sigma2;
		Brusselator_reaction::c[9] = Brusselator_reaction::DUm;
		Brusselator_reaction::c[10] = Brusselator_reaction::DVm;
		Brusselator_reaction::c[11] = Brusselator_reaction::DUc;
	}
	infile.close();
}