# include"Brusselator.h"
#include<time.h>
#include<sstream>
#include<iostream>
#include<fstream>
using namespace std;
long timee = time(NULL);
long timee1 = long(time(NULL));
double Brusselator_reaction::length = 0.1; // length stands for the size of compartment
//below are diffusion rates scaled with the length we set for each box(0.1)
double Brusselator_reaction::DVm = 0.012*1 / (Brusselator_reaction::length * Brusselator_reaction::length) ;
double Brusselator_reaction::DUm = 0.3*1 / (Brusselator_reaction::length * Brusselator_reaction::length);
double Brusselator_reaction::DUc = 0.3 / (Brusselator_reaction::length * Brusselator_reaction::length);
double Brusselator_reaction::V =  11 *Brusselator_reaction::length *120/400;      // change volume here
// V is the ratio of molecular number and concentration in one compartment, for V=1* length*120/400, it is around
//below are reaction coefficient beta, gamma, alpha, sigma, sigma1, epsilon(ep), beta1, sigma2
double Brusselator_reaction::beta = 1.5*1e-4/(pow(10 * Brusselator_reaction::length *120/400,2)); // we stop rescale reaction parameters for particle number here.
double Brusselator_reaction::gamma = 3.6;
double Brusselator_reaction::alpha =0.5;    // change here
double sigma_magnitude = 1;
double Brusselator_reaction::sigma = log(2)/50 *sigma_magnitude;
double Brusselator_reaction::sigma1 = log(2)/50 * sigma_magnitude;
double Brusselator_reaction::ep = 3 * Brusselator_reaction::sigma;
double Brusselator_reaction::beta1 = 0.1 *Brusselator_reaction::beta;     // change here
double Brusselator_reaction::sigma2 = 0.001 *Brusselator_reaction::sigma1;
// below are initial concentration of Um, Vm, Uc.
int Brusselator_reaction::Um0 = 123* Brusselator_reaction::V;
int Brusselator_reaction::Vm0 = 177* Brusselator_reaction::V;
int Brusselator_reaction::Uc0 = 100* Brusselator_reaction::V;
// below ar echemical reaction potential
double Brusselator_reaction::chemical_potential = log(Brusselator_reaction::beta / Brusselator_reaction::beta1  * Brusselator_reaction::gamma / Brusselator_reaction::alpha );
double Brusselator_reaction::c[M + 1] = { 0, beta, beta1, gamma, alpha, sigma, ep, sigma1, sigma2, DUm, DVm, DUc  };
// c[1] to c[8]: reaction coefficient. c[9] to c[11]: diffusion to next box
// For num of molecules see reference Murray 2017 Nature.Phys, Fig1.f around 10-30 molecules in each compartment
Brusselator_reaction & Brusselator_reaction::operator =(const Brusselator_reaction & s) {
	// copy function. I found I do not need it in our codes, but just attached it (History legacy...)
	if (this != &s) {
		int i;
		for (i = 0; i < N; i++) {
			this->S[i] = s.S[i];
		}
		for (i = 0; i < M + 1; i++) {
			this->a[i] = s.a[i];
		}
	};
	return *this;
	this->filename = s.filename;
	this->count = s.count;
};
Brusselator_reaction::~Brusselator_reaction() {
	// destructor
	delete[] S;
	delete[] a;
	if (output.is_open()) {
		output.close();
	}
}
Brusselator_reaction::Brusselator_reaction(int num1) { 
	// constructor
	long timee = time(NULL);
	S = new int[N]; // S[0]= Um, S[1]= Vm, S[2]= Uc
	a = new double[M + 1]; // reaction rate of different reactions
	count = 0; // count of number of reactions happened in our box
	num = num1; // number of index in the system
	S[0] = Um0;
	S[1] = Vm0;
	S[2] = Uc0;
}

void Brusselator_reaction::initialize_a() {
	// initialize reaction rate a
	// a[1] : Um + 2 Vm-> 3Vm,  a[2]: 3Vm-> Um+ 2 Vm, a[3]: Vm->Um,  a[4]: Um->Vm, a[5]: Um-> Uc, a[6]: Uc-> Um
	// a[7]: Vm-> Uc,  a[8]: Uc-> Vm
	a[1] = c[1] * S[0]* S[1] *S[1];
	a[2] = c[2] * S[1] *S[1] *S[1];
	a[3] = c[3] * S[1];
	a[4] = c[4] * S[0];
	a[5] = c[5] * S[0];
	a[6] = c[6] * S[2];
	a[7] = c[7] * S[1];
	a[8] = c[8] * S[2];
	// diffusion to the next box. a[9]: Um diffusion rate, a[10]: Vm reaction rate. a[11]: Uc reaction rate
	if (next != NULL) {
		a[9] = c[9] * abs(S[0] - next->S[0]);
		a[10] = c[10] * abs(S[1] - next->S[1]);
		a[11] = c[11] * abs(S[2] - next->S[2]);
	}
	else {
		a[9] = 0;
		a[10] = 0;
		a[11] = 0;
	}
	a[0] = 0;
	int i;
	// a[0]: total reaction rate for this reactions in this box.
	for (i = 1; i < M + 1; i++) {
		a[0] += a[i];
	}

}

void Brusselator_reaction::Gillespie_simulation(double &duration, double aa, double t, double &diffusion_W, double & reaction_W1, double & reaction_W2) {
	// Gillespie algorithm simulation in box. aa is random rate generated to decide the index of reaction we will perform. t is reaction time.
	// duration: duration of our simulation, used to check the speed of our simulation ( not used now, ignore it.)
	long timee = long(time(NULL));
	double random1;
	int miu1;
	double dW = 0;
	double bb1 = 0;
	int i;
	//double starttime, endtime;
	//starttime = clock();
	for (i = 1; i < M + 1; i++) {
		bb1 = bb1 + a[i];
		if (bb1 > aa) {
			miu1 = i;
			break;
		}
	}

	switch (miu1) {
	case 1:
		if (S[0] == 1 || S[1] == 0) {
			return;
		}
		S[0]--;
		S[1]++;
		dW = log(a[1] / a[2]);
		break;
	case 2:
		if (S[1] == 1 || S[0] == 0) {
			return;
		}
		S[0]++;
		S[1]--;
		dW = log(a[2] / a[1]);
		/*if (dW > 1e5) {
			printf("error");
		}*/
		break;
	case 3:
		if (S[1] == 1 || S[0] == 0) {
			return;
		}
		S[0]++;
		S[1]--;
		dW = log(a[3] / a[4]);
		break;
	case 4:
		if (S[0] == 1 || S[1] == 0) {
			return;
		}
		S[1]++;
		S[0]--;
		dW = log(a[4] / a[3]);
		break;
	case 5:
		if (S[0] == 1 || S[2] == 0) {
			return;
		}
		S[0]--;
		S[2]++;
		dW = log(a[5] / a[6]);
		break;
	case 6:
		if (S[2] == 1 || S[0] == 0) {
			return;
		}
		S[0]++;
		S[2]--;
		dW = log(a[6] / a[5]);
		break;
	case 7:
		if (S[1] == 1 || S[2]==0) {
			return;
		}
		S[1]--;
		S[2]++;
		dW = log(a[7] / a[8]);
		break;
	case 8:
		if (S[2] == 1 || S[1]==0) {
		return;
	}
		S[1]++;
		S[2]--;
		dW = log(a[8] / a[7]);
		break;
	case 9:
		if (S[0] > next->S[0]) {
			S[0]--;
			next->S[0]++;
		}
		else {
			S[0]++;
			next->S[0]--;
		}
		dW = abs((log(float(S[0]) / next->S[0])));
		break;
	case 10:
		if (S[1] > next->S[1]) {
			S[1]--;
			next->S[1]++;
		}
		else {
			S[1]++;
			next->S[1]--;
		}
		dW = abs((log(float(S[1]) / next->S[1])));
		break;
	case 11:
		if (S[2] > next->S[2]) {
			S[2]--;
			next->S[2]++;
		}
		else {
			S[2]++;
			next->S[2]--;
		}
		dW = abs((log(float(S[2]) / next->S[2])));
		break;
	}
	count++;
	// ����a
	a[1] = c[1] * S[0] * S[1] * S[1];
	a[2] = c[2] * S[1] * S[1] * S[1];
	a[3] = c[3] * S[1];
	a[4] = c[4] * S[0];
	a[5] = c[5] * S[0];
	a[6] = c[6] * S[2];
	a[7] = c[7] * S[1];
	a[8] = c[8] * S[2];
	// an extra step we have to take for using Gillespie algorithm to  simulate spatially extended system is to we have to update the change of 
	// diffusion rate in other box due to the change of species concentration
	if (next != NULL) {
		a[9] = c[9] * abs(S[0] - next->S[0]);
		a[10] = c[10] * abs(S[1] - next->S[1]);
		a[11] = c[11] * abs(S[2] - next->S[2]);
	}
	else {
		// last box.
		a[9] = 0;
		a[10] = 0;
		a[11] = 0;
	}
	if (num != 0) {
		// not first box, we have to updat previous box's diffusion rate and total reaction rate a[0]
		double previousrate = previous->a[9] + previous->a[10] + previous->a[11];
		previous->a[9] = c[9] * abs(previous->S[0] - S[0]);
		previous->a[10] = c[10] * abs(previous->S[1] - S[1]);
		previous->a[11] = c[11] * abs(previous->S[2] - S[2]);
		previous->a[0] = previous->a[0] - previousrate + previous->a[9] + previous->a[10] + previous->a[11];
	}
	if (miu1 == 9 || miu1 == 10 || miu1 == 11) {
		// species diffuse to next box, so we een have to update the next box's reaction rate
		next->a[1] = c[1] * next->S[0] * next->S[1] * next->S[1];
		next->a[2] = c[2] * next->S[1] * next->S[1] * next->S[1];
		next->a[3] = c[3] * next->S[1];
		next->a[4] = c[4] * next->S[0];
		next->a[5] = c[5] * next->S[0];
		next->a[6] = c[6] * next->S[2];
		next->a[7] = c[7] * next->S[1];
		next->a[8] = c[8] * next->S[2];
		if (next->next != NULL) {
			next->a[9] = c[9] * abs(next->S[0] - next->next->S[0]);
			next->a[10] = c[10] * abs(next->S[1] - next->next->S[1]);
			next->a[11] = c[11] * abs(next->S[2] - next->next->S[2]);
		}
		next->a[0] = 0;
		for (i = 1; i < M + 1; i++) {
			next->a[0] += next->a[i];
		}
	}
	a[0] = 0;
	for (i = 1; i < M + 1; i++) {
		a[0] += a[i];
	}
	// calculate energy dissipation
	dW = dW*1e-5;
	if (miu1 == 9 || miu1 == 10 || miu1==11) {
		diffusion_W += dW;
	}
	// reaction energy in small loop between Um and Vm
	else if (miu1==1 || miu1==2 || miu1==3 || miu1==4) {
		reaction_W1 += dW;
	}
	else {
		// reaction enregy between links of Vm, Uc, Um
		reaction_W2 += dW;
	}
}

