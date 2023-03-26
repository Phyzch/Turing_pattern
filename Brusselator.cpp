# include"Brusselator.h"
#include<time.h>
#include<sstream>
#include<iostream>
#include<fstream>
using namespace std;
long timee = time(NULL);
long timee1 = long(time(NULL));
// below are chemical reaction potential
double Brusselator_reaction::chemical_potential = log(Brusselator_reaction::beta_12 / Brusselator_reaction::beta_21  * Brusselator_reaction::k_21 / Brusselator_reaction::k_12 );
double Brusselator_reaction::c[M + 1] = {0, beta_12, beta_21, k_21, k_12, k_13, k_31, k_32, k_23, DX1, DX2, DX3  };
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
}
Brusselator_reaction::Brusselator_reaction(int num1) { 
	// constructor
	long timee = time(NULL);
	S = new int[N]; // S[0]= Um, S[1]= Vm, S[2]= Uc
	a = new double[M + 1]; // reaction rate of different reactions
	count = 0; // count of number of reactions happened in our box
	num = num1; // number of index in the system
	S[0] = X1_0;
	S[1] = X2_0;
	S[2] = X3_0;
}

void Brusselator_reaction::initialize_a() {
	// initialize reaction rate a
    // a[1] : X1 + 2 X2 -> 3X2, a[2] : 3X2 -> X1 + 2 X2, a[3]: X2 -> X1. a[4]: X1 -> X2. a[5]: X1 -> X3, a[6]: X3 -> X1.
    // a[7]: X2 -> X3,  a[8]: X3 -> X2
	a[1] = c[1] * S[0]* S[1] *S[1];
	a[2] = c[2] * S[1] *S[1] *S[1];
	a[3] = c[3] * S[1];
	a[4] = c[4] * S[0];
	a[5] = c[5] * S[0];
	a[6] = c[6] * S[2];
	a[7] = c[7] * S[1];
	a[8] = c[8] * S[2];
	// diffusion to the next box. a[9]: X1 diffusion rate, a[10]: X2 diffusion rate. a[11]: Uc diffusion rate
	if (next != NULL) {
		a[9] = c[9] * S[0];
		a[10] = c[10] * S[1];
		a[11] = c[11] * S[2];
		a[12]= c[9] * next->S[0];
		a[13]= c[10] * next->S[1];
		a[14]= c[11] * next->S[2];
	}
	else {
		a[9] = 0;
		a[10] = 0;
		a[11] = 0;
		a[12]=0;
		a[13]=0;
		a[14]=0;
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

	for (i = 1; i < M + 1; i++) {
		bb1 = bb1 + a[i];
		if (bb1 > aa) {
			miu1 = i;
			break;
		}
	}

	switch (miu1) {
        // S[0]: X1, S[1]: X2, S[2]: X3
	case 1:
        // X1 + 2X2 -> 3 X2
		if (S[0] == 1 || S[1] == 0) {
			return;
		}
		S[0]--;
		S[1]++;
		dW = log(a[1] / a[2]);
		break;
	case 2:
        // 3 X2 -> X1 + 2 X2
		if (S[1] == 1 || S[0] == 0) {
			return;
		}
		S[0]++;
		S[1]--;
		dW = log(a[2] / a[1]);
		break;
	case 3:
        // X2 -> X1
		if (S[1] == 1 || S[0] == 0) {
			return;
		}
		S[0]++;
		S[1]--;
		dW = log(a[3] / a[4]);
		break;
	case 4:
        // X1 -> X2
		if (S[0] == 1 || S[1] == 0) {
			return;
		}
		S[1]++;
		S[0]--;
		dW = log(a[4] / a[3]);
		break;
	case 5:
        // X1 -> X3
		if (S[0] == 1 || S[2] == 0) {
			return;
		}
		S[0]--;
		S[2]++;
		dW = log(a[5] / a[6]);
		break;
	case 6:
        // X3 -> X1
		if (S[2] == 1 || S[0] == 0) {
			return;
		}
		S[0]++;
		S[2]--;
		dW = log(a[6] / a[5]);
		break;
	case 7:
        // X3 -> X2
		if (S[1] == 1 || S[2]==0) {
			return;
		}
		S[1]++;
		S[2]--;
		dW = log(a[7] / a[8]);
		break;
	case 8:
        // X2 -> X3
		if (S[2] == 1 || S[1]==0) {
		return;
	}
		S[1]--;
		S[2]++;
		dW = log(a[8] / a[7]);
		break;
	case 9:
        // diffusion for X1 to next box
        S[0]--;
        next->S[0]++;
		dW = log(float(S[0]) / next->S[0]);
		break;
	case 10:
        S[1]--;
        next->S[1]++;
		dW = log(float(S[1]) / next->S[1]);
		break;
	case 11:
        S[2]--;
        next->S[2]++;
		dW = log(float(S[2]) / next->S[2]);
		break;
    case 12:
        S[0]++;
        next->S[0]--;
        dW=log(float(next->S[0])/S[0]);
        break;
    case 13:
        S[1]++;
        next->S[1]--;
        dW=log(float(next->S[1])/S[1]);
        break;
    case 14:
	    S[2]++;
	    next->S[2]--;
	    dW=log(float(next->S[2])/S[2]);
	    break;
	}
	count++;
	a[1] = c[1] * S[0] * S[1] * S[1];
	a[2] = c[2] * S[1] * S[1] * S[1];
	a[3] = c[3] * S[1];
	a[4] = c[4] * S[0];
	a[5] = c[5] * S[0];
	a[6] = c[6] * S[2];
	a[7] = c[7] * S[2];
	a[8] = c[8] * S[1];

	// an extra step we have to take for using Gillespie algorithm to  simulate spatially extended system is to we have to update the change of
	// diffusion rate in other box due to the change of species concentration
	if (next != NULL) {
        a[9] = c[9] * S[0];
        a[10] = c[10] * S[1];
        a[11] = c[11] * S[2];
        a[12]= c[9] * next->S[0];
        a[13]= c[10] * next->S[1];
        a[14]= c[11] * next->S[2];
	}
	else {
		// last box.
        a[9] = 0;
        a[10] = 0;
        a[11] = 0;
        a[12]=0;
        a[13]=0;
        a[14]=0;
	}
	if (num != 0) {
		// not first box, we have to update previous box's diffusion rate and total reaction rate a[0]
		double previousrate = previous->a[12] + previous->a[13] + previous->a[14];
		previous->a[12] = c[9] * S[0];
		previous->a[13] = c[10] * S[1];
		previous->a[14] = c[11] * S[2];
		previous->a[0] = previous->a[0] - previousrate + previous->a[12] + previous->a[13] + previous->a[14];
	}
	if (miu1 == 9 || miu1 == 10 || miu1 == 11|| miu1==12||miu1==13||miu1==14) {
		// species diffuse to next box, so we  have to update the next box's reaction rate
		next->a[1] = c[1] * next->S[0] * next->S[1] * next->S[1];
		next->a[2] = c[2] * next->S[1] * next->S[1] * next->S[1];
		next->a[3] = c[3] * next->S[1];
		next->a[4] = c[4] * next->S[0];
		next->a[5] = c[5] * next->S[0];
		next->a[6] = c[6] * next->S[2];
		next->a[7] = c[7] * next->S[1];
		next->a[8] = c[8] * next->S[2];
		if (next->next != NULL) {
			next->a[9] = c[9] * next->S[0] ;
			next->a[10] = c[10] *next->S[1];
			next->a[11] = c[11] * next->S[2];
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

	// calculate energy dissipation (scale it by 10^{-5}, therefore, the energy dissipation is in unit of 10^{-5} k_{B} T )
	dW = dW*1e-5;
	if (miu1 == 9 || miu1 == 10 || miu1==11) {
		diffusion_W += dW;
	}
	// reaction energy in small loop between X1 and X2 (Gamma is the ratio of product of rate)
	else if (miu1==1 || miu1==2 || miu1==3 || miu1==4) {
		reaction_W1 += dW;
	}
	else {
		// reaction enregy between links of X1, X2, X3 (Gamma' is the ratio of product of rate)
		reaction_W2 += dW;
	}
}

