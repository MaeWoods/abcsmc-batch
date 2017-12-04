
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///                                                                                                                                    ///
/// Evolutionary model ///  Model 1: insertions + deletions Model 2: insertions + deletions + translocations                           ///
///                         Model 3: insertions + deletions + selection Model 4: insertions + deletions + translocations + selection   ///
/// ------------------------------------------------------------------------------------------------------------                       ///
///                                                                                                                                    ///
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef simulation_H
#define simulation_H

#include <fstream>
#include <string>
#include <cstring>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <string.h>
#include <vector>
#include <random>
#include <sstream>
#include <ctime>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>

using namespace std;

// to build the program:
//   g++ SVModel.cpp -o SVModel -L/usr/local/lib/ -lgsl -I/usr/local/include
//
// to run the model:
//   ./SVModel Seed=1 Model=2


// order of parameters expected:  Npop, ngen, g_d, mu_i, p_tran, trp2, svtransgrad, threshold_g, SV_max, mu_k, gnmdou, maxchr, minchr

class simulation {

private:


  //constants
  int NChr, Npop, ngen;
  double lowB, chrlength, Curvemax, minchr, mu_b, minSV, Lchr_loss,Uchr_loss,Lchr_gain,Uchr_gain,Lp_tran,Up_tran,Lmu_ki,Umu_ki,Lmu_kd,Umu_kd,Lfitness,Ufitness,Lmaxchr,Umaxchr,LSV_mean,USV_mean,Lgnmdou,Ugnmdou;
  
  //parameters
  vector<double> All_mu_b, All_chr_loss, All_chr_gain, All_p_tran, All_fitness,  All_SV_mean, All_mu_ki, All_mu_kd, All_gnmdou, All_maxchr, All_minchr;

	vector<vector<vector<double> > > Mprev;
	vector<vector<vector<double> > > MCprev;
	vector<vector<vector<double> > > Cprev;
	vector<vector<vector<double> > > DivCprev;
	vector<vector<vector<double> > > CMixprev;	

	vector<vector<vector<double> > > NTprev;
	vector<vector<vector<double> > > NIprev;
	vector<vector<vector<double> > > NDprev;

	vector<vector<double> > rdelprev;
	vector<vector<double> > rinsprev;
	vector<vector<double> > rtransprev;
	vector<int> GDprev;
	vector<double> GDprobprev;
	vector<double> CSize;

	vector<vector<vector<double> > > M;
	vector<vector<vector<double> > > MC;
	vector<vector<vector<double> > > C;
	vector<vector<vector<double> > > DivC;
	vector<vector<vector<double> > > CMix;

	vector<vector<vector<double> > > NT;
	vector<vector<vector<double> > > NI;
	vector<vector<vector<double> > > ND;

	vector<vector<double> > rdel;
	vector<vector<double> > rins;
	vector<vector<double> > rtrans;
	vector<int> GD;
	vector<double> GDprob;


public:

simulation(void);

~simulation(void);

int read_params(const string& filename);

int read_constnts( const string& filename, int &NChr, int &Npop, int &ngen, double &lowB, double &chrlength, double &Curvemax, double &minchr, double &mu_b, double &minSV, double &Lchr_loss, double &Uchr_loss, double &Lchr_gain, double &Uchr_gain, double &Lp_tran, double &Up_tran, double &Lmu_ki, double &Umu_ki, double &Lmu_kd, double &Umu_kd, double &Lfitness, double &Ufitness, double &Lmaxchr, double &Umaxchr, double &LSV_mean, double &USV_mean, double &Lgnmdou, double &Ugnmdou );

double runiform(gsl_rng* r, double a, double b);

void simulate(vector<double> parameters, double epsilon, int position, int seed, int Rsim, gsl_rng* r);

};

#endif












