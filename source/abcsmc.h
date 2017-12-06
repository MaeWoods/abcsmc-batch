
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///                                                                                                                                    ///
/// Evolutionary model ///  Model 1: insertions + deletions Model 2: insertions + deletions + translocations                           ///
///                         Model 3: insertions + deletions + selection Model 4: insertions + deletions + translocations + selection   ///
/// ------------------------------------------------------------------------------------------------------------                       ///
///                                                                                                                                    ///
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef abcsmc_H
#define abcsmc_H
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

class abcsmc{

private:

int NChr;
int Npop;
int ngen;
double lowB;
int chrlength;
double Curvemax;
int minchr;
double mu_b;
double minSV;

double margins_prev;
double margins_curr;
vector<double> weights_prev;
vector<double> weights_curr;
vector<vector<double> > parameters;
vector<double> this_parameter;
vector<vector<double> > kernel;
vector<vector<double> > prior;
vector<double> distances;
int Nparticles;
int nparameters;
double epsilon;

public:

	
	// Default constructor
	
	abcsmc(void);
	
	// Destructor
	
	~abcsmc(void);
	
	double runiform(gsl_rng* r, double a, double b);
	void setup(void);
	int read_prior( const string& filename, vector<vector<double>> &prior);
	int set_prevs(int Job, vector<double> &weights_prev, vector<vector<double> > &kernel, double &margins_prev, int oldposition);
	int set_naccepted(int position, vector<double> &distances, vector<vector<double> > &parameters_prev);
	int read_position(int job);
	double getepsilon(int position, int Job);
	int run_ams(double final_epsilon, double alpha, int JobID, int npart, int nparam, int seed, int position);
	void compute_next_epsilon(double target_epsilon, double alpha, int position, int Job);
	void iterate_one_population(int position, double epsilon, int Job, int n_acc, gsl_rng* r, int position_prev, int naccep);
	void PrintOutRes(int Job, int position_prev);
	void modelMarginals();
	void getKernel(vector<vector<double> > &kernel);
	double perturbParticle(vector<double> &params, int thisp, gsl_rng* r);
	void computeParticleWeights(vector<double> &weights_curr);
	double getPdfParameterKernel(vector<double> params0, int j);
	void normalizeWeights(void);
	int sample_particle(gsl_rng* r);
	void sampleTheParameter(vector<double> &this_parameter, gsl_rng* r);
	void sampleTheParameterFromPrior(vector<double> &this_parameter);

};

#endif












