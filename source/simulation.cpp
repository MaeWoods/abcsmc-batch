
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///                                                                                                                                    ///
/// Evolutionary model ///  Model 1: insertions + deletions Model 2: insertions + deletions + translocations                           ///
///                         Model 3: insertions + deletions + selection Model 4: insertions + deletions + translocations + selection   ///
/// ------------------------------------------------------------------------------------------------------------                       ///
///                                                                                                                                    ///
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
#include <cstdlib>

#include "simulation.h"

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

double simulation::runiform(gsl_rng* r, double a, double b){
  double myrandom = a + (b-a)*gsl_rng_uniform (r);
 	
 	while(myrandom==0){
 	myrandom = a + (b-a)*gsl_rng_uniform (r);
 	}
 	return myrandom;
}

simulation::simulation(){

}

simulation::~simulation(){

}
  
void simulation::simulate(vector<double> this_parameters, double epsilon, int position, int seed, int Rsim, gsl_rng* r){

  int Nparticles = 1;
  int mod_choice = 4;
  
  string CFile = "ConstPars";

  string jobid = "";
  
  ////////////////////////////////////////////////////////////
  ///Generating parameters and random seed from bash script///
  ////////////////////////////////////////////////////////////


  //constants
  int NChr, Npop, ngen;
  double lowB, chrlength, Curvemax, minchr, mu_b, minSV, Lchr_loss,Uchr_loss,Lchr_gain,Uchr_gain,Lp_tran,Up_tran,Lmu_ki,Umu_ki,Lmu_kd,Umu_kd,Lfitness,Ufitness,Lmaxchr,Umaxchr,LSV_mean,USV_mean,Lgnmdou,Ugnmdou;
  
  //parameters
  vector<double> All_mu_b, All_chr_loss, All_chr_gain, All_p_tran, All_fitness,  All_SV_mean, All_mu_ki, All_mu_kd, All_gnmdou, All_maxchr, All_minchr;
  unsigned int nConstPars = 0;
  unsigned int nResimPars = 0;
  nConstPars = read_constnts(CFile, NChr, Npop, ngen, lowB, chrlength, Curvemax, minchr, mu_b, minSV,Lchr_loss,Uchr_loss,Lchr_gain,Uchr_gain,Lp_tran,Up_tran,Lmu_ki,Umu_ki,Lmu_kd,Umu_kd,Lfitness,Ufitness,Lmaxchr,Umaxchr,LSV_mean,USV_mean,Lgnmdou,Ugnmdou);
  std::cout << " NChr: " << NChr << " Npop: " << Npop << " ngen: " << ngen << " lowB: " << lowB << " chrlength: " << chrlength << " Curvemax: " << Curvemax << " minchr: " << minchr << " mu_b: " << mu_b << " minSV: " << minSV << std::endl;
  
  ////////////////////////////////////////////////////////////
  ///   Containers                                         ///
  ////////////////////////////////////////////////////////////
  
  	vector<vector<vector<double> > > Mprev;
	vector<vector<vector<double> > > MCprev;
	vector<vector<double> >  Cprev;
	vector<vector<double> >  DivCprev;
	vector<vector<double> >  CMixprev;	

	vector<vector<vector<double> > > NTprev;
	vector<vector<vector<double> > > NIprev;
	vector<vector<vector<double> > > NDprev;

	vector<vector<double> > rdelprev;
	vector<vector<double> > rinsprev;
	vector<vector<double> > rtransprev;
	
	vector<int> GDprev;
	vector<int> n_chr;
	vector<int> n_chrprev;
	vector<double> GDprobprev;
	vector<double> CSize(NChr);

	vector<vector<vector<double> > > M;
	vector<vector<vector<double> > > MC;
	vector<vector<double> >  C;
	vector<vector<double> >  DivC;
	vector<vector<double> >  CMix;

	vector<vector<vector<double> > > NT;
	vector<vector<vector<double> > > NI;
	vector<vector<vector<double> > > ND;

	vector<vector<double> > rdel;
	vector<vector<double> > rins;
	vector<vector<double> > rtrans;
	vector<int> GD;
	vector<double> GDprob;
	
	vector<int> Pcinplus;
	vector<int> Pcinminus;
	
	vector<int> Pcinplusprev;
	vector<int> Pcinminusprev;
       
    //parameters inferred
    double chr_gain = 0;
    double chr_loss = 0;
    double p_tran = 0;     
    double mu_ki = 0;
    double mu_kd = 0;
    double fitness = 0;
    double maxchr = 0; 
    double SV_mean = 0; 
    double gnmdou = 0; 
    
    int Pop_remaining = 1;
	vector<double> AllExp;
	int lenExp = 0;

    cout << "Rsim: " << Rsim << std::endl;
    if( Rsim == 0 ){
      ////////////////////////////////////////////////////////////
      ///   Sample from the prior                              ///
      ////////////////////////////////////////////////////////////
	  cout << "begin of sampling: " << std::endl;

      // probability of chromosome gain and loss
      chr_loss = runiform(r,Lchr_loss,Uchr_loss); //gsl_rng_uniform (r);
      chr_gain = runiform(r,Lchr_gain,Uchr_gain); //gsl_rng_uniform (r);

      // probability of translocation per generation
      p_tran = runiform(r,Lp_tran,Up_tran); //gsl_rng_uniform (r);

      // mutation rate of the chromosome for insertions
      mu_ki = runiform(r,Lmu_ki,Umu_ki); // gsl_rng_uniform (r);
      mu_kd = runiform(r,Lmu_kd,Umu_kd); // gsl_rng_uniform (r);
  
      // clim: the gradient of the fitness functions
      fitness = runiform(r,Lfitness,Ufitness);   //200 + (300-200)*gsl_rng_uniform (r);
      maxchr =  runiform(r,Lmaxchr,Umaxchr);
  
      // max and minimum size of SVs
      SV_mean = runiform(r,LSV_mean,USV_mean);
    
      // probability of genome doubling
      gnmdou = runiform(r,Lgnmdou,Ugnmdou);
      
      cout << "end of sampling: " << std::endl;
      
      }
      
      else{

      chr_gain =      this_parameters[0];
      chr_loss =     this_parameters[1];
      p_tran =      this_parameters[2];
      mu_ki =  this_parameters[3];
      mu_kd =  this_parameters[4];
      fitness = this_parameters[5];
      maxchr =      this_parameters[6];
      SV_mean =      this_parameters[7];
      gnmdou =      this_parameters[8];
      
    }
      
	//////////////////////////////////
	/////////Read in Exp//////////////
	//////////////////////////////////


 	 AllExp.resize(1);
  	string line;
 	 ifstream expfile ("data_Gastric.txt");
 	 if (expfile.is_open())
  	{
		std::cout << "0" << std::endl;
  
  	while (!expfile.eof( ))      //if not at end of file, continue reading numbers
     	{
     		string line;
    	 	getline (expfile,line,'\n');
     

 			double crint1 = stod(line);

   				 
   		if(lenExp==0){
   		AllExp[0] = crint1;
   		}
     	else{
     		AllExp.push_back(crint1);
     	}
     	lenExp += 1;
     	
     	}
     	
     	expfile.close();
     	
}
	
    
    
  
  	n_chr.resize(Npop);
  	n_chrprev.resize(Npop);
    GD.resize(Npop);
    GDprob.resize(Npop);
    GDprev.resize(Npop);
    GDprobprev.resize(Npop);
    M.resize(Npop);
    MC.resize(Npop);
    C.resize(Npop);
    DivC.resize(Npop);
    CMix.resize(Npop);
    
    //step containers
    Mprev.resize(Npop);
    MCprev.resize(Npop);
    Cprev.resize(Npop);
    DivCprev.resize(Npop);
    CMixprev.resize(Npop);
		
    NTprev.resize(Npop);
    NIprev.resize(Npop);
    NDprev.resize(Npop);
	
    rdelprev.resize(Npop);
    rinsprev.resize(Npop);
    rtransprev.resize(Npop);
	
    NT.resize(Npop);
    NI.resize(Npop);
    ND.resize(Npop);
	
    rdel.resize(Npop);
    rins.resize(Npop);
    rtrans.resize(Npop);    
    
    Pcinplus.resize(Npop);
    Pcinplusprev.resize(Npop);
    Pcinminus.resize(Npop);
    Pcinminusprev.resize(Npop);

	cout << "begin of initial: " << std::endl;
    //Initialize
    for (int i = 0; i < Npop; i++) {

	  n_chr[i] = NChr;
	  n_chrprev[i] = NChr;
      C[i].resize(NChr);
      DivC[i].resize(NChr);
      M[i].resize(NChr);
      MC[i].resize(NChr);
      
      CMix[i].resize(NChr);
      rdel[i].resize(NChr);
      rins[i].resize(NChr);
      rtrans[i].resize(NChr);
	
      NT[i].resize(NChr);
      NI[i].resize(NChr);
      ND[i].resize(NChr);
      
      //step containers
      n_chrprev[i] = NChr;
      Cprev[i].resize(NChr);
      CMixprev[i].resize(NChr);
      DivCprev[i].resize(NChr);
      Mprev[i].resize(NChr);
      MCprev[i].resize(NChr);
      rdelprev[i].resize(NChr);
      rinsprev[i].resize(NChr);
      rtransprev[i].resize(NChr);
	
      NTprev[i].resize(NChr);
      NIprev[i].resize(NChr);
      NDprev[i].resize(NChr);
	
      for(int j = 0; j < NChr; j++){
	
	NT[i][j].resize(NChr);
	NI[i][j].resize(NChr);
	ND[i][j].resize(NChr);
	NTprev[i][j].resize(NChr);
	NIprev[i][j].resize(NChr);
	NDprev[i][j].resize(NChr);
	
	M[i][j].resize(NChr);
	MC[i][j].resize(NChr);
	Mprev[i][j].resize(NChr);
	MCprev[i][j].resize(NChr);
	
	if(j==0){
	  C[i][j]=chrlength;
	  CMix[i][j]=chrlength;
	  CSize[j]=chrlength;
	  
	  Cprev[i][j]=C[i][j];
	CMixprev[i][j]=CMix[i][j];
	
	}
	else{
	  C[i][j]=chrlength/(j + j*1.5);
	  CMix[i][j]=chrlength/(j + j*1.5);
	  CSize[j]=chrlength/(j + j*1.5);
	  
	  Cprev[i][j]=C[i][j];
	CMixprev[i][j]=CMix[i][j];
	}
	DivC[i][j]=0;
	DivCprev[i][j]=DivC[i][j];
    
	for(int k = 0; k < NChr; k++){
	
		NT[i][j][k]=0;
		NI[i][j][k]=0;
		ND[i][j][k]=0;
		NTprev[i][j][k]=NT[i][j][k];
	  NIprev[i][j][k]=NI[i][j][k];
	  NDprev[i][j][k]=ND[i][j][k];
	  if(j==k){
	    M[i][j][k]=1;
	    Mprev[i][j][k]=M[i][j][k];
	    if(j==0){
	      MC[i][j][k]=chrlength;
	      MCprev[i][j][k]=MC[i][j][k];
	    }
	    else{
	      MC[i][j][k]=chrlength/(j + j*1.5);
	      MCprev[i][j][k]=MC[i][j][k];
	    }
	  }
	  else{
	    M[i][j][k]=0;
	    MC[i][j][k]=0;
		
	  Mprev[i][j][k]=M[i][j][k];
	  MCprev[i][j][k]=MC[i][j][k];
	  
	  }
	}
      }
    
    }
	
    for(int i=0; i<Npop; i++){
      GDprob[i] = gnmdou*(  runiform(r,0,1) );
      Pcinplus[i] = chr_gain*(  runiform(r,0,1) );
      Pcinplusprev[i] = Pcinplus[i];
      Pcinminus[i] = chr_loss*(  runiform(r,0,1) );
      Pcinminusprev[i] = Pcinminus[i];
      GDprobprev[i] = GDprob[i];
      GD[i] = 1;
      GDprev[i] = 1;
	
      for(int j = 0; j < NChr; j++){

	double qqi = ( runiform(r,0,1) );
	double p_ins = - log(qqi)*mu_ki;
	double qqd = ( runiform(r,0,1) );
	double p_del = - log(qqd)*mu_kd;
 
	  rdel[i][j] = mu_b*p_del;
	  rins[i][j] = mu_b*p_ins;	
	  rdelprev[i][j] = rdel[i][j];
	  rinsprev[i][j] = rins[i][j];
	  rtrans[i][j]= mu_b*p_tran*(   runiform(r,0,1) );
	  rtransprev[i][j] = rtrans[i][j];
	  
      }
      
    }
     	
cout << "end of initialization: " << std::endl;
    ////////////////////////////////////////////////////////////
    ///   Loop over the generations                          ///
    ////////////////////////////////////////////////////////////
    //std::cout << "Looping over generations" << std::endl;
    cout << "start of generations: " << std::endl;
    cout << "epsilon: " << epsilon << std::endl;

     for(int gg=1; gg<ngen; gg++){
      std::cout << "gg: " << gg << std::endl;
      //////////////////////////////////////////////////////
      /////////Mutate: insertions and deletions/////////////
      //////////////////////////////////////////////////////
	
    //elements for selection
    //std::cout << "\tSelection" << std::endl;
      vector<int> keep;
	
      keep.resize(Npop);
	
      int n_remaining = 0;
      int n_checklen = 0;

	//Genome doubling
    for(int i=0; i<Npop; i++){	
	double r_gd = runiform(r,0,1);
	if(r_gd<GDprob[i]){
		
	  //Insert vectors		
	  vector<vector<double> > MCins;
	  vector<vector<double> > Mins;
	  vector<double> CMixins;
	  MCins.resize(n_chrprev[i]);
	  Mins.resize(n_chrprev[i]);
	  CMixins.resize(n_chrprev[i]);
				
	  for(int j=0; j<n_chrprev[i]; j++){
	    MCins[j].resize(NChr);
	    Mins[j].resize(NChr);
				
	    for(int k=0; k<NChr; k++){
				
	      MCins[j][k]=MCprev[i][j][k];
	      Mins[j][k]=Mprev[i][j][k];
	      CMixins[j]=CMixprev[i][j];
				
	    }
	  }
			
	  MCprev[i].reserve(MCprev[i].size() + MCins.size());
	  MCprev[i].insert(MCprev[i].end(), MCins.begin(), MCins.end());
	  Mprev[i].reserve(Mprev[i].size() + Mins.size());
	  Mprev[i].insert(Mprev[i].end(), Mins.begin(), Mins.end());
	    
	  CMixprev[i].reserve(CMixprev[i].size() + CMixins.size());
	  CMixprev[i].insert(CMixprev[i].end(), CMixins.begin(), CMixins.end());
	  DivCprev[i][0] += Cprev[i][0];
	  DivCprev[i][1] += Cprev[i][1];
	  DivCprev[i][2] += Cprev[i][2];
				

	  for(int m=0; m<n_chrprev[i]; m++){
	  	for(int k=0; k<NChr; k++){
	    	NIprev[i][m][k] = NIprev[i][m][k] + NIprev[i][m][k];
	    	}
	  }

	  GDprev[i] += 1;
	  n_chr[i] = n_chrprev[i]*2;
	}
	
	double r_cd = runiform(r,0,1);
	double r_cl = runiform(r,0,1);
	
		//chromosome gain
		if(r_cd<Pcinplus[i]){

		//sample a chromosome
		int sample_chr = ceil(n_chrprev[i]*runiform(r,0,1))-1;

				//Insert vectors
		  vector<vector<double> > MCins;
		  vector<vector<double> > Mins;
		  vector<double> CMixins;
		  MCins.resize(NChr);
		  Mins.resize(NChr);
		  CMixins.resize(1);

				for(int k=0; k<NChr; k++){

				  MCins[0][k]=MCprev[i][sample_chr][k];
		  Mins[0][k]=Mprev[i][sample_chr][k];
		  CMixins[0]=CMixprev[i][sample_chr];

				}

		  MCprev[i].reserve(MCprev[i].size() + MCins.size());
		  MCprev[i].insert(MCprev[i].end(), MCins.begin(), MCins.end());
		  Mprev[i].reserve(Mprev[i].size() + Mins.size());
		  Mprev[i].insert(Mprev[i].end(), Mins.begin(), Mins.end());

		  CMixprev[i].reserve(CMixprev[i].size() + CMixins.size());
		  CMixprev[i].insert(CMixprev[i].end(), CMixins.begin(), CMixins.end());
		  DivCprev[i][0] += MCprev[i][sample_chr][0];
		  DivCprev[i][1] += MCprev[i][sample_chr][1];
		  DivCprev[i][2] += MCprev[i][sample_chr][2];

		  n_chrprev[i] = n_chrprev[i] + 1;

		  for(int m=0; m<NChr; m++){
				NIprev[i][sample_chr][m] = NIprev[i][sample_chr][m] + NIprev[i][sample_chr][m];
		  }

		}
		
		//chromosome loss
		if(r_cl<Pcinminus[i]){
		
		  //sample a chromosome
		  int sample_chr = ceil(n_chrprev[i]*runiform(r,0,1))-1;
		
		  DivCprev[i][0] -= MCprev[i][sample_chr][0];
		  DivCprev[i][1] -= MCprev[i][sample_chr][1];
		  DivCprev[i][2] -= MCprev[i][sample_chr][2];	
		  MCprev[i].erase(MCprev[i].begin() + sample_chr);	  
		  Mprev[i].erase(Mprev[i].begin() + sample_chr);
		  NIprev[i].erase(NIprev[i].begin() + sample_chr);
		  NDprev[i].erase(NDprev[i].begin() + sample_chr);
		  NTprev[i].erase(NTprev[i].begin() + sample_chr);
		  CMixprev[i].erase(CMixprev[i].begin() + sample_chr);	
		  n_chrprev[i] = n_chrprev[i] - 1;
	
		}
      
		//////////////////////////////////////////////////////
		/////////Mutate: insertions and deletions/////////////
		//////////////////////////////////////////////////////
		for(int g_d=0; g_d<n_chrprev[i]; g_d++){
			
			for(int k=0; k<NChr; k++){
		
			  int N_i = gsl_ran_poisson (r, MCprev[i][g_d][k]*rinsprev[i][k]);
			  int N_d = gsl_ran_poisson (r, MCprev[i][g_d][k]*rdelprev[i][k]);
					
			  for(int tot=0; tot<N_i; tot++){
				
			double qq = (   runiform(r,0,1) ); 
			double sz = - log(qq)*SV_mean; 
			if(sz<minSV){
			}
			else{
			DivCprev[i][k] += sz;
			Cprev[i][k] = Cprev[i][k] + sz;
			MCprev[i][g_d][k] = MCprev[i][g_d][k] + sz;
			NIprev[i][g_d][k] = NIprev[i][g_d][k] + 1;
			}
		  
			  }
				
			  for(int tot=0; tot<N_d; tot++){

			double qq = (   runiform(r,0,1) ); 
			double sz = - log(qq)*SV_mean; 
			if((MCprev[i][g_d][k]>=0)&&(sz>MCprev[i][g_d][k])){
		
			sz = MCprev[i][g_d][k];
		
			}
			if(sz<minSV){
			}
			else{
			DivCprev[i][k] -= sz;
			Cprev[i][k] = Cprev[i][k] - sz;
			MCprev[i][g_d][k] = MCprev[i][g_d][k] - sz;
			NDprev[i][g_d][k] = NDprev[i][g_d][k] + 1;
			}
		
			  }	      
		
			}	
		
		}
	
		/////////////////////////////
		/////////UpdateM/////////////
		/////////////////////////////
		for(int g_d=0; g_d<n_chrprev[i]; g_d++){
		
		  for(int j=0; j<NChr; j++){
			
			double sum_j = 0;
			
			for(int k=0; k<NChr; k++){
		
			double cprev_size = 0;
		
			cprev_size = Cprev[i][k];
		
			  //Avoid seg fault
			  if(cprev_size==0){
			  Mprev[i][g_d][k] = 0;
			  }
			  else{
			  Mprev[i][g_d][k] = MCprev[i][g_d][k]/cprev_size;
			  }
				
			  sum_j = Mprev[i][g_d][k]*Cprev[i][k] + sum_j;
							
			}
			
			CMixprev[i][g_d]=sum_j;

		  }
	  
		}
		
		//////////////////////////////////////////////////////
		/////////Mutate: translocations///////////////////////
		//////////////////////////////////////////////////////
		//std::cout << "\tTranslocations" << std::endl;
		  if((mod_choice==2)||(mod_choice==4)){
	  
	
		int jtick = 0;
		int jcounte = 0;
		for(int j=0; j<n_chrprev[i]; j++){
		
		  if(jtick==NChr){
			jcounte = jcounte + 1;
			jtick = 0;
		  }
			
		  for(int k=0; k<NChr; k++){
		  //translocate between C1 and C2//
			  int N_t = gsl_ran_poisson (r, MCprev[i][j][k]*rtransprev[i][k]);
			
			if( N_t > 0 ){
			
			  ///Mutate: second chomsome///
			  int j1tick = 0;
			  int j1counte = 0;
			  for(int j1=0; j1<n_chrprev[i]; j1++){
		  
			if(j1tick==NChr){
			  j1counte = j1counte + 1;
			  j1tick = 0;
			}
					
			if(j!=j1){
			
			  for(int k1=0; k1<NChr; k1++){
					
				if(k!=k1){
			  
				  int N_t1 = gsl_ran_poisson (r, MCprev[i][j1][k1]*rtransprev[i][k1]);
						
				  if( N_t1 > 0 ){
			
				double sum_j0 = 0;
				double sum_j1 = 0;
			
				for(int ka=0; ka<NChr; ka++){

				  double p = lowB + (1-lowB)*runiform(r,0,1);
				  double a1 = Mprev[i][j][ka];
				  double b1 = Mprev[i][j1][ka];
				  Mprev[i][j][ka] = p*(a1+b1);
				  Mprev[i][j1][ka] = (1-p)*(a1+b1);
				  sum_j0 += Mprev[i][j][ka]*Cprev[i][ka];
				  sum_j1 += Mprev[i][j1][ka]*Cprev[i][ka];
				
			
				}
			
				CMixprev[i][j]=sum_j0;
				CMixprev[i][j1]=sum_j1;
			
				NTprev[i][j][k] = NTprev[i][j][k] + 1;
				NTprev[i][j1][k1] = NTprev[i][j1][k1] + 1;
			
				  }
			
				}
			
			  }
					
			}
			
			
			j1tick = j1tick + 1;
			
			  }
			
			}
			
		  }
			
		  jtick = jtick + 1;
			
		}
	  
		  }
		  else{
	  
		  }
		  
		  
		//////////////////////////////////////////////////////
		/////////Update M & MC ///////////////////////////////
		//////////////////////////////////////////////////////
		
		for(int g_d=0; g_d<n_chrprev[i]; g_d++){
		
		  for(int j=0; j<NChr; j++){
			
			double sum_j = 0;
			
			for(int k=0; k<NChr; k++){
				
			  MCprev[i][g_d][k] = Mprev[i][g_d][k]*Cprev[i][k];
				
			  sum_j = Mprev[i][g_d][k]*Cprev[i][k] + sum_j;
							
			}
			
			CMixprev[i][g_d]=sum_j;
		
		  }
	  
		}
	
      
	

      //////////////////////////////////////////////////////
      /////////Selection////////////////////////////////////
      //////////////////////////////////////////////////////
	
     
	
	int check_len = 0;
	int check_lendiv = 0;
    keep[i] = 0;
		
	for(int g_d=0; g_d<n_chrprev[i]; g_d++){
			
	    double p1 = 1 - (Curvemax/(1+exp(-1*fitness*((CMixprev[i][g_d])-minchr))));
	    double p2 = (Curvemax/(1+exp(-1*fitness*((CMixprev[i][g_d])-maxchr))));
	      
	    double r1 =(   runiform(r,0,1) );
	    double r2 =(   runiform(r,0,1) );
	    
	    if((mod_choice==3)||(mod_choice==4)){
	    if(r1<p1){
			
	      check_len = 1;
	    }
	      
	    if(r2<p2){
	      check_len = 1;
	      
	    }
	    }
	    else{
	    check_len = 0;
	    }
	    
	    if(CMixprev[i][g_d]>100*maxchr){
	    
	    check_lendiv = 1;
	    
	    }
			
				
	  
			
	}
	  
	if(check_len==0){
		
	  keep[i] = 1;
	  n_remaining += 1;
	}
	if(check_lendiv==1){
	 n_checklen = 1;
	}
	
      }
      
      if((n_remaining==0)||(n_checklen==1)){
      
      ngen=gg;
      Pop_remaining = 0;

    //////////////////////////////////////////////////////
    /////////Print out parameters       //////////////////
    //////////////////////////////////////////////////////
    if(Rsim == 0){
      char Dat[100];
      char Dirr[100];
      int Dn;
      int Dr;
		Dr=sprintf(Dirr,"mkdir -p %d/err_params",position);

	  const int dir= system(Dirr);
      Dn=sprintf(Dat,"%d/err_params/%d.dat",position,seed);
   
      fstream myfileDatParam; 
      myfileDatParam.open(Dat,ios::out);

      if(myfileDatParam.is_open()){
	
	myfileDatParam << chr_loss << '\t';
	myfileDatParam << p_tran << '\t';
	myfileDatParam << fitness << '\t';
	myfileDatParam << SV_mean << '\t';
	myfileDatParam << mu_ki << '\t';
	myfileDatParam << mu_kd << '\t';
	myfileDatParam << gnmdou << '\t';
	myfileDatParam << maxchr << std::endl;
	
	myfileDatParam.close();
      }
      else{
	std::cerr << "Error: open of output parameter file unsuccessful: " << Dat << std::endl;

      }
    }
    
      
      
      }
      else{
	
      vector<int> sample_vec;
      sample_vec.resize(n_remaining);
      int counter = 0;
      for(int q=0; q<Npop; q++){
	
	if(keep[q]==1){
	
	  sample_vec[counter] = q;
	  counter += 1; 
	
	}
	
      }

      //////////////////////////////////////////////////////
      /////////Resample/////////////////////////////////////
      //////////////////////////////////////////////////////
      //std::cout << "\tResample" << std::endl;
      for (int i = 0; i < Npop; i++) {
	int newgd = 0;

	if(keep[i]==1){
	

	
	  if(n_chr[i]!=n_chrprev[i]){

	    MC[i].resize(MCprev[i].size());
	    M[i].resize(Mprev[i].size());
	    CMix[i].resize(CMixprev[i].size());

	    for(int j=0; j<M[i].size(); j++){
	
	      MC[i][j].resize(NChr);
	      M[i][j].resize(NChr);
	
	    }
	
	  }
	
	    for(int j = 0; j < NChr; j++){

	      C[i][j]=Cprev[i][j];
	      DivC[i][j]=DivCprev[i][j];
    
	      rdel[i][j] = rdelprev[i][j];
	      rins[i][j] = rinsprev[i][j];
	      rtrans[i][j] = rtransprev[i][j];
	      
	      for(int g_d=0; g_d<n_chrprev[i]; g_d++){
	      
	      CMix[i][g_d]=CMixprev[i][g_d];
    
	      NT[i][g_d][j] = NTprev[i][g_d][j];
	      NI[i][g_d][j] = NIprev[i][g_d][j];
	      ND[i][g_d][j] = NDprev[i][g_d][j];
	      
	      M[i][g_d][j]=Mprev[i][g_d][j];
		  MC[i][g_d][j]=MCprev[i][g_d][j];
	
	    }
	  }
	  newgd = GDprev[i];
	  GDprob[i] = GDprobprev[i];
	  Pcinplus[i] = Pcinplusprev[i];
	  Pcinminus[i] = Pcinminusprev[i];
	
	}
	else{

	  double rnew = runiform(r,0,1);
	  int position = floor(rnew*n_remaining);
	
	  if(n_chr[i]!=n_chrprev[sample_vec[position]]){
	
		MC[i].clear();
		M[i].clear();
		CMix[i].clear();
	    MC[i].resize(MCprev[sample_vec[position]].size());
	    M[i].resize(Mprev[sample_vec[position]].size());
	    CMix[i].resize(CMixprev[sample_vec[position]].size());
	
	    for(int j=0; j<M[i].size(); j++){
	
	      MC[i][j].resize(NChr);
	      M[i][j].resize(NChr);
	
	    }

	  }
	
	  for(int g_d=0; g_d<n_chrprev[sample_vec[position]]; g_d++){
	
	    for(int j = 0; j < NChr; j++){

	      C[i][j]=Cprev[sample_vec[position]][j];
	      DivC[i][j]=DivCprev[sample_vec[position]][j];
	      CMix[i][g_d]=CMixprev[sample_vec[position]][g_d];
    
	      rdel[i][j] = rdelprev[sample_vec[position]][j];
	      rins[i][j] = rinsprev[sample_vec[position]][j];
	      rtrans[i][j] = rtransprev[sample_vec[position]][j];
    
	      NT[i][g_d][j] = NTprev[sample_vec[position]][g_d][j];
	      NI[i][g_d][j] = NIprev[sample_vec[position]][g_d][j];
	      ND[i][g_d][j] = NDprev[sample_vec[position]][g_d][j];
    
	      M[i][g_d][j]=Mprev[sample_vec[position]][g_d][j];
		  MC[i][g_d][j]=MCprev[sample_vec[position]][g_d][j];
	
	
	    }

	  }


	  newgd = GDprev[sample_vec[position]];
	  GDprob[i] = GDprobprev[sample_vec[position]];
	  Pcinminus[i] = Pcinminusprev[sample_vec[position]];

	}
	
	GD[i] = newgd;
	
      }
	
      //////////////////////////////////////////////////////
      /////////Update vectors///////////////////////////////
      //////////////////////////////////////////////////////

      for (int i = 0; i < Npop; i++) {
	
	if(n_chrprev[i]!=n_chr[i]){

	  MCprev[i].resize(MC[i].size());
				
	  Mprev[i].resize(M[i].size());
			
	  CMixprev[i].resize(CMix[i].size());
	
	  for(int j=0; j<Mprev[i].size(); j++){
	
	    MCprev[i][j].resize(NChr);
	    Mprev[i][j].resize(NChr);
	
	  }	
	
	}
	
	for(int g_d=0; g_d<n_chrprev[i]; g_d++){

	  for(int j = 0; j < NChr; j++){	
	
	    Cprev[i][j]=C[i][j];
   
	    for(int k = 0; k < NChr; k++){   
     
	      Mprev[i][g_d][k]=M[i][g_d][k];
	      MCprev[i][g_d][k]=MC[i][g_d][k];
		
	    }

	    DivCprev[i][j]=DivC[i][j];
	    CMixprev[i][g_d]=CMix[i][g_d];
    	rdelprev[i][j] = rdel[i][j];
	    rinsprev[i][j] = rins[i][j];
	    rtransprev[i][j] = rtrans[i][j];
    
	    NTprev[i][g_d][j] = NT[i][g_d][j];
	    NIprev[i][g_d][j] = NI[i][g_d][j];
	    NDprev[i][g_d][j] = ND[i][g_d][j];
	
	  }
	
	}
	GDprev[i] = GD[i];
	Pcinminusprev[i] = Pcinminus[i];
	GDprobprev[i] = GDprob[i];
	
      }
      
      }
	
    } // Loop over generations

	if(Pop_remaining==0){
	
	}
	else{
      
    //////////////////////////////////////////////////////
    /////////Print out simulation results ////////////////
    //////////////////////////////////////////////////////
    vector<double > pertg(Npop*NChr);
    int gad = 0;
    for (int i = 0; i < Npop; i++) {
	for(int j = 0; j < NChr; j++){
	  
	  double log_signed_div = 0;
	
	  if(DivCprev[i][j]==0){
	    log_signed_div = 0;
	  }
	  else if(DivCprev[i][j]<0){
	    log_signed_div = -1*abs(DivCprev[i][j]);
	
	  }
	  else{
	    log_signed_div = DivCprev[i][j];
	  }
    
    pertg[gad] = log_signed_div/CSize[j];
    gad += 1;
    }
    }
     
     vector<double> Allmodsort;
     vector<double> AllExpsort;
    sort (pertg.begin(), pertg.begin() + gad); 	
    sort (AllExp.begin(), AllExp.begin() + lenExp);
    double sorted_dataexp[lenExp];
    double sorted_datamod[gad];
    size_t stride = 1;
    for(int k=0; k<lenExp; k++){
    
    sorted_dataexp[k] = AllExp[k];
    }
    for(int k=0; k<gad; k++){
    
    sorted_datamod[k] = pertg[k];
    }
    double del = 0;
    double myprob = 0;
    for(int b=0; b<90; b++){
    
    myprob = 0.05 + b*0.01;
    
    double qdatexp = gsl_stats_quantile_from_sorted_data (sorted_dataexp, stride, lenExp, myprob);
    
    double qdatmod = gsl_stats_quantile_from_sorted_data (sorted_datamod, stride, gad, myprob);
    del += (qdatexp-qdatmod)*(qdatexp-qdatmod);
      }
      
     if(del<epsilon){


char Datins[100];
    int Dnins;
      Dnins=sprintf(Datins,"results-resim-M%d-P%d-B%d-J%d.dat",position,seed,position,seed);
    fstream myfileIns; 
    myfileIns.open(Datins,ios::out);
 	
    if(myfileIns.is_open()){
 	
      myfileIns << "Num" << ", " << "CNum" <<  ", " << "CSize" ", " << "PDiv" << ", " << "CDiv" << ", " << "RelLen" ", " << "NTrans" << ", " << "NIns" << ", " << "NDel" << ", " << "Rins" << ", " << "Rdel" << ", " << "Rtrans" << '\n';
      
      for (int i = 0; i < Npop; i++) {
	for(int j = 0; j < NChr; j++){

	int NTpl = 0;
	int NIpl = 0;
	int NDpl = 0;
	
	for(int g_d=0; g_d<n_chrprev[i]; g_d++){
	
	NTpl += NT[i][g_d][j];
	NIpl += NI[i][g_d][j];
	NDpl += ND[i][g_d][j];

	}
	  
	  double log_signed_div = 0;
	
	  if(DivCprev[i][j]==0){
	    log_signed_div = 0;
	  }
	  else if(DivCprev[i][j]<0){
	    log_signed_div = -1*log10(abs(DivCprev[i][j]));
	
	  }
	  else{
	    log_signed_div = log10(DivCprev[i][j]);
	  }
    
	  myfileIns << i << ", "  << j << ", " << CSize[j] << ", " << log_signed_div/CSize[j] << ", " << log_signed_div << ", " << DivCprev[i][j]/CSize[j] << ", " << NTpl << ", " << NIpl << ", " << NDpl << ", " << rins[i][j] << ", " << rdel[i][j] << ", " << rtrans[i][j] << '\n';
	 
	 if(DivCprev[i][j]/CSize[j]<-2){
	 std::cout << "ERROR!!!  gd: " << GDprev[i] << " i: " << i << " j: " << j << " how much " << DivCprev[i][j]/CSize[j] << std::endl;
	 }
	 
	}
      }
      myfileIns.close();

    }
     else{
       std::cerr << "error: open of output results file unsuccessful: " << Datins << std::endl;
    }

    

	std::cout << "position: " << position << std::endl;

      char Dat[100];
      char Dirr[100];
      int Dn;
      int Dr;
      Dr=sprintf(Dirr,"mkdir -p %d/params_prev",position);
	const int dir= system(Dirr);


      Dn=sprintf(Dat,"%d/params_prev/%d.dat",position,seed);
   
      fstream myfileDatParam; 
      myfileDatParam.open(Dat,ios::out);

      if(myfileDatParam.is_open()){
	
	myfileDatParam << chr_loss << '\t';
	myfileDatParam << chr_gain << '\t';
	myfileDatParam << p_tran << '\t';
	myfileDatParam << mu_ki << '\t';
	myfileDatParam << mu_kd << '\t';
	myfileDatParam << fitness << '\t';
	myfileDatParam << maxchr << '\t';
	myfileDatParam << SV_mean << '\t';
	myfileDatParam << gnmdou << std::endl;
	myfileDatParam << del << std::endl;
	
	myfileDatParam.close();
      }
      else{
	std::cerr << "Error: open of output parameter file unsuccessful: " << Dat << std::endl;
	
      }
    
    
    }
    }
    
    ///////////////////////////
    // Clear containers ///////
    ///////////////////////////
    
    Mprev.clear();
	MCprev.clear();
	Cprev.clear();
	DivCprev.clear();
	CMixprev.clear();	

	NTprev.clear();
	NIprev.clear();
	NDprev.clear();

	rdelprev.clear();
	rinsprev.clear();
	rtransprev.clear();
	GDprev.clear();
	GDprobprev.clear();
	CSize.clear();
	Pcinplusprev.clear();
	Pcinminusprev.clear();

	M.clear();
	MC.clear();
	C.clear();
	DivC.clear();
	CMix.clear();

	NT.clear();
	NI.clear();
	ND.clear();

	rdel.clear();
	rins.clear();
	rtrans.clear();
	GD.clear();
	GDprob.clear();
	Pcinplus.clear();
	Pcinminus.clear();
      

  




  
}

// order of parameters expected:  NChr, Npop, ngen, lowB, chrlength, Curvemax, minchr, mu_b
// order of parameters expected:  NChr, Npop, ngen, lowB, chrlength, Curvemax, minchr, mu_b
int simulation::read_constnts( const string& filename, int &NChr, int &Npop, int &ngen, double &lowB, double &chrlength, double &Curvemax, double &minchr, double &mu_b, double &minSV, double &Lchr_loss, double &Uchr_loss, double &Lchr_gain, double &Uchr_gain, double &Lp_tran, double &Up_tran, double &Lmu_ki, double &Umu_ki, double &Lmu_kd, double &Umu_kd, double &Lfitness, double &Ufitness, double &Lmaxchr, double &Umaxchr, double &LSV_mean, double &USV_mean, double &Lgnmdou, double &Ugnmdou ){

  ifstream infile (filename.c_str());
  int counter = 0; 
  if (infile.is_open()){

    std::string line;
    while(!getline(infile,line).eof()){
      if(line.empty()) continue;

		std::vector<std::string> split;
      std::string buf; 
      stringstream ss(line);
      while (ss >> buf) split.push_back(buf);

		if(counter==0){
      
    
     NChr = atoi( split[0].c_str() ) ;
     Npop = atoi( split[1].c_str() ) ;
     ngen = atoi( split[2].c_str() ) ;
     lowB  = atof( split[3].c_str() ) ;
     chrlength = atof( split[4].c_str() ) ;
     Curvemax =  atof( split[5].c_str() ) ;
     minchr = atof( split[6].c_str() ) ;
     mu_b = atof( split[7].c_str() ) ;
     minSV = atof( split[8].c_str() ) ;
     
     }
     else{
      
      std::cout << "Line: " << split[0] << ' ' << split[1] << std::endl;
      Lchr_loss = atof( split[0].c_str() ) ; 
      Uchr_loss = atof( split[1].c_str() ) ; 
      Lchr_gain = atof( split[2].c_str() ) ;  
      Uchr_gain = atof( split[3].c_str() ) ;  
      Lp_tran = atof( split[4].c_str() ) ;  
      Up_tran = atof( split[5].c_str() ) ;  
      Lmu_ki = atof( split[6].c_str() ) ;  
      Umu_ki = atof( split[7].c_str() ) ;  
      Lmu_kd = atof( split[8].c_str() ) ;  
      Umu_kd = atof( split[9].c_str() ) ;  
      Lfitness = atof( split[10].c_str() ) ;  
      Ufitness = atof( split[11].c_str() ) ;  
      Lmaxchr = atof( split[12].c_str() ) ;  
      Umaxchr = atof( split[13].c_str() ) ;  
      LSV_mean = atof( split[14].c_str() ) ;  
      USV_mean = atof( split[15].c_str() ) ;  
      Lgnmdou = atof( split[16].c_str() ) ;  
      Ugnmdou = atof( split[17].c_str() ) ; 
     
     }
    
      counter++;
    }
      
  }else{
    std::cerr << "Error: open of constant parameters file unsuccessful: " <<  filename << std::endl;
  }
  
  return counter;
}














