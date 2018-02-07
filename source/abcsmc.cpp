
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
#include <dirent.h>

#include "simulation.h"
#include "abcsmc.h"


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

abcsmc::abcsmc(){

double margins_prev;
double margins_curr;
vector<double> weights_prev;
vector<double> weights_curr;
vector<double> distances;
vector<vector<double> > parameters;
vector<double> this_parameter;
vector<vector<double> > kernel;
vector<vector<double> > prior;
int Nparticles;
int nparameters;
	
}

abcsmc::~abcsmc(){
	
}

int abcsmc::read_prior( const string& filename, vector<vector<double>> &prior){

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
     
     }
     else{
      
      //std::cout << "Line: " << split[0] << ' ' << split[1] << std::endl;
      for(int k=0; k<nparameters; k++){
      prior[k][0] = atof( split[2*k].c_str() ) ; 
      prior[k][1] = atof( split[(2*k)+1].c_str() ) ; 
      std::cout << "k: " << k << " prior[k][0] " << prior[k][0] << std::endl;
      std::cout << "k: " << k << " prior[k][1] " << prior[k][1] << std::endl;
      }
     
     }
    
      counter++;
    }
      
  }else{
    std::cerr << "Error: open of prior parameters file unsuccessful: " <<  filename << std::endl;
  }
  
  return counter;
}


double abcsmc::runiform(gsl_rng* r, double a, double b){
  double myrandom = a + (b-a)*gsl_rng_uniform (r);
 	
 	while(myrandom==0){
 	myrandom = a + (b-a)*gsl_rng_uniform (r);
 	}
 	return myrandom;
}

int abcsmc::set_prevs(int Job, vector<double> &weights_prev, vector<vector<double> > &kernel, double &margins_prev, int oldposition){

	char Dat[100];
    int Dn;
    
    std::cout << "old position: " << oldposition << " Job: " << Job << std::endl;
	  //this file contains, position and epsilon
    Dn=sprintf(Dat,"%d/vars_curr-J%d.dat",oldposition, Job);
	std::string path = Dat;
	
	ifstream infile (path.c_str());
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
     //Npop
     margins_prev = atof( split[0].c_str() );
     //Curvemax
     
     for(int k=1; k<(1+Nparticles); k++){
	
		weights_prev[k-1] =  atof( split[k].c_str() );
	}
     
     }
     
     else{
     std::cout << "before kernel" << std::endl;
     	for(int k=0; k<nparameters; k++){
	 kernel[k][0]  =  atof( split[2*k].c_str() );
	kernel[k][1]  =  atof( split[2*k+1].c_str() );
	}
     
     }
    
      counter++;
    }
      
  }
  else{
    std::cerr << "Error: open of constant parameters file unsuccessful: " <<  path.c_str() << std::endl;
  }
}

int abcsmc::set_naccepted(int position, vector<double> &distances, vector<vector<double> > &parameters_prev){

	char Dat[100];
    int Dn;

	  //this file contains, position and epsilon
    Dn=sprintf(Dat,"%d/params_curr/",position);
	std::string path = Dat;
	std::cout << "path " << path << std::endl;
	
	DIR*    dir;
    dirent* pdir;
    std::vector<std::string> files;

    dir = opendir(path.c_str());

	int filecounter = 0;

    while ((pdir = readdir(dir))&&(filecounter<Nparticles+4)) {
        files.push_back(pdir->d_name);
std::cout << "path2 " << files[files.size()-1] << std::endl;        
std::cout << "path3 " << files.size()-1 << std::endl;

if(files[files.size()-1][(files[files.size()-1]).size()-1]=='t'){
	std::string ns = path + files[files.size()-1];

	char *cstr = new char[ns.length() + 1];
strcpy(cstr, ns.c_str());

        ifstream infile (cstr);

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
		if(filecounter<Nparticles){
      for(int k=0; k<nparameters; k++){
     parameters[filecounter][k] = atof( split[k].c_str() ) ;
     
     }
     }
     }
     else if(counter==1){
      
	std::cout << "Line1: " << std::endl;
      std::cout << "Line2: " << split[0] << std::endl;
      std::cout << "filecounter: " << filecounter << std::endl;
      if(filecounter<Nparticles){
      distances[filecounter] = atof( split[0].c_str() ) ; 
      }
     
     }
     else{}
    
      counter++;
    }
      
  }
  else{
    std::cerr << "Error: open of constant parameters file unsuccessful: " << std::endl;
  }
    
    delete [] cstr;
    filecounter += 1;

}
else{

}

    }
    
    int naccepted = files.size()-1;
    
    return(naccepted);

}

int abcsmc::read_position(int job){

int retpos;

char Dat[100];
    	int Dn;

	  //this file contains, position and epsilon
      Dn=sprintf(Dat,"position%d.dat",job);

ifstream infile (Dat);
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
			 //NChr
			 retpos = atoi( split[0].c_str() );
	     }
   
      counter++;
    }
      
  }
  else{
    std::cerr << "Error: open of constant position file unsuccessful: " <<  Dat << std::endl;
  }
return(retpos);

}

void abcsmc::setup(void){


weights_prev.resize(Nparticles);
weights_curr.resize(Nparticles);
parameters.resize(Nparticles);
distances.resize(Nparticles);
for(int k=0; k<Nparticles; k++){

parameters[k].resize(nparameters);

}
this_parameter.resize(nparameters);
kernel.resize(nparameters);
prior.resize(nparameters);
for(int k=0; k<nparameters; k++){

kernel[k].resize(2);
prior[k].resize(2);

}

}

double abcsmc::getepsilon(int position, int Job){

double retep;

char Dat[100];
    int Dn;

	  //this file contains, position and epsilon
    Dn=sprintf(Dat,"%d/epsilon_curr%d.dat",position, Job);
	std::string filename = Dat;

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
    		 retep = atof( split[0].c_str() ) ;
     	}
    
      counter++;
    }
      
  }
  else{
    std::cerr << "Error: open of epsilon file unsuccessful: " <<  filename << std::endl;
  }

return retep;

}

int abcsmc::run_ams(double final_epsilon, double alpha, int JobID, int npart, int nparam, int seed, int position){

		// initialise rng
  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  if( seed != 0 ){
    gsl_rng_set(r, seed );
  }else{
    gsl_rng_set (r, time(NULL) );
  }

int position_prev = 0;

		Nparticles = npart;
		nparameters = nparam;
		int naccepted = 0;
		setup();
		string ConstFile = "ConstPars";
		read_prior( ConstFile, prior);
		
		if(position==0){
		epsilon = 10;
		}
		else{

		position = read_position(JobID);
        
        naccepted = set_naccepted(position, distances, parameters);
		std::cout << "naccepted: " << naccepted << " Nparticles " << Nparticles << std::endl;
        
        position_prev = position;
        
        if(naccepted>Nparticles){
         position += 1;
        }
	else{
        set_naccepted(position-1, distances, parameters);
        }

        compute_next_epsilon(final_epsilon, alpha, position, JobID);
        epsilon = getepsilon(position, JobID);
        std::cout << "Got epsilon: " << epsilon << std::endl;
        
        }
        
		iterate_one_population(position, epsilon, JobID, naccepted,r, position_prev, naccepted);
		
	gsl_rng_free (r);

        return(0);
}

void abcsmc::compute_next_epsilon(double target_epsilon, double alpha, int position, int Job){

        vector<double> distance_values(1);
		distance_values[0] = distances[0];
        for(int i=1; i<Nparticles; i++){
                distance_values.push_back(distances[i]);
		std::cout << "Read in distance_values: " << distance_values[i] << std::endl;
		std::cout << "distances: " << distances[i] << std::endl;        

}
		
		std::cout << "Copied distance values: " << distance_values[0] << std::endl;
        // Important to remember that the initial sort on distance is done on the first distance value
        sort(distance_values.begin(), distance_values.begin()+distance_values.size());
        int ntar = floor( alpha * Nparticles );

        double new_epsilon = distance_values[ntar];
        double ret_epsilon = 0;
        
        // Set the next epsilon
            if(new_epsilon < epsilon){
                ret_epsilon = new_epsilon;
                }
            else{
                // This is an attempt to reduce epsilon even if the new and previous epsilon are equal
                ret_epsilon = 0.95*new_epsilon;
                }

        // See if we are finished
        int finished = 0;
     
            if((ret_epsilon < target_epsilon)|| ((ret_epsilon-target_epsilon) < 1e-6)){
                ret_epsilon = target_epsilon;
                finished = 1;
                }
            else{
         
                }
        

        
        char Dat[100];
    	int Dn;
    	
    	char Dirr[100];
      int Dr;
      Dr=sprintf(Dirr,"mkdir -p %d",position);
	const int dir= system(Dirr);

	  //this file contains epsilon
      Dn=sprintf(Dat,"%d/epsilon_curr%d.dat",position,Job);
   
      fstream my_vars_curr; 
      my_vars_curr.open(Dat,ios::out);

      if(my_vars_curr.is_open()){
	
			my_vars_curr << ret_epsilon << '\t';
			my_vars_curr << finished << std::endl;
			my_vars_curr.close();
      
      	}
      
      

}

void abcsmc::iterate_one_population(int position, double epsilon, int Job, int n_acc, gsl_rng* r, int position_prev, int naccepted){
 
 		simulation sim;
        int sampled = 0;
      
      //Print out the position of the next population 
      char pos[100];
      int ps;
      ps=sprintf(pos,"position%d.dat",Job);
      remove(pos);
   
      fstream my_vars_curr; 
      my_vars_curr.open(pos,ios::out);

      if(my_vars_curr.is_open()){
	
			my_vars_curr << position << std::endl;
			my_vars_curr.close();
      
      	}         		

    if(position>0){
		
        if(position>1){
		if(naccepted>Nparticles){
         		set_prevs(Job, weights_prev, kernel, margins_prev, position-1);  
       		 }
        	else{
        		set_prevs(Job, weights_prev, kernel, margins_prev, position-2);   
        	}   
			getKernel(kernel);
			std::cout << "kernel limits computed" << std::endl; 
            computeParticleWeights(weights_curr);
            std::cout << "weights calculated" << std::endl; 
            }
        else{
       	getKernel(kernel);
            for(int i=0; i<Nparticles; i++){
                weights_curr[i] = 1;
                }
            }
        normalizeWeights();
        
        //If developing for model selection
        modelMarginals();
        if(position==1){
        margins_prev = margins_curr;
        }
        
		PrintOutRes(Job, position);
        std::cout << "Printed weights, kernel" << std::endl;
        
        /////////////////////////////////
        // Prepare for next population //
        /////////////////////////////////
        int resim = 1;
		sampleTheParameter(this_parameter, r);		
		sim.simulate(this_parameter, epsilon, position, Job, resim, r);
        
        }
        
	else{
		int resim = 0;
		sampleTheParameterFromPrior(this_parameter);
		sim.simulate(this_parameter, epsilon, position, Job, resim, r);

	}        
         
}

void abcsmc::PrintOutRes(int Job, int pos){

	char Dat[100];
    int Dn;

	std::cout << "Printing position " << pos << " Job " << Job << std::endl;
	  //this file contains, position
      Dn=sprintf(Dat,"%d/vars_curr-J%d.dat",pos,Job);
   
      fstream my_vars_curr; 
      my_vars_curr.open(Dat,ios::out);

      if(my_vars_curr.is_open()){
	
	my_vars_curr << margins_curr << '\t';
	
	for(int k=0; k<(Nparticles-1); k++){
	my_vars_curr << weights_curr[k] << '\t';
	}
	
	my_vars_curr << weights_curr[Nparticles-1] << std::endl;
	
	for(int k=0; k<(nparameters-1); k++){
	my_vars_curr << kernel[k][0] << '\t';
	my_vars_curr << kernel[k][1] << '\t';
	}
	my_vars_curr << kernel[nparameters-1][0] << '\t';
	my_vars_curr << kernel[nparameters-1][1] << std::endl;
	
	my_vars_curr.close();
      }

}

void abcsmc::modelMarginals(){

   margins_curr = 0;
            for(int j=0; j<Nparticles; j++){   
                   margins_curr = 1;//margins_curr + weights_curr[j];
                   }
    
}

void abcsmc::getKernel(vector<vector<double> > &kernel){
  
		for(int param = 0; param<nparameters; param++){
		
			vector<double > pop_param(Nparticles);
			for(int q=0; q<Nparticles; q++){
			pop_param[q] = parameters[q][param];
			}
			double minimum=*min_element(pop_param.begin(),pop_param.end());
			double maximum=*max_element(pop_param.begin(),pop_param.end());
			double scale=(maximum-minimum);
			vector<double> a;
			a.resize(2);
			a[0] = -scale/2.0;
			a[1] = scale/2.0;

			kernel[param][0] = a[0];
			kernel[param][1] = a[1];
			std::cout << "kernel param: " << param << " kernel[param][0] " << kernel[param][0] << std::endl;
	  		std::cout << "kernel param: " << param << " kernel[param][1] " << kernel[param][1] << std::endl;
			}

}

// Here params refers to one particle
// The function changes params in place and returns the probability (which may be zero)
double abcsmc::perturbParticle(vector<double> &this_parameter, int thisp, gsl_rng* r){



		double prior_prob = 1;

		for(int n=0; n<nparameters; n++){
			
				this_parameter[n] = parameters[thisp][n] + runiform(r,kernel[n][0],kernel[n][1]);
		
		}

		//compute the likelihood
		prior_prob=1;
		for(int n=0; n<nparameters; n++){
			double x;			
			if ((this_parameter[n]>prior[n][1])||(this_parameter[n]<prior[n][0])){
                  x = 0.0;
                  }
        	else{
        		x = 1/(prior[n][1]-prior[n][0]);
        		}

			prior_prob = prior_prob*x;
			}

		return prior_prob;
}

void abcsmc::computeParticleWeights(vector<double> &weights_curr){

        for( int k=0; k<Nparticles; k++){

			vector <double> this_param(nparameters);
			for(int g=0; g<nparameters; g++){
            	this_param[g] = parameters[k][g];
            	std::cout << "this_param[g] " << this_param[g] << std::endl;
            }
            
            // calculate model prior probility 
            // particle prior probability
            double pprob = 1;
            for( int n=0; n<nparameters; n++){
                double x = 1.0;
                
                if ((this_param[n]>prior[n][1])||(this_param[n]<prior[n][0])){
                  x = 0.0;
                  }
        	else{
        		x = 1/(prior[n][1]-prior[n][0]);
        		}                              
                                                                             
                pprob = pprob*x;
            }

            double numer = pprob;
        
            double denom_m = 0;
            denom_m = margins_prev;
            
            double denom = 0;
            for(int j=0; j<Nparticles; j++){

                denom = denom + weights_prev[j] * getPdfParameterKernel(this_param, j);

			}
			std::cout << " numer " << numer << " denom_m " << denom_m << " denom " << denom << " margins_prev " << margins_prev << std::endl;		
			weights_curr[k] = numer/denom;
			
		}
}

double abcsmc::getPdfParameterKernel(vector<double> this_parameter, int j){

    double prob=1;
        for(int n=0; n<nparameters; n++){
			double kern; 
			if((parameters[j][n]>(this_parameter[n]+kernel[n][1]))||(parameters[j][n]<(this_parameter[n]+kernel[n][0]))){
				kern = 0.0;
			}
			else{
				kern = 1/((this_parameter[n]+kernel[n][1])-(this_parameter[n]+kernel[n][0]));
			}
				prob=prob*kern;
	    }
       
     return(prob);
 }
        
void abcsmc::normalizeWeights(void){
        double n = 0;
        for(int f=0; f<weights_curr.size(); f++){
          n += weights_curr[f];
         }
        for(int i = 0; i<Nparticles; i++){
            weights_curr[i] = weights_curr[i]/n;
		}
}

int abcsmc::sample_particle(gsl_rng* r){

    double u = runiform(r,0,margins_prev);
    double F = 0;
    int i;

    for(i=0; i<Nparticles; i++){
        
            F = F + weights_curr[i];
            

            if(F > u){
                break;
                }
                
            }

    return i;
    
}

void abcsmc::sampleTheParameterFromPrior(vector<double> &this_parameter){
 		
 		//void function as sampling done within simulation class
 	
			for(int n=0; n<nparameters; n++){
		
				this_parameter[n]=0;
			}
			
            
}

void abcsmc::sampleTheParameter(vector<double> &this_parameter, gsl_rng* r){

            double prior_prob = -1;            
            for(int i = 0; i<Nparticles; i++){
            std::cout << "weights_curr[i]" << weights_curr[i] << std::endl;
		}
            
            while(prior_prob <= 0){

				//sample putative particle from previous population
				int p = sample_particle(r);
				for(int nn=0; nn<nparameters; nn++){
					this_parameter[nn] = parameters[p][nn];
					}
					
				prior_prob = perturbParticle(this_parameter, p, r);
				
				if(prior_prob>0){
				
				for(int nn=0; nn<nparameters; nn++){
					std::cout << "this_parameter[nn]: " << this_parameter[nn] << " parameters[p][nn]: " << parameters[p][nn] << std::endl;
					}
				
				}						
	
			}            
}