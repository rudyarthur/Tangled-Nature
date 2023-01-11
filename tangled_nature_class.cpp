#include <complex>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <sys/time.h>
#include <stdio.h>
#include <bitset>
#include <stdlib.h>
#include <random>
#include <list>
#include <unordered_set>

#define L 20	  //length of the genome
#define N 1048576 //2^L

using namespace std;

/////////////////////
//Class (called Agents) to contain all the stuff about a species 
/////////////////////
class Agents {
public:
	int sa;				//species ID
	bitset<L> bin_sa;	//binary genome - sa in binary!
	int population;		//number of individuals of this species
	
	//initialise
	Agents(int sa_, int pop_) : sa(sa_), population(pop_) {
		bitset<L> tmp (sa);
		bin_sa = tmp;
	}
	
};

typedef list<Agents>::iterator Specie;

//given a list of species, check if species 'n' is in the list
Specie searchAgents(list<Agents> &species, int n) {
	for (Specie cur=species.begin(); cur != species.end(); ++cur){
		if(cur->sa == n) return cur;
	}
	return species.end();
}

///////////////////////////
//random number stuff
///////////////////////////




class Model {
public:
	
	/////////////////
	//parameters
	/////////////////
	int tgen = 0;				//number of generations elapsed
	double C = 100; 			//scaling parameter for interspecies interactions
	double sigma = 0;			//scaling parameter for species environment interactions, set to 0 for tandard TNM.
	double pkill = 0.2;			//probability of killing an individual
	double mu = 0.1; 			//resource abundance
	double pmut = 0.01; 		//mutation probability
	double A = 0; //2.0 + 4.0*pkill;	//Parameter controlling reproduction probability of fitness 0 individuals.
	int Npop_init = 500;		//initial total population
	int Npop = Npop_init;		//total population
	int max_gens = 100000;		//stop after this many generations
	double theta = 0.25;		//percentage of possible interspecies interactions which can are non-zero
	double mu_theta = 1.0;		//percentage of possible species-environment interactions which can be non-zero
	double nu = 0.000005;		//extra damping
	vector<bool> Jran1 = vector<bool>(N);				//helper list for implementing interspecies coupling matrix efficiently
	vector<double> Jran2 = vector<double>(N);			//helper list for implementing interspecies coupling matrix efficiently
	vector<double> Jran3 = vector<double>(N);			//helper list for implementing interspecies coupling matrix efficiently
	vector<bool> muran1 = vector<bool>(N);				//helper list for implementing species environment coupling matrix efficiently
	vector<double> muran2 = vector<double>(N);			//helper list for implementing species environment coupling matrix efficiently
	vector<double> muran3 = vector<double>(N);			//helper list for implementing species environment coupling matrix efficiently

	mt19937 mt_generator; // mt19937 is a standard //mersenne_twister_engine
	double rand_norm = 1.0/(double)(1.+mt_generator.max());
 
	list<Agents> species; //list of all extant species
	unordered_set<int> encountered; //set of all species encountered so far

	//random double between 0 and 1
	inline double mt_rand(){
		return (double)mt_generator()*rand_norm;
	} 
	//random double between -1 and 1
	inline double mt_rand_sym(){
		return -1.0 + 2.0*mt_generator()*rand_norm;
	} 
	//random int between 0 and N-1 inclusive
	inline int random_species(){
		return (int)( N*mt_rand() );
	}

	void copy_life_from(Model &m){
		Npop = m.Npop;
		species = m.species;
	}
	
	Model(int seed, int method=0) { 
		//seed the RNG
		default_random_engine generator(seed+123);
		mt_generator= mt19937(seed); 

		//set up the coupling matrices
		double oC = C;
		double omu = sqrt(mu);
		normal_distribution<double> distribution(0.0, 1.0);
		for(int i=0; i<N; i++){
			( mt_rand() < theta ) ? Jran1[i] = true : Jran1[i] = false;
			Jran2[i] = distribution(generator);
			Jran3[i] = oC*distribution(generator);
			
			( mt_rand() < mu_theta ) ? muran1[i] = true : muran1[i] = false;
			muran2[i] = distribution(generator);
			muran3[i] = sigma*distribution(generator);
		} 
		
		//Start off with Npop_init individuals of a single species
		if(method == 0){
			int rs = random_species();
			encountered.insert(rs);
			species.emplace_front(rs, Npop_init);
		} else if (method == 1){
			//start off with Npop_init species with population 1
			for(int i=0; i<Npop_init; ++i){
				int rs = random_species();
				encountered.insert(rs);
				species.emplace_front(rs, 1);
			}
		}
	}

	//does a lot of heavy lifting - returns 'fitness' of a species due to interspecies interactions
	inline double calc_HI(Specie elem){
		double sum = 0; 
		double musum = 0;
		for (Specie cur=species.begin(); cur != species.end(); ++cur){ //loop over species
			if(cur->bin_sa != elem->bin_sa){	//no self interaction
				bitset<L> bin_sab = cur->bin_sa^elem->bin_sa;		//trick to implement coupling matrix
				long int sab = bin_sab.to_ulong();				
				if( Jran1[sab] ){
					sum += Jran2[sab] * Jran3[cur->sa] * cur->population;	//inter species
				}
				if( muran1[sab] ){
					musum += muran2[sab] * muran3[cur->sa] * cur->population; //species environment
				}
			}
		}
		return sum/Npop - musum;
	}

	//returns E = total species environment interaction, only use this for printing
	inline double calc_E(Specie elem){
		double musum = 0;
		for (Specie cur=species.begin(); cur != species.end(); ++cur){
			if(cur->bin_sa != elem->bin_sa){
				bitset<L> bin_sab = cur->bin_sa^elem->bin_sa;
				long int sab = bin_sab.to_ulong();
				if( muran1[sab] ){
					musum += muran2[sab] * muran3[cur->sa] * cur->population;
				}
			}
		}
		return musum;
	}
	//returns F = interaction of species with every other species
	inline double calc_F(Specie elem){
		double sum = 0; 
		for (Specie cur=species.begin(); cur != species.end(); ++cur){
			if(cur->bin_sa != elem->bin_sa){
				bitset<L> bin_sab = cur->bin_sa^elem->bin_sa;
				long int sab = bin_sab.to_ulong();
				if( Jran1[sab] ){
					sum += Jran2[sab] * Jran3[cur->sa] * cur->population;
				}
			}
		}
		return sum;
	}

	//Total fitness, adds damping from carrying capacity.
	inline double calc_H(Specie elem){
		return calc_HI(elem) - mu*Npop - nu*Npop*Npop;;
	}

	//turn fitness into reproduction probability
	inline double poff(Specie elem){	return 1.0/(1.0 + exp(A-calc_H(elem))); }

	//generate offspring of species 'elem' with mutation
	inline void asexual(Specie elem){
		++Npop; 	//1 new individual
		bitset<L> bin_new;	//new individual genome
		for(int i=0; i<L; i++){
			if(mt_rand() < pmut){ bin_new[i] = !elem->bin_sa[i]; }	//each bit can mutate with probability pmut
			else{ bin_new[i] = elem->bin_sa[i]; }					
		}
		if( bin_new != elem->bin_sa ){ //if new species
			Specie tmpAgents = searchAgents(species, bin_new.to_ulong() ); //have we seen this species already?
			if(tmpAgents == species.end()){ //if new species not on list
				species.emplace_front(bin_new.to_ulong(), 1); //add to lsit
				encountered.insert(bin_new.to_ulong());
			} else { //if new species already on list
				++tmpAgents->population; //increase that species count by 1
			}	
		} else {  //no mutation
			++elem->population; //increase species count by 1
		}
	}
	
	//generate offspring of species 'elem' with mutation, both offspring can mutate
	inline void asexual2(Specie elem){
		++Npop; 	//1 new individual
		//make the new ones
		for(int offspring = 0; offspring<2; ++offspring){
			bitset<L> bin_new;	//new individual genome
			for(int i=0; i<L; i++){
				if(mt_rand() < pmut){ bin_new[i] = !elem->bin_sa[i]; }	//each bit can mutate with probability pmut
				else{ bin_new[i] = elem->bin_sa[i]; }					
			}
			if( bin_new != elem->bin_sa ){ //if new species
				Specie tmpAgents = searchAgents(species, bin_new.to_ulong() ); //have we seen this species already?
				if(tmpAgents == species.end()){ //if new species not on list
					species.emplace_front(bin_new.to_ulong(), 1); //add to lsit
					encountered.insert(bin_new.to_ulong());
				} else { //if new species already on list
					++tmpAgents->population; //increase that species count by 1
				}	
			} else {  //no mutation
				++elem->population; //increase species count by 1
			}
		}
		//kill the old one
		--elem->population;		//reduce species population by 1
		if(elem->population == 0){  //if species is now extinct, remove from list
			species.erase(elem);
		}		
	}

	//implement death
	inline Specie kill(){
		double rand = mt_rand();
		double sum = 0;
		for (Specie cur=species.begin(); cur != species.end(); ++cur){
			sum += cur->population;
			if( Npop*rand <= sum ){ //gives each *individual* an equal chance to be chosen
				if( mt_rand() < pkill ){	
					--cur->population;		//reduce species population by 1
					if(cur->population == 0){  //if species is now extinct, remove from list
						species.erase(cur);
					}
					--Npop;				//reduce total population
					return species.end();	
				} else {
					return cur; //don't kill an individual
				}
			}
		}
		cerr << "kill failed! Npop = " <<  Npop << endl;
		cerr << "rand " << rand << endl;
		exit(1);
	}

	//choose a random individual
	inline Specie choose(){
		double rand = mt_rand();
		int sum = 0;
		for (Specie cur=species.begin(); cur != species.end(); ++cur){
			sum += cur->population;
			if( Npop*rand <= sum ){
				return cur;
			}
		}
		
		cerr << "choose failed! Npop = " <<  Npop << endl;
		cerr << "rand " << rand << endl;
		exit(1);
	}

	///////////////////////////
	//Output stuff
	///////////////////////////
	//print basic global stats, called infrequently so doesn't need to be that efficient.
	inline void print_stats(ofstream &pop_file, double percent=0.05){ //percent = condition for being in the core
		
		double E = 0;
		double F = 0;
		int max = 0; 
		int diversity = 0;
		for (Specie cur=species.begin(); cur != species.end(); ++cur){
			if( cur->population >  max ){ max = cur->population; }
			++diversity;
			E += ( cur->population*calc_E(cur) / (double)Npop );
			F += ( cur->population*calc_F(cur) / (double)Npop );
		}

		int core_pop = 0;
		int core_size = 0;
		for (Specie cur=species.begin(); cur != species.end(); ++cur){
			if( cur->population > percent * (double) max ){ 
				core_pop += cur->population;
				++core_size;
			}
		}
							   
		pop_file << tgen << " " //generation number 
		<< Npop << " " //number of individuals
		<< diversity << " " //number of species
		<< encountered.size() << " " //number of species ever seen
		<< core_pop << " " //number ofindividuals in core 
		<< core_size << " " //number ofspecies in core
		<< E << " " //Life's effect on environment
		<< F //Life's effect on each other
		<< endl;			

	}

	//dump the interaction network out
	inline void print_species_network(ofstream &network_file){
		for (Specie elem=species.begin(); elem != species.end(); ++elem){
			for (Specie cur=species.begin(); cur != species.end(); ++cur){
				if(cur->bin_sa != elem->bin_sa){	
					network_file << elem->sa << " " << cur->sa << " ";
								
					bitset<L> bin_sab = cur->bin_sa^elem->bin_sa;
					long int sab = bin_sab.to_ulong();
						
					if( Jran1[sab] ){
						network_file << Jran2[sab] * Jran3[cur->sa] << " ";
					} else {
						network_file << "0 ";
					}
						
					if( muran1[sab] ){
						network_file << (muran2[sab] * muran3[cur->sa]) << " ";
					} else {
						network_file << "0 ";
					}
					network_file << cur->population << "\n";		
				}
			} 
		}	
	}

	//print species in core
	inline void print_core(ofstream &core_file, double percent=0.05){
		
		int max = 0; 
		for (Specie cur=species.begin(); cur != species.end(); ++cur){
			if( cur->population >  max ){ max = cur->population; }
		}
		core_file << tgen << " " << Npop;
		for (Specie cur=species.begin(); cur != species.end(); ++cur){
			if( cur->population > percent * (double) max ){ 
				core_file << " " << cur->sa;
			}
		}
		core_file << endl;			

	}
	
	string name_tag(int it){ 
		//put all the parameters in the filename
		string name_tag = "_seed"+ to_string(it);
		name_tag += "_mu"+ to_string(mu); 
		name_tag += "_nu"+ to_string(nu); 
		name_tag += "_sigma" + to_string(sigma); 
		name_tag += "_pmut"+ to_string(pmut); 
		name_tag +=".dat";
		return name_tag;
	}
		
		
	//run the model
	void run(ofstream &pop_file, int tmax=-1){
		if(tmax < 0){ tmax = max_gens; }

		int t = 0;
		double lgen = Npop/pkill; //how many birth/deaths per generation

		//start iteration
		do{
			
			Specie sa = kill(); if(Npop == 0){ break; } //choose an individual and kill with prob pkill
			if(sa == species.end()) sa = choose();	//if we killed the individual, choose another one
			if( mt_rand() < poff(sa) ){				//individual reproduces with probability poff
				asexual(sa);		//reproduce asexually
			}
			++t; //counter
		
			if(t >= lgen){	//generation is over
				t = 0; ++tgen; lgen = Npop/pkill; //recalculate generation length
				print_stats(pop_file);

				if(Npop > 1e5){		//annoying. Catches population explosions that make the simulation impossible to run.
					pop_file << "Large pop abort" << endl;
					break;
				}
				
			}
			
		} while(tgen < tmax); //stopping condition. run for max_gens

	}

};

int main(int argc, char *argv[]){

	
	if(argc < 3){ cerr << "usage: ./tangled_nature seed path\n e.g. ./tangled_nature 1 ./" << endl; exit(1); }

	int it = atoi(argv[1]);
	int seed = 123*it + 12345; //random seed
    string path = argv[2];
    
    /*{ //normal run
		Model t(seed);

		ofstream pop_file; 
		string output_filename = path + "single_run" + t.name_tag(it);
		pop_file.open (output_filename.c_str());
		
		t.run(pop_file, 1000);
		
		pop_file.close();
	}
	{ //do half the run, then the rest (note you have to put in the max time!)
		
		Model t(seed);

		ofstream pop_file; 
		string output_filename = path + "split_run" + t.name_tag(it);
		pop_file.open (output_filename.c_str());
		
		t.run(pop_file, 500);
		t.run(pop_file, 1000);
		
		pop_file.close();
		
	}*/
	/*{//copying the whole model
		
		Model t(seed);

		ofstream pop_file; 
		string output_filename = path + "copy_run" + t.name_tag(it);
		pop_file.open (output_filename.c_str());
		
		t.run(pop_file, 500);
		
		Model s(t); //initialise the model s with t
		//could change some parameters of the model here
		//e.g. s.mu /= 2
		s.run(pop_file, 1000);
		
		pop_file.close();
		
	}
	{//copying jsut the life
		Model t(seed);

		ofstream pop_file; 
		string output_filename = path + "copy2_run" + t.name_tag(it);
		pop_file.open (output_filename.c_str());
		
		t.run(pop_file, 500);
		
		Model s(seed); 
		s.copy_life_from(t); //This does not copy the random state, or the run parameters
		//can do that here (or not if you don't want to)
		s.tgen = t.tgen;
		s.mt_generator = t.mt_generator;
		
		s.run(pop_file, 1000);
		
		pop_file.close();
	}*/
	{ //example
		Model t(seed);		

		ofstream pop_file; 
		string output_filename = path + "param_run" + t.name_tag(it);
		pop_file.open (output_filename.c_str());


		//warm up
		t.run(pop_file, 500);
		
		t.mu /= 2;
		t.run(pop_file, 1000);

		t.mu *= 2;
		t.run(pop_file, 1500);
				
		pop_file.close();
		
	}
	return 0;
}


