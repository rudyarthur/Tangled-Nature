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
#include <unordered_map>

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
double totalPop(list<Agents> &species) {
	int sum = 0;
	for (Specie cur=species.begin(); cur != species.end(); ++cur){
		sum += cur->population;
	}
	return sum;
}

/////////////////
//parameters
/////////////////
int tgen = 0;				//number of generations elapsed
double C = 100; 			//scaling parameter for interspecies interactions
double sigma = 0.1;			//scaling parameter for species environment interactions, set to 0 for tandard TNM.
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
bool Jran1[N];				//helper list for implementing interspecies coupling matrix efficiently
double Jran2[N];			//helper list for implementing interspecies coupling matrix efficiently
double Jran3[N];			//helper list for implementing interspecies coupling matrix efficiently
bool muran1[N];				//helper list for implementing species environment coupling matrix efficiently
double muran2[N];			//helper list for implementing species environment coupling matrix efficiently
double muran3[N];			//helper list for implementing species environment coupling matrix efficiently


list<Agents> species; //list of all extant species
unordered_set<int> encountered; //set of all species encountered so far

///////////////////////////
//random number stuff
///////////////////////////
mt19937 mt_generator; // mt19937 is a standard //mersenne_twister_engine
double rand_norm = 1.0/(double)(1.+mt_generator.max());
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


//number of mutations to get from a to b
inline int hamming_distance( bitset<L> a, bitset<L> b){
	int sum = 0;
	for(int i=0;i<L;i++){ if(a[i] != b[i]) ++sum; }
	return sum;
}

//initialise everything
inline void init(int seed){

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
	int rs = random_species();
	encountered.insert(rs);
	//species.emplace_front(rs, Npop_init);
	//start off with Npop_init species with population 1
	Npop = 0;
	for(int i=0; i<20; ++i){
		int rs = random_species();
		encountered.insert(rs);
		int rp = max(1, (int) ( mt_rand()*100. ) );
		species.emplace_front(rs, rp);
		Npop += rp;
	}
}

//computes reproduction probability for the chosen species (sas)
inline unordered_map< int, double > poff(vector<Specie>& sas){
	
	unordered_map< int, double > offspring_prob;
	
	for (int i=0; i<sas.size(); ++i){ //loop over species
		
		if( offspring_prob.find( sas[i]->sa ) != offspring_prob.end() ){ continue; } //already calculated it
		
		Specie elem = sas[i];
		double sum = 0;
		double musum = 0;
		
		for (Specie cur=species.begin(); cur != species.end(); ++cur){
			
			if(cur != elem){ //no self
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
		offspring_prob[ sas[i]->sa ] = ( 1.0/(1.0 + exp(A - sum/Npop - musum - mu*Npop - nu*Npop*Npop)) );
		
	}
	return offspring_prob;


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


//generate offspring of species 'elem' with mutation
		asexual(sas, offspring_prob);

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
/*inline void asexual2(Specie elem){
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
}*/

//implement death
inline vector<Specie> kill(int num){ //argument is number of agents to advance

	double rand[num]; for(int i=0; i<num; ++i){ rand[i] = mt_rand(); }
	double sum = 0;
	vector<Specie> sas(num, species.end());
	
	for (Specie cur=species.begin(); cur != species.end(); ++cur){
		sum += cur->population;
		
		for(int i=0; i<num; ++i){ 
			
			if( sas[i] == species.end() && Npop*rand[i] <= sum ){ //gives each *individual* an equal chance to be chosen
			
				if( mt_rand() < pkill ){	//try to kill
					--cur->population;		//reduce species population by 1
					--Npop;					//reduce total population
					if(cur->population == 0){  //if species is now extinct, remove from list and abort
						cur = species.erase(cur); //cur now = the element after the one that was erased
						--cur; //decrement because we are about to increment in the loop counter!
						break;
					}
				} else {
					sas[i] = cur; //don't kill an individual
				}
			
			}
		}
		
		
	}
	return sas;
	
}

//choose a random individual
inline vector<Specie> choose(int num){
	double rand[num]; for(int i=0; i<num; ++i){ rand[i] = mt_rand(); }
	int sum = 0;
	vector<Specie> sas(num, species.end());
	
	for (Specie cur=species.begin(); cur != species.end(); ++cur){
		sum += cur->population;
		for(int i=0; i<num; ++i){
			if( sas[i] == species.end() && Npop*rand[i] <= sum ){
				sas[i] = cur;
			}
		}
	}
	return sas;

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


int main(int argc, char *argv[]){

	if(argc < 3){ cerr << "usage: tangled_nature seed path" << endl; exit(1); }

	int it = atoi(argv[1]);
	int seed = 123*it + 12345; //random seed
    string path = argv[2];
    int num = 1;
	if(argc >= 4){
		num = atoi(argv[3]);
	}
	cerr << "advance " << num << endl;

	init(seed);
	double t = 0;
	double lgen = Npop/pkill; //how many birth/deaths per generation

	//put all the parameters in the filename
	string name_tag = "_seed"+ to_string(it);
	name_tag += "_num"+ to_string(num); 
	name_tag += "_mu"+ to_string(mu); 
	name_tag += "_nu"+ to_string(nu); 
	name_tag += "_sigma" + to_string(sigma); 
	name_tag += "_pmut"+ to_string(pmut); 
	name_tag +=".dat";

	string name1 = path + "popplot"; name1 += name_tag;
	ofstream pop_file; pop_file.open (name1.c_str());

	//string name2 = path + "coreplot"; name2 += name_tag;
	//ofstream core_file; core_file.open (name2.c_str());

	//start iteration
	do{
		vector<Specie> sas = kill(num); if(Npop < num){ break; } //choose num individuals and kill with prob pkill
		//compute fitness for all extant individuals
		sas = choose(num); //if we killed someone, choose some new ones

		unordered_map<int, double> offspring_prob = poff(sas);
		asexual(sas, offspring_prob);
		
		cout << tgen << " " << t << "/" << lgen << " N = " << Npop << "==" << totalPop(species) << endl;
		t += num; //counter

		
		if(t >= lgen){	//generation is over
			t = 0; ++tgen; lgen = Npop/pkill; //recalculate generation length
			print_stats(pop_file);
			//print_core(core_file);

			if(Npop > 1e5){		//annoying. Catches population explosions that make the simulation impossible to run.
				//Fix this at some point...
				pop_file << "Large pop abort" << endl;
				break;
			}
			
		}
		
	}while(tgen < max_gens); //stopping condition. run for max_gens
	pop_file.close();
	//core_file.close();

	return 0;
}


