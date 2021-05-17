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
double sigma = 0;			//scaling parameter for species environment interactions, set to 0 for tandard TNM.
double pkill = 0.2;			//probability of killing an individual
double mu = 0.1; 			//resource abundance
double pmut = 0.01; 		//mutation probability
double A = 0; //2.0 + 4.0*pkill;	//Parameter controlling reproduction probability of fitness 0 individuals.
int Nactive_init = 500;		//total active population
int Nactive = Nactive_init;	//total population
int Ndormant = 0;			//total dormant population
int max_gens = 1000;		//stop after this many generations
double theta = 0.25;		//percentage of possible interspecies interactions which can are non-zero
double mu_theta = 1.0;		//percentage of possible species-environment interactions which can be non-zero
double nu = 0.00000;		//extra damping
bool Jran1[N];				//helper list for implementing interspecies coupling matrix efficiently
double Jran2[N];			//helper list for implementing interspecies coupling matrix efficiently
double Jran3[N];			//helper list for implementing interspecies coupling matrix efficiently
bool muran1[N];				//helper list for implementing species environment coupling matrix efficiently
double muran2[N];			//helper list for implementing species environment coupling matrix efficiently
double muran3[N];			//helper list for implementing species environment coupling matrix efficiently

double sleepiness[N];			//dormancy strategies
double wakiness[N];				//dormancy strategies
double S = 0;					//scaling parameter for sleepiness
double Sm = 1000;				//offset parameter for sleepiness
double eta = 0;					//less likely to die while dormant


list<Agents> active_species; //list of all extant species
list<Agents> dormant_species; //list of all extant species
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

		sleepiness[i] = S*distribution(generator) - Sm;		
		wakiness[i] = S*distribution(generator) - Sm;		
	} 
	
	//Start off with Nactive_init individuals of a single species
	//int rs = random_species();
	//encountered.insert(rs);
	//species.emplace_front(rs, Nactive_init);
	
	//start off with Nactive_init species with population 1
	for(int i=0; i<Nactive; ++i){
		int rs = random_species();
		encountered.insert(rs);
		active_species.emplace_front(rs, 1);
	}
}








//implement death
inline void kill(list<Agents> &species, vector<Specie> &sas, int &Npop, int num, double pk){ //argument is number of agents to advance

	bool chosen[num];
	double rand[num]; 
	for(int i=0; i<num; ++i){ chosen[i] = false; rand[i] = mt_rand(); }
	unordered_map<int, Specie> extinct;
			
	int ip = Npop;
	
	//cout << "START  KILL" << endl;		
	for (int i=0; i<sas.size(); ++i){
		Specie cur = sas[i];
		if( cur == species.end() ){ continue; }
		
		if( mt_rand() < pk ){	//try to kill
			//cout << "killed " << cur->sa << "," << cur->population << endl;
			if(cur->population > 0){ //we could theoreticall kill 2 members of a species of population 1...
				--cur->population;		//reduce species population by 1
				--Npop;					//reduce total population
				//cout << "killed (" << cur->sa << "," << cur->population << ")" << endl;
				if(cur->population == 0){ extinct[sas[i]->sa] = sas[i]; }
				
			}
		}
		
	}
	
	//very annoying, if the same species is selected twice but goes extinct causes memory errors. This fixes it.
	for(auto it = extinct.begin(); it != extinct.end(); ++it){ species.erase(it->second); }
		
		
}

//choose a random individual
inline void choose_num(int &num_active, int &num_dormant, int num){
	
	
	num_active = 0;
	num_dormant = 0;
	int Npop = Nactive+Ndormant;
	for (int i=0; i<num; ++i){
		if(Npop*mt_rand() <= Nactive ){
			++num_active;
		} else {
			++num_dormant;
		}
	}



}

//choose a random individual
inline vector<Specie> choose(list<Agents> &species, int Npop, int num){
	
	double rand[num]; 
	vector<Specie> sas(num, species.end());

	for(int i=0; i<num; ++i){ rand[i] = mt_rand(); }
	
	int sum = 0;
	for (Specie cur=species.begin(); cur != species.end(); ++cur){
		sum += cur->population;
		for(int i=0; i<num; ++i){
			if( sas[i] == species.end() && Npop*rand[i] <= sum ){
				sas[i] = cur;
			}
		}
	}

	
	
	/*cout << "choose returns ";
	for (int i=0; i<sas.size(); ++i){ //loop over species
		cout << "(" << sas[i]->sa << "," << sas[i]->population << ")";
	} cout << endl;*/
	
	return sas;

}

//computes reproduction probability for the chosen species (sas)
inline void calc_fitness(vector<Specie>& sas, unordered_map< int, double > &fitness){
	
	/*cout << "Computing interaction for ";
	for (int i=0; i<sas.size(); ++i){ //loop over species
		cout << "(" << sas[i]->sa << "," << sas[i]->population << ")";
	} cout << endl;*/
	
	for (int i=0; i<sas.size(); ++i){ //loop over species
		
		if( fitness.find( sas[i]->sa ) != fitness.end() ){ continue; } //dead or already calculated it
		
		Specie elem = sas[i];
		double sum = 0;
		double musum = 0;
		
		for (Specie cur=active_species.begin(); cur != active_species.end(); ++cur){  //only active species contribute?
			
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
		
		fitness[ sas[i]->sa ] = (sum/Nactive - musum - mu*Nactive - nu*Nactive*Nactive);
		
	}

	/*cout << "offspring prob for ";
	for (int i=0; i<sas.size(); ++i){ //loop over species
		cout << "(" << sas[i]->sa << "," << fitness[ sas[i]->sa ] << ")";
	} cout << endl;*/

}

inline double poff(double fitness){	return 1.0/(1.0 + exp(A-fitness)); }

//generate offspring of species with mutation
inline void asexual(vector<Specie> &sas, unordered_map<int, double> &fitness){

	for (int j=0; j<sas.size(); ++j){ //loop over species
		if( mt_rand() < poff(fitness[ sas[j]->sa ])  ){ //should we reproduce?
			
			//cout << "reproduce " << sas[j]->sa << " " << sas[j]->bin_sa.to_string() << endl;
			
			++Nactive; 	//1 new individual
			bitset<L> bin_new;	//new individual genome
			for(int i=0; i<L; i++){
				if(mt_rand() < pmut){ bin_new[i] = !sas[j]->bin_sa[i]; }	//each bit can mutate with probability pmut
				else{ bin_new[i] = sas[j]->bin_sa[i]; }					
			}
			if( bin_new != sas[j]->bin_sa ){ //if new species
				Specie tmpAgents = searchAgents(active_species, bin_new.to_ulong() ); //have we seen this species already? (TODO should use set for this)
				if(tmpAgents == active_species.end()){ //if new species not on list
					active_species.emplace_front(bin_new.to_ulong(), 1); //add to list
					encountered.insert(bin_new.to_ulong());
				} else { //if new species already on list
					++tmpAgents->population; //increase that species count by 1
				}	
			} else {  //no mutation
				++sas[j]->population; //increase species count by 1
			}
			
			//cout << "became " << bin_new.to_ulong() << " " << bin_new.to_string() << endl;

		}
	}
	
}

///////////////////////////
//DORMANCY/RESUSCITATION HAPPENS HERE
///////////////////////////
inline void sleep_wake(list<Agents> &from_species, list<Agents> &to_species, int &Nfrom, int &Nto, bool sleep, vector<Specie> &sas, unordered_map<int, double> &fitness){

	unordered_map<int, Specie> extinct;

	for (int i=0; i<sas.size(); ++i){ //loop over species
	
		if(sas[i]->population > 0){

			double sprob;
			if( sleep ){
				//TODO sprob = 1.0/(1.0 + exp( fitness[ sas[i]->sa ] - sleepiness[ sas[i]->sa ] ) ); 
				sprob = 0.01; 
			} else {
				//TODO sprob = 1 - 1.0/(1.0 + exp( fitness[ sas[i]->sa ] - wakiness[ sas[i]->sa ] ) ); 
				//TODO sprob = 1 - 1.0/(1.0 + exp( fitness[ sas[i]->sa ] - sleepiness[ sas[i]->sa ] ) ); 
				sprob = 0.01;
			}

			if( mt_rand() <  sprob ){	
				--sas[i]->population;		//reduce species population by 1
				if(sas[i]->population == 0){ extinct[sas[i]->sa] = sas[i]; }
				
				//increase species population by 1
				Specie tmpAgents = searchAgents(to_species, sas[i]->sa ); //have we seen this species already? 
				if(tmpAgents == to_species.end()){ //if new species not on list
					to_species.emplace_front(sas[i]->sa, 1); //add to list
				} else { //if species already on list
					++tmpAgents->population; //increase that species count by 1
				}	
						
				--Nfrom;					//reduce total active population
				++Nto;					//increase total dormant population				
				
			}
			
		}
		
	}
	
	for(auto it = extinct.begin(); it != extinct.end(); ++it){ from_species.erase(it->second); }
			
}

///////////////////////////
//Output stuff
///////////////////////////
//returns E = total species environment interaction, only use this for printing
inline double calc_E(list<Agents> &species, Specie elem){
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
inline double calc_F(list<Agents> &species, Specie elem){
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
//print basic global stats, called infrequently so doesn't need to be that efficient.
inline void print_stats(ofstream &pop_file, double percent=0.05){ //percent = condition for being in the core
	
	double E = 0;
	double F = 0;
	int max = 0; 
	int active_diversity = 0;
	int dormant_diversity = 0;
	int check = 0;
	for (Specie cur=active_species.begin(); cur != active_species.end(); ++cur){
		if( cur->population >  max ){ max = cur->population; }
		++active_diversity;
		E += ( cur->population*calc_E(active_species, cur) / (double)Nactive );
		F += ( cur->population*calc_F(active_species, cur) / (double)Nactive );
	}
	for (Specie cur=dormant_species.begin(); cur != dormant_species.end(); ++cur){
		++dormant_diversity;
	}
	
	int core_pop = 0;
	int core_size = 0;
	for (Specie cur=active_species.begin(); cur != active_species.end(); ++cur){
		if( cur->population > percent * (double) max ){ 
			core_pop += cur->population;
			++core_size;
		}
	}
	                       
	pop_file << tgen << " " //generation number 
	<< Nactive << " " //number of active individuals
	<< Ndormant << " " //number of dormant individuals
	<< active_diversity << " " //number of species
	<< dormant_diversity << " " //number of species
	<< encountered.size() << " " //number of species ever seen
	<< core_pop << " " //number ofindividuals in core 
	<< core_size << " " //number ofspecies in core
	<< E << " " //Life's effect on environment
	<< F //Life's effect on each other
	<< endl;			

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
	double lgen = Nactive/pkill; //how many birth/deaths per generation

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
		int num_active, num_dormant;
		choose_num(num_active, num_dormant, num);		
		
		//choose some individuals
		vector<Specie> active_individuals = (num_active > 0) ? choose(active_species, Nactive, num_active) : vector<Specie>{};
		vector<Specie> dormant_individuals = (num_dormant > 0) ? choose(dormant_species, Ndormant, num_dormant) : vector<Specie>{};

		//cout << tgen << " " << t << " " << Nactive << " " << Ndormant << " " << active_individuals.size() << " " << dormant_individuals.size() << endl;

		//try to kill them
		kill(active_species, active_individuals, Nactive, num, pkill);
		kill(dormant_species, dormant_individuals, Ndormant, num, eta*pkill);
		if(Nactive + Ndormant == 0){ break; }

		//choose some more individuals
		choose_num(num_active, num_dormant, num);			
		active_individuals = (num_active > 0) ? choose(active_species, Nactive, num_active) : vector<Specie>{};
		dormant_individuals = (num_dormant > 0) ? choose(dormant_species, Ndormant, num_dormant) : vector<Specie>{};

		//compute fitness for all extant individuals
		unordered_map< int, double > fitness; //the map thing saves time by not calculating the interaction for the same species twice
		calc_fitness(active_individuals, fitness);
		calc_fitness(dormant_individuals, fitness);
 
		//reproduce the active guys
		asexual(active_individuals, fitness);

		//(active -> dormant)
		if(Nactive > 0 && num_active > 0){
			sleep_wake(active_species, dormant_species, Nactive, Ndormant, true, active_individuals, fitness);
		}

		//(dormant -> active)
		if(Ndormant > 0 && num_dormant > 0){
			sleep_wake(dormant_species, active_species, Ndormant, Nactive, false, dormant_individuals, fitness);
		}


		t += num; //counter
		
		if(t >= lgen){	//generation is over
			t = 0; ++tgen; lgen = Nactive/pkill; //TODO, this is probably wrong

			print_stats(pop_file);

			if(Nactive > 1e5){		//Catches population explosions that make the simulation impossible to run.
				pop_file << "Large pop abort" << endl;
				break;
			}
			
		}

	}while(tgen < max_gens); //stopping condition. run for max_gens

	pop_file.close();

	return 0;
}


