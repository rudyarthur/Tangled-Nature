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

////////////////////
//Node = species
////////////////////
class Node {
public:
	int sa;				//species ID
	bitset<L> bin_sa;	//binary genome - sa in binary!
	int population;		//number of individuals of this species
	
	//initialise
	Node(int sa_, int pop_) : sa(sa_), population(pop_) {
		bitset<L> tmp (sa);
		bin_sa = tmp;
	}
	
};
//given a list of species, check if species 'n' is in the list
list<Node>::iterator searchNode(list<Node> &species, int n) {
	for (list<Node>::iterator cur=species.begin(); cur != species.end(); ++cur){
		if(cur->sa == n) return cur;
	}
	return species.end();
}

/////////////////
//parameters
/////////////////
int tgen = 0;				//number of generations
double C = 100; 			//scaling parameter for interspecies interactions
double sigma = 0;			//scaling parameter for species environment interactions
double pkill = 0.2;			//probability of killing an individual
double mu = 0.1; 			//resource abundance
double pmut = 0.01; 		//mutation probability
double A = 0; //2.0 + 4.0*pkill;	//Parameter controlling reproduction probability of fitness 0 individuals.
int Npop_init = 500;		//initial total population
int Npop = Npop_init;		//total population
int max_gens = 100000;		//stop after this many generations
double theta = 0.25;		//percentage of possible interspecies interactions which can are non-zero
double mu_theta = 1.0;		//percentage of possible species-environment interactions which can be non-zero
bool Jran1[N];				//helper list for implementing interspecies coupling matrix efficiently
double Jran2[N];			//helper list for implementing interspecies coupling matrix efficiently
double Jran3[N];			//helper list for implementing interspecies coupling matrix efficiently
bool muran1[N];				//helper list for implementing species environment coupling matrix efficiently
double muran2[N];			//helper list for implementing species environment coupling matrix efficiently
double muran3[N];			//helper list for implementing species environment coupling matrix efficiently


list<Node> species; //list of all extant species
unordered_set<int> encountered;

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

	default_random_engine generator(seed+123);
	mt_generator= mt19937(seed); 

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
	
	int rs = random_species();
	encountered.insert(rs);
	species.emplace_front(rs, Npop_init);
}

//does a lot of heavy lifting - returns 'fitness' of a species
inline double calc_HI(list<Node>::iterator elem){
	double sum = 0; 
	double musum = 0;
	for (list<Node>::iterator cur=species.begin(); cur != species.end(); ++cur){ //loop over species
		if(cur->bin_sa != elem->bin_sa){	//no self interaction
			bitset<L> bin_sab = cur->bin_sa^elem->bin_sa;		//trick to implement coupling matrix
			long int sab = bin_sab.to_ulong();				
			if( Jran1[sab] ){
				sum += Jran2[sab] * Jran3[cur->sa] * cur->population;	//inter species
			}
			if( muran1[sab] ){
				musum += muran2[sab] * muran3[cur->sa] * cur->population * cur->population; //species environment
			}
		}
	}
	return sum - musum;
}




//returns E = total species environment interaction
inline double calc_E(list<Node>::iterator elem){
	double musum = 0;
	for (list<Node>::iterator cur=species.begin(); cur != species.end(); ++cur){
		if(cur->bin_sa != elem->bin_sa){
			bitset<L> bin_sab = cur->bin_sa^elem->bin_sa;
			long int sab = bin_sab.to_ulong();
			if( muran1[sab] ){
				musum += muran2[sab] * muran3[cur->sa] * cur->population * cur->population;
			}
		}
	}
	return musum;
}
//returns F = interaction of species with every other species
inline double calc_F(list<Node>::iterator elem){
	double sum = 0; 
	for (list<Node>::iterator cur=species.begin(); cur != species.end(); ++cur){
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

//global resource constraints
inline double calc_H(list<Node>::iterator elem){
	return calc_HI(elem)/Npop - mu*Npop;
}

//reproduction probability	
inline double poff(list<Node>::iterator elem){	return 1.0/(1.0 + exp(A-calc_H(elem))); }

//generate offspring of species 'elem' with mutation
inline void asexual(list<Node>::iterator elem){
	++Npop; 	//1 new individual
	bitset<L> bin_new;	//new individual genome
	for(int i=0; i<L; i++){
		if(mt_rand() < pmut){ bin_new[i] = !elem->bin_sa[i]; }	//each bit can mutate with probability pmut
		else{ bin_new[i] = elem->bin_sa[i]; }					
	}
	if( bin_new != elem->bin_sa ){ //if new species
		list<Node>::iterator tmpNode = searchNode(species, bin_new.to_ulong() ); //have we seen this species already?
		if(tmpNode == species.end()){ //if new species not on list
			species.emplace_front(bin_new.to_ulong(), 1); //add to lsit
			encountered.insert(bin_new.to_ulong());
		} else { //if new species already on list
			++tmpNode->population; //increase that species count by 1
		}	
	} else {  //no mutation
		++elem->population; //increase species count by 1
	}
}

//implement death
inline list<Node>::iterator kill(){
	double rand = mt_rand();
	double sum = 0;
	for (list<Node>::iterator cur=species.begin(); cur != species.end(); ++cur){
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
inline list<Node>::iterator choose(){
	double rand = mt_rand();
	int sum = 0;
	for (list<Node>::iterator cur=species.begin(); cur != species.end(); ++cur){
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
	for (list<Node>::iterator cur=species.begin(); cur != species.end(); ++cur){
		if( cur->population >  max ){ max = cur->population; }
		++diversity;
		E += ( cur->population*calc_E(cur) / (double)Npop );
		F += ( cur->population*calc_F(cur) / (double)Npop );
	}

	int core_pop = 0;
	int core_size = 0;
	for (list<Node>::iterator cur=species.begin(); cur != species.end(); ++cur){
		if( cur->population > percent * (double) max ){ 
			core_pop += cur->population;
			++core_size;
		}
	}
	//generation number   number of individuals    number of species    individuals in core     species in core     effect on environment    effect on each other     external resources
	pop_file << tgen << " " << Npop << " " << diversity << " " << encountered.size() << " " << core_pop << " " << core_size << " " << E << " " << F << " " << mu << endl;			

}

//dump the interaction network out
inline void print_species_network(ofstream &network_file){
	for (list<Node>::iterator elem=species.begin(); elem != species.end(); ++elem){
		for (list<Node>::iterator cur=species.begin(); cur != species.end(); ++cur){
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

//print species growth rates
inline void print_growth(ofstream &growth_file, double time, double percent=0.05){
	
	int max = 0; 
	for (list<Node>::iterator cur=species.begin(); cur != species.end(); ++cur){
		if( cur->population >  max ){ max = cur->population; }
	}
	
	for (list<Node>::iterator elem=species.begin(); elem != species.end(); ++elem){
		if( elem->population > percent * (double) max ){
			double ri = calc_F(elem) / (double)Npop;
			double ei = calc_E(elem) / (double)Npop;
			growth_file << setprecision(16) << time << " " << elem->sa << " " << ri << " " << ei << " " << elem->population << endl;				
		}
	}
	
}

//print species in core
inline void print_core(ofstream &core_file, double percent=0.05){
	
	int max = 0; 
	for (list<Node>::iterator cur=species.begin(); cur != species.end(); ++cur){
		if( cur->population >  max ){ max = cur->population; }
	}
	core_file << tgen << " " << Npop;
	for (list<Node>::iterator cur=species.begin(); cur != species.end(); ++cur){
		if( cur->population > percent * (double) max ){ 
			core_file << " " << cur->sa;
		}
	}
	core_file << endl;			

}

//print core core interaction matrix
inline void print_core_core(string path, string name_tag, double percent=0.05){
	

	int max = 0; 
	for (list<Node>::iterator cur=species.begin(); cur != species.end(); ++cur){
		if( cur->population >  max ){ max = cur->population; }
	}

	string name2 = path + "corecore";   name2 += name_tag; ofstream Jcc_file; Jcc_file.open (name2.c_str());
	string name3 = path + "corecloud";  name3 += name_tag; ofstream Jcl_file; Jcl_file.open (name3.c_str());
	
	for (list<Node>::iterator elem=species.begin(); elem != species.end(); ++elem){
		
		if( elem->population > percent * (double) max ){ //elem is core
			for (list<Node>::iterator cur=species.begin(); cur != species.end(); ++cur){
				if(cur->bin_sa != elem->bin_sa){	
					if( cur->population > percent * (double) max ){ //cur is core
						Jcc_file << elem->sa << " " << cur->sa << " ";
							
						bitset<L> bin_sab = cur->bin_sa^elem->bin_sa;
						long int sab = bin_sab.to_ulong();
						
						if( Jran1[sab] ){
							Jcc_file << Jran2[sab] * Jran3[cur->sa] << " ";
						} else {
							Jcc_file << "0 ";
						}
						
						if( muran1[sab] ){
							Jcc_file << (muran2[sab] * muran3[cur->sa]) << " ";
						} else {
							Jcc_file << "0 ";
						}
						Jcc_file << cur->population << "\n";
						
						
					} else { //cur is cloud
						Jcl_file << elem->sa << " " << cur->sa << " ";
							
						bitset<L> bin_sab = cur->bin_sa^elem->bin_sa;
						long int sab = bin_sab.to_ulong();
						
						if( Jran1[sab] ){
							Jcl_file << Jran2[sab] * Jran3[cur->sa] << " ";
						} else {
							Jcl_file << "0 ";
						}
						
						if( muran1[sab] ){
							Jcl_file << (muran2[sab] * muran3[cur->sa]) << " ";
						} else {
							Jcl_file << "0 ";
						}
						Jcl_file << cur->population << "\n";

					}
				}
			}
		} 
	}
	
	Jcc_file.close();	
	Jcl_file.close();	

}



int main(int argc, char *argv[]){

	if(argc < 2){ cerr << "usage: tangled_nature seed " << endl; }

	int it = atoi(argv[1]);
	int seed = 123*it + 12345; //random seed
    string path = argv[2];

	init(seed);
	int t = 0;
	double lgen = Npop/pkill; //how many birth/deaths per generation

	//put all the parameters in the filename
	string name_tag = "_seed"+ to_string(it);
	name_tag += "_C"+ to_string(C); 
	name_tag += "_A"+ to_string(A); 
	name_tag += "_mu"+ to_string(mu); 
	name_tag += "_theta" + to_string(theta); 
	name_tag += "_mutheta" + to_string(mu_theta); 
	name_tag += "_pmut"+ to_string(pmut); 
	name_tag += "_L" + to_string(L);
	name_tag +=".dat";

	string name1 = path + "popplot"; name1 += name_tag;
	ofstream pop_file; pop_file.open (name1.c_str());

	//string name2 = path + "coreplot"; name2 += name_tag;
	//ofstream core_file; core_file.open (name2.c_str());

	//start iteration
	do{
		list<Node>::iterator sa = kill(); if(Npop == 0){ break; } //choose an individual and kill with prob pkill
		if(sa == species.end()) sa = choose();	//if we killed the individual, choose another one
		if( mt_rand() < poff(sa) ){				//individual reproduces with probability poff
			asexual(sa);		//reproduce asexually
		}
		++t; //counter

		
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


