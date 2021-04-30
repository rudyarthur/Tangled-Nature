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
//Class (called Node) to contain all the stuff about a species 
////////////////////
class Node {
public:
	int sa;				//species ID
	bitset<L> bin_sa;	//binary genome - sa in binary!
	int active_population;		//number of active individuals of this species
	int dormant_population;		//number of dormant individuals of this species
	
	//initialise
	Node(int sa_, int pop_){
		sa = sa_;
		active_population = pop_;
		bitset<L> tmp (sa);
		bin_sa = tmp;
		dormant_population = 0;
	}
	
};
//given a list of species, check if species 'n' is in the list
list<Node>::iterator searchNode(list<Node> &species, int n) {
	for (list<Node>::iterator cur=species.begin(); cur != species.end(); ++cur){
		if(cur->sa == n) return cur;
	}
	return species.end();
}
//total number of dormant individuals
int totalDormant(list<Node> &species) {
	int sum = 0;
	int count = 0;
	for (list<Node>::iterator cur=species.begin(); cur != species.end(); ++cur){
		sum += cur->dormant_population;
		if( cur->dormant_population > 0){ ++count; }		
	}
	return sum;
}
//total number of active individuals
int totalActive(list<Node> &species) {
	int sum = 0;
	int count = 0;	
	for (list<Node>::iterator cur=species.begin(); cur != species.end(); ++cur){
		sum += cur->active_population;
		if( cur->active_population > 0){ ++count; }
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
double eta = 1;				//less likely to die while dormant
double mu = 0.1; 			//resource abundance
double pmut = 0.01; 		//mutation probability
double A = 0; //2.0 + 4.0*pkill;	//Parameter controlling reproduction probability of fitness 0 individuals.
int Npop_init = 500;		//initial total population
int Nactive = Npop_init;		//total active population
int Ndormant = 0;				//total dormant population
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

double sleepiness[N];			//dormancy strategies
double S = 0;					//scaling parameter for sleepiness
double Sm = 1000;					//offset parameter for sleepiness


list<Node> species; //list of all extant species
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
	} 
	
	//Start off with Npop_init individuals of a single species
	int rs = random_species();
	encountered.insert(rs);
	species.emplace_front(rs, Npop_init);
	//start off with Npop_init species with population 1
	/*for(int i=0; i<Npop_init; ++i){
		int rs = random_species();
		encountered.insert(rs);
		species.emplace_front(rs, 1);
	}*/
}

//does a lot of heavy lifting - returns 'fitness' of a species due to interspecies interactions
inline pair<double, double> calc_HI(list<Node>::iterator elem){
	double sum = 0; 
	double musum = 0;
	for (list<Node>::iterator cur=species.begin(); cur != species.end(); ++cur){ //loop over species
		if(cur->bin_sa != elem->bin_sa){	//no self interaction
			bitset<L> bin_sab = cur->bin_sa^elem->bin_sa;		//trick to implement coupling matrix
			long int sab = bin_sab.to_ulong();				
			if( Jran1[sab] ){
				sum += Jran2[sab] * Jran3[cur->sa] * cur->active_population;	//inter species
			}
			if( muran1[sab] ){
				musum += muran2[sab] * muran3[cur->sa] * cur->active_population; //species environment
			}
		}
	}
	//return sum/Nactive - musum;
	return make_pair(sum, musum);
	
}




//returns E = total species environment interaction, only use this for printing
inline double calc_E(list<Node>::iterator elem){
	double musum = 0;
	for (list<Node>::iterator cur=species.begin(); cur != species.end(); ++cur){
		if(cur->bin_sa != elem->bin_sa){
			bitset<L> bin_sab = cur->bin_sa^elem->bin_sa;
			long int sab = bin_sab.to_ulong();
			if( muran1[sab] ){
				musum += muran2[sab] * muran3[cur->sa] * cur->active_population;
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
				sum += Jran2[sab] * Jran3[cur->sa] * cur->active_population;
			}
		}
	}
	return sum;
}

//Total fitness, adds damping from carrying capacity.
//inline double calc_H(list<Node>::iterator elem){return calc_HI(elem) - mu*Nactive - nu*Nactive*Nactive;}

//turn fitness into reproduction probability
inline double poff(double sum, double musum){	
	return 1.0/(1.0 + exp(A - (sum/Nactive - musum - mu*Nactive - nu*Nactive*Nactive) )  ); 
}

//generate offspring of species 'elem' with mutation
inline void asexual(list<Node>::iterator elem){
	++Nactive; 	//1 new individual
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
			++tmpNode->active_population; //increase that species count by 1
		}	
	} else {  //no mutation
		++elem->active_population; //increase species count by 1
	}
}
//generate offspring of species 'elem' with mutation, both offspring can mutate
/*inline void asexual2(list<Node>::iterator elem){
	++Nactive; 	//1 new individual
	//make the new ones
	for(int offspring = 0; offspring<2; ++offspring){
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
				++tmpNode->active_population; //increase that species count by 1
			}	
		} else {  //no mutation
			++elem->active_population; //increase species count by 1
		}
	}
	//kill the old one
	--elem->active_population;		//reduce species active_population by 1
	if(elem->active_population == 0){  //if species is now extinct, remove from list
		species.erase(elem);
	}		
}*/


//////////////////////////////////////////////
//HERE IS ALL THE NEW STUFF FOR DORMANCY!
//////////////////////////////////////////////

//kill active
inline list<Node>::iterator kill_active(list<Node>::iterator cur){

	if(Nactive == 0){ return species.end(); } //this shouldn't happen!

	if( mt_rand() < pkill ){	
		--cur->active_population;		//reduce species population by 1
		if(cur->active_population + cur->dormant_population == 0){  //if species is now extinct, remove from list
			species.erase(cur);
		}
		--Nactive;				//reduce total population
		return species.end();	
	} else {
		return cur; //don't kill an individual
	}

	//never see this
	cerr << "kill_active failed! Nactive = " <<  Nactive << " " << totalActive(species) << endl;
	cerr << "rand " << rand << endl;
	exit(1);
}
inline list<Node>::iterator kill_dormant(list<Node>::iterator cur){
	
	if(Ndormant == 0){ return species.end(); }
	
	if( mt_rand() < eta*pkill ){	
		--cur->dormant_population;		//reduce species population by 1
		if(cur->active_population + cur->dormant_population == 0){  //if species is now extinct, remove from list
			species.erase(cur);
		}
		--Ndormant;				//reduce total population
		return species.end();	
	} else {
		return cur; //don't kill an individual
	}

	cerr << "kill_dormant failed! Nactive,Ndormant  = " <<  Nactive << " " << Ndormant << endl;
	cerr << "rand " << rand << endl;
	exit(1);
}

//choose an individual
inline pair< list<Node>::iterator, bool > choose(){ //id, active/dormant

	if(Nactive+Ndormant == 0){ return make_pair( species.end(), true ); }
	
	int Npop = Nactive+Ndormant;
	double rand = mt_rand();
	int sum = 0;
	//I think this is right?
	for (list<Node>::iterator cur=species.begin(); cur != species.end(); ++cur){
		sum += cur->active_population;
		if( Npop*rand <= sum ){
			return make_pair( cur, true );
		}
		sum += cur->dormant_population;
		if( Npop*rand <= sum ){
			return make_pair( cur, false );
		}		
	}
	
	cerr << "choose active failed! Nactive, Ndormant = " <<  Nactive << "," << Ndormant << endl;
	cerr << "rand " << rand << endl;
	exit(1);
}
//choose a random active individual
inline list<Node>::iterator choose_active(){

	if(Nactive == 0){ return species.end(); }
	
	double rand = mt_rand();
	int sum = 0;
	for (list<Node>::iterator cur=species.begin(); cur != species.end(); ++cur){
		sum += cur->active_population;
		if( Nactive*rand <= sum ){
			return cur;
		}
	}
	cerr << "choose active failed! Nactive = " <<  Nactive << endl;
	cerr << "rand " << rand << endl;
	exit(1);
}
//choose a random dormant individual
inline list<Node>::iterator choose_dormant(){

	if(Ndormant == 0){ return species.end(); }
	
	double rand = mt_rand();
	int sum = 0;
	for (list<Node>::iterator cur=species.begin(); cur != species.end(); ++cur){
		sum += cur->dormant_population;
		if( Ndormant*rand <= sum ){
			return cur;
		}
	}
	cerr << "choose dorm failed! Ndormant = " <<  Ndormant << endl;
	cerr << "rand " << rand << endl;
	exit(1);
}

///////////////////////////
//DORMANCY ACTUALLY HAPPEN HERE
///////////////////////////
//implement dormancy
//takes species and 'conditions'
//MAKE SURE I GOT ALL THE SIGNS RIGHT!
inline list<Node>::iterator go_dormant(list<Node>::iterator cur, double f){

	if(Nactive == 0){ return species.end(); }

	
	double sprob = 1.0/(1.0 + exp( f-sleepiness[ cur->sa ] ) ); 

	//cout << "F = " << f << " si = " << sprob << " sleepiness " << sleepiness[ cur->sa ] << endl;
	if( mt_rand() <  sprob ){	

		--cur->active_population;		//reduce species active population by 1
		++cur->dormant_population;		//increase species dormant population by 1
		--Nactive;					//reduce total active population
		++Ndormant;					//increase total dormant population
		return species.end();	
	} else {
		return cur; //don't go to sleep
	}
	
	cerr << "go_dormant failed! Nactive, Ndormant = " <<  Nactive << " " << Ndormant << endl;
	cerr << "rand " << rand << endl;
	exit(1);
}
//implement resuscitation
inline list<Node>::iterator go_active(list<Node>::iterator cur, double f){
	
	if(Ndormant == 0){ return species.end(); }

	double sprob = 1.0/(1.0 + exp( f-sleepiness[ cur->sa ] ) ); 

	
	if( mt_rand() <  1-sprob ){	
		++cur->active_population;		//increase species active population by 1
		--cur->dormant_population;		//decrease species dormant population by 1
		++Nactive;					//increase total active population
		--Ndormant;				//decrease total dormant population
		return species.end();	
	} else {
		return cur; //don't wake up
	}

	cerr << "go_active failed! Nactive, Ndormant = " <<  Nactive << " " << Ndormant << endl;
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
		if( cur->active_population >  max ){ max = cur->active_population; }
		++diversity;
		E += ( cur->active_population*calc_E(cur) / (double)Nactive );
		F += ( cur->active_population*calc_F(cur) / (double)Nactive );
	}

	int core_pop = 0;
	int core_size = 0;
	for (list<Node>::iterator cur=species.begin(); cur != species.end(); ++cur){
		if( cur->active_population > percent * (double) max ){ 
			core_pop += cur->active_population;
			++core_size;
		}
	}
	                       
	pop_file << tgen << " " //generation number 
	<< Nactive << " " //number of active individuals
	<< Ndormant << " " //number of dormant individuals
	<< diversity << " " //number of species (dormant or active)
	<< encountered.size() << " " //number of species ever seen
	<< core_pop << " " //number ofindividuals in core 
	<< core_size << " " //number ofspecies in core
	<< E << " " //Life's effect on environment
	<< F //Life's effect on each other
	<< endl;			

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
				network_file << cur->active_population << "\n";		
			}
		} 
	}	
}

//print species in core
inline void print_core(ofstream &core_file, double percent=0.05){
	
	int max = 0; 
	for (list<Node>::iterator cur=species.begin(); cur != species.end(); ++cur){
		if( cur->active_population >  max ){ max = cur->active_population; }
	}
	core_file << tgen << " " << Nactive;
	for (list<Node>::iterator cur=species.begin(); cur != species.end(); ++cur){
		if( cur->active_population > percent * (double) max ){ 
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

	init(seed);
	int t = 0;
	double lgen = Nactive/pkill; //how many birth/deaths per generation

	//put all the parameters in the filename
	string name_tag = "_seed"+ to_string(it);
	name_tag += "_mu"+ to_string(mu); 
	name_tag += "_nu"+ to_string(nu); 
	name_tag += "_sigma" + to_string(sigma); 
	name_tag += "_pmut"+ to_string(pmut); 
	name_tag += "_eta"+ to_string(eta); 
	name_tag += "_S"+ to_string(S); 
	name_tag += "_Sm"+ to_string(Sm); 
	name_tag +=".dat";

	string name1 = path + "popplot"; name1 += name_tag;
	ofstream pop_file; pop_file.open (name1.c_str());

	//string name2 = path + "coreplot"; name2 += name_tag;
	//ofstream core_file; core_file.open (name2.c_str());

	//start iteration
	list<Node>::iterator sa;
	bool active;
	pair<double, double> ss;	
	
	do{
		/////////////////////
		//Individual update//
		/////////////////////
		
		//choose an individual of either type
		pair< list<Node>::iterator, bool > choice = choose();
		sa = choice.first; active = choice.second;

		if( active ){
			sa = kill_active(sa); //try to kill individual
			if(sa != species.end()){ //if I didn't kill it
				ss = calc_HI(sa);		//calculate its environment
				sa = go_dormant(sa, (ss.first/Nactive - ss.second - mu*Nactive - nu*Nactive*Nactive)  ); 	//put it to sleep if conditions are bad
			}
			if(sa != species.end()){ //if is isn't asleep
				if( mt_rand() < poff(ss.first, ss.second) ){	//can it reproduce?			
					asexual(sa);		//reproduce asexually
					//asexual2(sa);		//reproduce asexually, copy twice
				}
			}
			++t; //counter
		} else {
			sa = kill_dormant(sa);
			if(sa != species.end()){ //if it is alive
				ss = calc_HI(sa);		//calculate its environment
				sa = go_active(sa, (ss.first/Nactive - ss.second - mu*Nactive - nu*Nactive*Nactive)  ); 	//wake up if conditions are good
			}
			++t; //should dormant update ticker?
		}
		
		/////////////////////
		//Parallel update//
		/////////////////////
		/*sa = choose_active();
		sa = kill_active(sa); //try to kill individual
		if(sa != species.end()){ //if I didn't kill it
			ss = calc_HI(sa);		//calculate its environment
			sa = go_dormant(sa, (ss.first/Nactive - ss.second - mu*Nactive - nu*Nactive*Nactive)  ); 	//put it to sleep if conditions are bad
		}
		if(sa != species.end()){ //if is isn't asleep
			if( mt_rand() < poff(ss.first, ss.second) ){	//can it reproduce?			
				asexual(sa);		//reproduce asexually
				//asexual2(sa);		//reproduce asexually, copy twice
			}
		}
		sa = choose_dormant();
		sa = kill_dormant(sa);
		if(sa != species.end()){ //if it is alive
			ss = calc_HI(sa);		//calculate its environment
			sa = go_active(sa, (ss.first/Nactive - ss.second - mu*Nactive - nu*Nactive*Nactive)  ); 	//wake up if conditions are good
		}
		++t; */
		
		
		if(t >= lgen){	//generation is over
			t = 0; ++tgen; lgen = Nactive/pkill; //recalculate generation length. Should this be Nactive+Ndormant?
			print_stats(pop_file);
			//print_core(core_file);

			if(Nactive > 1e5){		//annoying. Catches population explosions that make the simulation impossible to run.
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


