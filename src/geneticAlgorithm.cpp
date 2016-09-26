/**
*	Larissa Fernandes Leijôto <larissa.leijoto@gmail.com>
*	
* 	This file is part of genetic algorithm developed to look for matches
* 	in different proteins.
*/
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <time.h>
#include <cmath>
#include <assert.h>
#include <string>
#include <vector>
#include <map>
#include <iostream>     // std::cout
#include <algorithm>    // std::shuffle
#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock
#include <cfloat>
#include <limits>
#include "util.h"
#include "geneticAlgorithm.h"
#include "database.h"

using namespace std;


#define dbg 0
#define time_program 0
#define time_ga 0
#define saveFile 0


/*
 * Sanity check for genetic_algorithm().
 */
#define SANITY_CHECK()         \
	assert(popsize > 0);       \
	assert(ngen > 0);          \
	assert(g->m_rate >= 0.0);  \
	assert(g->m_rate <= 1.0);  \
	assert(g->c_rate >= 0.0);  \
	assert(g->c_rate <= 1.0);  \
	assert(g->e_rate > 0.0);   \
	assert(g->e_rate < 1.0);   \

/*============================================================================*
 *                              genetic Operators                             *
 *============================================================================*/
/**
 * @brief Creates a Individual.
 * 
 * @returns A new Individual.
 */
 
Individual *individual_create(void)
{
	// Individual *g = (Individual *)smalloc(sizeof(Individual));
	Individual *g = new Individual;

    g->aminoacids_chain = (int **)smalloc(nproteins*sizeof(int *));

    for(int i = 0; i < nproteins; i++)
		 g->aminoacids_chain[i] = (int *)smalloc(database.minaminoacids * sizeof(int));
	
	// g->objectives = (double *)smalloc(numberOfObjectives*sizeof(double));
		
	g->rank = numeric_limits<int>::max();
	g->distance = 0.0;
	g->dom_count = 0;

	return (g);
}
	

/**
 * @brief Destroys a Individual.
 */
 /*
void individual_destroy(Individual *individual)
{
	cout<<"TESTE 1"<<endl;
	for(int i = 0; i < nproteins;i++ )
		free(individual->aminoacids_chain[i]);
		
	cout<<"TESTE 2"<<endl;
	free(individual->objectives);
	cout<<"TESTE 3"<<endl;
	free(individual->aminoacids_chain);
	cout<<"TESTE 4"<<endl;
	free(individual);
}
*/


int recoverBLosumScore(std::string AA1, std::string AA2)
{
	
	return -1;
}

	
int recoverClesumScore(std::string AA1, std::string AA2)
{
	
	return -1;
	
}

/**
 * @brief creates a random Individual.
 * 
 * @returns A random Individual.
 */
Individual *individual_random(void)
{	
	int i = 0, k = 0;
	vector <int> aminoacids_list;
	
	#if(dbg>0)
	fprintf(stderr, "info: Individual random\n");
	#endif
	
	Individual *g = individual_create();
	
    /* create Individual.*/
    for (i = 0; i < nproteins; i++)
    {	
		for( k = 0; k < database.naminoacids[i]; k++)
			aminoacids_list.push_back(k);		
		
		 std::random_shuffle(aminoacids_list.begin(), aminoacids_list.end());
			
		for( k = 0; k < database.minaminoacids;k++)
		{	
			g->aminoacids_chain[i][k] = aminoacids_list.back();
			aminoacids_list.pop_back();
		}
		aminoacids_list.clear();		
	} 

	return (g);
}


/**
 * @brief Asserts if a Individual has a feature.
 * 
 * @returns One if the Individual has the feature, and zero otherwise.
 */
int has_feature(int *g, int feature, int n)
{	
	int resp = 0;

	for (int i = 0; i < n; i++)
	{
		if (g[i] == feature)
		{
			resp = 1;
			i = n;
		}
	}
		return resp;
}

/**
* @brief Funtion used to make a copy of a Individual..
* 
* @param Result individual offspring
* @param Individual1 First Individual.
* @param Individual2 Second Individual.
*/
/*
Individual *copy_individual(Individual *individual)
{
	Individual *g = individual_create();

	for(int wprotein = 0; wprotein < nproteins;wprotein++)
		memcpy(g->aminoacids_chain[wprotein],individual->aminoacids_chain[wprotein], database.minaminoacids*sizeof(int));

	memcpy(g->objectives, individual->objectives, numberOfObjectives*sizeof(double));

	g->rank = individual->rank;
	g->distance = individual->distance;
		
	return g;
}
*/

/**
 * @brief Crossovers two amino acids chains.
 * 
 * @param Result individual offspring
 * @param Individual1 First Individual.
 * @param Individual2 Second Individual.
 */
void aminoacids_crossover(Individual *offspring, Individual *individual1, Individual *individual2)
{	
	int point1 = 0;    		 	 						/* First crossover point.   */
	int point2 = 0;     	 		 					/* Second crsossover point. */
	int nbegin, nmiddle, nend; 	 						/* Size of gene parts.      */
	int i, j; 						 					/* Loop index.				*/
	int *begin, *middle, *end;		 					/* Gene parts.              */
	int **map = (int **)smalloc(2*sizeof(int *));
	bool search;

	/* Sanity check. */
	assert(n >= 0);
		
	/* Generate crossover points. */
	int tmp;
	do
	{
		point1 = rand()%database.minaminoacids;
		point2 = rand()%database.minaminoacids;
		
	} while (((point1<=1) || (abs(point1-point2)<=1) || (point2<=1)|| ((database.minaminoacids-point2)<=1)));

	if (point1 > point2)
	{
		tmp = point1;
		point1 = point2;
		point2 = tmp;
	}
	/* Size of gene parts. */
	nbegin = point1 - 1 ;
	nmiddle = (point2 - 1) - nbegin;
	nend = database.minaminoacids - (nbegin+nmiddle);

	for(int wprotein = 0; wprotein < nproteins; wprotein++)
	{		
		/* Gene parts. */
		begin  = (int *)smalloc(nbegin*sizeof(int));
		middle = (int *)smalloc(nmiddle*sizeof(int));
		end    = (int *)smalloc(nend*sizeof(int));
		
		map[0] = (int *)smalloc(nmiddle*sizeof(int));
		map[1] = (int *)smalloc(nmiddle*sizeof(int));
		
		memcpy(map[0],individual1->aminoacids_chain[wprotein]+ nbegin , nmiddle*sizeof(int));
		memcpy(map[1],individual2->aminoacids_chain[wprotein] + nbegin , nmiddle*sizeof(int));
	
		
		memcpy(begin, individual1->aminoacids_chain[wprotein], nbegin *sizeof(int));
		memcpy(middle,individual2->aminoacids_chain[wprotein] + nbegin , nmiddle*sizeof(int));
		memcpy(end,   individual1->aminoacids_chain[wprotein] + nbegin  + nmiddle , nend*sizeof(int));
		
		do
		{
			search = false;
			
			for ( i = 0; i < nmiddle; i++)
			{
				for ( j = 0; j < nbegin; j++)
				{
					if (map[1][i] == begin[j])
					{					
						begin[j] = map[0][i];
						search = true;
						break;				
					}
				}			
			
				for ( j = 0; j < nend; j++)
				{
					if (map[1][i] == end[j])
					{		
						end[j] = map[0][i];				
						search = true;
						break;
					}
				}
			}				
		}while(search);			
		
				
		memcpy(offspring->aminoacids_chain[wprotein], begin, nbegin*sizeof(int)); 
		memcpy(offspring->aminoacids_chain[wprotein] + nbegin, middle, nmiddle*sizeof(int));
		memcpy(offspring->aminoacids_chain[wprotein] + nbegin + nmiddle, end, nend*sizeof(int));

		/* House keeping. */
		free(begin);
		free(middle);
		free(end);
		free(map[0]);
		free(map[1]);
	}
	free(map);	
}


/**
 * @brief Crossovers two Individuals.
 * 
 * @param Individual1 First Individual.
 * @param Individual2 Second Individual.
 */
void individual_crossover(Individual *offspring, Individual *mom, Individual *dad)
{
	#if(dbg>0)
	fprintf(stderr, "info: Individual crossover...\n");
	#endif
	
	aminoacids_crossover(offspring, mom,  dad);
}

/**
 * @brief Mutates a Individual.
 * 
 * @returns The mutated Individual.
 */
void individual_mutation(Individual *individual)
{	
	#if(dbg>0)
	fprintf(stderr, "info: Individual mutation...\n");
	#endif	
	
	int i, j;       		/* Mutation point. */
	int aminoacid; 	/* Feature.        */
	int tmp;
	
	for(int wprotein = 0; wprotein < nproteins; wprotein++)
	{
		if(database.naminoacids[wprotein] != database.minaminoacids)
		{
			i = rand()%database.minaminoacids;
			do
			{
				aminoacid = rand()%database.naminoacids[wprotein];
			} while (has_feature(individual->aminoacids_chain[wprotein], aminoacid, database.minaminoacids));

			individual->aminoacids_chain[wprotein][i] = aminoacid;
		}
		else
		{
			i =  rand()%database.minaminoacids;
			j =  rand()%database.minaminoacids;
			
			tmp = individual->aminoacids_chain[wprotein][i];
			individual->aminoacids_chain[wprotein][i] = individual->aminoacids_chain[wprotein][j];
			individual->aminoacids_chain[wprotein][j] = tmp;
			
		}
 
	}
}

double distanceCalculate(Individual *i1)
{
	double x, y, z, dist = 0;
	
	for(int wAminoacid = 0; wAminoacid < database.minaminoacids; wAminoacid++)
	{
		
		int first_aminoacid  = i1->aminoacids_chain[ 0 ][ wAminoacid ];
		int second_aminoacid = i1->aminoacids_chain[ 1 ][ wAminoacid ];
		
		x = database.data_aminoacids[0][first_aminoacid]->listAtoms[1]->x - database.data_aminoacids[1][second_aminoacid]->listAtoms[1]->x;
		y = database.data_aminoacids[0][first_aminoacid]->listAtoms[1]->y - database.data_aminoacids[1][second_aminoacid]->listAtoms[1]->y;
		z = database.data_aminoacids[0][first_aminoacid]->listAtoms[1]->z - database.data_aminoacids[1][second_aminoacid]->listAtoms[1]->z;

		dist += pow(x,2)+ pow(y,2) + pow(z,2); /* Calculating distance by euclidean formula. */
	}             

    return dist;
}


/**
 * @brief Evaluates the fitness of a Individual.
 */
void individual_evaluate(Individual *ind)
{
	
	#if(dbg>0)
	fprintf(stderr, "info: Individual evaluate...\n");
	#endif	
	double matches_aminoacids = 0.000001; /*	Quantity of matches in the amino acids chain.	*/
	

		for(int wProtein1 = 0; wProtein1 < nproteins; wProtein1++)
		{
			for(int wProtein2 = 0; wProtein2 < wProtein1; wProtein2++)
			{	
				for(int wAminoacid = 0; wAminoacid < database.minaminoacids; wAminoacid++)
				{
					int first_aminoacid  = ind->aminoacids_chain[ wProtein1 ][ wAminoacid ];
					int second_aminoacid = ind->aminoacids_chain[ wProtein2 ][ wAminoacid ];
										
					if (strcmp( (database.data_aminoacids[wProtein1][ first_aminoacid ]->aminoacids_aa).c_str(),
						(database.data_aminoacids[wProtein2][ second_aminoacid ]->aminoacids_aa).c_str()) == 0)
					{
						matches_aminoacids++;
					}
						
				}
		
			}
		}
			
	// Initializing objectives
	ind->objectives.push_back(1/matches_aminoacids);
	ind->objectives.push_back(1/(float)
	(1 + (distanceCalculate(ind)/(float)(pow(1.24*pow(((database.minaminoacids -15) -1.8), 1.0/3.0), 2)))));

}



/*============================================================================*
 *                              genetic Definitions                           *
 *============================================================================*/
/*
*  Configuration genome.
*/
struct genome problem = 
{
		0.05,           /* Mutation rate.    */
		0.80,           /* Crossover rate.   */
		0.01,           /* Elitism rate.     */
		1.00           	/* Replacement rate. */
};

/*
 * genetic algorithm for proteins's matches.
 */
void protein_matches(int popsize, int ngen)
{
	geneticAlgorithm(&problem, popsize, ngen);
} 

/*============================================================================*
 *                              genetic utilities                             *
 *============================================================================*/
void print_individual(Individual *individual)
{
	fprintf(stderr,"\n");
		
		for(int j = 0; j < nproteins; j++)
		{		
			for(int k = 0; k < database.minaminoacids;k++)
			{
				//fprintf(stderr,"%d ",individual->aminoacids_chain[j][k]);
				fprintf(stderr, "%s ",((database.data_aminoacids[j][individual->aminoacids_chain[j][k]])->aminoacids_aa).c_str());
			}
			fprintf(stderr,"\n");
		}
		
		fprintf(stderr,"\n");
		for(int i = 0; i < numberOfObjectives; i++)
			fprintf(stderr,"%f ",individual->objectives[i]);
		
		fprintf(stderr,"%f ",individual->distance);	
  		fprintf(stderr,"%d ",individual->rank);
		fprintf(stderr,"\n");
	
}

bool Individual::dominates(Individual *individual)
{
	for (int i = 0; i < numberOfObjectives; i++) 
	{
		if (this->objectives[i] > individual->objectives[i])
			return false;
	}
	return true;
}

vector<vector<Individual *>> fastNondominatedSort(vector<Individual *> &population) 
{
	/* Discover Pareto fronts in P, based on non-domination criterion. */
	vector<vector<Individual*>> fronts(1);

	// Building initial frontier
	for (Individual *p1 : population)
	{
		p1->dom_count = 0;
		p1->dom_set.clear();

		for (Individual *p2 : population)
		{
			
			if (p1 == p2) continue;

			if (p1->dominates(p2)) p1->dom_set.push_back(p2);
			else if (p2->dominates(p1)) p1->dom_count += 1;

		}

		if (p1->dom_count == 0)
		{
			p1->rank = 0;
			fronts[0].push_back(p1);
		}

	}

	int curr = 0;
	
	do {

		const auto& currFront(fronts[curr]);
		std::vector<Individual *> nextFront;

		for (Individual* p1 : currFront) 
		{
			for (Individual* p2 : p1->dom_set) 
			{
				p2->dom_count -= 1;
				
				if (p2->dom_count == 0) 
				{
					p2->rank = (curr+1);
					nextFront.push_back(p2);
				}

			}
		}

		curr += 1;
		
		if(not nextFront.empty())
			fronts.push_back(nextFront);
			
	} while (curr < fronts.size());

	return fronts;

}

void crowdingDistanceAssignment(vector<Individual *> &front)
{
	/* Assign a crowding distance for each solution in the front. */
	for (int i = 0; i < (signed)front.size(); i++) 
		front.at(i)->distance = 0;	

	for (int i = 0; i < numberOfObjectives; i++) 
	{
		
		sortObjective(front, i);

		double rge = front[0]->objectives[i] - front.back()->objectives[i];

		front[0]->distance = numeric_limits<double>::infinity();
		front.back()->distance = numeric_limits<double>::infinity();

		if (rge == 0.0) continue;

		for (int j = 1; j < (signed)(front.size() - 1); j++) 
			front.at(j)->distance += ((front.at(j + 1)->objectives[i] - front.at(j - 1)->objectives[i])/rge);

	}

}

void sortObjective(vector<Individual *> &front, int objective_index)
{
	for (int i = (signed)(front.size() - 1); i >= 0; i--) 
	{
		for (int j = 1; j < i; j++) 
		{
			if (front[j-1]->objectives[objective_index] > front[j]->objectives[objective_index]) 
			{
				Individual *temp = front.at(j - 1);
				front.at(j - 1) = front.at(j);
				front.at(j) = temp;
			}
		}
	}
}

void sortCrowding(vector<Individual *> &P) 
{
	for (int i = (signed)(P.size() - 1); i >= 0; i--) 
	{
		for (int j = 1; j < i; j++)
		{
			if (crowdedComparison(P.at(j - 1), P.at(j)) < 0) 
			{
				Individual *temp = P.at(j - 1);
				P.at(j - 1) = P.at(j);
				P.at(j) = temp;
			}
		}
	}
}

int crowdedComparison(const Individual *s1, const Individual *s2) 
{
	if (s1->rank < s2->rank)
	 	return 1;
	else if (s1->rank > s2->rank) 
		return -1;
	else if (s1->distance > s2->distance) 
		return 1;
	else if (s1->distance < s2->distance) 
		return -1;
	else 
		return 0;
	
}

vector<Individual *> reproduce(genome *g, vector<Individual *> &population)
{
	vector<Individual *> children;
	srand (time(NULL));
	/* Make new population Q, offspring of P. */
	while (children.size() < population.size()) 
	{
		if (((double) rand() / RAND_MAX) < g->c_rate)
		{
			Individual *childSolution = individual_create();
			int r1 = rand()%popsize;
			int r2 = rand()%popsize;
			
			individual_crossover(childSolution, population[r1], population[r2]);
			
			if (((double) rand() / RAND_MAX) < g->m_rate) 
				individual_mutation(childSolution);
		
			children.push_back(childSolution);
		}
	}
		
	return children;
}

Individual *better( Individual *one, Individual *two)
{
	if(one->distance != numeric_limits<double>::infinity() and one->rank == two->rank)
	{
		if(one->distance > two->distance)
			return one;
		else
			return two;
	}
	else
	{
		if(one->rank < two->rank)
			return one;
		else
			return two;
	}

}

vector<Individual *> select_parents(vector<vector<Individual *>> &fronts)
{	

	for (auto& front : fronts) {
		assert(not front.empty() and "empty front! =/");
		crowdingDistanceAssignment(front);
	}

	int last_front = 0;
	vector<Individual *> offspring;

	for (auto& front : fronts) 
	{

		if ((offspring.size() + front.size()) > popsize)
			break;

		offspring.insert(offspring.end(), front.begin(), front.end());
		++last_front;
		
	}

	int remaining = popsize - offspring.size();
	if (remaining > 0 and last_front < fronts.size()) 
	{		

		assert(last_front < front.size() and "Só Jesus na causa =/");

		auto& front = fronts[last_front];
		sortCrowding(front);	

		offspring.insert(offspring.end(), front.begin(), front.begin() + remaining);

	}
	
	
	return offspring;

}
/*============================================================================*
 *                              genetic Algoritmh                             *
 *============================================================================*/
void geneticAlgorithm(genome *g, int popsize, int ngen)
{

	srand (time(NULL));
	fprintf(stderr,"popsize: %d\n", popsize);
	fprintf(stderr,"ngen: %d\n", ngen);
	fprintf(stderr,"Crossover rate: %.2f\n", g->c_rate);
	fprintf(stderr,"Elitism rate: %.2f\n", g->e_rate);
	fprintf(stderr,"Repalcement rate: %.2f\n", g->r_rate);
	fprintf(stderr,"Mutation rate: %.2f\n", g->m_rate);
	
	vector<Individual *> population;					/* Current population.   */
	vector<Individual *> selected;						/* Population selected.  */
	vector<Individual *> children;						/* Children created.	 */
	vector<Individual *> parents;						/**/
	
	SANITY_CHECK();
	
		/* Build initial population. */
		for(int i = 0 ; i < popsize; i++)		
			population.push_back(individual_random());

		for(int i = 0 ; i < popsize; i++)
			individual_evaluate(population[i]);
			
		// ---------------------------------------------------------------
		fastNondominatedSort(population);

		while(selected.size() < population.size())
			selected.push_back(better(population[rand()%population.size()],population[rand()%population.size()]));

		children = reproduce(g, selected);
				
		for(int i = 0 ; i < popsize; i++)
			individual_evaluate(children[i]);
					
		for(int count = 0; count < ngen; count++)
		{	
			fprintf(stderr,"generation: %d\n", count);

			vector<Individual *> Union;
			/* Union between Population and Children to contruct Union. */
			Union.reserve(population.size() + children.size());
			Union.insert(Union.end(), population.begin(), population.end());
			Union.insert(Union.end(), children.begin(), children.end());

			vector<vector<Individual*>> fronts = fastNondominatedSort(Union);
			parents = select_parents(fronts);

			while(selected.size() < parents.size())
				selected.push_back(better(parents[rand()%parents.size()], parents[rand()%parents.size()]));

			population = children;
			children = reproduce(g, selected);

			for(int i = 0 ; i < (signed)children.size(); i++)
				individual_evaluate(children[i]);

		}	

		vector<Individual *> Union;
		Union.reserve(population.size() + children.size());
		Union.insert(Union.end(), population.begin(), population.end());
		Union.insert(Union.end(), children.begin(), children.end());
		
		vector<vector<Individual*>>fronts = fastNondominatedSort(Union);

		parents = select_parents(fronts);

		for(int i = 0 ; i < (signed)parents.size(); i++)
			print_individual(parents[i]);

}
