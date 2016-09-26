/**
*	Larissa Fernandes Leijôto <larissa.leijoto@gmail.com>
*	
* 	This file is part of genetic algorithm developed to look for matches
* 	in different proteins.
*/

#ifndef _GENETIC_ALGORITHM_H_
#define _GENETIC_ALGORITHM_H_
#include <vector>

	/**
	 *  Gene.
	 */
	struct Individual
	{
		int **aminoacids_chain;    		 /* Chains of amino acids.  */
		int rank;
		double distance;
		// double *objectives;
		std::vector<double> objectives;
			
		int dom_count;
		std::vector<Individual *> dom_set;

		bool dominates(Individual *individual);
	
	};

	/*
	 * Genome.
	 */
	struct genome
	{
		/* Attributes. */
		double m_rate; /* Mutation rate.    */
		double c_rate; /* Crossover rate.   */
		double e_rate; /* Elitism rate.     */
		double r_rate; /* Replacement rate. */
	};
	
	/* 
	* Creates a Individual.
	*/
	Individual *individual_create(void);
		
	/*
	 * Destroys an individual.
	 */
	void individual_destroy(Individual *individual);
	
	/*
	 * Mates two individuals.
	 */
	void individual_crossover(Individual *offspring, Individual *mom, Individual *dad);
	
	/*
	 * Mates two amino acids chain.
	 */
	void aminoacids_crossover(Individual *offspring,  Individual *mom, Individual *dad);
	
		/*
	 * Mutates a individual.
	 */
	void individual_mutation(Individual *individual);
	
	/*
	 * Generates a random individual.
	 */
	Individual *individual_random(void);

	/**
	*  One if the Individual has the feature, and zero otherwise.
	*/
	int has_feature(int *g, int feature, int n);
	
	/*
	 * Genetic algorithm.
	 */
	void geneticAlgorithm(genome *g, int popsize, int ngen);

	/**
	 * Call genetic algorithm for protein matches
	 */
	void protein_matches(int popsize, int ngen);
	
	/**
	 * Prints the individual and his respective amino acids and atoms.
	 */
	void print_individual(Individual *individual);
	
	/**
	 * Evaluates each individual 
	 */
	void individual_evaluate(std::vector<struct Individual> &population);
			
	/**
	 * Funtion used to make a copy of a Individual.
	 */
	Individual *copy_individual(Individual *individual);
	
	int recoverBLosumScore(std::string AA1, std::string AA2);
	
	int recoverClesumScore(std::string AA1, std::string AA2);
	
	std::vector<std::vector<Individual *>> fastNondominatedSort(std::vector<Individual *> &R);
	
	void crowdingDistanceAssignment(std::vector<Individual *> &front);
	
	void sortObjective(std::vector<Individual *> &front, int objective_indexi);
	
	void sortCrowding(std::vector<Individual *> &P);
	
	int crowdedComparison(const Individual *s1, const Individual *s2);
	
	std::vector<Individual *> reproduce(genome *g, std::vector<Individual *> &P);
	
	double distanceCalculate(Individual *individual);

	Individual *better(Individual *one, Individual *two);
	
	std::vector<Individual *> select_parents(std::vector<std::vector<Individual *>> &fronts);


	/* Global parameters. */
	extern int gen;       			 /* Number of generations elapsed. */
	extern int popsize;   			 /* Population size.               */
	extern genome *g; 					 /* Genome.                        */
	extern int numberOfObjectives;


#endif /* _GENETIC_ALGORITHM_H_ */
