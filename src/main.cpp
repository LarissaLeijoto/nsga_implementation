/**
*	Larissa Fernandes Leijôto <larissa.leijoto@gmail.com>
*	
* 	This file is part of genetic algorithm developed to look for matches
* 	in different proteins.
* 
* 	
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <time.h>
#include <omp.h>
#include <time.h>
#include "database.h"
#include "util.h"

using namespace std;

	
#define dbg  0
#define time_program 0

/* Program parameters. */
static const char **filenames = NULL; /* Name of input files.               */
int nproteins = 0;               	  /* Number of proteins (input files).  */
int popsize = 0;		         	  /* Population size.                   */
int ngen = 0;          		 		  /* Number of generations.             */
int numberOfObjectives = 0;			  /* Number of external individuals.	*/
/**
 * @brief Database.
 */
struct database database;


/**
 * @brief Prints program usage and exits.
 * 
 * @details Prints program usage and exits.
 */
static void usage(void)
{
	printf("Usage: genetic algorithm <popsize> <ngen> <numberOfObjectivies> <protein files>");
	exit(EXIT_SUCCESS);
}

/**
 * @brief Reads command line arguments.
 * 
 * @details Reads and parses command line arguments pointed to by @p argv.
 * 
 * @todo Read and parse command line arguments.
 */
static void readargs(int argc, char **argv)
{
	
	/* Missing arguments. */
	if (argc < 4)
		usage();
	
	popsize = atoi(argv[1]);
	ngen = atoi(argv[2]);
	numberOfObjectives = atoi(argv[3]);
	
	
	/* Count the number of proteins. */
	for (int i = 4; argv[i] != NULL; i++)
		nproteins++;
		
	filenames = (const char **)smalloc(nproteins*sizeof(char *));
	
	/* Extract protein files. */
	for (int i = 0; i < nproteins; i++)
		filenames[i] = argv[4 + i];
	
	/* Assert program parameters. */
	 if (popsize == 0)
		error("invalid population size");
	else if (ngen == 0)
		error("invalid number of generations");
	else if(numberOfObjectives < 1)
		error("invalid number of objectives");
}

int main(int argc, char **argv)
{
	clock_t start_time = 0;
	start_time = clock();
	srand( (unsigned)time(NULL) );
	
	//omp_set_num_threads(4);

	readargs(argc, argv);
	
	#if(dbg>0)
	fprintf(stderr, "info: parsing database! Ok \n");
	#endif
	/* Parse database in order to determine the largest number of amino acids among all proteins.*/
	database_parse(filenames, nproteins);
	
	#if(dbg>0)
	fprintf(stderr, "info: reading database! Ok\n");
	#endif
	/* Read database. */
	database_read(filenames, nproteins);
	
	//readMatrix();
	
	#if(dbg>0)
	fprintf(stderr, "info: proteins matches...\n");
	#endif	
	protein_matches(popsize,ngen);
	
	//print_base();
	
	/* House keeping. */
	//database_destroy();
	
	double time_in_seconds = (clock() - start_time) / (double)CLOCKS_PER_SEC;
	fprintf(stderr, "program time main: %.2f\n",time_in_seconds );
	
	
	free(filenames);
	
	return (EXIT_SUCCESS);
}
