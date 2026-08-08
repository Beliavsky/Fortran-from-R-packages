#include <stdio.h>
#include <stdlib.h>


#define FALSE 0
#define TRUE 1

#ifndef _MULTI_MCGA_H_
#define _MULTI_MCGA_H_ 1
#endif

#ifndef _MCGATYPES_H_
#define _MCGATYPES_H_ 1
#include "mcgatypes.h"
#endif


/*
 * Holds information for a Chromosome.
 * Note that the variable genes is type of double* and it is
 * created dynamically.
 * */
struct MultiChromosome
{
	int ch_size;
	double *genes;
	double *cost;
	int numfunc;
	double rank;
};


/*
Holds information for a Population.
*/
struct MultiMcga
{
    int pop_size;	                        /* Population size */
    int ch_size;    	                    /* Chromosome size */
    struct MultiChromosome *chromosomes;    /* Array of chromosomes */
    double mutate_prob;                     /* Mutation probability */
    double cross_prob;                	    /* Cross over probability */
    int cross_type;                    	    /* Cross over type */
    int elitism;                            /* Number of elitisim */
    int max_iter;                           /* Number of maximum iterations */
    int current_iter;                       /* Current iteration while the process is being performed */
    double *(*cost_function)(struct MultiChromosome*);    /* Pointer to the function subject to be minimized. */
    int numfunc;			    /* Number of objective functions */
};

/*
Constructs a new 'struct MultiMcga' with given parameters.
popsize     : Population size
chsize      : Chromosome size
crossprob   : Cross over probability
mutateprob  : Mutation probability
crosstype   : Cross over type
elitism     : Number of chromosomes to be passed to next generation
cost_function:  Function pointer subject to be minimized. A cost function has to be defined as

            double *a_simple_cost_function (Chromosome *c)

and the genes of chromosome c is accessible via field 'genes'. A sample_function is given in the main.c file.
*/
struct MultiMcga *multi_mcga_create (int popsize, int chsize, double crossprob, double mutateprob, int crosstype, int elitism, int numfunc, double *(*cost_function)(struct MultiChromosome*));

/*
Clones the structure mcga.
Note that only a new structure is created but chromosomes are not copied into new one.
*/
struct MultiMcga *multi_mcga_clone (struct MultiMcga *mcga);

/*
Frees the memory that used by struct mcga.
*/
void multi_mcga_free (struct MultiMcga *mcga);

/*
Dumps the given chromosome as a byte array.
Standard input is used as output device
*/
void multi_mcga_dump_chromosome_chars(struct MultiChromosome *c);

/*
Dumps the given chromosome as real values.
Standard input is used as output device
*/
void multi_mcga_dump_chromosome (struct MultiChromosome *c, int idx);

/*
Dumps the given population as real values.
Standard input is used as output device
*/
void multi_mcga_dump (struct MultiMcga *mcga);

/*
Randomizes the given population with the distribution of uniform(min,max)
*/
void multi_mcga_randomize(struct MultiMcga *mcga, double min, double max);

/*
Sets the cost of 'idx'th chromosome of population mcga
*/
void multi_mcga_set_cost (struct MultiMcga *mcga, int idx, double *val);

/*
Getter method for the cost of 'idx'th chromosome of population mcga
*/
double *multi_mcga_get_cost (struct MultiMcga *mcga, int idx);

/*
Calculate the cost values using the function pointer that given in 'struct mcga'
This operation is done for entire population.
*/
void multi_mcga_calculate_fitness(struct MultiMcga* mcga);

/*
Mutation operation in mcga. Mutates a chromosome with probability 'prob'
If u < prob then the given chromosome 'c' is mutated, where u is a random value
from Uniform(0,1). Mutation is done using the byte representation of the chromosome 'c'
in the memory.
*/
void multi_mcga_mutate(struct MultiChromosome *c,  double prob);

/*
Performs one point cross over operation on the given chromosomes 'c1' and 'c2'.
This operation is done using the byte representations of the chromosomes. The logic
underlying of this operation is same as the classical cross over operation.
*/
void multi_mcga_onepoint_crossover (struct MultiChromosome *c1, struct MultiChromosome *c2);

/*
Performs uniform cross over operation on the given chromosomes 'c1' and 'c2'.
This operation is done using the byte representations of the chromosomes. The logic
underlying of this operation is same as the classical cross over operation.
*/
void multi_mcga_uniform_crossover (struct MultiChromosome *c1, struct MultiChromosome *c2);

/*
Comparetor function for sorting operations.
*/
int multi_mcga_sorter(const void *c1, const void *c2);

/*
Sorts the given population by the cost values of chromosomes in ascending order.
*/
void multi_mcga_sortpopulation(struct MultiMcga *mcga);

/*
 * Sorts population by the ranks of chromosomes using non-dominated sorting algorithm
*/
void multi_mcga_calculate_ranks(struct MultiMcga *mcga);

/*
In each generation of mcga, selection operation is done in a temporary population.
This method copies temporary population to current population.
*/
void multi_mcga_swap_population(struct MultiMcga *source, struct MultiMcga *target);

/*
Performs tournament selection. In each iteration two chromosomes are selected. The the best one is labeled as 'Parent 1'.
Then two chromosomes are selected again and the best one is labeled as 'Parent 2'. Then these parents are crossed over
and mutated with respective probabilities.
*/
void multi_mcga_tournament_selection (struct MultiMcga *source, struct MultiMcga *target);

/*
Combines those procedures in order to make MCGA compact.
mcga: Main population
temppop: Temporary population
maxiter: Number of maximum generations
min and max: Parameters of a uniform distribution for randomized initial population.
In the file main.c, there is a simple usage of mcga_start().
*/
void multi_mcga_start(struct MultiMcga *mcga, struct MultiMcga *temppop, int maxiter, double min, double max);

