#include "multi_mcga.h"
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

struct MultiMcga *multi_mcga_create (int popsize, int chsize, double crossprob, double mutateprob, int crosstype, int elitism, int numfunc, double *(*cost_function)(struct MultiChromosome*)){
    int i;
    struct MultiMcga *mcga = malloc (sizeof(struct MultiMcga));
    mcga->pop_size=popsize;
    mcga->ch_size = chsize;
    mcga->cross_type=crosstype;
    mcga->cross_prob = crossprob;
    mcga->mutate_prob = mutateprob;
    mcga->elitism = elitism;
    mcga->numfunc = numfunc;

    mcga->chromosomes = malloc(sizeof(struct MultiChromosome) * popsize);
    

    for (i=0; i<popsize; i++)
    {
        mcga->chromosomes[i].ch_size = chsize;
        mcga->chromosomes[i].genes = malloc (sizeof(double) * chsize);
        mcga->chromosomes[i].cost = malloc(sizeof(double) * numfunc);
	mcga->chromosomes[i].numfunc = numfunc;
	mcga->chromosomes[i].rank = 0.0;
    }

    mcga->cost_function = cost_function;
    return(mcga);
}

struct MultiMcga *multi_mcga_clone (struct MultiMcga *mcga)
{
    struct MultiMcga *mcga2 = multi_mcga_create(mcga->pop_size, mcga->ch_size, mcga->cross_prob, mcga->mutate_prob, mcga->cross_type, mcga->elitism, mcga->numfunc, mcga->cost_function);
    return(mcga2);
}


void multi_mcga_free (struct MultiMcga *mcga){
	int i;
    for (i=0;i<mcga->pop_size;i++){
		free(mcga->chromosomes[i].genes);
		//free(mcga->chromosomes[i].cost);
    }
	free(mcga);
}


void multi_mcga_dump_chromosome_chars(struct MultiChromosome *c)
{
    int i;
    unsigned char *ch = (unsigned char*)c->genes;
    for (i=0; i< ((c->ch_size) * sizeof(double)); i++)
    {
        Rprintf("%d ", *ch);
        ch++;
    }
    Rprintf("\n");
}

void multi_mcga_dump_chromosome (struct MultiChromosome *c, int idx)
{
    int j;
	Rprintf("%d) [R%d]",idx,(int)c->rank);
	for (j=0;j<c->numfunc;j++){
		Rprintf("[%f] ", c->cost[j]);
	}
    for(j=0; j<c->ch_size; j++)
    {
        Rprintf("%.5f ",c->genes[j]);
    }
    Rprintf("%s", "\n");
}


void multi_mcga_dump (struct MultiMcga *mcga)
{
    int i;
    for (i =0; i<mcga->pop_size; i++)
    {
        multi_mcga_dump_chromosome(&mcga->chromosomes[i], i);
	Rprintf("\n");
    }
}


void multi_mcga_randomize(struct MultiMcga *mcga, double min, double max)
{
    int i,j,h;
    for (i =0; i<mcga->pop_size; i++)
    {
        for (j=0; j<mcga->ch_size; j++)
        {
            //mcga->chromosomes[i].genes[j] = min + ((double)rand() / RAND_MAX) * (max-min);
            mcga->chromosomes[i].genes[j] = runif(min,max);
            for (h=0;h<mcga->numfunc;h++){
				mcga->chromosomes[i].cost[h] = 0.0;
			}
        }
    }
}


void multi_mcga_set_cost (struct MultiMcga *mcga, int idx, double *val)
{
	int i;
	for (i = 0; i< mcga->chromosomes[idx].numfunc; i++){
	    mcga->chromosomes[idx].cost[i] = val[i];
	}
}


double *multi_mcga_get_cost (struct MultiMcga *mcga, int idx)
{
    return(mcga->chromosomes[idx].cost);
}


void multi_mcga_calculate_fitness(struct MultiMcga* mcga)
{
    int i;
    for (i =0; i<mcga->pop_size; i++)
    {
	mcga->chromosomes[i].cost = mcga->cost_function(&mcga->chromosomes[i]);
    }
}

void multi_mcga_mutate(struct MultiChromosome *c,  double prob)
{
    int i;
    unsigned char *pch = (unsigned char *) c->genes;
    for (i = 0 ; i< (sizeof(double)* c->ch_size); i++)
    {
        /* will mutate current byte? */
        //if ( ((double)rand() / RAND_MAX)<prob)
        if ( runif(0,1) < prob)
        {
            /* Mutate upwards or downwards ?*/
            //int amount = (int)(((double)rand() / RAND_MAX) * 255);
            //if ( ((double)rand() / RAND_MAX)< 0.5)
            if ( runif(0,1) < 0.5)
            {
                (*pch)+= 1; //amount;
            }
            else
            {
                (*pch)-= 1; //amount;
            }
        }
        pch++;
    }
}


void multi_mcga_onepoint_crossover (struct MultiChromosome *c1, struct MultiChromosome *c2)
{
    unsigned char *pch1 = (unsigned char*) c1->genes;
    unsigned char *pch2 = (unsigned char*) c2->genes;
    unsigned char tmpch;
    //int cutpoint = (int)(((double)rand() / RAND_MAX) * c1->ch_size);
    int cutpoint = (int)(runif(0,1) * c1->ch_size);
    int i;
    pch1 += cutpoint;
    pch2 += cutpoint;
    for (i = cutpoint; i < sizeof(double) * c1->ch_size ; i++)
    {
        tmpch = *pch1;
        *pch1 = *pch2;
        *pch2 = tmpch;
        pch1++;
        pch2++;
    }
}


void multi_mcga_uniform_crossover (struct MultiChromosome *c1, struct MultiChromosome *c2)
{
    unsigned char *pch1 = (unsigned char*) c1->genes;
    unsigned char *pch2 = (unsigned char*) c2->genes;
    unsigned char tmpch;
    int i;
    for (i = 0; i < sizeof(double) * c1->ch_size ; i++)
    {
        //if( ((double) rand() / RAND_MAX) < 0.5)
        if( runif(0,1) < 0.5)
        {
            tmpch = *pch1;
            *pch1 = *pch2;
            *pch2 = tmpch;
        }
        pch1++;
        pch2++;
    }
}


int multi_mcga_sorter(const void *c1, const void *c2)
{
    struct MultiChromosome *cc1=(struct MultiChromosome*)c1;
    struct MultiChromosome *cc2=(struct MultiChromosome*)c2;
   
    if(cc1->rank >= cc2->rank)
    {
        return(-1);
    }
    else
    {
        return(1);
    }
}

void multi_mcga_calculate_ranks(struct MultiMcga *mcga){
	int i,j,h;
	int rank=0;
	for (i = 0; i < mcga->pop_size; i++){
		rank=0;
		for (j = 0; j < mcga->pop_size; j++){
			for (h=0; h < mcga->numfunc; h++){
				if(mcga->chromosomes[i].cost[h] < mcga->chromosomes[j].cost[h]){
					rank++;
					break;
				}
			}
		}
	mcga->chromosomes[i].rank = rank;
	}	
}

void multi_mcga_sortpopulation(struct MultiMcga *mcga)
{
    qsort(mcga->chromosomes, mcga->pop_size, sizeof(struct MultiChromosome), multi_mcga_sorter);
}

void multi_mcga_fill_chromosome (struct MultiMcga *mcga, struct MultiChromosome *c, int idx)
{
    int i;
    if(c->ch_size != mcga->ch_size) return;
    for (i=0; i<c->ch_size; i++)
    {
        mcga->chromosomes[idx].genes[i] = c->genes[i];
    }
    mcga->chromosomes[idx].ch_size = c->ch_size;
    mcga->chromosomes[idx].rank = c->rank;
    for (i=0; i< c->numfunc; i++){
	    mcga->chromosomes[idx].cost[i] = c->cost[i];
    }
}

void multi_mcga_swap_population(struct MultiMcga *source, struct MultiMcga *target)
{
    int i;
    for (i=0; i<source->pop_size; i++)
    {
        multi_mcga_fill_chromosome(source, &target->chromosomes[i], i);
    }
}

void multi_mcga_tournament_selection (struct MultiMcga *source, struct MultiMcga *target)
{
    int totalselected = 0;
    int idx1, idx2, idx3, idx4;
    int winner1, winner2;
    struct MultiChromosome *c1;
    struct MultiChromosome *c2;
    int last_filled1 = 0, last_filled2 =0;

    if(source->elitism > 0)
    {
        int i;
        multi_mcga_sortpopulation(source);
        for (i =0; i<source->elitism; i++)
        {
            multi_mcga_fill_chromosome(target, &source->chromosomes[i], i);
            totalselected++;
        }
    }

    while(totalselected < source->pop_size)
    {
rebuildIndexes1:
        //idx1 = (int)(((double)rand() / RAND_MAX) * source->pop_size);
        //idx2 = (int)(((double)rand() / RAND_MAX) * source->pop_size);
        idx1 = (int)(runif(0,1) * source->pop_size);
        idx2 = (int)(runif(0,1) * source->pop_size);
        if(idx1 == idx2) goto rebuildIndexes1;
rebuildIndexes2:
        //idx3 = (int)(((double)rand() / RAND_MAX) * source->pop_size);
        //idx4 = (int)(((double)rand() / RAND_MAX) * source->pop_size);
        idx3 = (int)(runif(0,1) * source->pop_size);
        idx4 = (int)(runif(0,1) * source->pop_size);

        if(idx3 == idx4) goto rebuildIndexes2;

        if(source->chromosomes[idx1].rank > source->chromosomes[idx2].rank )
        {
            winner1 = idx1;
        }
        else
        {
            winner1 = idx2;
        }
        if(source->chromosomes[idx3].rank < source->chromosomes[idx4].rank )
        {
            winner2 = idx3;
        }
        else
        {
            winner2 = idx4;
        }


        c1 = &source->chromosomes[winner1];
        c2 = &source->chromosomes[winner2];
        last_filled1 = totalselected;
        multi_mcga_fill_chromosome(target, c1, totalselected);
        totalselected++;
        if(totalselected >= source->pop_size) break;
        multi_mcga_fill_chromosome(target, c2, totalselected);
        last_filled2 = totalselected;
        totalselected++;

        /*
        Crossing over
        */
        //if( ((double) rand() / RAND_MAX) < source->cross_prob)
        if( runif(0,1) < source->cross_prob)
        {
            if(source->cross_type == UNIFORM_CROSS_OVER)
            {
                multi_mcga_uniform_crossover(&target->chromosomes[last_filled1], &target->chromosomes[last_filled2]);
            }
        }

        /*
        Mutating
        */
        multi_mcga_mutate(&target->chromosomes[last_filled1], source->mutate_prob);
        multi_mcga_mutate(&target->chromosomes[last_filled2], source->mutate_prob);
    }
    multi_mcga_swap_population(source, target);
}


void multi_mcga_start(struct MultiMcga *mcga, struct MultiMcga *temppop, int maxiter, double min, double max)
{
    mcga->max_iter = maxiter;

    multi_mcga_randomize(mcga, min, max);
    multi_mcga_calculate_fitness(mcga);
    for (mcga->current_iter = 0; mcga->current_iter < mcga->max_iter; mcga->current_iter++)
    {
        multi_mcga_tournament_selection(mcga, temppop);
        multi_mcga_calculate_fitness(mcga);
    }

}

