#include "mcga.h"
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

struct Mcga *mcga_create (int popsize, int chsize, double crossprob, double mutateprob, int crosstype, int elitism, double (*cost_function)(struct Chromosome*))
{
    int i;

    struct Mcga *mcga = malloc (sizeof(struct Mcga));

	 /*
	 We are removing 'assert' lines for version 2.14
    assert(mcga!=NULL);
    assert(crossprob>=0.0 && crossprob<=1.0);
    assert(mutateprob>=0.0 && mutateprob<=1.0);
    assert(elitism>=0 && elitism<popsize);
    assert(popsize>=2);
    assert(chsize>=1);
	 */
    mcga->pop_size=popsize;
    mcga->ch_size = chsize;
    mcga->cross_type=crosstype;
    mcga->cross_prob = crossprob;
    mcga->mutate_prob = mutateprob;
    mcga->elitism = elitism;

    mcga->chromosomes = malloc(sizeof(struct Chromosome) * popsize);
    /* assert(mcga->chromosomes != NULL); */

    for (i=0; i<popsize; i++)
    {
        mcga->chromosomes[i].ch_size = chsize;
        mcga->chromosomes[i].genes = malloc (sizeof(double) * chsize);
        /* assert(mcga->chromosomes[i].genes != NULL); */
        mcga->chromosomes[i].cost = 0.0;
    }

    mcga->cost_function = cost_function;
    return(mcga);
}

struct Mcga *mcga_clone (struct Mcga *mcga)
{
    struct Mcga *mcga2 = mcga_create(mcga->pop_size, mcga->ch_size, mcga->cross_prob, mcga->mutate_prob, mcga->cross_type, mcga->elitism, mcga->cost_function);
    return(mcga2);
}


void mcga_free (struct Mcga *mcga)
{
    free(mcga->chromosomes->genes);
    free(mcga->chromosomes);
    free(mcga);
}


void mcga_dump_chromosome_chars(struct Chromosome *c)
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

void mcga_dump_chromosome (struct Chromosome *c, int idx)
{
    int j;
    Rprintf("%d)[%.3f] ",idx,c->cost);
    for(j=0; j<c->ch_size; j++)
    {
        Rprintf("%.5f ",c->genes[j]);
    }
    Rprintf("%s", "\n");
}


void mcga_dump (struct Mcga *mcga)
{
    int i;
    for (i =0; i<mcga->pop_size; i++)
    {
        mcga_dump_chromosome(&mcga->chromosomes[i], i);
    }
}


void mcga_randomize(struct Mcga *mcga, double min, double max)
{
    int i,j;
    for (i =0; i<mcga->pop_size; i++)
    {
        for (j=0; j<mcga->ch_size; j++)
        {
            //mcga->chromosomes[i].genes[j] = min + ((double)rand() / RAND_MAX) * (max-min);
            mcga->chromosomes[i].genes[j] = runif(min, max);
            mcga->chromosomes[i].cost = 0.0;
        }
    }
}


void mcga_set_cost (struct Mcga *mcga, int idx, double val)
{
    mcga->chromosomes[idx].cost = val;
}


double mcga_get_cost (struct Mcga *mcga, int idx)
{
    return(mcga->chromosomes[idx].cost);
}


void mcga_calculate_fitness(struct Mcga* mcga)
{
    int i;
    for (i =0; i<mcga->pop_size; i++)
    {
        /* assert(mcga->chromosomes[i].genes!=NULL); */
        mcga->chromosomes[i].cost = mcga->cost_function(&mcga->chromosomes[i]);
    }
}

void mcga_mutate(struct Chromosome *c,  double prob)
{
    int i;
    unsigned char *pch = (unsigned char *) c->genes;
    for (i = 0 ; i< (sizeof(double)* c->ch_size); i++)
    {
        /* will mutate current byte? */
        //if ( ((double)rand() / RAND_MAX)<prob)
        if ( runif(0.0,1.0) < prob)
        {
            /* Mutate upwards or downwards ?*/
            //int amount = (int)(((double)rand() / RAND_MAX) * 255);
            //if ( ((double)rand() / RAND_MAX)< 0.5)
            if ( runif(0.0,1.0) < 0.5)
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


void mcga_onepoint_crossover (struct Chromosome *c1, struct Chromosome *c2)
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


void mcga_uniform_crossover (struct Chromosome *c1, struct Chromosome *c2)
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

int mcga_sorter(const void *c1, const void *c2)
{
    struct Chromosome *cc1=(struct Chromosome*)c1;
    struct Chromosome *cc2=(struct Chromosome*)c2;
    if(cc1->cost <= cc2->cost)
    {
        return(-1);
    }
    else
    {
        return(1);
    }
}




void mcga_sortpopulation(struct Mcga *mcga)
{
    qsort(mcga->chromosomes, mcga->pop_size, sizeof(struct Chromosome), mcga_sorter);
}

void mcga_fill_chromosome (struct Mcga *mcga, struct Chromosome *c, int idx)
{
    int i;
	if(c->ch_size != mcga->ch_size) return;
    for (i=0; i<c->ch_size; i++)
    {
        mcga->chromosomes[idx].genes[i] = c->genes[i];
    }
    mcga->chromosomes[idx].ch_size = c->ch_size;
    mcga->chromosomes[idx].cost = c->cost;
}

void mcga_swap_population(struct Mcga *source, struct Mcga *target)
{
    int i;
    for (i=0; i<source->pop_size; i++)
    {
        mcga_fill_chromosome(source, &target->chromosomes[i], i);
    }
}

void mcga_tournament_selection (struct Mcga *source, struct Mcga *target)
{
    int totalselected = 0;
    int idx1, idx2, idx3, idx4;
    int winner1, winner2;
    struct Chromosome *c1;
    struct Chromosome *c2;
    int last_filled1 = 0, last_filled2 =0;

    if(source->elitism > 0)
    {
        int i;
        mcga_sortpopulation(source);
        for (i =0; i<source->elitism; i++)
        {
            mcga_fill_chromosome(target, &source->chromosomes[i], i);
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

        if(source->chromosomes[idx1].cost < source->chromosomes[idx2].cost )
        {
            winner1 = idx1;
        }
        else
        {
            winner1 = idx2;
        }
        if(source->chromosomes[idx3].cost < source->chromosomes[idx4].cost )
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
        mcga_fill_chromosome(target, c1, totalselected);
        totalselected++;
        if(totalselected >= source->pop_size) break;
        mcga_fill_chromosome(target, c2, totalselected);
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
                mcga_uniform_crossover(&target->chromosomes[last_filled1], &target->chromosomes[last_filled2]);
            }
        }

        /*
        Mutating
        */
        mcga_mutate(&target->chromosomes[last_filled1], source->mutate_prob);
        mcga_mutate(&target->chromosomes[last_filled2], source->mutate_prob);
    }
    mcga_swap_population(source, target);
}


void mcga_start(struct Mcga *mcga, struct Mcga *temppop, int maxiter, double min, double max)
{
    /*
	 assert(mcga != NULL);
    assert(temppop != NULL);
    */
    mcga->max_iter = maxiter;

    mcga_randomize(mcga, min, max);
    mcga_calculate_fitness(mcga);
    for (mcga->current_iter = 0; mcga->current_iter < mcga->max_iter; mcga->current_iter++)
    {
        mcga_tournament_selection(mcga, temppop);
        mcga_calculate_fitness(mcga);
    }

}

