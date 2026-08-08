#include <R.h>
#include <Rdefines.h>
#include <stdio.h>
#include "mcga.h"
#include "multi_mcga.h"

SEXP
mcga (SEXP popsize, SEXP chsize, SEXP crossp, SEXP mutatep, SEXP elitism,
      SEXP min, SEXP max, SEXP maxiter, SEXP par, SEXP best, SEXP evalFunc,
      SEXP population, SEXP costs, SEXP env)
{
  int popsize1 = NUMERIC_POINTER (popsize)[0];
  int chsize1 = NUMERIC_POINTER (chsize)[0];
  int elitism1 = NUMERIC_POINTER (elitism)[0];
  double crossp1 = NUMERIC_POINTER (crossp)[0];
  double mutatep1 = NUMERIC_POINTER (mutatep)[0];
  double min1 = NUMERIC_POINTER (min)[0];
  double max1 = NUMERIC_POINTER (max)[0];
  int maxiter1 = NUMERIC_POINTER (maxiter)[0];
  double currentfitness;

  SEXP fn, sexp_fvec;
  struct Mcga *mcga =
    mcga_create (popsize1, chsize1, crossp1, mutatep1, UNIFORM_CROSS_OVER,
		 elitism1, NULL);
  struct Mcga *mcga2 =
    mcga_create (popsize1, chsize1, crossp1, mutatep1, UNIFORM_CROSS_OVER,
		 elitism1, NULL);
  mcga_randomize (mcga, min1, max1);

  for (int i = 0; i < maxiter1; i++)
    {
      mcga_tournament_selection (mcga, mcga2);
      for (int j = 0; j < popsize1; j++)
	{
	  for (int h = 0; h < chsize1; h++)
	    {
	      NUMERIC_POINTER (par)[h] = mcga->chromosomes[j].genes[h];
	    }
	  fn = PROTECT(lang2 (evalFunc, par));
	  sexp_fvec = PROTECT(eval (fn, env));
	  currentfitness = NUMERIC_POINTER (sexp_fvec)[0];
	  mcga->chromosomes[j].cost = currentfitness;
	  UNPROTECT(2);
	}
    }

  mcga_sortpopulation (mcga);

  for (int i = 0; i < popsize1; i++)
    {
      for (int j = 0; j < chsize1; j++)
	{
	  NUMERIC_POINTER (population)[i * chsize1 + j] =
	    mcga->chromosomes[i].genes[j];
	}
      NUMERIC_POINTER (costs)[i] = mcga->chromosomes[i].cost;
    }

	mcga_free(mcga);
	mcga_free(mcga2);

  return (popsize);
}








SEXP
multi_mcga (SEXP popsize, SEXP chsize, SEXP crossp, SEXP mutatep, SEXP elitism,
      SEXP min, SEXP max, SEXP maxiter, SEXP par, SEXP best, SEXP evalFunc, SEXP numfunc,
      SEXP population, SEXP costs, SEXP ranks, SEXP env)
{
  int popsize1 = NUMERIC_POINTER (popsize)[0];
  int chsize1 = NUMERIC_POINTER (chsize)[0];
  int elitism1 = NUMERIC_POINTER (elitism)[0];
  double crossp1 = NUMERIC_POINTER (crossp)[0];
  double mutatep1 = NUMERIC_POINTER (mutatep)[0];
  double min1 = NUMERIC_POINTER (min)[0];
  double max1 = NUMERIC_POINTER (max)[0];
  int maxiter1 = NUMERIC_POINTER (maxiter)[0];
  int numfunc1 = NUMERIC_POINTER (numfunc)[0];

  SEXP fn, sexp_fvec;
  struct MultiMcga *mcga =
    multi_mcga_create (popsize1, chsize1, crossp1, mutatep1, UNIFORM_CROSS_OVER,
		 elitism1, numfunc1, NULL);
  struct MultiMcga *mcga2 =
    multi_mcga_create (popsize1, chsize1, crossp1, mutatep1, UNIFORM_CROSS_OVER,
		 elitism1, numfunc1, NULL);
  multi_mcga_randomize (mcga, min1, max1);

  for (int i = 0; i < maxiter1; i++)
    {
      multi_mcga_tournament_selection (mcga, mcga2);
      for (int j = 0; j < popsize1; j++)
	{
	  for (int h = 0; h < chsize1; h++)
	    {
	      NUMERIC_POINTER (par)[h] = mcga->chromosomes[j].genes[h];
	    }
	  fn = PROTECT(lang2 (evalFunc, par));
	  sexp_fvec = PROTECT(eval (fn, env));
	  /* currentfitness = NUMERIC_POINTER (sexp_fvec); */
	  /* mcga->chromosomes[j].cost = currentfitness; */
	 for (int h=0; h<numfunc1; h++){
		mcga->chromosomes[j].cost[h]=NUMERIC_POINTER(sexp_fvec)[h];
	  }
	  UNPROTECT(2);
	}
      multi_mcga_calculate_ranks(mcga);
    }

  multi_mcga_calculate_ranks(mcga);
  multi_mcga_sortpopulation (mcga);
  for (int i = 0; i < popsize1; i++)
    {
      for (int j = 0; j < chsize1; j++)
	{
	  NUMERIC_POINTER (population)[i * chsize1 + j] =
	    mcga->chromosomes[i].genes[j];
	}
      for (int h=0; h<numfunc1; h++){
	      NUMERIC_POINTER (costs)[i * numfunc1 + h] = mcga->chromosomes[i].cost[h];
      }
      NUMERIC_POINTER (ranks)[i] = mcga->chromosomes[i].rank;
    }

	multi_mcga_free(mcga);
	multi_mcga_free(mcga2);

  return (popsize);
}


