#include <gep_internal.h>
#include <utils.h>

#include <math.h>
#include <omp.h>

/*----------------------------------------------------------------------------*/
//!
static void gep_CalculateRoulette(gep_population_t *population)
{
  int population_size = population->params->population_size;

  if (population->params->use_selection_probability_density)
  {
    double all_diffs = 0;
    for (int i = 0; i < population_size; ++i)
    {
      for (int j = 0; j < population_size; ++j)
      {
        if (i != j)
        {
          double diff = population->fitnesses[i] - population->fitnesses[j];
          diff = (diff < 0) ? (-diff) : (diff);
          all_diffs += diff;
        }
      }
    }

    for (int k = 0; k < population_size; ++k)
    {
      double sum_diff = 0;
      for (int i = 0; i < population_size; ++i)
      {
        double diff = population->fitnesses[k] - population->fitnesses[i];
        diff = (diff < 0) ? (-diff) : (diff);
        sum_diff += diff;
      }
      population->probabilities[k] = sum_diff / all_diffs;
    }
  }
  else
  {
    double sum = 0;

    for (int i = 0; i < population_size; ++i)
    {
      if (population->individuals[i].viable)
      {
        sum += population->fitnesses[i] / 100.0;
      }
    }

    for (int i = 0; i < population_size; ++i)
    {
      if (population->individuals[i].viable)
      {
        population->probabilities[i] = (population->fitnesses[i] / 100.0) / sum;
      }
      else
      {
        population->probabilities[i] = 0;
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
//!
static gep_individual_t* gep_SelectFromRoulette(gep_population_t *population, gep_individual_t *forbidden)
{
  int i = 0;
  double sum = population->probabilities[i];
  double probe = (rand() % 100000) / 100000.0;

  if (1 == population->viable_count)
  {
    fprintf(stderr, "Impossible! Can't select from roulette!\n");
    return forbidden;
  }

  while (i < (population->params->population_size - 1) && sum < probe)
  {
    if (population->individuals[i].viable)
    {
      sum += population->probabilities[i + 1];
    }
    ++i;
  }

  if ((i >= population->params->population_size) || ((population->individuals + i) == forbidden))
  {
    int step = (rand() % (population->viable_count - 1)) + 1;
    for (int k = 0; k < step; ++k)
    {
      do
      {
        ++i;
        if (population->params->population_size <= i)
        {
          i = 0;
        }
      } while (!population->individuals[i].viable);
    }
  }

  if (i < 0 || i >= population->params->population_size)
  {
    GEP_ABORT("Broken selection");
  }

  return population->individuals + i;
}

/*----------------------------------------------------------------------------*/
//!
void GEP_RecalculatePopulation(gep_problem_t *problem, gep_population_t *population)
{
  double best_fitness = 0;
  double worst_fitness = GEP_MAX_FITNESS;
  double sum_fitness = 0;

  int viable_count = 0;
  gep_params_t *params = population->params;

  gep_individual_t *ind_best  = &population->individuals[0];
  gep_individual_t *ind_worst = &population->individuals[0];

  #pragma omp parallel for shared(population, params, best_fitness, worst_fitness, sum_fitness, viable_count, ind_best, ind_worst)
  for (int i = 0; i < params->population_size; ++i)
  {
    gep_individual_t *ind = &population->individuals[i];

    population->fitnesses[i] = 0;
    double fitness = 0;

    switch (params->fitness_type)
    {
      case GEP_FITNESS__MSE:
        fitness = GEP_fitness_MSE(problem, ind, params);
        break;

      case GEP_FITNESS__PARTIAL_MSE:
        fitness = GEP_fitness_MSE_Partial(problem, ind, params, population->mean_fitness);
        break;

      case GEP_FITNESS__R_SQUARED:
        fitness = GEP_fitness_Rsquared(problem, ind, params);
        break;

      default:
        GEP_ABORT("Invalid fitness type");
    }

    ind->viable = !isnan(fitness) && !isinf(fitness) && (fitness > 0) && (fitness <= GEP_MAX_FITNESS);

    if (ind->viable)
    {
      #pragma omp critical
      {
        population->fitnesses[i] = fitness;

        ++viable_count;
        sum_fitness += fitness;
        if (fitness > best_fitness)
        {
          best_fitness = fitness;
          ind_best = ind;
        }

        if (fitness < worst_fitness)
        {
          worst_fitness = fitness;
          ind_worst = ind;
        }
      } /* omp critical */
    }
  }

  population->best  = ind_best;
  population->worst = ind_worst;

  if (0 == viable_count)
  {
    GEP_ABORT("No viable individuals left");
  }
  else if (1 == viable_count)
  {
    fprintf(stderr, "Just one viable!\n");
  }
  else if ((population->params->population_size / 4) > viable_count)
  {
    fprintf(stderr, "Less than half of viable: %d%%!\n", (100 * viable_count) / population->params->population_size);
  }

  population->best_fitness = best_fitness;
  population->mean_fitness = sum_fitness / viable_count;
  population->viable_count = viable_count;

  if (!params->use_tournament)
  {
    gep_CalculateRoulette(population);
  }
}

/*----------------------------------------------------------------------------*/
//!
void GEP_Iterate(gep_problem_t *problem, gep_population_t *population)
{
  gep_params_t *params = population->params;

  if (NULL == population->best)
  {
    GEP_ABORT("Invalid argument: not calculated population");
  }

  //Elitism
  GEP_CopyChromosome(&population->temp_individuals[0], params, population->best, params);

  int times_selected[params->population_size];
  for (int i = 0; i < params->population_size; ++i)
  {
    times_selected[i] = 0;
  }
  times_selected[0] = 1;

  if (params->use_tournament)
  {
    for (int i = 1; i < params->population_size; ++i)
    {
      int j1 = rand() % params->population_size; 
      int j2 = rand() % params->population_size;
      while (j1 == j2)
      {
        j2 = rand() % params->population_size;
      }
      double f_1 = population->fitnesses[j1];
      double f_2 = population->fitnesses[j2];

      gep_individual_t *src_ind = (f_1 > f_2) ? &population->individuals[j1] : &population->individuals[j2];
      gep_individual_t *dst_ind = &population->temp_individuals[i];

      ++times_selected[src_ind - population->individuals];

      GEP_CopyChromosome(dst_ind, params, src_ind, params);
      gep_Replicate(dst_ind, params);
    }
  }
  else
  {
    for (int i = 1; i < params->population_size; ++i)
    {
      gep_individual_t *src_ind = gep_SelectFromRoulette(population, NULL);
      gep_individual_t *dst_ind = &population->temp_individuals[i];

      ++times_selected[src_ind - population->individuals];

      GEP_CopyChromosome(dst_ind, params, src_ind, params);
      gep_Replicate(dst_ind, params);
    }
  }

  /* Recombination */
  for (int i = 1; i < params->population_size; ++i)    
  {
    double recombination_probe = normal_rand();
    if (recombination_probe < params->global_recombination_rate)
    {
      int j = rand() % (params->population_size - 1) + 1;
      while (j == i)
      {
        j = rand() % (params->population_size - 1) + 1;
      }

      gep_individual_t *ind1 = &population->temp_individuals[i];
      gep_individual_t *ind2 = &population->temp_individuals[j];

      gep_Recombinate(ind1, ind2, params);
    }
  }

  gep_individual_t *individuals_swap = population->individuals;
  population->individuals            = population->temp_individuals;
  population->temp_individuals       = individuals_swap;

  GEP_RecalculatePopulation(problem, population);

  if (params->use_replace_worst)
  {
    int replaced_count = 0;
    for (int i = 0; i < params->population_size; ++i)
    {
      if ((population->fitnesses[i] < (population->best_fitness * 0.5)) || !population->individuals[i].viable)
      {
        ++replaced_count;
        gep_CreateRandomTree(&population->individuals[i], params);
      }
    }
    GEP_RecalculatePopulation(problem, population);
  }
}
