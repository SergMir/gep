#include <gep.h>
#include <gep_ops.h>
#include <gep_internal.h>
#include <utils.h>

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <time.h>


/*----------------------------------------------------------------------------*/
//!
struct gep_ctx_s
{
  gep_problem_t      *problem;
  gep_population_t    population;
  gep_population_t    aux_population;

  int flag_first_init;
  int flag_reinit;
  volatile int flag_stop;

  int tree_depth_max;
  int tree_depth_step;
  int tree_depth_iters;
  int tree_depth_min;
  int tree_depth_current;

  int generations_count;
  int current_generation;

  double best_fitness;
  gep_individual_t best_individual;
  gep_params_t max_params;

  char output_samples_filename[256];
};


/*----------------------------------------------------------------------------*/
//!
static void recalcParamSize(gep_params_t *params, int tree_depth)
{
  params->gene_tail_idx = 0;
  int nodes_on_depth = 1;
  for (int i = 1; i < tree_depth; ++i)
  {
    params->gene_tail_idx += nodes_on_depth;
    nodes_on_depth *= GEP_MAX_ARGS;
  }
  params->gene_size       = params->gene_tail_idx + nodes_on_depth;
  params->chromosome_size = params->gene_size * params->genes_count;
  params->max_depth       = tree_depth;

  params->mutation_rate = ((double)params->mutations_per_chromosome) / ((double)params->chromosome_size);
}

/*----------------------------------------------------------------------------*/
//!
static gep_ctx_t* gep_CreateDefaultContext(void)
{
  gep_ctx_t *ctx = malloc(sizeof(gep_ctx_t));

  memset(ctx, 0, sizeof(*ctx));

  ctx->problem            = NULL;
  
  ctx->population.params  = malloc(sizeof(*ctx->population.params));
  gep_params_t *params    = ctx->population.params;
  memset(params, 0, sizeof(*params));
  if (params->use_additional_population)
  {
    ctx->aux_population.params = params;
  }

  params->mutation_rate                 = 0.08;
  params->constants_mutation_rate       = 0.3;
  params->constants_mutation_intensity  = 0.1;
  params->inversion_rate                = 0.1;
  params->IS_transposition_rate         = 0.1;
  params->RIS_transposition_rate        = 0.1;
  params->genes_transposition_rate      = 0.1;
  params->global_recombination_rate     = 0.7;
  params->one_point_recombination_rate  = 0.3;
  params->two_points_recombination_rate = 0.3;
  params->subtrees_operation            = 0; /* Index of first operation on op_set */

  ctx->flag_first_init    = 1;
  ctx->flag_reinit        = 1;
  ctx->flag_stop          = 0;

  ctx->tree_depth_max     = 1;
  ctx->tree_depth_step    = 1;
  ctx->tree_depth_iters   = 1;
  ctx->tree_depth_min     = 3;
  ctx->tree_depth_current = 3;

  ctx->generations_count  = 0;
  ctx->current_generation = 0;

  ctx->best_fitness       = 0;

  ctx->output_samples_filename[0] = '\0';

  return ctx;
}


/*----------------------------------------------------------------------------*/
//!
gep_ctx_t* GEP_CreateContext(const gep_short_params_t *short_params, gep_stat_params_t *stat_params)
{
  srand(time(NULL));

  gep_ctx_t    *ctx          = gep_CreateDefaultContext();
  gep_params_t *params       = ctx->population.params;
  ctx->problem               = GEP_ReadProblemFromFile(short_params->input_samples_filename);
  strncpy(ctx->output_samples_filename, short_params->output_samples_filename, sizeof(ctx->output_samples_filename));

  ctx->generations_count     = short_params->generations_count;
  params->population_size    = short_params->population_size;
  params->genes_count        = short_params->genes_count;
  params->input_dimensions_count
                             = ctx->problem->input_dimensions_count;
  params->ops_set            = GEP_GetOpSet(short_params->ops_preset);
  params->max_depth          = short_params->tree_depth;
  params->mutations_per_chromosome
                             = short_params->mutations_per_chromosome;
  params->fitness_type       = short_params->fitness_type;
  params->mse_coefficient    = short_params->mse_coefficient;

  params->use_tournament     = short_params->use_tournament;
  params->use_replace_worst  = short_params->use_replace_worst;
  params->use_incremental_evolution
                             = short_params->use_incremental_evolution;
  params->use_selection_probability_density
                             = short_params->use_selection_probability_density;
  params->use_additional_population
                             = short_params->use_additional_population;
  params->use_dynamic_constants
                             = short_params->use_dynamic_constants;

  ctx->tree_depth_max        = short_params->tree_depth;
  if (!params->use_incremental_evolution)
  {
    ctx->tree_depth_current  = ctx->tree_depth_min = short_params->tree_depth;
  }
  ctx->tree_depth_step       = (ctx->tree_depth_max - ctx->tree_depth_min) / 10;
  if (0 == ctx->tree_depth_step)
  {
    ctx->tree_depth_step     = 1;
  }
  ctx->tree_depth_iters      = (ctx->tree_depth_max - ctx->tree_depth_min) / ctx->tree_depth_step + 1;

  memcpy(&ctx->max_params, params, sizeof(ctx->max_params));
  recalcParamSize(&ctx->max_params, ctx->tree_depth_max);

  ctx->best_individual.chromosome      = malloc(ctx->max_params.chromosome_size * NODE_SIZE);
  ctx->best_individual.head_lengths    = malloc(ctx->max_params.genes_count * sizeof(ctx->best_individual.head_lengths[0]));
  ctx->best_individual.working_lengths = malloc(ctx->max_params.genes_count * sizeof(ctx->best_individual.working_lengths[0]));
  ctx->best_individual.viable          = 0;

  stat_params->dimensions_count        = params->input_dimensions_count;
  stat_params->input_samples_count     = ctx->problem->samples_count;
  stat_params->total_generations_count = ctx->generations_count * ctx->tree_depth_iters;

  return ctx;
}


/*----------------------------------------------------------------------------*/
//!
void GEP_Run(gep_ctx_t *ctx, gep_stat_callback_t callback)
{
  gep_params_t     *params     =  ctx->population.params;
  gep_population_t *population = &ctx->population;
  gep_problem_t    *problem    =  ctx->problem;

  callback(0, 0, 0, 0, problem->outputs, NULL);

  double *predicted_points = malloc(problem->samples_count * sizeof(predicted_points[0]));

  while (!ctx->flag_stop)
  {
    if (ctx->flag_reinit)
    {
      recalcParamSize(params, ctx->tree_depth_current);

      GEP_InitPopulation(population);
      if (!ctx->flag_first_init)
      {
        GEP_CopyChromosome(population->best, params, &ctx->best_individual, &ctx->max_params);
      }
      GEP_RecalculatePopulation(problem, population);

      if (params->use_additional_population)
      {
        ctx->aux_population.params = params;
        GEP_InitPopulation(&ctx->aux_population);
        GEP_RecalculatePopulation(problem, &ctx->aux_population);
      }

      ctx->flag_first_init = 0;
      ctx->current_generation = 0;
      ctx->flag_reinit = 0;
    }

    GEP_Iterate(problem, population);

    if (params->use_additional_population)
    {
      GEP_Iterate(problem, &ctx->aux_population);

      if (population->best_fitness < ctx->aux_population.best_fitness)
      {
        GEP_CopyChromosome(population->worst, params, ctx->aux_population.best, params);
        GEP_RecalculatePopulation(problem, population);
      }
    }

    GEP_CopyChromosome(&ctx->best_individual, &ctx->max_params, population->best, params);

    int flag_update_points = 0;
    if (population->best_fitness > ctx->best_fitness)
    {
      ctx->best_fitness = population->best_fitness;
      flag_update_points = 1;
    }

    double *inputs = problem->inputs;
    for (int i = 0; i < problem->samples_count; ++i)
    {
      predicted_points[i] = GEP_CalculateET(population->best, params, inputs);

      inputs += problem->input_dimensions_count;
    }

    callback(
      population->best_fitness,
      population->mean_fitness,
      GEP_MSE(problem, population->best, params),
      ctx->current_generation + ctx->generations_count * (ctx->tree_depth_iters - 1),
      NULL,
      flag_update_points ? predicted_points : NULL);

    if (++ctx->current_generation >= ctx->generations_count)
    {
      if ((ctx->tree_depth_current + ctx->tree_depth_step) > ctx->tree_depth_max)
      {
        ctx->flag_stop = 1;
      }
      else
      {
        GEP_DeinitPopulation(population);
        if (params->use_additional_population)
        {
          GEP_DeinitPopulation(&ctx->aux_population);
        }

        ctx->tree_depth_current += ctx->tree_depth_step;
        ctx->flag_reinit = 1;
      }
    }
  }

  free(predicted_points);

  fflush(stdout);
  if (0 != strcmp("", ctx->output_samples_filename))
  {
    double points[problem->samples_count];
    double *inputs = problem->inputs;

    for (int i = 0; i < problem->samples_count; ++i)
    {
      points[i] = GEP_CalculateET(population->best, params, inputs);

      inputs  += problem->input_dimensions_count;
    }
    GEP_WriteProblemToFile(problem, points, ctx->output_samples_filename);
  }

  //GEP_CalculateSize(population.best, params);
}


/*----------------------------------------------------------------------------*/
//!
void GEP_Stop(gep_ctx_t *ctx)
{
  ctx->flag_stop = 1;
}


/*----------------------------------------------------------------------------*/
//!
void GEP_DestroyContext(gep_ctx_t *ctx)
{
  gep_params_t *params =  ctx->population.params;

  free(ctx->best_individual.chromosome);
  free(ctx->best_individual.head_lengths);
  free(ctx->best_individual.working_lengths);

  free(params->ops_set->ops);
  free(params->ops_set);

  GEP_DeinitPopulation(&ctx->population);
  if (params->use_additional_population)
  {
    GEP_DeinitPopulation(&ctx->aux_population);
  }

  free(params);

  GEP_DeallocateProblem(ctx->problem);
  free(ctx);
}
