#ifndef _H_GEP_INTERNAL_
#define _H_GEP_INTERNAL_

#include <gep.h>
#include <gep_ops.h>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <memory.h>

#define GEP_DBG_PREFIX "GEP "

#define NODE_SIZE sizeof(gep_node_t)


#define GEP_ABORT(msg) \
  fprintf(stderr, GEP_DBG_PREFIX "[%s:%d] Internal error: %s!\n", __FUNCTION__, __LINE__, msg); \
  abort();

#define GEP_COND_ABORT(cond, msg) \
  if (!(cond)) \
  { \
    if (NULL != msg) \
    { \
      GEP_ABORT(msg); \
    } \
    break; \
  }

/*----------------------------------------------------------------------------*/
//!
typedef struct
{
  int     samples_count;
  int     input_dimensions_count;
  int    *shuffle;
  double *inputs;
  int    *inputs_sizes;
  double *inputs_mins;
  double *inputs_maxs;
  double *outputs;
  double  output_min;
  double  output_scale;
  double *MSEs;
} gep_problem_t;


/*----------------------------------------------------------------------------*/
//!
typedef enum
{
  GEP_NODE__FUNCTION = 0,
  GEP_NODE__INPUT,        /* GEP_NODE__TERMINAL */
  GEP_NODE__CONSTANT,     /* GEP_NODE__TERMINAL */
  GEP_NODE_MAX
} gep_node_type_t;


/*----------------------------------------------------------------------------*/
//!
typedef struct gep_node_s
{
  gep_node_type_t node_type;

  union
  {
    union
    {
      unsigned input_index;
      double value;
    } terminal;

    struct
    {
      int    operation_type;
      struct gep_node_s *left;
      struct gep_node_s *right;
      struct gep_node_s *middle;
    } function;
  } node;

  int calculated;
  double value;
} gep_node_t;


/*----------------------------------------------------------------------------*/
//!
typedef struct
{
  gep_node_t *chromosome;
  int        *head_lengths;
  int        *working_lengths;
  int         tree_size;
  int         et_valid;
  int         viable;
} gep_individual_t;


/*----------------------------------------------------------------------------*/
//!
typedef struct
{
  int     population_size;
  int     gene_tail_idx;            //< Index of first tail element
  int     gene_size;
  int     genes_count;
  int     chromosome_size;          //< genes_count * gene_size
  int     max_depth;
  int     input_dimensions_count;
  int     subtrees_operation;
  gep_op_set_t
         *ops_set;
  gep_coding_type_t
          coding_type;
  gep_fitness_type_t
          fitness_type;
  int     mutations_per_chromosome;
  double  mse_coefficient;

  int     use_tournament;
  int     use_replace_worst;
  int     use_incremental_evolution;
  int     use_selection_probability_density;
  int     use_additional_population;
  int     use_dynamic_constants;

  double  mutation_rate;
  double  constants_mutation_rate;
  double  constants_mutation_intensity;
  double  inversion_rate;
  double  IS_transposition_rate;
  double  RIS_transposition_rate;
  double  genes_transposition_rate;
  double  global_recombination_rate;
  double  one_point_recombination_rate;
  double  two_points_recombination_rate;
} gep_params_t;


/*----------------------------------------------------------------------------*/
//!
typedef struct
{
  gep_params_t     *params;
  gep_individual_t *individuals;
  gep_individual_t *temp_individuals;

  //Output values
  double *fitnesses;
  double *probabilities;
  gep_individual_t *best;
  gep_individual_t *worst;
  double best_fitness;
  double mean_fitness;
  int viable_count;
} gep_population_t;


double gep_GenerateConstant(void);

void gep_CreateRandomTree(gep_individual_t *individual, const gep_params_t *params);

void gep_CreateRandomTerminal(gep_node_t *node, const gep_params_t *params);
void gep_CreateRandomNode    (gep_node_t *node, const gep_params_t *params);

void gep_Replicate  (gep_individual_t *individual,                   const gep_params_t *params);
void gep_Recombinate(gep_individual_t *ind1, gep_individual_t *ind2, const gep_params_t *params);

double GEP_MSE                (gep_problem_t *problem, gep_individual_t *individual, const gep_params_t *params);
double GEP_fitness_MSE        (gep_problem_t *problem, gep_individual_t *individual, const gep_params_t *params);
double GEP_fitness_MSE_Partial(gep_problem_t *problem, gep_individual_t *individual, const gep_params_t *params, double margin);
double GEP_fitness_Rsquared   (gep_problem_t *problem, gep_individual_t *individual, const gep_params_t *params);

int    GEP_CopyChromosome(gep_individual_t *dst,        gep_params_t *dst_params, gep_individual_t *src, gep_params_t *src_params);
double GEP_CalculateET   (gep_individual_t *individual, const gep_params_t *params, const double inputs[]);
void   GEP_PrintTree     (gep_individual_t *individual, gep_params_t *params);
int    GEP_CalculateSize (gep_individual_t *individual, gep_params_t *params);

int  GEP_InitPopulation  (gep_population_t *population);
void GEP_DeinitPopulation(gep_population_t *population);

void GEP_RecalculatePopulation(gep_problem_t *problem, gep_population_t *population);
void GEP_Iterate              (gep_problem_t *problem, gep_population_t *population);

gep_problem_t* GEP_ReadProblemFromFile(const char *filename);
void           GEP_WriteProblemToFile (const gep_problem_t *problem, const double *points, const char *filename);
void           GEP_DeallocateProblem(gep_problem_t *problem);

#endif /* _H_GEP_INTERNAL_ */
