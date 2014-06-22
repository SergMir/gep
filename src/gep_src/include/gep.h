#ifndef _H_GEP_
#define _H_GEP_

#include <stdio.h>

#define GEP_MAX_FITNESS      1000.0

//! Preset of functional set
/**
 * "short":  '+', '-', '*', '/'
 * "middle": "short" + '-' (unary), 'sqrt'. 'exp'
 * "big":    "middle" + 'pow', 'ln', 'sin', 'cos'
 * "full":   "big" + 'log2' + 'log10' + 'tan' + 'sinh' + 'cosh' + 'tanh'
 */
typedef enum
{
  GEP_OPLIST__SHORT = 0,
  GEP_OPLIST__MIDDLE,
  GEP_OPLIST__BIG,
  GEP_OPLIST__FULL,
  GEP_OPLIST_LAST
} gep_ops_preset_t;


//! Type of tree coding
/**
 * "breadth-first": traditional coding, used by C.Ferreira
 * "prefix": depth-first traverse of tree
 * "overlapped": depth-first traverse with re-usage of nodes
 */
typedef enum
{
  GEP_CODING__BREADTH_FIRST = 0,
  GEP_CODING__PREFIX,
  GEP_CODING__OVERLAPPED,
  GEP_CODING_LAST
} gep_coding_type_t;


//! Type of fitness function
/**
 * "mse": mean squared error-based fitness
 * "partial_mse": quick implementation of mse, with threshold for interrupting calculation
 * "r_squared": fitness function based on R^2 coefficient of determination
 */
typedef enum
{
  GEP_FITNESS__MSE = 0,
  GEP_FITNESS__PARTIAL_MSE,
  GEP_FITNESS__R_SQUARED,
  GEP_FITNESS_LAST
} gep_fitness_type_t;

/*----------------------------------------------------------------------------*/
//!
typedef struct
{
  int                population_size;
  int                tree_depth;
  int                genes_count;
  int                generations_count;
  int                mutations_per_chromosome;
  gep_ops_preset_t   ops_preset;
  gep_coding_type_t  coding_type;
  gep_fitness_type_t fitness_type;
  double             mse_coefficient;
  int                use_tournament;
  int                use_replace_worst;
  int                use_incremental_evolution;
  int                use_selection_probability_density;
  int                use_additional_population;
  int                use_dynamic_constants;

  char               input_samples_filename[256];
  char               output_samples_filename[256];
} gep_short_params_t;


/*----------------------------------------------------------------------------*/
//!
typedef struct
{
  int dimensions_count;
  int input_samples_count;
  int total_generations_count;
} gep_stat_params_t;


typedef struct gep_ctx_s gep_ctx_t;

typedef void (*gep_stat_callback_t)(double best_fitness, double avg_fitness, double best_mse, int generation, double *target_points, double *predicted_points);

gep_ctx_t* GEP_CreateContext(const gep_short_params_t *short_params, gep_stat_params_t *stat_params);

void GEP_Run           (gep_ctx_t *ctx, gep_stat_callback_t callback);
void GEP_Stop          (gep_ctx_t *ctx);
void GEP_DestroyContext(gep_ctx_t *ctx);

#endif /* _H_GEP_ */
