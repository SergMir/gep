#include <gep_internal.h>
#include <utils.h>

#include <math.h>
#include <stdio.h>

static const double gold_constants[]   = {0.1316, 0.2128, 0.3441, 0.5571, 0.9015, 1.4588, 2.3605, 3.8195, 6.1804, 10.0007};
static const int    gold_constants_cnt = sizeof(gold_constants) / sizeof(gold_constants[0]);


/*----------------------------------------------------------------------------*/
//!
double gep_GenerateConstant(void)
{
    return 2.0 * (normal_rand() - 0.5) * gold_constants[rand() % gold_constants_cnt];
}


/*----------------------------------------------------------------------------*/
//!
int GEP_InitPopulation(gep_population_t *population)
{
  int status = 1;

  do
  {
    int    pop_size                  = population->params->population_size;
    int    chromosome_size           = population->params->chromosome_size;
    
    size_t fitnesses_malloc_size     = pop_size * sizeof(population->fitnesses[0]);
    size_t probabilities_malloc_size = pop_size * sizeof(population->probabilities[0]);
    size_t individuals_malloc_size   = pop_size * sizeof(population->individuals[0]);
    population->fitnesses            = malloc(fitnesses_malloc_size);
    population->probabilities        = malloc(probabilities_malloc_size);
    population->individuals          = malloc(individuals_malloc_size);
    population->temp_individuals     = malloc(individuals_malloc_size);

    GEP_COND_ABORT(NULL != population->fitnesses,        "Can't malloc fitnesses");
    GEP_COND_ABORT(NULL != population->probabilities,    "Can't malloc probabilities");
    GEP_COND_ABORT(NULL != population->individuals,      "Can't malloc individuals");
    GEP_COND_ABORT(NULL != population->temp_individuals, "Can't malloc temp_individuals");

    size_t all_chromosomes_malloc_size              = pop_size * chromosome_size * NODE_SIZE;
    population->individuals[0].chromosome           = malloc(all_chromosomes_malloc_size);
    population->temp_individuals[0].chromosome      = malloc(all_chromosomes_malloc_size);

    GEP_COND_ABORT(NULL != population->individuals[0].chromosome,      "Can't malloc chromosomes");
    GEP_COND_ABORT(NULL != population->temp_individuals[0].chromosome, "Can't malloc chromosomes");

    size_t all_working_lengths_malloc_size          = pop_size * population->params->genes_count * sizeof(population->individuals[0].working_lengths[0]);
    population->individuals[0].working_lengths      = malloc(all_working_lengths_malloc_size);
    population->temp_individuals[0].working_lengths = malloc(all_working_lengths_malloc_size);
    population->individuals[0].head_lengths         = malloc(all_working_lengths_malloc_size);
    population->temp_individuals[0].head_lengths    = malloc(all_working_lengths_malloc_size);

    GEP_COND_ABORT(NULL != population->individuals[0].working_lengths,      "Can't malloc working_lengths");
    GEP_COND_ABORT(NULL != population->temp_individuals[0].working_lengths, "Can't malloc working_lengths");
    GEP_COND_ABORT(NULL != population->individuals[0].head_lengths,         "Can't malloc head_lengths");
    GEP_COND_ABORT(NULL != population->temp_individuals[0].head_lengths,    "Can't malloc head_lengths");

    for (int i = 1; i < pop_size; ++i)
    {
      population->individuals[i].chromosome      = population->individuals[i - 1].chromosome      + chromosome_size;
      population->temp_individuals[i].chromosome = population->temp_individuals[i - 1].chromosome + chromosome_size;

      population->individuals[i].working_lengths      = population->individuals[i - 1].working_lengths      + population->params->genes_count;
      population->temp_individuals[i].working_lengths = population->temp_individuals[i - 1].working_lengths + population->params->genes_count;

      population->individuals[i].head_lengths      = population->individuals[i - 1].head_lengths      + population->params->genes_count;
      population->temp_individuals[i].head_lengths = population->temp_individuals[i - 1].head_lengths + population->params->genes_count;
    }

    for (int i = 0; i < pop_size; ++i)
    {
      population->individuals[i].et_valid             = 0;
      population->temp_individuals[i].et_valid        = 0;

      gep_CreateRandomTree(&population->individuals[i], population->params);
    }

    population->best         = &population->individuals[0];
    population->worst        = &population->individuals[1];
    population->best_fitness = GEP_MAX_FITNESS;
    population->mean_fitness = GEP_MAX_FITNESS;

    status = 0;

  } while (0);

  return status;
}


/*----------------------------------------------------------------------------*/
//!
void GEP_DeinitPopulation(gep_population_t *population)
{
  free(population->fitnesses);
  free(population->probabilities);
  free(population->individuals[0].chromosome);
  free(population->individuals[0].working_lengths);
  free(population->individuals[0].head_lengths);
  free(population->temp_individuals[0].chromosome);
  free(population->temp_individuals[0].working_lengths);
  free(population->temp_individuals[0].head_lengths);
  free(population->individuals);
  free(population->temp_individuals);
}


/*----------------------------------------------------------------------------*/
//!
static gep_problem_t* gep_StartAllocateProblem(int input_dimensions_count, int samples_count)
{
  gep_problem_t *problem = malloc(sizeof(gep_problem_t));
  problem->input_dimensions_count
                         = input_dimensions_count;
  problem->samples_count = samples_count;
  problem->inputs        = malloc(sizeof(double) * input_dimensions_count * samples_count);
  problem->outputs       = malloc(sizeof(double) * samples_count);
  problem->MSEs          = malloc(sizeof(double) * samples_count);
  problem->shuffle       = malloc(sizeof(int)    * samples_count);
  problem->inputs_mins   = malloc(sizeof(double) * input_dimensions_count);
  problem->inputs_maxs   = malloc(sizeof(double) * input_dimensions_count);
  problem->inputs_sizes  = malloc(sizeof(int)    * input_dimensions_count);

  for (int i = 0; i < problem->samples_count; ++i)
  {
    problem->shuffle[i] = i;
  }
  for (int i = problem->samples_count - 1; i > 1 ; --i)
  {
    int swap_i = rand() % i;
    int swap = problem->shuffle[i];
    problem->shuffle[i] = problem->shuffle[swap_i];
    problem->shuffle[swap_i] = swap;
  }

  return problem;
}


/*----------------------------------------------------------------------------*/
//!
void GEP_DeallocateProblem(gep_problem_t *problem)
{
  free(problem->inputs);
  free(problem->outputs);
  free(problem->MSEs);
  free(problem->shuffle);
  free(problem->inputs_mins);
  free(problem->inputs_maxs);
  free(problem->inputs_sizes);
  free(problem);
}


/*----------------------------------------------------------------------------*/
//!
static void gep_PrintExBranch(gep_node_t *node, gep_op_set_t *ops_set, int *overflow_idx)
{
  if (NULL == node)
  {
    GEP_ABORT("Null node");
  }

  switch (node->node_type)
  {
    case GEP_NODE__FUNCTION:
    {
      switch (ops_set->ops[node->node.function.operation_type]->arguments_count)
      {
        case 1:
          fprintf(stderr, "(%s(", ops_set->ops[node->node.function.operation_type]->op_name);
          gep_PrintExBranch(node->node.function.left, ops_set, overflow_idx);
          fprintf(stderr, ")");
          break;

        case 2:
          fprintf(stderr, "(");
          gep_PrintExBranch(node->node.function.left,  ops_set, overflow_idx);
          fprintf(stderr, " %s ", ops_set->ops[node->node.function.operation_type]->op_name);
          gep_PrintExBranch(node->node.function.right, ops_set, overflow_idx);
          break;

        case 3:
          fprintf(stderr, "%s(", ops_set->ops[node->node.function.operation_type]->op_name);
          gep_PrintExBranch(node->node.function.left,   ops_set, overflow_idx);
          fprintf(stderr, ", ");
          gep_PrintExBranch(node->node.function.middle, ops_set, overflow_idx);
          fprintf(stderr, ", ");
          gep_PrintExBranch(node->node.function.right,  ops_set, overflow_idx);
          break;

        default:
          GEP_ABORT("Bad arguments count");
          break;
      }
      fprintf(stderr, ")");
      break;
    }

    case GEP_NODE__INPUT:
      fprintf(stderr, "x%u", node->node.terminal.input_index);
      break;

    case GEP_NODE__CONSTANT:
      fprintf(stderr, "%f", node->node.terminal.value);
      break;

    default:
      GEP_ABORT("Bad node type");
      break;
  }
}


/*----------------------------------------------------------------------------*/
//!
void GEP_PrintTree(gep_individual_t *individual, gep_params_t *params)
{
#if 1
  fprintf(stderr, "Genome: ");
  for (int gene_idx = 0; gene_idx < params->genes_count; ++gene_idx)
  {
    fprintf(stderr, "G%d: ", gene_idx);
    for (int i = 0; i < params->gene_size; ++i)
    {
      gep_node_t *node = &individual->chromosome[gene_idx * params->gene_size + i];
      switch (node->node_type)
      {
        case GEP_NODE__FUNCTION:
          fprintf(stderr, "%s ", params->ops_set->ops[node->node.function.operation_type]->op_name);
          break;

        case GEP_NODE__INPUT:
          fprintf(stderr, "x%u ", node->node.terminal.input_index);
          break;

        case GEP_NODE__CONSTANT:
          fprintf(stderr, "%f ", node->node.terminal.value);
          break;

        default:
          GEP_ABORT("Bad node type");
          break;
      }
    }
    fprintf(stderr, "; ");
  }
  fprintf(stderr, "\n");
#endif
  fprintf(stderr, "Expression: ");
  for (int gene_idx = 0; gene_idx < params->genes_count; ++gene_idx)
  {
    if (0 != gene_idx)
    {
      fprintf(stderr, "%s", params->ops_set->ops[params->subtrees_operation]->op_name);
    }
    int overflow_idx = 0;
    gep_PrintExBranch(individual->chromosome + gene_idx * params->gene_size, params->ops_set, &overflow_idx);
  }
  fprintf(stderr, "\n");
}


/*----------------------------------------------------------------------------*/
//!
static void gep_CalculateSizeExBranch(gep_node_t *node, gep_op_set_t *ops_set, int *ops_cnt, int *const_cnt, int *inp_cnt)
{
  if (NULL == node)
  {
    return;
  }

  switch (node->node_type)
  {
    case GEP_NODE__FUNCTION:
      {
        switch (ops_set->ops[node->node.function.operation_type]->arguments_count)
        {
          case 1:
            gep_CalculateSizeExBranch(node->node.function.left, ops_set, ops_cnt, const_cnt, inp_cnt);
            break;

          case 2:
            gep_CalculateSizeExBranch(node->node.function.left,  ops_set, ops_cnt, const_cnt, inp_cnt);
            gep_CalculateSizeExBranch(node->node.function.right, ops_set, ops_cnt, const_cnt, inp_cnt);
            break;

          case 3:
            gep_CalculateSizeExBranch(node->node.function.left,   ops_set, ops_cnt, const_cnt, inp_cnt);
            gep_CalculateSizeExBranch(node->node.function.middle, ops_set, ops_cnt, const_cnt, inp_cnt);
            gep_CalculateSizeExBranch(node->node.function.right,  ops_set, ops_cnt, const_cnt, inp_cnt);
            break;

          default:
            GEP_ABORT("Bad arguments count");
            break;
        }
      }
      *ops_cnt += 1;
      break;

    case GEP_NODE__INPUT:
      *inp_cnt += 1;
      break;

    case GEP_NODE__CONSTANT:
      *const_cnt += 1;
      break;

    default:
      GEP_ABORT("Bad node type");
      break;
  }
}

/*----------------------------------------------------------------------------*/
//!
int GEP_CalculateSize(gep_individual_t *individual, gep_params_t *params)
{
  int ops_cnt = 0, const_cnt = 0, inp_cnt = 0;

  for (int gene_idx = 0; gene_idx < params->genes_count; ++gene_idx)
  {
    gep_CalculateSizeExBranch(individual->chromosome + gene_idx * params->gene_size, params->ops_set, &ops_cnt, &const_cnt, &inp_cnt);
  }

  return (ops_cnt * (3 + 2) + const_cnt * (32 + 2) + inp_cnt * (2 + 2)) / 8;
}

/*----------------------------------------------------------------------------*/
//!
int GEP_CopyChromosome(gep_individual_t *dst, gep_params_t *dst_params, gep_individual_t *src, gep_params_t *src_params)
{
  int ret = -1;

  if (dst_params->gene_size >= src_params->gene_size && dst_params->genes_count == src_params->genes_count)
  {
    int genes_count   =  dst_params->genes_count;
    int dst_gene_size =  dst_params->gene_size;
    int src_gene_size =  src_params->gene_size;
    int min_gene_size = (dst_gene_size < src_gene_size) ? dst_gene_size : src_gene_size;
    int cp_size       =  min_gene_size * NODE_SIZE;

    for (int i = 0; i < genes_count; ++i)
    {
      dst->head_lengths[i]    = src->head_lengths[i];
      dst->working_lengths[i] = src->working_lengths[i];
      memcpy(dst->chromosome + i * dst_gene_size, src->chromosome + i * src_gene_size, cp_size);
    }

    dst->et_valid = 0;
    ret = 0;
  }
  else
  {
    GEP_ABORT("Dst ind is smaller than src");
  }

  return ret;
}

/*----------------------------------------------------------------------------*/
//!
gep_problem_t* GEP_ReadProblemFromFile(const char *filename)
{
  FILE *fid = fopen(filename, "r");
  int input_dimensions_count = 0, samples_count = 1;

  if (1 != fscanf(fid, "%d", &input_dimensions_count))
  {
    fprintf(stderr, "[%s : %d] Error reading file %s", __FILE__, __LINE__, filename);
  }

  int inputs_sizes[input_dimensions_count];

  for (int i = 0; i < input_dimensions_count; ++i)
  {
    if (1 != fscanf(fid, "%d", &inputs_sizes[i]))
    {
      fprintf(stderr, "[%s : %d] Error reading file %s", __FILE__, __LINE__, filename);
    }
    samples_count *= inputs_sizes[i];
  }

  gep_problem_t *problem = gep_StartAllocateProblem(input_dimensions_count, samples_count);

  for (int i = 0; i < input_dimensions_count; ++i)
  {
    if (2 != fscanf(fid, "%lf %lf", &problem->inputs_mins[i], &problem->inputs_maxs[i]))
    {
      fprintf(stderr, "[%s : %d] Error reading file %s", __FILE__, __LINE__, filename);
    }
  }

  problem->output_min = 0;
  if (1 != fscanf(fid, "%lf", &problem->output_scale))
  {
    fprintf(stderr, "[%s : %d] Error reading file %s", __FILE__, __LINE__, filename);
  }

  memcpy(problem->inputs_sizes, inputs_sizes, sizeof(int) * input_dimensions_count);

  for (int i = 0; i < samples_count; ++i)
  {
    if (1 != fscanf(fid, "%lf", &problem->outputs[i]))
    {
      fprintf(stderr, "[%s : %d] Error reading file %s\n", __FILE__, __LINE__, filename);
    }
    problem->outputs[i] /= problem->output_scale;
  }

  fclose(fid);

  //printf("Scanned input file: dims_count = %d; samples_count = %d\n", input_dimensions_count, samples_count);

  double dx = (problem->inputs_maxs[0] - problem->inputs_mins[0]) / (inputs_sizes[0] - 1);

  if (1 == input_dimensions_count)
  {
    for (int x = 0; x < inputs_sizes[0]; ++x)
    {
      double x_val = problem->inputs_mins[0] + dx * x;
      problem->inputs[x] = x_val;
    }
  }
  else if (2 == input_dimensions_count)
  {
    double dy = (problem->inputs_maxs[1] - problem->inputs_mins[1]) / (inputs_sizes[1] - 1);

    for (int y = 0; y < inputs_sizes[1]; ++y)
    {
      double y_val = problem->inputs_mins[0] + dy * y;
      for (int x = 0; x < inputs_sizes[0]; ++x)
      {
        double x_val = problem->inputs_mins[0] + dx * x;
        int out_idx = 2 * (y * problem->inputs_sizes[0] + x);

        problem->inputs[out_idx]     = y_val;
        problem->inputs[out_idx + 1] = x_val;
      }
    }
  }
  else if (3 == input_dimensions_count)
  {
    double dy = (problem->inputs_maxs[1] - problem->inputs_mins[1]) / (inputs_sizes[1] - 1);
    double dz = (problem->inputs_maxs[2] - problem->inputs_mins[2]) / (inputs_sizes[2] - 1);

    for (int y = 0; y < inputs_sizes[1]; ++y)
    {
      double y_val = problem->inputs_mins[0] + dy * y;

      for (int x = 0; x < inputs_sizes[0]; ++x)
      {
        double x_val = problem->inputs_mins[0] + dx * x;

        for (int z = 0; z < inputs_sizes[2]; ++z)
        {
          double z_val = problem->inputs_mins[2] + dz * z;
          int out_idx = 3 * (y * problem->inputs_sizes[0] * problem->inputs_sizes[2] + x * problem->inputs_sizes[2] + z);

          problem->inputs[out_idx]     = y_val;
          problem->inputs[out_idx + 1] = x_val;
          problem->inputs[out_idx + 2] = z_val;
        }
      }
    }
  }

  return problem;
}


/*----------------------------------------------------------------------------*/
//!
void GEP_WriteProblemToFile(const gep_problem_t *problem, const double *points, const char *filename)
{
  FILE *fid = fopen(filename, "w");
  int input_dimensions_count = problem->input_dimensions_count, samples_count = problem->samples_count;

  fprintf(fid, "%d\n", input_dimensions_count);

  if (1 == input_dimensions_count)
  {
    int inp_size = samples_count;
    fprintf(fid, "%d\n", inp_size);
    fprintf(fid, "%lf %lf\n", problem->inputs_mins[0], problem->inputs_maxs[0]);
    fprintf(fid, "%lf\n", problem->output_scale);

    for (int i = 0; i < samples_count; ++i)
    {
      fprintf(fid, "%lf ", points[i]);
    }
    fprintf(fid, "\n");
  }
  else if (2 == input_dimensions_count)
  {
    fprintf(fid, "%d %d\n", problem->inputs_sizes[0], problem->inputs_sizes[1]);
    fprintf(fid, "%lf %lf\n", problem->inputs_mins[0], problem->inputs_maxs[0]);
    fprintf(fid, "%lf %lf\n", problem->inputs_mins[1], problem->inputs_maxs[1]);
    fprintf(fid, "%lf\n", problem->output_scale);

    for (int y = 0; y < problem->inputs_sizes[1]; ++y)
    {
      for (int x = 0; x < problem->inputs_sizes[0]; ++x)
      {
        int z_idx = y * problem->inputs_sizes[0] + x;
        fprintf(fid, "%lf ", points[z_idx] * problem->output_scale);
      }
      fprintf(fid, "\n");
    }
  }
  else if (3 == input_dimensions_count)
  {
    fprintf(fid, "%d %d %d\n", problem->inputs_sizes[0], problem->inputs_sizes[1], problem->inputs_sizes[2]);
    fprintf(fid, "%lf %lf\n", problem->inputs_mins[0], problem->inputs_maxs[0]);
    fprintf(fid, "%lf %lf\n", problem->inputs_mins[1], problem->inputs_maxs[1]);
    fprintf(fid, "%lf %lf\n", problem->inputs_mins[2], problem->inputs_maxs[2]);
    fprintf(fid, "%lf\n", problem->output_scale);

    for (int y = 0; y < problem->inputs_sizes[1]; ++y)
    {
      for (int x = 0; x < problem->inputs_sizes[0]; ++x)
      {
        for (int z = 0; z < problem->inputs_sizes[2]; ++z)
        {
          int out_idx = y * problem->inputs_sizes[0] * problem->inputs_sizes[2] + x * problem->inputs_sizes[2] + z;
          fprintf(fid, "%lf ", points[out_idx] * problem->output_scale);
        }
        fprintf(fid, "   ");
      }
      fprintf(fid, "\n");
    }
  }

  fclose(fid);
}
