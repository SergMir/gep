#include <gep_internal.h>

#include <math.h>

/*----------------------------------------------------------------------------*/
//!
double GEP_MSE(gep_problem_t *problem, gep_individual_t *individual, const gep_params_t *params)
{
  double  mse = 0;
  double *inputs = problem->inputs;
  double *output = problem->outputs;

  for (int i = 0; i < problem->samples_count; ++i)
  {
    double et = GEP_CalculateET(individual, params, inputs);
    double diff = et - *output;

    mse += diff * diff;

    inputs  += problem->input_dimensions_count;
    output  += 1;
  }

  mse /= problem->samples_count;
  mse  = sqrt(mse);

  return mse;
}


/*----------------------------------------------------------------------------*/
//!
double GEP_fitness_MSE(gep_problem_t *problem, gep_individual_t *individual, const gep_params_t *params)
{
  double mse = GEP_MSE(problem, individual, params);

  return GEP_MAX_FITNESS / (1.0 + mse * params->mse_coefficient);
}


/*----------------------------------------------------------------------------*/
//!
double GEP_fitness_MSE_Partial(gep_problem_t *problem, gep_individual_t *individual, const gep_params_t *params, double threshold)
{
  double mse          = 0;
  int    portion_idx  = 0;
  int    portion_size = 10;
  double penalty      = 1;

  threshold *= 0.7;

  for (int i = 0; i < problem->samples_count; ++i)
  {
    double *inputs = &problem->inputs[problem->shuffle[i] * problem->input_dimensions_count];
    double  output =  problem->outputs[problem->shuffle[i]];

    double et = GEP_CalculateET(individual, params, inputs);
    double diff = et - output;

    mse += diff * diff;

    if (++portion_idx == portion_size)
    {
      double test_mse = sqrt(mse / (i + 1.0));
      double current_fit = GEP_MAX_FITNESS / (1.0 + test_mse * params->mse_coefficient);
      if (current_fit < threshold)
      {
        penalty = 1.0 / (problem->samples_count / (((double)i) + 1.0));
        break;
      }
      portion_idx = 0;
    }
  }

  mse = sqrt(mse / problem->samples_count);

  return penalty * GEP_MAX_FITNESS / (1.0 + mse * params->mse_coefficient);
}


/*----------------------------------------------------------------------------*/
//!
double GEP_fitness_Rsquared(gep_problem_t *problem, gep_individual_t *individual, const gep_params_t *params)
{
  double sT    = 0;
  double sT_sq = 0;
  double sP    = 0;
  double sP_sq = 0;
  double sTP   = 0;
  double n = problem->samples_count;

  int viable = 1;

  double *inputs = problem->inputs;
  double *output = problem->outputs;

  for (int i = 0; i < problem->samples_count; ++i)
  {
    double P = GEP_CalculateET(individual, params, inputs);
    double T = *output;

    if (isnan(P) || isinf(P))
    {
      viable = 0;
      break;
    }

    sTP   += T * P;
    sT    += T;
    sT_sq += T * T;
    sP    += P;
    sP_sq += P * P;

    inputs  += problem->input_dimensions_count;
    output  += 1;
  }

  double fitness = NAN;

  if (viable)
  {
    double R = (n * sTP - sT * sP) / sqrt(  (n * sT_sq - sT * sT)  * (n * sP_sq - sP * sP) );
    if (!isnan(R) && !isinf(R))
    {
      fitness = GEP_MAX_FITNESS * R * R;
    }
  }

  return fitness;
}
