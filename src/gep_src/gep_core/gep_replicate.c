#include <gep_internal.h>
#include <utils.h>

/*----------------------------------------------------------------------------*/
//!
void gep_Mutate(gep_individual_t *individual, const gep_params_t *params)
{
  for (int j = 0; j < params->genes_count; ++j)
  {
    for (int i = 0; i < params->gene_size; ++i)
    {
      gep_node_t *node = &individual->chromosome[j * params->gene_size + i];
      double mut_rate = params->mutation_rate;
      if (normal_rand() < mut_rate)
      {
        if (GEP_CODING__BREADTH_FIRST == params->coding_type)
        {
          if (i < individual->head_lengths[j])
          {
            gep_CreateRandomNode(node, params);
          }
          else
          {
            gep_CreateRandomTerminal(node, params);
          }
        }
        else
        {
          gep_CreateRandomNode(node, params);
        }
      }
    }
  }
  individual->et_valid = 0;
}


/*----------------------------------------------------------------------------*/
//!
void gep_ConstantsMutate(gep_individual_t *individual, const gep_params_t *params)
{
  if (params->use_dynamic_constants)
  {
    for (int i = 0; i < params->chromosome_size; ++i)
    {
      if (GEP_NODE__CONSTANT == individual->chromosome[i].node_type)
      {
        double *value   = &individual->chromosome[i].node.terminal.value;
        double const_V  = gep_GenerateConstant();
        int    const_OP = rand() % 4;

        switch (const_OP)
        {
          case 0:
            *value += const_V;
            break;

          case 1:
            *value -= const_V;
            break;

          case 2:
            *value *= const_V;
            break;

          case 3:
            *value /= const_V;
            break;
        }
      }
    }
  }
  else
  {
    for (int i = 0; i < params->chromosome_size; ++i)
    {
      if (GEP_NODE__CONSTANT == individual->chromosome[i].node_type)
      {
        double mutation = ((rand() % 100) - 50) * 0.1 * params->constants_mutation_intensity;
        individual->chromosome[i].node.terminal.value *= (1.0 + mutation);
      }
    }
  }

  individual->et_valid = 0;
}


/*----------------------------------------------------------------------------*/
//!
void gep_Transpose_Inversion(gep_individual_t *individual, const gep_params_t *params)
{
  int gene_idx  = rand() % params->genes_count;
  int start_idx = rand() % params->gene_size;
  int len = rand() % (params->gene_size - start_idx);

  if (GEP_CODING__BREADTH_FIRST == params->coding_type)
  {
    int tail_idx = individual->head_lengths[gene_idx];
    start_idx = rand() % tail_idx;
    len = 0;

    if (start_idx < tail_idx)
    {
      len = rand() % (tail_idx - start_idx);
    }
    else
    {
      len = rand() % (params->gene_size - start_idx);
    }
  }

  gep_node_t *start = individual->chromosome + gene_idx * params->gene_size + start_idx;

  for (int i = 0; i < (len / 2); ++i)
  {
    gep_node_t *node_left  = start + i;
    gep_node_t *node_right = start + len - i - 1;
    gep_node_t node_swap;
    memcpy(&node_swap,  node_left,  NODE_SIZE);
    memcpy(node_left,   node_right, NODE_SIZE);
    memcpy(node_right,  &node_swap, NODE_SIZE);
  }

  individual->et_valid = 0;
}


/*----------------------------------------------------------------------------*/
//!
void gep_Transpose_IS(gep_individual_t *individual, const gep_params_t *params)
{
  int gene_idx       =  rand() %  params->genes_count;
  gep_node_t *start = individual->chromosome + gene_idx * params->gene_size;

  int trans_idx  = (rand() % (params->gene_size - 1)) + 1;
  int trans_len  =  rand() % (params->gene_size - trans_idx);
  int target_idx =  rand() % trans_idx;

  if (target_idx > 0 && target_idx != trans_idx)
  {
    gep_node_t trans[trans_len];
    memcpy(trans,             start + trans_idx,   NODE_SIZE * trans_len);

    int shift_len  = params->gene_size - (target_idx + trans_len);
    if (GEP_CODING__BREADTH_FIRST == params->coding_type)
    {
      shift_len  = individual->head_lengths[gene_idx] - (target_idx + trans_len);
    }

    if (shift_len > 0)
    {
      gep_node_t shift[shift_len];
      memcpy(shift, start + target_idx, NODE_SIZE * shift_len);
      memcpy(start + target_idx + trans_len, shift, NODE_SIZE * shift_len);
    }
    memcpy(start + target_idx, trans, NODE_SIZE * trans_len);
  }

  individual->et_valid = 0;
}


/*----------------------------------------------------------------------------*/
//!
void gep_Transpose_RIS(gep_individual_t *individual, const gep_params_t *params)
{
  int gene_idx       =  rand() %  params->genes_count;
  gep_node_t *start = individual->chromosome + gene_idx * params->gene_size;

  int tail_idx = individual->head_lengths[gene_idx];
  int trans_idx  = (rand() % (tail_idx - 1)) + 1;

  while (GEP_NODE__FUNCTION != (start + trans_idx)->node_type
    && trans_idx < tail_idx)
  {
    ++trans_idx;
  }

  if (GEP_NODE__FUNCTION == (start + trans_idx)->node_type)
  {
    int trans_len  =  rand() % (params->gene_size - trans_idx);
    gep_node_t trans[trans_len];
    memcpy(trans,             start + trans_idx,   NODE_SIZE * trans_len);

    int shift_len  = tail_idx - trans_len;

    if (shift_len > 0)
    {
      gep_node_t shift[shift_len];
      memcpy(shift,             start, NODE_SIZE * shift_len);
      memcpy(start + trans_len, shift, NODE_SIZE * shift_len);
    }

    memcpy(start, trans, NODE_SIZE * trans_len);
  }

  individual->et_valid = 0;
}


/*----------------------------------------------------------------------------*/
//!
void gep_Transpose_Genes(gep_individual_t *individual, const gep_params_t *params)
{
  if (params->genes_count > 1)
  {
    int gene = (rand() % (params->genes_count - 1)) + 1;
    gep_node_t transponon[params->gene_size];
    memcpy(transponon, individual->chromosome + (gene * params->gene_size), NODE_SIZE * params->gene_size);

    for (int i = 0; i < gene; ++i)
    {
      gep_node_t *src_ptr = individual->chromosome + (i * params->gene_size);
      gep_node_t *dst_ptr = individual->chromosome + ((i + 1) * params->gene_size);
      memcpy(dst_ptr, src_ptr, NODE_SIZE * params->gene_size);
    }
    memcpy(individual->chromosome, transponon, NODE_SIZE * params->gene_size);
  }

  individual->et_valid = 0;
}


/*----------------------------------------------------------------------------*/
//!
void gep_Cross_1Point(gep_individual_t *ind1, gep_individual_t *ind2, const gep_params_t *params)
{
  int point;
  do
  {
    point = rand() % (params->chromosome_size - 1);
  } while (0 == point);

  int size = point * NODE_SIZE;
  gep_node_t buffer[point];

  memcpy(buffer, ind2->chromosome, size);
  memcpy(ind2->chromosome, ind1->chromosome, size);
  memcpy(ind1->chromosome, buffer, size);

  ind1->et_valid = 0;
  ind2->et_valid = 0;
}


/*----------------------------------------------------------------------------*/
//!
void gep_Cross_2Points(gep_individual_t *ind1, gep_individual_t *ind2, const gep_params_t *params)
{
  int point1;
  int point2;

  do
  {
    point1 = rand() % (params->chromosome_size - 1);
  } while (point1 == 0);

  do
  {
    point2 = rand() % (params->chromosome_size - 1);
  } while (point2 == 0 || point2 == point1);

  if (point1 > point2)
  {
    int swap = point1;
    point1 = point2;
    point2 = swap;
  }
  int size = (point2 - point1) * NODE_SIZE;
  gep_node_t buffer[point2 - point1];

  memcpy(buffer, ind2->chromosome + point1, size);
  memcpy(ind2->chromosome + point1, ind1->chromosome + point1, size);
  memcpy(ind1->chromosome + point1, buffer, size);

  ind1->et_valid = 0;
  ind2->et_valid = 0;
}


/*----------------------------------------------------------------------------*/
//!
void gep_Cross_Genes(gep_individual_t *ind1, gep_individual_t *ind2, const gep_params_t *params)
{
  int gene1_idx = rand() % params->genes_count;
  int gene2_idx = rand() % params->genes_count;
  int point1 = gene1_idx * params->gene_size;
  int point2 = gene2_idx * params->gene_size;

  int size = params->gene_size * NODE_SIZE;
  gep_node_t buffer[params->gene_size];

  memcpy(buffer, ind2->chromosome + point2, size);
  memcpy(ind2->chromosome + point2, ind1->chromosome + point1, size);
  memcpy(ind1->chromosome + point1, buffer, size);

  ind1->et_valid = 0;
  ind2->et_valid = 0;
}


/*----------------------------------------------------------------------------*/
//!
void gep_Replicate(gep_individual_t *individual, const gep_params_t *params)
{
  gep_Mutate(individual, params);

  double constants_mutation_probe = normal_rand();
  if (constants_mutation_probe < params->constants_mutation_rate)
  {
    gep_ConstantsMutate(individual, params);
  }

  double inversion_probe = normal_rand();
  if (inversion_probe < params->inversion_rate)
  {
    gep_Transpose_Inversion(individual, params);
  }

  double IS_probe = normal_rand();
  if (IS_probe < params->IS_transposition_rate)
  {
    gep_Transpose_IS(individual, params);
  }

  if (GEP_CODING__BREADTH_FIRST == params->coding_type)
  {
    double RIS_probe = normal_rand();
    if (RIS_probe < params->RIS_transposition_rate)
    {
      gep_Transpose_RIS(individual, params);
    }
  }

  double genes_probe = normal_rand();
  if (genes_probe < params->genes_transposition_rate)
  {
    gep_Transpose_Genes(individual, params);
  }
}


/*----------------------------------------------------------------------------*/
//!
void gep_Recombinate(gep_individual_t *ind1, gep_individual_t *ind2, const gep_params_t *params)
{
  double recombination_probe = normal_rand();

  if (recombination_probe < params->one_point_recombination_rate)
  {
    gep_Cross_1Point(ind1, ind2, params);
  }
  else if (recombination_probe < (params->one_point_recombination_rate + params->two_points_recombination_rate))
  {
    gep_Cross_2Points(ind1, ind2, params);
  }
  else
  {
    gep_Cross_Genes(ind1, ind2, params);
  }
}
