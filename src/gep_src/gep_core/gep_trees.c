#include <gep_internal.h>
#include <utils.h>

#include <math.h>


#define TREE_SIZE_LIMIT (3 * params->chromosome_size)


static gep_node_t gep_overlapped_overflow_nodes[GEP_MAX_ARGS] =
{
  { .node_type = GEP_NODE__INPUT,    .node.terminal.input_index = 0 },
  { .node_type = GEP_NODE__CONSTANT, .node.terminal.value       = 1.0 },
  { .node_type = GEP_NODE__CONSTANT, .node.terminal.value       = 0 }
};


/*----------------------------------------------------------------------------*/
//!
static double gep_FunctionProbability(int level, int max_level)
{
  double x = ((double)level - 1.0) / (((double)max_level) - 1.0);
  return sqrt(1.0 - x*x);
}


/*----------------------------------------------------------------------------*/
//!
void gep_CreateRandomTerminal(gep_node_t *node, const gep_params_t *params)
{
  int choise = rand() & (1 + params->input_dimensions_count);
  if (choise < params->input_dimensions_count)
  {
    node->node_type = GEP_NODE__INPUT;
    node->node.terminal.input_index = choise;
  }
  else
  {
    node->node_type = GEP_NODE__CONSTANT;
    node->node.terminal.value = gep_GenerateConstant();
  }
}


/*----------------------------------------------------------------------------*/
//!
void gep_CreateRandomNode(gep_node_t *node, const gep_params_t *params)
{
  int choise = rand() & (params->ops_set->ops_count + params->input_dimensions_count + 1);

  if (choise < params->ops_set->ops_count)
  {
    node->node_type = GEP_NODE__FUNCTION;
    node->node.function.operation_type = choise;
  }
  else if (choise < (params->ops_set->ops_count + params->input_dimensions_count))
  {
    node->node_type = GEP_NODE__INPUT;
    node->node.terminal.input_index = choise - params->ops_set->ops_count;
  }
  else
  {
    node->node_type = GEP_NODE__CONSTANT;
    node->node.terminal.value = gep_GenerateConstant();
  }
}


/*----------------------------------------------------------------------------*/
//!
static void gep_CreateRandomTreeBranch(gep_node_t *gene, int i, int level, int *used, const gep_params_t *params)
{
  double function_probability = gep_FunctionProbability(level, params->max_depth);
  int gene_size = params->gene_size;

  if (i >= gene_size)
  {
    return;
  }

  if (1 == used[i])
  {
    return;
  }

  used[i] = 1;

  gep_node_t *node = &gene[i];

  double f_probe = normal_rand();

  if (f_probe < function_probability)
  {
    node->node_type = GEP_NODE__FUNCTION;
    node->node.function.operation_type = rand() % params->ops_set->ops_count;

    int arguments_count = params->ops_set->ops[node->node.function.operation_type]->arguments_count;

    node->node.function.left   = NULL;
    node->node.function.middle = NULL;
    node->node.function.right  = NULL;

    if (GEP_CODING__OVERLAPPED == params->coding_type)
    {
      switch (arguments_count)
      {
        case 1:
          node->node.function.left   = ((i + 1) < gene_size) ? (&gene[i + 1]) : (&gep_overlapped_overflow_nodes[0]);
          gep_CreateRandomTreeBranch(gene, i + 1, level + 1, used, params);
          break;

        case 2:
          node->node.function.left   = ((i + 1) < gene_size) ? (&gene[i + 1]) : (&gep_overlapped_overflow_nodes[i + 1 - gene_size]);
          node->node.function.right  = ((i + 2) < gene_size) ? (&gene[i + 2]) : (&gep_overlapped_overflow_nodes[i + 2 - gene_size]);
          gep_CreateRandomTreeBranch(gene, i + 1, level + 1, used, params);
          gep_CreateRandomTreeBranch(gene, i + 2, level + 1, used, params);
          break;

        case 3:
          node->node.function.left   = ((i + 1) < gene_size) ? (&gene[i + 1]) : (&gep_overlapped_overflow_nodes[i + 1 - gene_size]);
          node->node.function.middle = ((i + 2) < gene_size) ? (&gene[i + 2]) : (&gep_overlapped_overflow_nodes[i + 2 - gene_size]);
          node->node.function.right  = ((i + 3) < gene_size) ? (&gene[i + 3]) : (&gep_overlapped_overflow_nodes[i + 3 - gene_size]);
          gep_CreateRandomTreeBranch(gene, i + 1, level + 1, used, params);
          gep_CreateRandomTreeBranch(gene, i + 2, level + 1, used, params);
          gep_CreateRandomTreeBranch(gene, i + 3, level + 1, used, params);
          break;

        default:
          GEP_ABORT("Bad arguments count");
          break;
      }
    }
    else if (GEP_CODING__PREFIX == params->coding_type)
    {
      int i_left = i + 1;
      while (i_left < gene_size && 0 != used[i_left])
      {
        ++i_left;
      }

      if (i_left < gene_size)
      {
        node->node.function.left = &gene[i_left];
        gep_CreateRandomTreeBranch(gene, i_left, level + 1, used, params);
      }
      else
      {
        node->node.function.left = &gep_overlapped_overflow_nodes[0];
      }

      if (2 == arguments_count)
      {
        int i_right = i_left + 1;
        while (i_right < gene_size && 0 != used[i_right])
        {
          ++i_right;
        }

        if (i_right < gene_size)
        {
          node->node.function.right = &gene[i_right];
          gep_CreateRandomTreeBranch(gene, i_right, level + 1, used, params);
        }
        else
        {
          node->node.function.right = (i_left < gene_size) ? (&gep_overlapped_overflow_nodes[1]) : (&gep_overlapped_overflow_nodes[0]);
        }
      }
      else if (3 == arguments_count)
      {
        int i_middle = i_left + 1;
        while (i_middle < gene_size && 0 != used[i_middle])
        {
          ++i_middle;
        }

        if (i_middle < gene_size)
        {
          node->node.function.middle = &gene[i_middle];
          gep_CreateRandomTreeBranch(gene, i_middle, level + 1, used, params);
        }
        else
        {
          node->node.function.middle = (i_left < gene_size) ? (&gep_overlapped_overflow_nodes[1]) : (&gep_overlapped_overflow_nodes[0]);
        }

        int i_right = i_middle + 1;
        while (i_right < gene_size && 0 != used[i_right])
        {
          ++i_right;
        }

        if (i_right < gene_size)
        {
          node->node.function.right = &gene[i_right];
          gep_CreateRandomTreeBranch(gene, i_right, level + 1, used, params);
        }
        else
        {
          if (i_left < gene_size)
          {
            node->node.function.right = &gep_overlapped_overflow_nodes[2];
          }
          else if (i_middle < gene_size)
          {
            node->node.function.right = &gep_overlapped_overflow_nodes[1];
          }
          else
          {
            node->node.function.right = &gep_overlapped_overflow_nodes[0];
          }
        }
      }
      else if (1 != arguments_count)
      {
        GEP_ABORT("Bad arguments count");
      }
    }
    else
    {
      GEP_ABORT("Invalid coding type");
    }
  }
  else
  {
    gep_CreateRandomTerminal(node, params);
  }
}


/*----------------------------------------------------------------------------*/
//!
void gep_CreateRandomTree(gep_individual_t *individual, const gep_params_t *params)
{
  if (GEP_CODING__OVERLAPPED == params->coding_type || GEP_CODING__PREFIX == params->coding_type)
  {
    for (int gene_idx = 0; gene_idx < params->genes_count; ++gene_idx)
    {
      int used[params->gene_size + 3];

      individual->head_lengths[gene_idx] = rand() % (params->gene_size - 1) + 1;

      used[0] = 0;
      for (int i = 1; i < params->gene_size + 3; ++i)
      {
        used[i] = 0;
      }

      gep_CreateRandomTreeBranch(&individual->chromosome[gene_idx * params->gene_size], 0, 1, used, params);

      for (int i = 0; i < params->gene_size; ++i)
      {
        if (0 == used[i])
        {
          gep_node_t *node = &individual->chromosome[gene_idx * params->gene_size + i];
          gep_CreateRandomNode(node, params);
        }
      }
    }
  }
  else if (GEP_CODING__BREADTH_FIRST == params->coding_type)
  {
    for (int gene_idx = 0; gene_idx < params->genes_count; ++gene_idx)
    {
      int level_size = 1;
      int i = 0;

      individual->head_lengths[gene_idx] = params->gene_tail_idx;

      for (int depth = 1; depth <= params->max_depth; ++depth)
      {
        double function_probability = gep_FunctionProbability(depth, params->max_depth);
        int i_start = i;
        for (; i < (i_start + level_size); ++i)
        {
          gep_node_t *node = &individual->chromosome[gene_idx * params->gene_size + i];
          if (normal_rand() < function_probability)
          {
            node->node_type = GEP_NODE__FUNCTION;
            node->node.function.operation_type = rand() % params->ops_set->ops_count;
          }
          else
          {
            gep_CreateRandomTerminal(node, params);
          }
        }
        level_size *= GEP_MAX_ARGS;
      }  
    }
  }
  else
  {
    GEP_ABORT("Invalid coding type");
  }

  individual->et_valid = 0;
}


/*----------------------------------------------------------------------------*/
//!
static int gep_CalculateETBranchSize(const gep_node_t *node, const gep_params_t *params)
{
  int ret = 1;

  if (NULL == node)
  {
    GEP_ABORT("Null node");
  }

  switch (node->node_type)
  {
    case GEP_NODE__FUNCTION:
    {
      int arguments_count = params->ops_set->ops[node->node.function.operation_type]->arguments_count;
      ret = gep_CalculateETBranchSize(node->node.function.left, params);

      if (ret < TREE_SIZE_LIMIT)
      {
        switch (arguments_count)
        {
          case 1:
            break;
          case 2:
            ret += gep_CalculateETBranchSize(node->node.function.right,  params);
            break;
          case 3:
            ret += gep_CalculateETBranchSize(node->node.function.middle, params);
            if (ret < TREE_SIZE_LIMIT)
            {
              ret += gep_CalculateETBranchSize(node->node.function.right,  params);
            }
            break;
          default:
            GEP_ABORT("Bad arguments count");
            break;
        }
      }
      break;
    }
    case GEP_NODE__INPUT:
    {
      ret = 1;
      break;
    }
    case GEP_NODE__CONSTANT:
    {
      ret = 1;
      break;
    }
    default:
      GEP_ABORT("Bad node type");
      break;
  }

  return ret;
}


/*----------------------------------------------------------------------------*/
//!
static int gep_RebuildET(gep_individual_t *individual, const gep_params_t *params)
{
  gep_node_t *chromosome = individual->chromosome;

  if (individual->et_valid)
  {
    return 0;
  }

  for (int gene_idx = 0; gene_idx < params->genes_count; ++gene_idx)
  {
    int i = gene_idx * params->gene_size;
    int i_stop = i + params->gene_size;
    int used[params->gene_size + 3];

    used[0] = 1;
    for (int j = 1; j < params->gene_size; ++j)
    {
      used[j] = 0;
    }

    if (GEP_CODING__OVERLAPPED == params->coding_type)
    {
      for (int j = i; j < i_stop; ++j)
      {
        gep_node_t *dst_node = &chromosome[j];
        if (GEP_NODE__FUNCTION == dst_node->node_type)
        {
          int arguments_count = params->ops_set->ops[dst_node->node.function.operation_type]->arguments_count;

          switch (arguments_count)
          {
            case 1:
              dst_node->node.function.left   = ((j + 1) < i_stop) ? (dst_node + 1) : (&gep_overlapped_overflow_nodes[j + 1 - i_stop]);
              dst_node->node.function.middle = NULL;
              dst_node->node.function.right  = NULL;

              if (used[j - i])
              {
                used[j - i + 1] = 1;
              }
              break;

            case 2:
              dst_node->node.function.left   = ((j + 1) < i_stop) ? (dst_node + 1) : (&gep_overlapped_overflow_nodes[j + 1 - i_stop]);
              dst_node->node.function.middle = NULL;
              dst_node->node.function.right  = ((j + 2) < i_stop) ? (dst_node + 2) : (&gep_overlapped_overflow_nodes[j + 2 - i_stop]);

              if (used[j - i])
              {
                used[j - i + 1] = 1;
                used[j - i + 2] = 1;
              }
              break;

            case 3:
              dst_node->node.function.left   = ((j + 1) < i_stop) ? (dst_node + 1) : (&gep_overlapped_overflow_nodes[j + 1 - i_stop]);
              dst_node->node.function.middle = ((j + 2) < i_stop) ? (dst_node + 2) : (&gep_overlapped_overflow_nodes[j + 2 - i_stop]);
              dst_node->node.function.right  = ((j + 3) < i_stop) ? (dst_node + 3) : (&gep_overlapped_overflow_nodes[j + 3 - i_stop]);

              if (used[j - i])
              {
                used[j - i + 1] = 1;
                used[j - i + 2] = 1;
                used[j - i + 3] = 1;
              }
              break;

            default:
              GEP_ABORT("Bad arguments count");
              break;
          }
        }
      }
    }
    else if (GEP_CODING__PREFIX == params->coding_type)
    {
      for (int j = i_stop - 1; j >= i; --j)
      {
        gep_node_t *dst_node = &chromosome[j];
        if (GEP_NODE__FUNCTION == dst_node->node_type)
        {
          int arguments_count = params->ops_set->ops[dst_node->node.function.operation_type]->arguments_count;
          int arguments_left = arguments_count;

          dst_node->node.function.left   = NULL;
          dst_node->node.function.middle = NULL;
          dst_node->node.function.right  = NULL;

          for (int z = j + 1; z < i_stop - 1; ++z)
          {
            if (0 == used[z - i])
            {
              used[z - i] = 1;
              gep_node_t *ptr_node = &chromosome[z];

              if (1 == arguments_count)
              {
                dst_node->node.function.left = ptr_node;
              }
              else if (2 == arguments_count)
              {
                if (2 == arguments_left)
                {
                  dst_node->node.function.left = ptr_node;
                }
                else if (1 == arguments_left)
                {
                  dst_node->node.function.right = ptr_node;
                }
              }
              else if (3 == arguments_count)
              {
                if (3 == arguments_left)
                {
                  dst_node->node.function.left = ptr_node;
                }
                else if (2 == arguments_left)
                {
                  dst_node->node.function.middle = ptr_node;
                }
                else if (1 == arguments_left)
                {
                  dst_node->node.function.right = ptr_node;
                }
              }

              --arguments_left;
              if (0 == arguments_left)
              {
                break;
              }
            } //used
          } //for z

          if (0 != arguments_left)
          {
            switch (arguments_count)
            {
              case 1:
                dst_node->node.function.left = &gep_overlapped_overflow_nodes[0];
                break;

              case 2:
                if (2 == arguments_left)
                {
                  dst_node->node.function.left  = &gep_overlapped_overflow_nodes[0];
                  dst_node->node.function.right = &gep_overlapped_overflow_nodes[1];
                }
                else if (1 == arguments_left)
                {
                  dst_node->node.function.right = &gep_overlapped_overflow_nodes[0];
                }
                break;

              case 3:
                if (3 == arguments_left)
                {
                  dst_node->node.function.left   = &gep_overlapped_overflow_nodes[0];
                  dst_node->node.function.middle = &gep_overlapped_overflow_nodes[1];
                  dst_node->node.function.right  = &gep_overlapped_overflow_nodes[2];
                }
                else if (2 == arguments_left)
                {
                  dst_node->node.function.middle = &gep_overlapped_overflow_nodes[0];
                  dst_node->node.function.right  = &gep_overlapped_overflow_nodes[1];
                }
                else if (1 == arguments_left)
                {
                  dst_node->node.function.right  = &gep_overlapped_overflow_nodes[0];
                }
                break;
            }
          }
        } //GEP_NODE__FUNCTION
      }
    }
    else if (GEP_CODING__BREADTH_FIRST == params->coding_type)
    {
      int j = i + 1;

      do
      {
        gep_node_t *dst_node = &chromosome[i];
        gep_node_t *ptr_node = &chromosome[j];
        if (GEP_NODE__FUNCTION == dst_node->node_type)
        {
          int arguments_count = params->ops_set->ops[dst_node->node.function.operation_type]->arguments_count;

          switch (arguments_count)
          {
            case 1:
              dst_node->node.function.left   = ptr_node;
              dst_node->node.function.middle = NULL;
              dst_node->node.function.right  = NULL;

              used[j - gene_idx * params->gene_size] = 1;
              break;

            case 2:
              dst_node->node.function.left   = ptr_node;
              dst_node->node.function.middle = NULL;
              dst_node->node.function.right  = ptr_node + 1;

              used[j - gene_idx * params->gene_size]     = 1;
              used[j - gene_idx * params->gene_size + 1] = 1;
              break;

            case 3:
              dst_node->node.function.left   = ptr_node;
              dst_node->node.function.middle = ptr_node + 1;
              dst_node->node.function.right  = ptr_node + 2;

              used[j - gene_idx * params->gene_size]     = 1;
              used[j - gene_idx * params->gene_size + 1] = 1;
              used[j - gene_idx * params->gene_size + 2] = 1;
              break;

            default:
              GEP_ABORT("Bad arguments count");
              break;
          }

          j += arguments_count;
        }
        ++i;
      } while (i != j);

      if (i > i_stop || j > i_stop)
      {
        GEP_ABORT("Decoding overflow");
      }
    }
    else
    {
      GEP_ABORT("Invalid coding type");
    }

    int used_cnt = 0;
    for (int j = 0; j < params->gene_size; ++j)
    {
      if (used[j])
      {
        ++used_cnt;
      }
      else
      {
        break;
      }
    }

    individual->working_lengths[gene_idx] = used_cnt;
  }
  individual->et_valid = 1;

  individual->tree_size = (1 == params->genes_count) ? 0 : 1;

  for (int gene_idx = 0; gene_idx < params->genes_count; ++gene_idx)
  {
    int i = gene_idx * params->gene_size;
    individual->tree_size += gep_CalculateETBranchSize(&chromosome[i], params);
  }

  return 0;
}


/*----------------------------------------------------------------------------*/
//!
static double gep_CalculateETBranch(gep_node_t *node, const double inputs[], int *viable, const gep_params_t *params)
{
  double ret = 0;
  
  if (NULL == node)
  {
    GEP_ABORT("Null node");
  }

  if (node->calculated)
  {
    return node->value;
  }

  switch (node->node_type)
  {
    case GEP_NODE__FUNCTION:
    {
      double arg1 = gep_CalculateETBranch(node->node.function.left, inputs, viable, params);
      double arg2 = 0;
      double arg3 = 0;

      if (*viable)
      {
        switch (params->ops_set->ops[node->node.function.operation_type]->arguments_count)
        {
          case 1:
            break;

          case 2:
            arg2 = gep_CalculateETBranch(node->node.function.right, inputs, viable, params);
            break;

          case 3:
            arg2 = gep_CalculateETBranch(node->node.function.middle, inputs, viable, params);
            if (*viable)
            {
              arg3 = gep_CalculateETBranch(node->node.function.right, inputs, viable, params);
            }

          default:
            GEP_ABORT("Bad arguments count");
            break;
        }
      }

      if (*viable)
      {
        ret = params->ops_set->ops[node->node.function.operation_type]->op_function(arg1, arg2, arg3);
        *viable = !isnan(ret) && !isinf(ret);
      }
      break;
    }

    case GEP_NODE__INPUT:
    {
      ret = inputs[node->node.terminal.input_index];
      break;
    }

    case GEP_NODE__CONSTANT:
    {
      ret = node->node.terminal.value;
      break;
    }

    default:
      GEP_ABORT("Bad node type");
      break;
  }

  node->value = ret;
  node->calculated = 1;

  return ret;
}


/*----------------------------------------------------------------------------*/
//!
double GEP_CalculateET(gep_individual_t *individual, const gep_params_t *params, const double inputs[])
{
  double ret = 0;

  if (!individual->et_valid)
  {
    int rebuild_status = gep_RebuildET(individual, params);
    if (0 != rebuild_status)
    {
      GEP_ABORT("Can't rebuild tree");
    }
  }

  if (individual->tree_size > TREE_SIZE_LIMIT)
  {
    individual->viable = 0;
    return 0;
  }

  for (int i = 0; i < params->chromosome_size; ++i)
  {
    individual->chromosome[i].calculated = 0;
  }

  individual->viable = 1;
  for (int gene_idx = 0; gene_idx < params->genes_count; ++gene_idx)
  {
    int node_idx = params->gene_size * gene_idx;
    ret = params->ops_set->ops[params->subtrees_operation]->op_function(ret,  gep_CalculateETBranch(&individual->chromosome[node_idx], inputs, &individual->viable, params), 0);
    if (!individual->viable)
    {
      break;
    }
  }

  return ret;
}
