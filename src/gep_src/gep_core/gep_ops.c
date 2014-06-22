#include <gep_ops.h>
#include <gep_internal.h>

#include <math.h>

/*----------------------------------------------------------------------------*/
//!
typedef enum
{
  GEP_OP__SUM = 0,
  GEP_OP__SUB,
  GEP_OP__MUL,
  GEP_OP__DIV,
  GEP_OP__INV,
  GEP_OP__SQR,
  GEP_OP__SQRT,
  GEP_OP__CUBE,
  GEP_OP__EXP,
  GEP_OP__POW,
  GEP_OP__LN,
  GEP_OP__LOG2,
  GEP_OP__LOG10,
  GEP_OP__SIN,
  GEP_OP__COS,
  GEP_OP__TAN,
  GEP_OP__SINH,
  GEP_OP__COSH,
  GEP_OP__TANH,
  GEP_OP_LAST,
} gep_ops_t;

/*----------------------------------------------------------------------------*/
//!
typedef struct
{
  gep_ops_t *ops_list;
  int        ops_count;
} gep_ops_list_t;

/*----------------------------------------------------------------------------*/
//!
static gep_ops_t ops_list_short[] =
{
  GEP_OP__SUM,
  GEP_OP__SUB,
  GEP_OP__MUL,
  GEP_OP__DIV
};

/*----------------------------------------------------------------------------*/
//!
static gep_ops_t ops_list_middle[] =
{
  GEP_OP__SUM,
  GEP_OP__SUB,
  GEP_OP__MUL,
  GEP_OP__DIV,
  GEP_OP__INV,
  GEP_OP__SQRT,
  GEP_OP__EXP,
};

/*----------------------------------------------------------------------------*/
//!
static gep_ops_t ops_list_big[] =
{
  GEP_OP__SUM,
  GEP_OP__SUB,
  GEP_OP__MUL,
  GEP_OP__DIV,
  GEP_OP__INV,
  GEP_OP__SQRT,
  GEP_OP__EXP,
  GEP_OP__POW,
  GEP_OP__LN,
  GEP_OP__SIN,
  GEP_OP__COS
};

/*----------------------------------------------------------------------------*/
//!
static gep_ops_t ops_list_full[] =
{
  GEP_OP__SUM,
  GEP_OP__SUB,
  GEP_OP__MUL,
  GEP_OP__DIV,
  GEP_OP__INV,
  GEP_OP__SQR,
  GEP_OP__SQRT,
  GEP_OP__CUBE,
  GEP_OP__EXP,
  GEP_OP__POW,
  GEP_OP__LN,
  GEP_OP__LOG2,
  GEP_OP__LOG10,
  GEP_OP__SIN,
  GEP_OP__COS,
  GEP_OP__TAN,
  GEP_OP__SINH,
  GEP_OP__COSH,
  GEP_OP__TANH,
};

/*----------------------------------------------------------------------------*/
//!
static gep_ops_list_t ops_lists[GEP_OPLIST_LAST] =
{
  {ops_list_short,  sizeof(ops_list_short)  / sizeof(ops_list_short[0])},
  {ops_list_middle, sizeof(ops_list_middle) / sizeof(ops_list_middle[0])},
  {ops_list_big,    sizeof(ops_list_big)    / sizeof(ops_list_big[0])},
  {ops_list_full,   sizeof(ops_list_full)   / sizeof(ops_list_full[0])}
};

/*----------------------------------------------------------------------------*/
//!
static double gep_function_sum(double arg1, double arg2, double arg3)
{
  arg3 = arg3; /* Anti-warning */
  return arg1 + arg2;
}

/*----------------------------------------------------------------------------*/
//!
static double gep_function_sub(double arg1, double arg2, double arg3)
{
  arg3 = arg3; /* Anti-warning */
  return arg1 - arg2;
}

/*----------------------------------------------------------------------------*/
//!
static double gep_function_mul(double arg1, double arg2, double arg3)
{
  arg3 = arg3; /* Anti-warning */
  return arg1 * arg2;
}

/*----------------------------------------------------------------------------*/
//!
static double gep_function_div(double arg1, double arg2, double arg3)
{
  arg3 = arg3; /* Anti-warning */
  return arg1 / arg2;
}

/*----------------------------------------------------------------------------*/
//!
static double gep_function_inv(double arg1, double arg2, double arg3)
{
  arg2 = arg2; /* Anti-warning */
  arg3 = arg3; /* Anti-warning */
  return -arg1;
}

/*----------------------------------------------------------------------------*/
//!
static double gep_function_sqr(double arg1, double arg2, double arg3)
{
  arg2 = arg2; /* Anti-warning */
  arg3 = arg3; /* Anti-warning */
  return arg1 * arg1;
}

/*----------------------------------------------------------------------------*/
//!
static double gep_function_sqrt(double arg1, double arg2, double arg3)
{
  arg2 = arg2; /* Anti-warning */
  arg3 = arg3; /* Anti-warning */
  return sqrt(arg1);
}

/*----------------------------------------------------------------------------*/
//!
static double gep_function_cube(double arg1, double arg2, double arg3)
{
  arg2 = arg2; /* Anti-warning */
  arg3 = arg3; /* Anti-warning */
  return arg1 * arg1 * arg1;
}

/*----------------------------------------------------------------------------*/
//!
static double gep_function_exp(double arg1, double arg2, double arg3)
{
  arg2 = arg2; /* Anti-warning */
  arg3 = arg3; /* Anti-warning */
  double ret = exp(arg1);
  return ret;
}

/*----------------------------------------------------------------------------*/
//!
static double gep_function_pow(double arg1, double arg2, double arg3)
{
  arg3 = arg3; /* Anti-warning */
  double ret = pow(arg1, arg2);
  return ret;
}

/*----------------------------------------------------------------------------*/
//!
static double gep_function_ln(double arg1, double arg2, double arg3)
{
  arg2 = arg2; /* Anti-warning */
  arg3 = arg3; /* Anti-warning */
  double ret = log(arg1);
  return ret;
}

/*----------------------------------------------------------------------------*/
//!
static double gep_function_log2(double arg1, double arg2, double arg3)
{
  arg2 = arg2; /* Anti-warning */
  arg3 = arg3; /* Anti-warning */
  double ret = log2(arg1);
  return ret;
}

/*----------------------------------------------------------------------------*/
//!
static double gep_function_log10(double arg1, double arg2, double arg3)
{
  arg2 = arg2; /* Anti-warning */
  arg3 = arg3; /* Anti-warning */
  double ret = log10(arg1);
  return ret;
}

/*----------------------------------------------------------------------------*/
//!
static double gep_function_sin(double arg1, double arg2, double arg3)
{
  arg2 = arg2; /* Anti-warning */
  arg3 = arg3; /* Anti-warning */
  double ret = sin(arg1);
  return ret;
}

/*----------------------------------------------------------------------------*/
//!
static double gep_function_cos(double arg1, double arg2, double arg3)
{
  arg2 = arg2; /* Anti-warning */
  arg3 = arg3; /* Anti-warning */
  double ret = cos(arg1);
  return ret;
}

/*----------------------------------------------------------------------------*/
//!
static double gep_function_tan(double arg1, double arg2, double arg3)
{
  arg2 = arg2; /* Anti-warning */
  arg3 = arg3; /* Anti-warning */
  double ret = tan(arg1);
  return ret;
}

/*----------------------------------------------------------------------------*/
//!
static double gep_function_sinh(double arg1, double arg2, double arg3)
{
  arg2 = arg2; /* Anti-warning */
  arg3 = arg3; /* Anti-warning */
  double ret = sinh(arg1);
  return ret;
}

/*----------------------------------------------------------------------------*/
//!
static double gep_function_cosh(double arg1, double arg2, double arg3)
{
  arg2 = arg2; /* Anti-warning */
  arg3 = arg3; /* Anti-warning */
  double ret = cosh(arg1);
  return ret;
}

/*----------------------------------------------------------------------------*/
//!
static double gep_function_tanh(double arg1, double arg2, double arg3)
{
  arg2 = arg2; /* Anti-warning */
  arg3 = arg3; /* Anti-warning */
  double ret = tanh(arg1);
  return ret;
}

/*----------------------------------------------------------------------------*/
//!
static gep_op_t gep_all_ops[GEP_OP_LAST] =
{
  { gep_function_sum,    2, "+"     },
  { gep_function_sub,    2, "-"     },
  { gep_function_mul,    2, "*"     },
  { gep_function_div,    2, "/"     },
  { gep_function_inv,    1, "-"     },
  { gep_function_sqr,    1, "sqr"   },
  { gep_function_sqrt,   1, "sqrt"  },
  { gep_function_cube,   1, "cube"  },
  { gep_function_exp,    1, "exp"   },
  { gep_function_pow,    2, "^"     },
  { gep_function_ln,     1, "ln"    },
  { gep_function_log2,   1, "log2"  },
  { gep_function_log10,  1, "log10" },
  { gep_function_sin,    1, "sin"   },
  { gep_function_cos,    1, "cos"   },
  { gep_function_tan,    1, "tan"   },
  { gep_function_sinh,   1, "sinh"  },
  { gep_function_cosh,   1, "cosh"  },
  { gep_function_tanh,   1, "tanh"  },
};

/*----------------------------------------------------------------------------*/
//!
gep_op_set_t* GEP_GetOpSet(gep_ops_preset_t preset)
{
  gep_op_set_t* op_set = malloc(sizeof(gep_op_set_t));
  op_set->ops_count    = ops_lists[preset].ops_count;
  op_set->ops          = malloc(op_set->ops_count * sizeof(gep_op_t*));
  

  for (int i = 0; i < op_set->ops_count; ++i)
  {
    op_set->ops[i] = &gep_all_ops[ops_lists[preset].ops_list[i]];
  }

  return op_set;
}
