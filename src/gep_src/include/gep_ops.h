#ifndef _H_GEP_OPS_
#define _H_GEP_OPS_

#include <gep.h>

#define GEP_MAX_OPS  50
#define GEP_MAX_ARGS 3

/*----------------------------------------------------------------------------*/
//!
typedef double (*gep_op_f_t)(double arg1, double arg2, double arg3);


/*----------------------------------------------------------------------------*/
//!
typedef struct
{
  gep_op_f_t  op_function;
  int         arguments_count;
  const char *op_name;
} gep_op_t;


/*----------------------------------------------------------------------------*/
//!
typedef struct
{
  gep_op_t** ops;
  int        ops_count;
} gep_op_set_t;


gep_op_set_t* GEP_GetOpSet(gep_ops_preset_t preset);

#endif /* _H_GEP_OPS_ */
