#include <stdio.h>
#include <math.h>

#define DIM_ZISE 32

/*----------------------------------------------------------------------------*/
//!
typedef double (*genexpr_function)(double arg1, double arg2, double arg3);

/*----------------------------------------------------------------------------*/
//!
typedef struct
{
	genexpr_function  f;
	int               dimensions;
	double            begin;
	double            end;
	double            f_max;
	const char       *name;
} tests_fun_t;

/*----------------------------------------------------------------------------*/
//!
static double fun_2D_sin(double x, double unused_0, double unused_1)
{
	unused_0 = unused_0;
	unused_1 = unused_1;

	return sin(x);
}

/*----------------------------------------------------------------------------*/
//!
static double fun_2D_ferr_1(double x, double unused_0, double unused_1)
{
	unused_0 = unused_0;
	unused_1 = unused_1;

	return x*x*x*x + x*x*x + x*x + x;
}

/*----------------------------------------------------------------------------*/
//!
static double fun_3D_OGR(double x, double y, double unused)
{
	unused = unused;

	double z1 = 0.5  * exp(- (pow(9.0 * x - 2.0, 2.0) + pow(9.0 * y - 2.0, 2.0)) / 4.0);
	double z2 = 0.75 * exp(-  pow(9.0 * x + 1.0, 2.0) / 49.0 - pow(9.0 * y + 1.0, 2.0) / 10.0);
	double z3 = 0.5  * exp(- (pow(9.0 * x - 7.0, 2.0) + pow(9.0 * y - 3.0, 2.0)) / 4.0);
	double z4 = 0.2  * exp(-  pow(9.0 * x - 4.0, 2.0) - pow(9.0 * y - 7.0, 2.0));

	return z1 + z2 + z3 - z4;
}

/*----------------------------------------------------------------------------*/
//!
static double fun_3D_OGR2(double x, double y, double unused)
{
	//return 0.6 * x + 0.3 * y + exp(1);
	return 5 * x + y * (y - 4) + x * y;
}

/*----------------------------------------------------------------------------*/
//!
static double fun_3D_gabor(double x, double y, double unused)
{
	return (3.1416 / 2.0) * exp(-2.0 * (x * x + y * y)) * cos(2 * 3.1416 * (x + y)) + exp(1);
}

/*----------------------------------------------------------------------------*/
//!
static double fun_3D_Rosenbrock(double x, double y, double unused)
{
	return (1 - x) * (1 - x) + 100.0 * (y - x*x) * (y - x*x);
}

/*----------------------------------------------------------------------------*/
//!
static tests_fun_t function_table[] =
{
	{fun_2D_ferr_1,     2,    1,   5,  1000, "test_ferr_1"},
	{fun_2D_sin,        2,   -5, 4.6,     1, "test_sin"},
	{fun_3D_OGR,        3,   -2,   2,   0.8, "test_OGR1"},
	{fun_3D_OGR2,       3, -100, 100, 50000, "test_OGR2"},
	{fun_3D_gabor,      3,    0,   1,    10, "test_gabor"},
	{fun_3D_Rosenbrock, 3,   -3,   3, 10000, "test_Rosenbrock"},
};

/*----------------------------------------------------------------------------*/
//!
static void generate_test(tests_fun_t *test)
{
	FILE *fid = fopen(test->name, "w");

	genexpr_function f = test->f;
	int dims_count     = test->dimensions - 1;
	int dim_size       = DIM_ZISE;
	double begin       = test->begin;
	double end         = test->end;
	double f_max       = test->f_max;

	double delta = (end - begin) / ((double)(dim_size - 1));

	fprintf(fid, "%d\n",  dims_count);
	if (2 == dims_count)
	{
		fprintf(fid, "%d %d\n",   dim_size, dim_size);
		fprintf(fid, "%lf %lf\n", begin,    end);
		fprintf(fid, "%lf %lf\n", begin,    end);
	}
	else if (1 == dims_count)
	{
		fprintf(fid, "%d\n",      dim_size);
		fprintf(fid, "%lf %lf\n", begin, end);
	}
	fprintf(fid, "1.0\n");

	if (2 == dims_count)
	{
		for (double y = begin; y < end; y += delta)
		{
			for (double x = begin; x < end; x += delta)
			{
				fprintf(fid, "%lf ", f(x, y, 0) / f_max);
			}
			fprintf(fid, "\n");
		}
	}
	else if (1 == dims_count)
	{
		for (double x = begin; x < end; x += delta)
		{
			fprintf(fid, "%lf ", f(x, 0, 0) / f_max);
		}
		fprintf(fid, "\n");
	}

	fclose(fid);
}

/*----------------------------------------------------------------------------*/
//!
int main(void)
{
	int tests_count = sizeof(function_table) / sizeof(tests_fun_t);

	for (int i = 0; i < tests_count; ++i)
	{
		generate_test(&function_table[i]);
	}

	return 0;
}
