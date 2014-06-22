#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/*----------------------------------------------------------------------------*/
//!
int main(int argc, char *argv[])
{
  if (argc != 5)
  {
    fprintf(stderr, "Invalid usage\n");
    return -1;
  }

  int status = 0;

  FILE *fid_in_1 = fopen(argv[2], "r");
  FILE *fid_in_2 = fopen(argv[3], "r");
  FILE *fid_out  = fopen(argv[4], "w");

  int inputs_count = 0, samples_count = 1;
  double scale;

  status |= (1 != fscanf(fid_in_1, "%d", &inputs_count));
  status |= (1 != fscanf(fid_in_2, "%d", &inputs_count));
  fprintf(fid_out, "%d\n", inputs_count);

  int inputs_sizes[inputs_count];

  for (int i = 0; i < inputs_count; ++i)
  {
    int inp_size = 0;
    status |= (1 != fscanf(fid_in_1, "%d", &inp_size));
    status |= (1 != fscanf(fid_in_2, "%d", &inp_size));
    inputs_sizes[i] = inp_size;  
    fprintf(fid_out, "%d ", inp_size);
    samples_count *= inp_size;
  }
  status |= (1 != fscanf(fid_in_1, "%lf", &scale));
  status |= (1 != fscanf(fid_in_2, "%lf", &scale));

  double *points_1   = malloc(sizeof(double) * samples_count);
  double *points_2   = malloc(sizeof(double) * samples_count);
  double *points_out = malloc(sizeof(double) * samples_count);
  double out_min = 100000, out_max = -out_min;

  double mse = 0;

  for (int i = 0; i < samples_count; ++i)
  {
    status |= (1 != fscanf(fid_in_1, "%lf", &points_1[i]));
    status |= (1 != fscanf(fid_in_2, "%lf", &points_2[i]));

    points_out[i] = ('+' == argv[1][0]) ? (points_1[i] + points_2[i]) : (points_1[i] - points_2[i]);
    if (points_out[i] > out_max)
    {
      out_max = points_out[i];
    }
    if (points_out[i] < out_min)
    {
      out_min = points_out[i];
    }
  }

  //fprintf(fid_out, "\n%lf\n", out_max - out_min);
  fprintf(fid_out, "\n1.0\n");

  if (1 == inputs_count)
  {
    for (int i = 0; i < samples_count; ++i)
    {
      double diff = ('+' == argv[1][0]) ? (points_1[i] + points_2[i]) : (points_1[i] - points_2[i]);
      fprintf(fid_out, "%lf ", diff);

      mse += diff * diff;
    }
    fprintf(fid_out, "\n");
  }
  else if (2 == inputs_count)
  {
    for (int y = 0; y < inputs_sizes[1]; ++y)
    {
      for (int x = 0; x < inputs_sizes[0]; ++x)
      {
        int z_idx = y * inputs_sizes[0] + x;
        fprintf(fid_out, "%lf ", points_out[z_idx]);

        double diff = points_1[z_idx] - points_2[z_idx];
        mse += diff * diff;
      }
      fprintf(fid_out, "\n");
    }
  }
  else if (3 == inputs_count)
  {
    for (int y = 0; y < inputs_sizes[1]; ++y)
    {
      for (int x = 0; x < inputs_sizes[0]; ++x)
      {
        for (int z = 0; z < inputs_sizes[2]; ++z)
        {
          int out_idx = y * inputs_sizes[0] * inputs_sizes[2] + x * inputs_sizes[2] + z;
          fprintf(fid_out, "%lf ", points_out[out_idx]);

          double diff = points_1[out_idx] - points_2[out_idx];
          mse += diff * diff;
        }
        fprintf(fid_out, "   ");
      }
      fprintf(fid_out, "\n");
    }
  }

  mse /= samples_count;
  mse = sqrt(mse);

  #if 0
    if ('-' == argv[1][0])
    {
      double psnr = 10 * log10( scale * scale / mse);
      printf("%s mse=%lf (psnr=%lf)\n", argv[4], mse, psnr);
    }
  #endif

  fclose(fid_in_1);
  fclose(fid_in_2);
  fclose(fid_out);

  return status;
}