#include <utils.h>

#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

/*----------------------------------------------------------------------------*/
//!
double normal_rand(void)
{
  return (rand() % 100000) / 100000.0;
}


/*----------------------------------------------------------------------------*/
//!
unsigned long getTimeNs(void)
{
  struct timeval  tv;
  gettimeofday(&tv, NULL);
  unsigned long ret = tv.tv_sec * 1000000 + tv.tv_usec;
  return ret;
}


/*----------------------------------------------------------------------------*/
//!
double getTimeDiffMs(unsigned long time_start, unsigned long time_stop)
{
  double ret = (time_stop - time_start) / 1000;
  return ret;
}
