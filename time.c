#include <stdlib.h>
#include <sys/times.h>
#include <unistd.h>


static clock_t ticks_per_sec=1;
static clock_t start_time=0;

void init_run_time(void)
{
	ticks_per_sec=sysconf(_SC_CLK_TCK);
	start_time=times(NULL);
}

double get_run_time(void)
{
	clock_t curr_time=times(NULL);
	return (curr_time-start_time)*1.0/ticks_per_sec;
}


