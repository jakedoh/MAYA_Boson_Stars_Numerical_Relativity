
#include <stdio.h> 
#include <string.h> 
#include <sys/resource.h>
#include <unistd.h>
#include <malloc.h>
#include <limits.h>

#include "cctk.h" 
#include "cctk_Arguments.h" 
#include "cctk_Parameters.h" 

// a "safe" downcast to int
static CCTK_INT int_limit(size_t x)
{
  if(((CCTK_INT)x) != x) // kludge to detect if CCTK_INT can hold x
    return INT_MAX;
  else
    return (CCTK_INT)x;
}

size_t get_rss()
{
  size_t size; //       total program size
  size_t resident;//   resident set size
  size_t share;//      shared pages
  size_t text;//       text (code)
  size_t lib;//        library
  size_t data;//       data/stack
  size_t dt;//         dirty pages (unused in Linux 2.6)

  unsigned int page_size = (unsigned int)sysconf(_SC_PAGESIZE);

  char buf[30];
  snprintf(buf, 30, "/proc/%u/statm", (unsigned)getpid());
  FILE* pf = fopen(buf, "r");
  if (pf) {
    fscanf(pf, "%zu %zu %zu %zu %zu %zu", &size, &resident, &share, &text, &lib, &data);
  }
  fclose(pf);
  return resident * page_size;
}

size_t get_majflt()
{
  int pid;
  char exe[256];
  char  state;
  unsigned long int majflt;

  unsigned int page_size = (unsigned int)sysconf(_SC_PAGESIZE);

  char buf[30];
  snprintf(buf, 30, "/proc/%u/stat", (unsigned)getpid());
  FILE* pf = fopen(buf, "r");
  if (pf) {
    fscanf(pf, "%d %s %c %*d %*d %*d %*d %*d %*lu %*lu %*lu %lu", &pid, exe, &state, &majflt);
  }
  fclose(pf);
  return majflt * page_size;
}

void get_node_stats(size_t *total, size_t *freemem, size_t *active, size_t *inactive, size_t *swaptotal, 
                    size_t *swapfree)
{
  FILE *f = fopen("/proc/meminfo", "r");
  const int buf_len = 100;
  char buffer[buf_len];

  if (f == 0)
  {
    printf("Could not open /proc/meminfo\n");
    return;
  }

  while(!feof(f))
  {
    fgets(buffer, buf_len, f);
    if (!feof(f))
    {
      char key[100];
      char unit[100];
      size_t val, factor;
      int fields;

      val = 0.;factor = 1.;
      fields = sscanf(buffer, "%s %zu %s", key, &val, unit);
      if(fields == 3)
      {
        if(strcmp(unit, "kB") == 0)
          factor = 1000UL;
        else if(strcmp(unit, "MB") == 0)
          factor = 1000000UL;
        else if(strcmp(unit, "GB") == 0)
          factor = 1000000000UL;
        else
          CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
	    "Could not parse /proc/meminfo line '%s', unknown unit '%s'\n", buffer, unit);
      }
      else if(fields == 2) 
      {
        factor = 1.;
      }
      else
      {
        factor = 1.;
        CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
	  "Could not parse /proc/meminfo line '%s', read %d fields\n", buffer, fields);
      }
      val *= factor; // apply unit

      if (strcmp(key, "MemTotal:") == 0)
      {
        if (total) *total = val;
      }
      else if (strcmp(key, "MemFree:") == 0)
      {
        if (freemem) *freemem = val;
      }
      else if (strcmp(key, "Active:") == 0)
      {
        if (active) *active = val;
      }
      else if (strcmp(key, "Inactive:") == 0)
      {
        if (inactive) *inactive = val;
      }
      else if (strcmp(key, "SwapTotal:") == 0)
      {
        if (swaptotal) *swaptotal = val;
      }
      else if (strcmp(key, "SwapFree:") == 0)
      {
        if (swapfree) *swapfree = val;
      }
    }
  }
  fclose(f);
}

FILE *open_output_file(char *leaf, int *first_time, char *header)
{
  const int buf_size = 1024;
  const char *mode = *first_time ? "w" : "a";
  CCTK_STRING *out_dir  = (CCTK_STRING *) CCTK_ParameterGet("out_dir", "IOUtil", NULL);
  char output_name[buf_size];

  if (*out_dir == 0)
  {
    CCTK_WARN(1,"Parameter IOUtil::out_dir not found");
    return 0;
  }

  if (strlen(*out_dir)+1 > buf_size)
  {
    CCTK_WARN(1,"Parameter IOUtil::out_dir is too long");
    return 0;
  }

  sprintf(output_name, "%s/%s", *out_dir, leaf);
  FILE *fp = fopen(output_name, mode);
  
  if (fp == 0)
  {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
		"Could not open output file %s\n", output_name);
    return 0;
  }

  if (*first_time)
    fputs(header, fp);
  
  *first_time = 0;

  return fp;
}

/* Output a file listing the memory statistics for each processor  */
void RunStats_OutputAllMemory(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if (cctk_iteration % output_mem_every != 0)
    return;

  // Open output file
  static int first_time = 1;
  const int leaf_size = 100;
  char leaf[leaf_size];
  sprintf(leaf, "MemStats%.3d.asc", CCTK_MyProc(cctkGH));
  FILE *f = open_output_file(leaf, &first_time, 
          "# 1:cctk_time 2:maxrss 3:majflt 4:arena 5:ordblks 6:hblks 7:hblkhd 8:uordblks 9:fordblks 10:keepcost\n");

  if (f)
  {
    // if Carpet supported INT8 out of the box, then these would be %zu...
    fprintf(f, "%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", 
            cctk_time, *maxrss, *majflt, *arena, *ordblks, *hblks, *hblkhd, 
            *uordblks, *fordblks, *keepcost);
    fclose(f);
  }
}

/* Output a file listing the memory statistics for the system  */
void RunStats_OutputSystemMemory(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if (cctk_iteration % output_mem_every != 0)
    return;

  // Open output file
  static int first_time = 1;
  const int leaf_size = 100;
  char leaf[leaf_size];
  sprintf(leaf, "MemStatsSystem%.3d.asc", CCTK_MyProc(cctkGH));
  FILE *f = open_output_file(leaf, &first_time, 
          "# 1:cctk_time 2:active 3:swapfree 4:total 5:freemem 6:inactive 7:swaptotal\n");

  // Work out memory usage

  if (f)
  {
    size_t total = 0, freemem = 0, active = 0, inactive = 0, swaptotal = 0, swapfree = 0;
    get_node_stats(&total, &freemem, &active, &inactive, &swaptotal, 
                   &swapfree);
    fprintf(f, "%f\t%zu\t%zu\t%zu\t%zu\t%zu\t%zu\n", 
            cctk_time, active, swapfree, total, freemem, inactive, swaptotal);
    fclose(f);
  }
}


void RunStats_GetMemoryStatistics(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  int mb = 1024*1024;
  int kb = 1024;
  size_t maxrssL, majfltL;
  static int maxrss_warned = 0, majflt_warned = 0;

  maxrssL = get_rss();
  majfltL = get_majflt();
  if(maxrssL > INT_MAX && !maxrss_warned)
  {
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
		"maxxrss (%zu) exceeds INT_MAX (%d), limiting.\n", maxrssL, INT_MAX);
    maxrss_warned = 1;
  }
  if(majfltL > INT_MAX && !majflt_warned)
  {
    CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
		"majflt (%zu) exceeds INT_MAX (%d), limiting.\n", majfltL, INT_MAX);
    majflt_warned = 1;
  }

  // all but maxrss and majflt should fit in an int according to 'info mallinfo'
  // arena is fishy here since it says "total size of allocated memory in bytes"
  *maxrss = int_limit(maxrssL);
  *majflt = int_limit(majfltL);
  *arena = mallinfo().arena;
  *ordblks = mallinfo().ordblks;
  *hblks = mallinfo().hblks;
  *hblkhd = mallinfo().hblkhd;
  *uordblks = mallinfo().uordblks;
  *fordblks = mallinfo().fordblks;
  *keepcost = mallinfo().keepcost;

  *maxrss_mb =  int_limit(maxrssL / mb);
  *majflt_mb = int_limit(majfltL / mb);
  *arena_mb = *arena / mb;
  // the following ones do not seem to make that much sense to me since they
  // are counters of chunks, not byte sizes
  *ordblks_mb = *ordblks / mb;
  *hblks_mb = *hblks / mb;
  *hblkhd_mb = *hblkhd / mb;
  *uordblks_mb = *uordblks / mb;
  *fordblks_mb = *fordblks / mb;
  *keepcost_mb = *keepcost / mb;

  *maxrss_kb = int_limit(maxrssL / kb);
  *majflt_kb = int_limit(majfltL / kb);
  *arena_kb = *arena / kb;
  // the following ones do not seem to make that much sense to me since they
  // are counters of chunks, not byte sizes
  *ordblks_kb = *ordblks / kb;
  *hblks_kb = *hblks / kb;
  *hblkhd_kb = *hblkhd / kb;
  *uordblks_kb = *uordblks / kb;
  *fordblks_kb = *fordblks / kb;
  *keepcost_kb = *keepcost / kb;
}

