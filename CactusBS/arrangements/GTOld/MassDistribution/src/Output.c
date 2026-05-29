#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "MassDistribution.h"

void MassDistribution_Output(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  const CCTK_REAL NANVAL = atof("nan");
  const char *fmode;  /* file mode: append if already written */
  char *filename;
  char basename[1024];
  char format_str_real[256];
  int vol, bin;
  size_t len_written, sz_filename;
  CCTK_REAL bin_width, bin_value, bin_label;
  CCTK_REAL total_mass, mass_below_threshold, mass_above_threshold;
  enum bin_scaling_type scaling;
  FILE *file;

  const CCTK_INT myproc= CCTK_MyProc(cctkGH);

  /* output only on one processor */
  if(myproc != 0) {
    return;
  }

  if (verbose >= 2) {
    CCTK_VInfo(CCTK_THORNSTRING, "writing output");
  }

  for (vol = 0 ; vol < intvolumes ; vol++)
  {
    /* check if we have to write output for this volume */
    /* NOTE: we also produce output (NAN) for invalid histograms to keep the
     * counting of datasets in sync */
    if(!HistogramActive[vol]) {
      continue;
    }

    /* filename */
    sprintf(basename, outfile_basename, vol);

    sz_filename = strlen (out_dir) + strlen (basename) + 2;
    filename = (char *) malloc (sz_filename);
    assert(filename);
    len_written = snprintf (filename, sz_filename, "%s/%s", out_dir, basename);
    assert(len_written > 0);
    assert(len_written < sz_filename);

    /* open file */
    fmode = file_created[vol] ? "a" : "w"; 
    file = fopen (filename, fmode);
    if (!file) {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Could not open 1d output file '%s'",
                  filename);
    }

    /* write header on startup */
    if (!file_created[vol]) {
      fprintf(file,"# MassDistribution mass totals\n");
      fprintf(file,"# gnuplot column index:\n");
      if (CCTK_IsFunctionAliased ("UniqueBuildID")) {
        fprintf(file, "# Build ID: %s\n", (const char *)UniqueBuildID(cctkGH));
      }
      if (CCTK_IsFunctionAliased ("UniqueSimulationID")) {
        fprintf(file, "# Simulation ID: %s\n", (const char *)UniqueSimulationID(cctkGH));
      }
      if (CCTK_IsFunctionAliased ("UniqueRunID")) {
        fprintf(file, "# Run ID: %s\n", (const char *)UniqueRunID(cctkGH));
      }
      fprintf(file,"# 1:it 2:t 3:bin 4:%s 5:Delta_%s/Delta_%s\n", bin_variable[vol], mass_variable[vol], bin_variable[vol]);
      file_created[vol]=1;
    }

    /* format string */
    len_written = snprintf (format_str_real, sizeof(format_str_real),
            "%%d\t%s\t%%d\t%s\t%s\n",out_format,out_format,out_format);
    assert(len_written > 0);
    assert(len_written < sizeof(format_str_real));

    /* write header */
    scaling = MassDistribution_resolve_keyword_parameter("bin_scaling", bin_scaling[vol], NUM_BIN_SCALINGS, bin_scalings);
    switch (scaling) /* size of bin for header */
    {
      case bin_scaling_linear:
        bin_width = (HistogramMaximum[vol] - HistogramMinimum[vol])/number_of_bins;
        break;
      case bin_scaling_logarithmic:
        bin_width = log10(HistogramMaximum[vol]/HistogramMinimum[vol])/number_of_bins;
        break;
      default:
        assert(0 == "Internal error");
        break;
    }
    total_mass = MassPerSpecificEnergy[vol*TOTAL_NUMBER_OF_BINS + TOTAL_MASS_BIN];
    mass_below_threshold = MassPerSpecificEnergy[vol*TOTAL_NUMBER_OF_BINS + BELOW_MASS_BIN];
    mass_above_threshold = MassPerSpecificEnergy[vol*TOTAL_NUMBER_OF_BINS + ABOVE_MASS_BIN];
    fprintf(file, "# Time = %f min = %g max = %g width = %g total = %g below_threshold = %g above_threshold = %g\n", 
            cctk_time, HistogramMinimum[vol], HistogramMaximum[vol],
            bin_width, total_mass, mass_below_threshold, mass_above_threshold);

    /* write data */
    for (bin = 0 ; bin < number_of_bins ; bin++)
    {
      if (divide_by_bin_width[vol]) {
        bin_width = HistogramLabels[vol*NUMBER_OF_LABELS + bin + 1] - HistogramLabels[vol*NUMBER_OF_LABELS + bin];
        bin_value = MassPerSpecificEnergy[vol*TOTAL_NUMBER_OF_BINS + bin] / bin_width;
      } else {
        bin_value = MassPerSpecificEnergy[vol*TOTAL_NUMBER_OF_BINS + bin];
      }
      bin_label = HistogramLabels[vol*NUMBER_OF_LABELS + bin];

      fprintf(file, format_str_real, cctk_iteration, cctk_time, bin,
              bin_label, HistogramValid[vol] ?  bin_value : NANVAL);
    }
    fputs("\n\n", file);

    fclose(file);
    free(filename);
  }
}
