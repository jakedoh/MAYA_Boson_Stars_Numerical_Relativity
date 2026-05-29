#include <stdlib.h>

#include "cctk.h"
#include "cctk_FortranString.h"
#include "CactusBase/IOUtil/src/ioutil_AdvertisedFiles.h"


static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusEinstein_Extract_Advertise_c)

void CCTK_FCALL CCTK_FNAME(Extract_Advertise)
     (cGH **GH, ONE_FORTSTRING_ARG);

void CCTK_FCALL CCTK_FNAME(Extract_Advertise)
     (cGH **GH, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE(file)
  /* advertise new files for downloading */
  ioAdvertisedFileDesc advertisedFile;

  advertisedFile.slice = "";
  advertisedFile.thorn = CCTK_THORNSTRING;
  advertisedFile.varname = "Waveform";
  advertisedFile.description = "Extract files";
  advertisedFile.mimetype = "application/x-graph";

  IOUtil_AdvertiseFile (*GH, file, &advertisedFile);
  free(file);
}

