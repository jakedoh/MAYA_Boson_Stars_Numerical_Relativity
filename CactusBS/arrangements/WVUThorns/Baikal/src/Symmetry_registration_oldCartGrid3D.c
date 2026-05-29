
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

void Baikal_Symmetry_registration(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Stores gridfunction parity across x=0, y=0, and z=0 planes, respectively
  int sym[3];

  // Next register parities for each gridfunction based on its name
  //    (to ensure this algorithm is robust, gridfunctions with integers
  //     in their base names are forbidden in NRPy+).

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      sym[0] *= -1;
      sym[0] *= -1;
      SetCartSymVN(cctkGH, sym, "Baikal::aDD00GF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      sym[0] *= -1;
      sym[1] *= -1;
      SetCartSymVN(cctkGH, sym, "Baikal::aDD01GF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      sym[0] *= -1;
      sym[2] *= -1;
      SetCartSymVN(cctkGH, sym, "Baikal::aDD02GF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      sym[1] *= -1;
      sym[1] *= -1;
      SetCartSymVN(cctkGH, sym, "Baikal::aDD11GF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      sym[1] *= -1;
      sym[2] *= -1;
      SetCartSymVN(cctkGH, sym, "Baikal::aDD12GF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      sym[2] *= -1;
      sym[2] *= -1;
      SetCartSymVN(cctkGH, sym, "Baikal::aDD22GF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      // (this gridfunction is a scalar -- no need to change default sym[]'s!)
      SetCartSymVN(cctkGH, sym, "Baikal::alphaGF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      sym[0] = -1;
      SetCartSymVN(cctkGH, sym, "Baikal::betU0GF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      sym[1] = -1;
      SetCartSymVN(cctkGH, sym, "Baikal::betU1GF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      sym[2] = -1;
      SetCartSymVN(cctkGH, sym, "Baikal::betU2GF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      // (this gridfunction is a scalar -- no need to change default sym[]'s!)
      SetCartSymVN(cctkGH, sym, "Baikal::cfGF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      sym[0] *= -1;
      sym[0] *= -1;
      SetCartSymVN(cctkGH, sym, "Baikal::hDD00GF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      sym[0] *= -1;
      sym[1] *= -1;
      SetCartSymVN(cctkGH, sym, "Baikal::hDD01GF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      sym[0] *= -1;
      sym[2] *= -1;
      SetCartSymVN(cctkGH, sym, "Baikal::hDD02GF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      sym[1] *= -1;
      sym[1] *= -1;
      SetCartSymVN(cctkGH, sym, "Baikal::hDD11GF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      sym[1] *= -1;
      sym[2] *= -1;
      SetCartSymVN(cctkGH, sym, "Baikal::hDD12GF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      sym[2] *= -1;
      sym[2] *= -1;
      SetCartSymVN(cctkGH, sym, "Baikal::hDD22GF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      sym[0] = -1;
      SetCartSymVN(cctkGH, sym, "Baikal::lambdaU0GF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      sym[1] = -1;
      SetCartSymVN(cctkGH, sym, "Baikal::lambdaU1GF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      sym[2] = -1;
      SetCartSymVN(cctkGH, sym, "Baikal::lambdaU2GF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      // (this gridfunction is a scalar -- no need to change default sym[]'s!)
      SetCartSymVN(cctkGH, sym, "Baikal::trKGF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      sym[0] = -1;
      SetCartSymVN(cctkGH, sym, "Baikal::vetU0GF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      sym[1] = -1;
      SetCartSymVN(cctkGH, sym, "Baikal::vetU1GF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      sym[2] = -1;
      SetCartSymVN(cctkGH, sym, "Baikal::vetU2GF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      sym[0] *= -1;
      sym[0] *= -1;
      SetCartSymVN(cctkGH, sym, "Baikal::RbarDD00GF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      sym[0] *= -1;
      sym[1] *= -1;
      SetCartSymVN(cctkGH, sym, "Baikal::RbarDD01GF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      sym[0] *= -1;
      sym[2] *= -1;
      SetCartSymVN(cctkGH, sym, "Baikal::RbarDD02GF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      sym[1] *= -1;
      sym[1] *= -1;
      SetCartSymVN(cctkGH, sym, "Baikal::RbarDD11GF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      sym[1] *= -1;
      sym[2] *= -1;
      SetCartSymVN(cctkGH, sym, "Baikal::RbarDD12GF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      sym[2] *= -1;
      sym[2] *= -1;
      SetCartSymVN(cctkGH, sym, "Baikal::RbarDD22GF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      SetCartSymVN(cctkGH, sym, "Baikal::T4UU00GF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      sym[0] *= -1;
      SetCartSymVN(cctkGH, sym, "Baikal::T4UU01GF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      sym[1] *= -1;
      SetCartSymVN(cctkGH, sym, "Baikal::T4UU02GF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      sym[2] *= -1;
      SetCartSymVN(cctkGH, sym, "Baikal::T4UU03GF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      sym[0] *= -1;
      sym[0] *= -1;
      SetCartSymVN(cctkGH, sym, "Baikal::T4UU11GF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      sym[0] *= -1;
      sym[1] *= -1;
      SetCartSymVN(cctkGH, sym, "Baikal::T4UU12GF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      sym[0] *= -1;
      sym[2] *= -1;
      SetCartSymVN(cctkGH, sym, "Baikal::T4UU13GF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      sym[1] *= -1;
      sym[1] *= -1;
      SetCartSymVN(cctkGH, sym, "Baikal::T4UU22GF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      sym[1] *= -1;
      sym[2] *= -1;
      SetCartSymVN(cctkGH, sym, "Baikal::T4UU23GF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      sym[2] *= -1;
      sym[2] *= -1;
      SetCartSymVN(cctkGH, sym, "Baikal::T4UU33GF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      // (this gridfunction is a scalar -- no need to change default sym[]'s!)
      SetCartSymVN(cctkGH, sym, "Baikal::HGF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      sym[0] = -1;
      SetCartSymVN(cctkGH, sym, "Baikal::MU0GF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      sym[1] = -1;
      SetCartSymVN(cctkGH, sym, "Baikal::MU1GF");

      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
      sym[2] = -1;
      SetCartSymVN(cctkGH, sym, "Baikal::MU2GF");
}
