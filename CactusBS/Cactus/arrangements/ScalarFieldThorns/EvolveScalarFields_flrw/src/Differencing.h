#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd1(u) ((-KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0))*p1o2dx)
#else
#  define PDstandard2nd1(u) (PDstandard2nd1_impl(u,p1o2dx,cdj,cdk))
static CCTK_REAL PDstandard2nd1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard2nd1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0))*p1o2dx;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd2(u) ((-KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0))*p1o2dy)
#else
#  define PDstandard2nd2(u) (PDstandard2nd2_impl(u,p1o2dy,cdj,cdk))
static CCTK_REAL PDstandard2nd2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard2nd2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0))*p1o2dy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd3(u) ((-KRANC_GFOFFSET3D(u,0,0,-1) + KRANC_GFOFFSET3D(u,0,0,1))*p1o2dz)
#else
#  define PDstandard2nd3(u) (PDstandard2nd3_impl(u,p1o2dz,cdj,cdk))
static CCTK_REAL PDstandard2nd3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard2nd3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard2nd2_impl(u, p1o2dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd11(u) ((-2*KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0))*p1odx2)
#else
#  define PDstandard2nd11(u) (PDstandard2nd11_impl(u,p1odx2,cdj,cdk))
static CCTK_REAL PDstandard2nd11_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odx2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard2nd11_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odx2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-2*KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0))*p1odx2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd22(u) ((-2*KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0))*p1ody2)
#else
#  define PDstandard2nd22(u) (PDstandard2nd22_impl(u,p1ody2,cdj,cdk))
static CCTK_REAL PDstandard2nd22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1ody2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard2nd22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1ody2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-2*KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0))*p1ody2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd33(u) ((-2*KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,0,0,-1) + KRANC_GFOFFSET3D(u,0,0,1))*p1odz2)
#else
#  define PDstandard2nd33(u) (PDstandard2nd33_impl(u,p1odz2,cdj,cdk))
static CCTK_REAL PDstandard2nd33_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odz2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard2nd33_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odz2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard2nd22_impl(u, p1odz2, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd12(u) ((KRANC_GFOFFSET3D(u,-1,-1,0) - KRANC_GFOFFSET3D(u,-1,1,0) - KRANC_GFOFFSET3D(u,1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0))*p1o4dxdy)
#else
#  define PDstandard2nd12(u) (PDstandard2nd12_impl(u,p1o4dxdy,cdj,cdk))
static CCTK_REAL PDstandard2nd12_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard2nd12_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (KRANC_GFOFFSET3D(u,-1,-1,0) - KRANC_GFOFFSET3D(u,-1,1,0) - KRANC_GFOFFSET3D(u,1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0))*p1o4dxdy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd13(u) ((KRANC_GFOFFSET3D(u,-1,0,-1) - KRANC_GFOFFSET3D(u,-1,0,1) - KRANC_GFOFFSET3D(u,1,0,-1) + KRANC_GFOFFSET3D(u,1,0,1))*p1o4dxdz)
#else
#  define PDstandard2nd13(u) (PDstandard2nd13_impl(u,p1o4dxdz,cdj,cdk))
static CCTK_REAL PDstandard2nd13_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard2nd13_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard2nd12_impl(u, p1o4dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd21(u) ((KRANC_GFOFFSET3D(u,-1,-1,0) - KRANC_GFOFFSET3D(u,-1,1,0) - KRANC_GFOFFSET3D(u,1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0))*p1o4dxdy)
#else
#  define PDstandard2nd21(u) (PDstandard2nd21_impl(u,p1o4dxdy,cdj,cdk))
static CCTK_REAL PDstandard2nd21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard2nd21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard2nd12_impl(u, p1o4dxdy, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd23(u) ((KRANC_GFOFFSET3D(u,0,-1,-1) - KRANC_GFOFFSET3D(u,0,-1,1) - KRANC_GFOFFSET3D(u,0,1,-1) + KRANC_GFOFFSET3D(u,0,1,1))*p1o4dydz)
#else
#  define PDstandard2nd23(u) (PDstandard2nd23_impl(u,p1o4dydz,cdj,cdk))
static CCTK_REAL PDstandard2nd23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard2nd23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (KRANC_GFOFFSET3D(u,0,-1,-1) - KRANC_GFOFFSET3D(u,0,-1,1) - KRANC_GFOFFSET3D(u,0,1,-1) + KRANC_GFOFFSET3D(u,0,1,1))*p1o4dydz;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd31(u) ((KRANC_GFOFFSET3D(u,-1,0,-1) - KRANC_GFOFFSET3D(u,-1,0,1) - KRANC_GFOFFSET3D(u,1,0,-1) + KRANC_GFOFFSET3D(u,1,0,1))*p1o4dxdz)
#else
#  define PDstandard2nd31(u) (PDstandard2nd31_impl(u,p1o4dxdz,cdj,cdk))
static CCTK_REAL PDstandard2nd31_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard2nd31_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard2nd12_impl(u, p1o4dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard2nd32(u) ((KRANC_GFOFFSET3D(u,0,-1,-1) - KRANC_GFOFFSET3D(u,0,-1,1) - KRANC_GFOFFSET3D(u,0,1,-1) + KRANC_GFOFFSET3D(u,0,1,1))*p1o4dydz)
#else
#  define PDstandard2nd32(u) (PDstandard2nd32_impl(u,p1o4dydz,cdj,cdk))
static CCTK_REAL PDstandard2nd32_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard2nd32_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard2nd23_impl(u, p1o4dydz, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard4th1(u) ((-8*KRANC_GFOFFSET3D(u,-1,0,0) + 8*KRANC_GFOFFSET3D(u,1,0,0) + KRANC_GFOFFSET3D(u,-2,0,0) - KRANC_GFOFFSET3D(u,2,0,0))*p1o12dx)
#else
#  define PDstandard4th1(u) (PDstandard4th1_impl(u,p1o12dx,cdj,cdk))
static CCTK_REAL PDstandard4th1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o12dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard4th1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o12dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-8*KRANC_GFOFFSET3D(u,-1,0,0) + 8*KRANC_GFOFFSET3D(u,1,0,0) + KRANC_GFOFFSET3D(u,-2,0,0) - KRANC_GFOFFSET3D(u,2,0,0))*p1o12dx;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard4th2(u) ((-8*KRANC_GFOFFSET3D(u,0,-1,0) + 8*KRANC_GFOFFSET3D(u,0,1,0) + KRANC_GFOFFSET3D(u,0,-2,0) - KRANC_GFOFFSET3D(u,0,2,0))*p1o12dy)
#else
#  define PDstandard4th2(u) (PDstandard4th2_impl(u,p1o12dy,cdj,cdk))
static CCTK_REAL PDstandard4th2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o12dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard4th2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o12dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-8*KRANC_GFOFFSET3D(u,0,-1,0) + 8*KRANC_GFOFFSET3D(u,0,1,0) + KRANC_GFOFFSET3D(u,0,-2,0) - KRANC_GFOFFSET3D(u,0,2,0))*p1o12dy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard4th3(u) ((-8*KRANC_GFOFFSET3D(u,0,0,-1) + 8*KRANC_GFOFFSET3D(u,0,0,1) + KRANC_GFOFFSET3D(u,0,0,-2) - KRANC_GFOFFSET3D(u,0,0,2))*p1o12dz)
#else
#  define PDstandard4th3(u) (PDstandard4th3_impl(u,p1o12dz,cdj,cdk))
static CCTK_REAL PDstandard4th3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o12dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard4th3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o12dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard4th2_impl(u, p1o12dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard4th11(u) ((30*KRANC_GFOFFSET3D(u,0,0,0) - 16*(KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0)) + KRANC_GFOFFSET3D(u,-2,0,0) + KRANC_GFOFFSET3D(u,2,0,0))*pm1o12dx2)
#else
#  define PDstandard4th11(u) (PDstandard4th11_impl(u,pm1o12dx2,cdj,cdk))
static CCTK_REAL PDstandard4th11_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o12dx2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard4th11_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o12dx2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (30*KRANC_GFOFFSET3D(u,0,0,0) - 16*(KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0)) + KRANC_GFOFFSET3D(u,-2,0,0) + KRANC_GFOFFSET3D(u,2,0,0))*pm1o12dx2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard4th22(u) ((30*KRANC_GFOFFSET3D(u,0,0,0) - 16*(KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0)) + KRANC_GFOFFSET3D(u,0,-2,0) + KRANC_GFOFFSET3D(u,0,2,0))*pm1o12dy2)
#else
#  define PDstandard4th22(u) (PDstandard4th22_impl(u,pm1o12dy2,cdj,cdk))
static CCTK_REAL PDstandard4th22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o12dy2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard4th22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o12dy2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (30*KRANC_GFOFFSET3D(u,0,0,0) - 16*(KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0)) + KRANC_GFOFFSET3D(u,0,-2,0) + KRANC_GFOFFSET3D(u,0,2,0))*pm1o12dy2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard4th33(u) ((30*KRANC_GFOFFSET3D(u,0,0,0) - 16*(KRANC_GFOFFSET3D(u,0,0,-1) + KRANC_GFOFFSET3D(u,0,0,1)) + KRANC_GFOFFSET3D(u,0,0,-2) + KRANC_GFOFFSET3D(u,0,0,2))*pm1o12dz2)
#else
#  define PDstandard4th33(u) (PDstandard4th33_impl(u,pm1o12dz2,cdj,cdk))
static CCTK_REAL PDstandard4th33_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o12dz2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard4th33_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o12dz2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard4th22_impl(u, pm1o12dz2, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard4th12(u) ((-64*(KRANC_GFOFFSET3D(u,-1,1,0) + KRANC_GFOFFSET3D(u,1,-1,0)) + 64*(KRANC_GFOFFSET3D(u,-1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0)) + 8*(KRANC_GFOFFSET3D(u,-1,2,0) + KRANC_GFOFFSET3D(u,1,-2,0) + KRANC_GFOFFSET3D(u,-2,1,0) + KRANC_GFOFFSET3D(u,2,-1,0)) - 8*(KRANC_GFOFFSET3D(u,-1,-2,0) + KRANC_GFOFFSET3D(u,1,2,0) + KRANC_GFOFFSET3D(u,-2,-1,0) + KRANC_GFOFFSET3D(u,2,1,0)) + KRANC_GFOFFSET3D(u,-2,-2,0) - KRANC_GFOFFSET3D(u,-2,2,0) - KRANC_GFOFFSET3D(u,2,-2,0) + KRANC_GFOFFSET3D(u,2,2,0))*p1o144dxdy)
#else
#  define PDstandard4th12(u) (PDstandard4th12_impl(u,p1o144dxdy,cdj,cdk))
static CCTK_REAL PDstandard4th12_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard4th12_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-64*(KRANC_GFOFFSET3D(u,-1,1,0) + KRANC_GFOFFSET3D(u,1,-1,0)) + 64*(KRANC_GFOFFSET3D(u,-1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0)) + 8*(KRANC_GFOFFSET3D(u,-1,2,0) + KRANC_GFOFFSET3D(u,1,-2,0) + KRANC_GFOFFSET3D(u,-2,1,0) + KRANC_GFOFFSET3D(u,2,-1,0)) - 8*(KRANC_GFOFFSET3D(u,-1,-2,0) + KRANC_GFOFFSET3D(u,1,2,0) + KRANC_GFOFFSET3D(u,-2,-1,0) + KRANC_GFOFFSET3D(u,2,1,0)) + KRANC_GFOFFSET3D(u,-2,-2,0) - KRANC_GFOFFSET3D(u,-2,2,0) - KRANC_GFOFFSET3D(u,2,-2,0) + KRANC_GFOFFSET3D(u,2,2,0))*p1o144dxdy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard4th13(u) ((-64*(KRANC_GFOFFSET3D(u,-1,0,1) + KRANC_GFOFFSET3D(u,1,0,-1)) + 64*(KRANC_GFOFFSET3D(u,-1,0,-1) + KRANC_GFOFFSET3D(u,1,0,1)) + 8*(KRANC_GFOFFSET3D(u,-1,0,2) + KRANC_GFOFFSET3D(u,1,0,-2) + KRANC_GFOFFSET3D(u,-2,0,1) + KRANC_GFOFFSET3D(u,2,0,-1)) - 8*(KRANC_GFOFFSET3D(u,-1,0,-2) + KRANC_GFOFFSET3D(u,1,0,2) + KRANC_GFOFFSET3D(u,-2,0,-1) + KRANC_GFOFFSET3D(u,2,0,1)) + KRANC_GFOFFSET3D(u,-2,0,-2) - KRANC_GFOFFSET3D(u,-2,0,2) - KRANC_GFOFFSET3D(u,2,0,-2) + KRANC_GFOFFSET3D(u,2,0,2))*p1o144dxdz)
#else
#  define PDstandard4th13(u) (PDstandard4th13_impl(u,p1o144dxdz,cdj,cdk))
static CCTK_REAL PDstandard4th13_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard4th13_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard4th12_impl(u, p1o144dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard4th21(u) ((-64*(KRANC_GFOFFSET3D(u,-1,1,0) + KRANC_GFOFFSET3D(u,1,-1,0)) + 64*(KRANC_GFOFFSET3D(u,-1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0)) + 8*(KRANC_GFOFFSET3D(u,-1,2,0) + KRANC_GFOFFSET3D(u,1,-2,0) + KRANC_GFOFFSET3D(u,-2,1,0) + KRANC_GFOFFSET3D(u,2,-1,0)) - 8*(KRANC_GFOFFSET3D(u,-1,-2,0) + KRANC_GFOFFSET3D(u,1,2,0) + KRANC_GFOFFSET3D(u,-2,-1,0) + KRANC_GFOFFSET3D(u,2,1,0)) + KRANC_GFOFFSET3D(u,-2,-2,0) - KRANC_GFOFFSET3D(u,-2,2,0) - KRANC_GFOFFSET3D(u,2,-2,0) + KRANC_GFOFFSET3D(u,2,2,0))*p1o144dxdy)
#else
#  define PDstandard4th21(u) (PDstandard4th21_impl(u,p1o144dxdy,cdj,cdk))
static CCTK_REAL PDstandard4th21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard4th21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard4th12_impl(u, p1o144dxdy, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard4th23(u) ((-64*(KRANC_GFOFFSET3D(u,0,-1,1) + KRANC_GFOFFSET3D(u,0,1,-1)) + 64*(KRANC_GFOFFSET3D(u,0,-1,-1) + KRANC_GFOFFSET3D(u,0,1,1)) + 8*(KRANC_GFOFFSET3D(u,0,-1,2) + KRANC_GFOFFSET3D(u,0,1,-2) + KRANC_GFOFFSET3D(u,0,-2,1) + KRANC_GFOFFSET3D(u,0,2,-1)) - 8*(KRANC_GFOFFSET3D(u,0,-1,-2) + KRANC_GFOFFSET3D(u,0,1,2) + KRANC_GFOFFSET3D(u,0,-2,-1) + KRANC_GFOFFSET3D(u,0,2,1)) + KRANC_GFOFFSET3D(u,0,-2,-2) - KRANC_GFOFFSET3D(u,0,-2,2) - KRANC_GFOFFSET3D(u,0,2,-2) + KRANC_GFOFFSET3D(u,0,2,2))*p1o144dydz)
#else
#  define PDstandard4th23(u) (PDstandard4th23_impl(u,p1o144dydz,cdj,cdk))
static CCTK_REAL PDstandard4th23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard4th23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-64*(KRANC_GFOFFSET3D(u,0,-1,1) + KRANC_GFOFFSET3D(u,0,1,-1)) + 64*(KRANC_GFOFFSET3D(u,0,-1,-1) + KRANC_GFOFFSET3D(u,0,1,1)) + 8*(KRANC_GFOFFSET3D(u,0,-1,2) + KRANC_GFOFFSET3D(u,0,1,-2) + KRANC_GFOFFSET3D(u,0,-2,1) + KRANC_GFOFFSET3D(u,0,2,-1)) - 8*(KRANC_GFOFFSET3D(u,0,-1,-2) + KRANC_GFOFFSET3D(u,0,1,2) + KRANC_GFOFFSET3D(u,0,-2,-1) + KRANC_GFOFFSET3D(u,0,2,1)) + KRANC_GFOFFSET3D(u,0,-2,-2) - KRANC_GFOFFSET3D(u,0,-2,2) - KRANC_GFOFFSET3D(u,0,2,-2) + KRANC_GFOFFSET3D(u,0,2,2))*p1o144dydz;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard4th31(u) ((-64*(KRANC_GFOFFSET3D(u,-1,0,1) + KRANC_GFOFFSET3D(u,1,0,-1)) + 64*(KRANC_GFOFFSET3D(u,-1,0,-1) + KRANC_GFOFFSET3D(u,1,0,1)) + 8*(KRANC_GFOFFSET3D(u,-1,0,2) + KRANC_GFOFFSET3D(u,1,0,-2) + KRANC_GFOFFSET3D(u,-2,0,1) + KRANC_GFOFFSET3D(u,2,0,-1)) - 8*(KRANC_GFOFFSET3D(u,-1,0,-2) + KRANC_GFOFFSET3D(u,1,0,2) + KRANC_GFOFFSET3D(u,-2,0,-1) + KRANC_GFOFFSET3D(u,2,0,1)) + KRANC_GFOFFSET3D(u,-2,0,-2) - KRANC_GFOFFSET3D(u,-2,0,2) - KRANC_GFOFFSET3D(u,2,0,-2) + KRANC_GFOFFSET3D(u,2,0,2))*p1o144dxdz)
#else
#  define PDstandard4th31(u) (PDstandard4th31_impl(u,p1o144dxdz,cdj,cdk))
static CCTK_REAL PDstandard4th31_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard4th31_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard4th12_impl(u, p1o144dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard4th32(u) ((-64*(KRANC_GFOFFSET3D(u,0,-1,1) + KRANC_GFOFFSET3D(u,0,1,-1)) + 64*(KRANC_GFOFFSET3D(u,0,-1,-1) + KRANC_GFOFFSET3D(u,0,1,1)) + 8*(KRANC_GFOFFSET3D(u,0,-1,2) + KRANC_GFOFFSET3D(u,0,1,-2) + KRANC_GFOFFSET3D(u,0,-2,1) + KRANC_GFOFFSET3D(u,0,2,-1)) - 8*(KRANC_GFOFFSET3D(u,0,-1,-2) + KRANC_GFOFFSET3D(u,0,1,2) + KRANC_GFOFFSET3D(u,0,-2,-1) + KRANC_GFOFFSET3D(u,0,2,1)) + KRANC_GFOFFSET3D(u,0,-2,-2) - KRANC_GFOFFSET3D(u,0,-2,2) - KRANC_GFOFFSET3D(u,0,2,-2) + KRANC_GFOFFSET3D(u,0,2,2))*p1o144dydz)
#else
#  define PDstandard4th32(u) (PDstandard4th32_impl(u,p1o144dydz,cdj,cdk))
static CCTK_REAL PDstandard4th32_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard4th32_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard4th23_impl(u, p1o144dydz, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard6th1(u) ((-45*KRANC_GFOFFSET3D(u,-1,0,0) + 45*KRANC_GFOFFSET3D(u,1,0,0) + 9*KRANC_GFOFFSET3D(u,-2,0,0) - 9*KRANC_GFOFFSET3D(u,2,0,0) - KRANC_GFOFFSET3D(u,-3,0,0) + KRANC_GFOFFSET3D(u,3,0,0))*p1o60dx)
#else
#  define PDstandard6th1(u) (PDstandard6th1_impl(u,p1o60dx,cdj,cdk))
static CCTK_REAL PDstandard6th1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o60dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard6th1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o60dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-45*KRANC_GFOFFSET3D(u,-1,0,0) + 45*KRANC_GFOFFSET3D(u,1,0,0) + 9*KRANC_GFOFFSET3D(u,-2,0,0) - 9*KRANC_GFOFFSET3D(u,2,0,0) - KRANC_GFOFFSET3D(u,-3,0,0) + KRANC_GFOFFSET3D(u,3,0,0))*p1o60dx;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard6th2(u) ((-45*KRANC_GFOFFSET3D(u,0,-1,0) + 45*KRANC_GFOFFSET3D(u,0,1,0) + 9*KRANC_GFOFFSET3D(u,0,-2,0) - 9*KRANC_GFOFFSET3D(u,0,2,0) - KRANC_GFOFFSET3D(u,0,-3,0) + KRANC_GFOFFSET3D(u,0,3,0))*p1o60dy)
#else
#  define PDstandard6th2(u) (PDstandard6th2_impl(u,p1o60dy,cdj,cdk))
static CCTK_REAL PDstandard6th2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o60dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard6th2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o60dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-45*KRANC_GFOFFSET3D(u,0,-1,0) + 45*KRANC_GFOFFSET3D(u,0,1,0) + 9*KRANC_GFOFFSET3D(u,0,-2,0) - 9*KRANC_GFOFFSET3D(u,0,2,0) - KRANC_GFOFFSET3D(u,0,-3,0) + KRANC_GFOFFSET3D(u,0,3,0))*p1o60dy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard6th3(u) ((-45*KRANC_GFOFFSET3D(u,0,0,-1) + 45*KRANC_GFOFFSET3D(u,0,0,1) + 9*KRANC_GFOFFSET3D(u,0,0,-2) - 9*KRANC_GFOFFSET3D(u,0,0,2) - KRANC_GFOFFSET3D(u,0,0,-3) + KRANC_GFOFFSET3D(u,0,0,3))*p1o60dz)
#else
#  define PDstandard6th3(u) (PDstandard6th3_impl(u,p1o60dz,cdj,cdk))
static CCTK_REAL PDstandard6th3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o60dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard6th3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o60dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard6th2_impl(u, p1o60dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard6th11(u) ((-490*KRANC_GFOFFSET3D(u,0,0,0) + 270*(KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0)) - 27*(KRANC_GFOFFSET3D(u,-2,0,0) + KRANC_GFOFFSET3D(u,2,0,0)) + 2*(KRANC_GFOFFSET3D(u,-3,0,0) + KRANC_GFOFFSET3D(u,3,0,0)))*p1o180dx2)
#else
#  define PDstandard6th11(u) (PDstandard6th11_impl(u,p1o180dx2,cdj,cdk))
static CCTK_REAL PDstandard6th11_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o180dx2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard6th11_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o180dx2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-490*KRANC_GFOFFSET3D(u,0,0,0) + 270*(KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0)) - 27*(KRANC_GFOFFSET3D(u,-2,0,0) + KRANC_GFOFFSET3D(u,2,0,0)) + 2*(KRANC_GFOFFSET3D(u,-3,0,0) + KRANC_GFOFFSET3D(u,3,0,0)))*p1o180dx2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard6th22(u) ((-490*KRANC_GFOFFSET3D(u,0,0,0) + 270*(KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0)) - 27*(KRANC_GFOFFSET3D(u,0,-2,0) + KRANC_GFOFFSET3D(u,0,2,0)) + 2*(KRANC_GFOFFSET3D(u,0,-3,0) + KRANC_GFOFFSET3D(u,0,3,0)))*p1o180dy2)
#else
#  define PDstandard6th22(u) (PDstandard6th22_impl(u,p1o180dy2,cdj,cdk))
static CCTK_REAL PDstandard6th22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o180dy2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard6th22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o180dy2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-490*KRANC_GFOFFSET3D(u,0,0,0) + 270*(KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0)) - 27*(KRANC_GFOFFSET3D(u,0,-2,0) + KRANC_GFOFFSET3D(u,0,2,0)) + 2*(KRANC_GFOFFSET3D(u,0,-3,0) + KRANC_GFOFFSET3D(u,0,3,0)))*p1o180dy2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard6th33(u) ((-490*KRANC_GFOFFSET3D(u,0,0,0) + 270*(KRANC_GFOFFSET3D(u,0,0,-1) + KRANC_GFOFFSET3D(u,0,0,1)) - 27*(KRANC_GFOFFSET3D(u,0,0,-2) + KRANC_GFOFFSET3D(u,0,0,2)) + 2*(KRANC_GFOFFSET3D(u,0,0,-3) + KRANC_GFOFFSET3D(u,0,0,3)))*p1o180dz2)
#else
#  define PDstandard6th33(u) (PDstandard6th33_impl(u,p1o180dz2,cdj,cdk))
static CCTK_REAL PDstandard6th33_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o180dz2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard6th33_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o180dz2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard6th22_impl(u, p1o180dz2, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard6th12(u) ((-2025*(KRANC_GFOFFSET3D(u,-1,1,0) + KRANC_GFOFFSET3D(u,1,-1,0)) + 2025*(KRANC_GFOFFSET3D(u,-1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0)) + 405*(KRANC_GFOFFSET3D(u,-1,2,0) + KRANC_GFOFFSET3D(u,1,-2,0) + KRANC_GFOFFSET3D(u,-2,1,0) + KRANC_GFOFFSET3D(u,2,-1,0)) - 405*(KRANC_GFOFFSET3D(u,-1,-2,0) + KRANC_GFOFFSET3D(u,1,2,0) + KRANC_GFOFFSET3D(u,-2,-1,0) + KRANC_GFOFFSET3D(u,2,1,0)) - 81*(KRANC_GFOFFSET3D(u,-2,2,0) + KRANC_GFOFFSET3D(u,2,-2,0)) + 81*(KRANC_GFOFFSET3D(u,-2,-2,0) + KRANC_GFOFFSET3D(u,2,2,0)) - 45*(KRANC_GFOFFSET3D(u,-1,3,0) + KRANC_GFOFFSET3D(u,1,-3,0) + KRANC_GFOFFSET3D(u,-3,1,0) + KRANC_GFOFFSET3D(u,3,-1,0)) + 45*(KRANC_GFOFFSET3D(u,-1,-3,0) + KRANC_GFOFFSET3D(u,1,3,0) + KRANC_GFOFFSET3D(u,-3,-1,0) + KRANC_GFOFFSET3D(u,3,1,0)) + 9*(KRANC_GFOFFSET3D(u,-2,3,0) + KRANC_GFOFFSET3D(u,2,-3,0) + KRANC_GFOFFSET3D(u,-3,2,0) + KRANC_GFOFFSET3D(u,3,-2,0)) - 9*(KRANC_GFOFFSET3D(u,-2,-3,0) + KRANC_GFOFFSET3D(u,2,3,0) + KRANC_GFOFFSET3D(u,-3,-2,0) + KRANC_GFOFFSET3D(u,3,2,0)) + KRANC_GFOFFSET3D(u,-3,-3,0) - KRANC_GFOFFSET3D(u,-3,3,0) - KRANC_GFOFFSET3D(u,3,-3,0) + KRANC_GFOFFSET3D(u,3,3,0))*p1o3600dxdy)
#else
#  define PDstandard6th12(u) (PDstandard6th12_impl(u,p1o3600dxdy,cdj,cdk))
static CCTK_REAL PDstandard6th12_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o3600dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard6th12_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o3600dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-2025*(KRANC_GFOFFSET3D(u,-1,1,0) + KRANC_GFOFFSET3D(u,1,-1,0)) + 2025*(KRANC_GFOFFSET3D(u,-1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0)) + 405*(KRANC_GFOFFSET3D(u,-1,2,0) + KRANC_GFOFFSET3D(u,1,-2,0) + KRANC_GFOFFSET3D(u,-2,1,0) + KRANC_GFOFFSET3D(u,2,-1,0)) - 405*(KRANC_GFOFFSET3D(u,-1,-2,0) + KRANC_GFOFFSET3D(u,1,2,0) + KRANC_GFOFFSET3D(u,-2,-1,0) + KRANC_GFOFFSET3D(u,2,1,0)) - 81*(KRANC_GFOFFSET3D(u,-2,2,0) + KRANC_GFOFFSET3D(u,2,-2,0)) + 81*(KRANC_GFOFFSET3D(u,-2,-2,0) + KRANC_GFOFFSET3D(u,2,2,0)) - 45*(KRANC_GFOFFSET3D(u,-1,3,0) + KRANC_GFOFFSET3D(u,1,-3,0) + KRANC_GFOFFSET3D(u,-3,1,0) + KRANC_GFOFFSET3D(u,3,-1,0)) + 45*(KRANC_GFOFFSET3D(u,-1,-3,0) + KRANC_GFOFFSET3D(u,1,3,0) + KRANC_GFOFFSET3D(u,-3,-1,0) + KRANC_GFOFFSET3D(u,3,1,0)) + 9*(KRANC_GFOFFSET3D(u,-2,3,0) + KRANC_GFOFFSET3D(u,2,-3,0) + KRANC_GFOFFSET3D(u,-3,2,0) + KRANC_GFOFFSET3D(u,3,-2,0)) - 9*(KRANC_GFOFFSET3D(u,-2,-3,0) + KRANC_GFOFFSET3D(u,2,3,0) + KRANC_GFOFFSET3D(u,-3,-2,0) + KRANC_GFOFFSET3D(u,3,2,0)) + KRANC_GFOFFSET3D(u,-3,-3,0) - KRANC_GFOFFSET3D(u,-3,3,0) - KRANC_GFOFFSET3D(u,3,-3,0) + KRANC_GFOFFSET3D(u,3,3,0))*p1o3600dxdy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard6th13(u) ((-2025*(KRANC_GFOFFSET3D(u,-1,0,1) + KRANC_GFOFFSET3D(u,1,0,-1)) + 2025*(KRANC_GFOFFSET3D(u,-1,0,-1) + KRANC_GFOFFSET3D(u,1,0,1)) + 405*(KRANC_GFOFFSET3D(u,-1,0,2) + KRANC_GFOFFSET3D(u,1,0,-2) + KRANC_GFOFFSET3D(u,-2,0,1) + KRANC_GFOFFSET3D(u,2,0,-1)) - 405*(KRANC_GFOFFSET3D(u,-1,0,-2) + KRANC_GFOFFSET3D(u,1,0,2) + KRANC_GFOFFSET3D(u,-2,0,-1) + KRANC_GFOFFSET3D(u,2,0,1)) - 81*(KRANC_GFOFFSET3D(u,-2,0,2) + KRANC_GFOFFSET3D(u,2,0,-2)) + 81*(KRANC_GFOFFSET3D(u,-2,0,-2) + KRANC_GFOFFSET3D(u,2,0,2)) - 45*(KRANC_GFOFFSET3D(u,-1,0,3) + KRANC_GFOFFSET3D(u,1,0,-3) + KRANC_GFOFFSET3D(u,-3,0,1) + KRANC_GFOFFSET3D(u,3,0,-1)) + 45*(KRANC_GFOFFSET3D(u,-1,0,-3) + KRANC_GFOFFSET3D(u,1,0,3) + KRANC_GFOFFSET3D(u,-3,0,-1) + KRANC_GFOFFSET3D(u,3,0,1)) + 9*(KRANC_GFOFFSET3D(u,-2,0,3) + KRANC_GFOFFSET3D(u,2,0,-3) + KRANC_GFOFFSET3D(u,-3,0,2) + KRANC_GFOFFSET3D(u,3,0,-2)) - 9*(KRANC_GFOFFSET3D(u,-2,0,-3) + KRANC_GFOFFSET3D(u,2,0,3) + KRANC_GFOFFSET3D(u,-3,0,-2) + KRANC_GFOFFSET3D(u,3,0,2)) + KRANC_GFOFFSET3D(u,-3,0,-3) - KRANC_GFOFFSET3D(u,-3,0,3) - KRANC_GFOFFSET3D(u,3,0,-3) + KRANC_GFOFFSET3D(u,3,0,3))*p1o3600dxdz)
#else
#  define PDstandard6th13(u) (PDstandard6th13_impl(u,p1o3600dxdz,cdj,cdk))
static CCTK_REAL PDstandard6th13_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o3600dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard6th13_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o3600dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard6th12_impl(u, p1o3600dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard6th21(u) ((-2025*(KRANC_GFOFFSET3D(u,-1,1,0) + KRANC_GFOFFSET3D(u,1,-1,0)) + 2025*(KRANC_GFOFFSET3D(u,-1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0)) + 405*(KRANC_GFOFFSET3D(u,-1,2,0) + KRANC_GFOFFSET3D(u,1,-2,0) + KRANC_GFOFFSET3D(u,-2,1,0) + KRANC_GFOFFSET3D(u,2,-1,0)) - 405*(KRANC_GFOFFSET3D(u,-1,-2,0) + KRANC_GFOFFSET3D(u,1,2,0) + KRANC_GFOFFSET3D(u,-2,-1,0) + KRANC_GFOFFSET3D(u,2,1,0)) - 81*(KRANC_GFOFFSET3D(u,-2,2,0) + KRANC_GFOFFSET3D(u,2,-2,0)) + 81*(KRANC_GFOFFSET3D(u,-2,-2,0) + KRANC_GFOFFSET3D(u,2,2,0)) - 45*(KRANC_GFOFFSET3D(u,-1,3,0) + KRANC_GFOFFSET3D(u,1,-3,0) + KRANC_GFOFFSET3D(u,-3,1,0) + KRANC_GFOFFSET3D(u,3,-1,0)) + 45*(KRANC_GFOFFSET3D(u,-1,-3,0) + KRANC_GFOFFSET3D(u,1,3,0) + KRANC_GFOFFSET3D(u,-3,-1,0) + KRANC_GFOFFSET3D(u,3,1,0)) + 9*(KRANC_GFOFFSET3D(u,-2,3,0) + KRANC_GFOFFSET3D(u,2,-3,0) + KRANC_GFOFFSET3D(u,-3,2,0) + KRANC_GFOFFSET3D(u,3,-2,0)) - 9*(KRANC_GFOFFSET3D(u,-2,-3,0) + KRANC_GFOFFSET3D(u,2,3,0) + KRANC_GFOFFSET3D(u,-3,-2,0) + KRANC_GFOFFSET3D(u,3,2,0)) + KRANC_GFOFFSET3D(u,-3,-3,0) - KRANC_GFOFFSET3D(u,-3,3,0) - KRANC_GFOFFSET3D(u,3,-3,0) + KRANC_GFOFFSET3D(u,3,3,0))*p1o3600dxdy)
#else
#  define PDstandard6th21(u) (PDstandard6th21_impl(u,p1o3600dxdy,cdj,cdk))
static CCTK_REAL PDstandard6th21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o3600dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard6th21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o3600dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard6th12_impl(u, p1o3600dxdy, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard6th23(u) ((-2025*(KRANC_GFOFFSET3D(u,0,-1,1) + KRANC_GFOFFSET3D(u,0,1,-1)) + 2025*(KRANC_GFOFFSET3D(u,0,-1,-1) + KRANC_GFOFFSET3D(u,0,1,1)) + 405*(KRANC_GFOFFSET3D(u,0,-1,2) + KRANC_GFOFFSET3D(u,0,1,-2) + KRANC_GFOFFSET3D(u,0,-2,1) + KRANC_GFOFFSET3D(u,0,2,-1)) - 405*(KRANC_GFOFFSET3D(u,0,-1,-2) + KRANC_GFOFFSET3D(u,0,1,2) + KRANC_GFOFFSET3D(u,0,-2,-1) + KRANC_GFOFFSET3D(u,0,2,1)) - 81*(KRANC_GFOFFSET3D(u,0,-2,2) + KRANC_GFOFFSET3D(u,0,2,-2)) + 81*(KRANC_GFOFFSET3D(u,0,-2,-2) + KRANC_GFOFFSET3D(u,0,2,2)) - 45*(KRANC_GFOFFSET3D(u,0,-1,3) + KRANC_GFOFFSET3D(u,0,1,-3) + KRANC_GFOFFSET3D(u,0,-3,1) + KRANC_GFOFFSET3D(u,0,3,-1)) + 45*(KRANC_GFOFFSET3D(u,0,-1,-3) + KRANC_GFOFFSET3D(u,0,1,3) + KRANC_GFOFFSET3D(u,0,-3,-1) + KRANC_GFOFFSET3D(u,0,3,1)) + 9*(KRANC_GFOFFSET3D(u,0,-2,3) + KRANC_GFOFFSET3D(u,0,2,-3) + KRANC_GFOFFSET3D(u,0,-3,2) + KRANC_GFOFFSET3D(u,0,3,-2)) - 9*(KRANC_GFOFFSET3D(u,0,-2,-3) + KRANC_GFOFFSET3D(u,0,2,3) + KRANC_GFOFFSET3D(u,0,-3,-2) + KRANC_GFOFFSET3D(u,0,3,2)) + KRANC_GFOFFSET3D(u,0,-3,-3) - KRANC_GFOFFSET3D(u,0,-3,3) - KRANC_GFOFFSET3D(u,0,3,-3) + KRANC_GFOFFSET3D(u,0,3,3))*p1o3600dydz)
#else
#  define PDstandard6th23(u) (PDstandard6th23_impl(u,p1o3600dydz,cdj,cdk))
static CCTK_REAL PDstandard6th23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o3600dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard6th23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o3600dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-2025*(KRANC_GFOFFSET3D(u,0,-1,1) + KRANC_GFOFFSET3D(u,0,1,-1)) + 2025*(KRANC_GFOFFSET3D(u,0,-1,-1) + KRANC_GFOFFSET3D(u,0,1,1)) + 405*(KRANC_GFOFFSET3D(u,0,-1,2) + KRANC_GFOFFSET3D(u,0,1,-2) + KRANC_GFOFFSET3D(u,0,-2,1) + KRANC_GFOFFSET3D(u,0,2,-1)) - 405*(KRANC_GFOFFSET3D(u,0,-1,-2) + KRANC_GFOFFSET3D(u,0,1,2) + KRANC_GFOFFSET3D(u,0,-2,-1) + KRANC_GFOFFSET3D(u,0,2,1)) - 81*(KRANC_GFOFFSET3D(u,0,-2,2) + KRANC_GFOFFSET3D(u,0,2,-2)) + 81*(KRANC_GFOFFSET3D(u,0,-2,-2) + KRANC_GFOFFSET3D(u,0,2,2)) - 45*(KRANC_GFOFFSET3D(u,0,-1,3) + KRANC_GFOFFSET3D(u,0,1,-3) + KRANC_GFOFFSET3D(u,0,-3,1) + KRANC_GFOFFSET3D(u,0,3,-1)) + 45*(KRANC_GFOFFSET3D(u,0,-1,-3) + KRANC_GFOFFSET3D(u,0,1,3) + KRANC_GFOFFSET3D(u,0,-3,-1) + KRANC_GFOFFSET3D(u,0,3,1)) + 9*(KRANC_GFOFFSET3D(u,0,-2,3) + KRANC_GFOFFSET3D(u,0,2,-3) + KRANC_GFOFFSET3D(u,0,-3,2) + KRANC_GFOFFSET3D(u,0,3,-2)) - 9*(KRANC_GFOFFSET3D(u,0,-2,-3) + KRANC_GFOFFSET3D(u,0,2,3) + KRANC_GFOFFSET3D(u,0,-3,-2) + KRANC_GFOFFSET3D(u,0,3,2)) + KRANC_GFOFFSET3D(u,0,-3,-3) - KRANC_GFOFFSET3D(u,0,-3,3) - KRANC_GFOFFSET3D(u,0,3,-3) + KRANC_GFOFFSET3D(u,0,3,3))*p1o3600dydz;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard6th31(u) ((-2025*(KRANC_GFOFFSET3D(u,-1,0,1) + KRANC_GFOFFSET3D(u,1,0,-1)) + 2025*(KRANC_GFOFFSET3D(u,-1,0,-1) + KRANC_GFOFFSET3D(u,1,0,1)) + 405*(KRANC_GFOFFSET3D(u,-1,0,2) + KRANC_GFOFFSET3D(u,1,0,-2) + KRANC_GFOFFSET3D(u,-2,0,1) + KRANC_GFOFFSET3D(u,2,0,-1)) - 405*(KRANC_GFOFFSET3D(u,-1,0,-2) + KRANC_GFOFFSET3D(u,1,0,2) + KRANC_GFOFFSET3D(u,-2,0,-1) + KRANC_GFOFFSET3D(u,2,0,1)) - 81*(KRANC_GFOFFSET3D(u,-2,0,2) + KRANC_GFOFFSET3D(u,2,0,-2)) + 81*(KRANC_GFOFFSET3D(u,-2,0,-2) + KRANC_GFOFFSET3D(u,2,0,2)) - 45*(KRANC_GFOFFSET3D(u,-1,0,3) + KRANC_GFOFFSET3D(u,1,0,-3) + KRANC_GFOFFSET3D(u,-3,0,1) + KRANC_GFOFFSET3D(u,3,0,-1)) + 45*(KRANC_GFOFFSET3D(u,-1,0,-3) + KRANC_GFOFFSET3D(u,1,0,3) + KRANC_GFOFFSET3D(u,-3,0,-1) + KRANC_GFOFFSET3D(u,3,0,1)) + 9*(KRANC_GFOFFSET3D(u,-2,0,3) + KRANC_GFOFFSET3D(u,2,0,-3) + KRANC_GFOFFSET3D(u,-3,0,2) + KRANC_GFOFFSET3D(u,3,0,-2)) - 9*(KRANC_GFOFFSET3D(u,-2,0,-3) + KRANC_GFOFFSET3D(u,2,0,3) + KRANC_GFOFFSET3D(u,-3,0,-2) + KRANC_GFOFFSET3D(u,3,0,2)) + KRANC_GFOFFSET3D(u,-3,0,-3) - KRANC_GFOFFSET3D(u,-3,0,3) - KRANC_GFOFFSET3D(u,3,0,-3) + KRANC_GFOFFSET3D(u,3,0,3))*p1o3600dxdz)
#else
#  define PDstandard6th31(u) (PDstandard6th31_impl(u,p1o3600dxdz,cdj,cdk))
static CCTK_REAL PDstandard6th31_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o3600dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard6th31_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o3600dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard6th12_impl(u, p1o3600dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandard6th32(u) ((-2025*(KRANC_GFOFFSET3D(u,0,-1,1) + KRANC_GFOFFSET3D(u,0,1,-1)) + 2025*(KRANC_GFOFFSET3D(u,0,-1,-1) + KRANC_GFOFFSET3D(u,0,1,1)) + 405*(KRANC_GFOFFSET3D(u,0,-1,2) + KRANC_GFOFFSET3D(u,0,1,-2) + KRANC_GFOFFSET3D(u,0,-2,1) + KRANC_GFOFFSET3D(u,0,2,-1)) - 405*(KRANC_GFOFFSET3D(u,0,-1,-2) + KRANC_GFOFFSET3D(u,0,1,2) + KRANC_GFOFFSET3D(u,0,-2,-1) + KRANC_GFOFFSET3D(u,0,2,1)) - 81*(KRANC_GFOFFSET3D(u,0,-2,2) + KRANC_GFOFFSET3D(u,0,2,-2)) + 81*(KRANC_GFOFFSET3D(u,0,-2,-2) + KRANC_GFOFFSET3D(u,0,2,2)) - 45*(KRANC_GFOFFSET3D(u,0,-1,3) + KRANC_GFOFFSET3D(u,0,1,-3) + KRANC_GFOFFSET3D(u,0,-3,1) + KRANC_GFOFFSET3D(u,0,3,-1)) + 45*(KRANC_GFOFFSET3D(u,0,-1,-3) + KRANC_GFOFFSET3D(u,0,1,3) + KRANC_GFOFFSET3D(u,0,-3,-1) + KRANC_GFOFFSET3D(u,0,3,1)) + 9*(KRANC_GFOFFSET3D(u,0,-2,3) + KRANC_GFOFFSET3D(u,0,2,-3) + KRANC_GFOFFSET3D(u,0,-3,2) + KRANC_GFOFFSET3D(u,0,3,-2)) - 9*(KRANC_GFOFFSET3D(u,0,-2,-3) + KRANC_GFOFFSET3D(u,0,2,3) + KRANC_GFOFFSET3D(u,0,-3,-2) + KRANC_GFOFFSET3D(u,0,3,2)) + KRANC_GFOFFSET3D(u,0,-3,-3) - KRANC_GFOFFSET3D(u,0,-3,3) - KRANC_GFOFFSET3D(u,0,3,-3) + KRANC_GFOFFSET3D(u,0,3,3))*p1o3600dydz)
#else
#  define PDstandard6th32(u) (PDstandard6th32_impl(u,p1o3600dydz,cdj,cdk))
static CCTK_REAL PDstandard6th32_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o3600dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandard6th32_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o3600dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandard6th23_impl(u, p1o3600dydz, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDplus1(u) ((-KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,1,0,0))*p1odx)
#else
#  define PDplus1(u) (PDplus1_impl(u,p1odx,cdj,cdk))
static CCTK_REAL PDplus1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDplus1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,1,0,0))*p1odx;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDplus2(u) ((-KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,0,1,0))*p1ody)
#else
#  define PDplus2(u) (PDplus2_impl(u,p1ody,cdj,cdk))
static CCTK_REAL PDplus2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1ody, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDplus2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1ody, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,0,1,0))*p1ody;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDplus3(u) ((-KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,0,0,1))*p1odz)
#else
#  define PDplus3(u) (PDplus3_impl(u,p1odz,cdj,cdk))
static CCTK_REAL PDplus3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDplus3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDplus2_impl(u, p1odz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDminus1(u) ((KRANC_GFOFFSET3D(u,0,0,0) - KRANC_GFOFFSET3D(u,-1,0,0))*p1odx)
#else
#  define PDminus1(u) (PDminus1_impl(u,p1odx,cdj,cdk))
static CCTK_REAL PDminus1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDminus1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (KRANC_GFOFFSET3D(u,0,0,0) - KRANC_GFOFFSET3D(u,-1,0,0))*p1odx;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDminus2(u) ((KRANC_GFOFFSET3D(u,0,0,0) - KRANC_GFOFFSET3D(u,0,-1,0))*p1ody)
#else
#  define PDminus2(u) (PDminus2_impl(u,p1ody,cdj,cdk))
static CCTK_REAL PDminus2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1ody, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDminus2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1ody, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (KRANC_GFOFFSET3D(u,0,0,0) - KRANC_GFOFFSET3D(u,0,-1,0))*p1ody;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDminus3(u) ((KRANC_GFOFFSET3D(u,0,0,0) - KRANC_GFOFFSET3D(u,0,0,-1))*p1odz)
#else
#  define PDminus3(u) (PDminus3_impl(u,p1odz,cdj,cdk))
static CCTK_REAL PDminus3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDminus3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDminus2_impl(u, p1odz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDplus11(u) ((KRANC_GFOFFSET3D(u,0,0,0) - 2*KRANC_GFOFFSET3D(u,1,0,0) + KRANC_GFOFFSET3D(u,2,0,0))*p1odx2)
#else
#  define PDplus11(u) (PDplus11_impl(u,p1odx2,cdj,cdk))
static CCTK_REAL PDplus11_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odx2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDplus11_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odx2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (KRANC_GFOFFSET3D(u,0,0,0) - 2*KRANC_GFOFFSET3D(u,1,0,0) + KRANC_GFOFFSET3D(u,2,0,0))*p1odx2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDplus12(u) ((KRANC_GFOFFSET3D(u,0,0,0) - KRANC_GFOFFSET3D(u,0,1,0) - KRANC_GFOFFSET3D(u,1,0,0) + KRANC_GFOFFSET3D(u,1,1,0))*p1odxdy)
#else
#  define PDplus12(u) (PDplus12_impl(u,p1odxdy,cdj,cdk))
static CCTK_REAL PDplus12_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDplus12_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (KRANC_GFOFFSET3D(u,0,0,0) - KRANC_GFOFFSET3D(u,0,1,0) - KRANC_GFOFFSET3D(u,1,0,0) + KRANC_GFOFFSET3D(u,1,1,0))*p1odxdy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDplus13(u) ((KRANC_GFOFFSET3D(u,0,0,0) - KRANC_GFOFFSET3D(u,0,0,1) - KRANC_GFOFFSET3D(u,1,0,0) + KRANC_GFOFFSET3D(u,1,0,1))*p1odxdz)
#else
#  define PDplus13(u) (PDplus13_impl(u,p1odxdz,cdj,cdk))
static CCTK_REAL PDplus13_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDplus13_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDplus12_impl(u, p1odxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDplus21(u) ((KRANC_GFOFFSET3D(u,0,0,0) - KRANC_GFOFFSET3D(u,0,1,0) - KRANC_GFOFFSET3D(u,1,0,0) + KRANC_GFOFFSET3D(u,1,1,0))*p1odxdy)
#else
#  define PDplus21(u) (PDplus21_impl(u,p1odxdy,cdj,cdk))
static CCTK_REAL PDplus21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDplus21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDplus12_impl(u, p1odxdy, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDplus22(u) ((KRANC_GFOFFSET3D(u,0,0,0) - 2*KRANC_GFOFFSET3D(u,0,1,0) + KRANC_GFOFFSET3D(u,0,2,0))*p1ody2)
#else
#  define PDplus22(u) (PDplus22_impl(u,p1ody2,cdj,cdk))
static CCTK_REAL PDplus22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1ody2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDplus22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1ody2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (KRANC_GFOFFSET3D(u,0,0,0) - 2*KRANC_GFOFFSET3D(u,0,1,0) + KRANC_GFOFFSET3D(u,0,2,0))*p1ody2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDplus23(u) ((KRANC_GFOFFSET3D(u,0,0,0) - KRANC_GFOFFSET3D(u,0,0,1) - KRANC_GFOFFSET3D(u,0,1,0) + KRANC_GFOFFSET3D(u,0,1,1))*p1odydz)
#else
#  define PDplus23(u) (PDplus23_impl(u,p1odydz,cdj,cdk))
static CCTK_REAL PDplus23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDplus23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (KRANC_GFOFFSET3D(u,0,0,0) - KRANC_GFOFFSET3D(u,0,0,1) - KRANC_GFOFFSET3D(u,0,1,0) + KRANC_GFOFFSET3D(u,0,1,1))*p1odydz;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDplus31(u) ((KRANC_GFOFFSET3D(u,0,0,0) - KRANC_GFOFFSET3D(u,0,0,1) - KRANC_GFOFFSET3D(u,1,0,0) + KRANC_GFOFFSET3D(u,1,0,1))*p1odxdz)
#else
#  define PDplus31(u) (PDplus31_impl(u,p1odxdz,cdj,cdk))
static CCTK_REAL PDplus31_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDplus31_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDplus12_impl(u, p1odxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDplus32(u) ((KRANC_GFOFFSET3D(u,0,0,0) - KRANC_GFOFFSET3D(u,0,0,1) - KRANC_GFOFFSET3D(u,0,1,0) + KRANC_GFOFFSET3D(u,0,1,1))*p1odydz)
#else
#  define PDplus32(u) (PDplus32_impl(u,p1odydz,cdj,cdk))
static CCTK_REAL PDplus32_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDplus32_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDplus23_impl(u, p1odydz, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDplus33(u) ((KRANC_GFOFFSET3D(u,0,0,0) - 2*KRANC_GFOFFSET3D(u,0,0,1) + KRANC_GFOFFSET3D(u,0,0,2))*p1odz2)
#else
#  define PDplus33(u) (PDplus33_impl(u,p1odz2,cdj,cdk))
static CCTK_REAL PDplus33_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odz2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDplus33_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odz2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDplus22_impl(u, p1odz2, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesided2nd1(u) ((3*KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,2*dir1,0,0) - 4*KRANC_GFOFFSET3D(u,dir1,0,0))*pm1o2dx*dir1)
#else
#  define PDonesided2nd1(u) (PDonesided2nd1_impl(u,pm1o2dx,cdj,cdk,dir1,dir2,dir3))
static CCTK_REAL PDonesided2nd1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o2dx, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDonesided2nd1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o2dx, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (3*KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,2*dir1,0,0) - 4*KRANC_GFOFFSET3D(u,dir1,0,0))*pm1o2dx*dir1;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesided2nd2(u) ((3*KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,0,2*dir2,0) - 4*KRANC_GFOFFSET3D(u,0,dir2,0))*pm1o2dy*dir2)
#else
#  define PDonesided2nd2(u) (PDonesided2nd2_impl(u,pm1o2dy,cdj,cdk,dir1,dir2,dir3))
static CCTK_REAL PDonesided2nd2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o2dy, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDonesided2nd2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o2dy, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (3*KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,0,2*dir2,0) - 4*KRANC_GFOFFSET3D(u,0,dir2,0))*pm1o2dy*dir2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesided2nd3(u) ((3*KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,0,0,2*dir3) - 4*KRANC_GFOFFSET3D(u,0,0,dir3))*pm1o2dz*dir3)
#else
#  define PDonesided2nd3(u) (PDonesided2nd3_impl(u,pm1o2dz,cdj,cdk,dir1,dir2,dir3))
static CCTK_REAL PDonesided2nd3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o2dz, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDonesided2nd3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o2dz, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDonesided2nd2_impl(u, pm1o2dz, cdk, cdj, dir1, dir2, dir3);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDlopsided4th1(u) ((-10*KRANC_GFOFFSET3D(u,0,0,0) - 6*KRANC_GFOFFSET3D(u,2*dir1,0,0) + KRANC_GFOFFSET3D(u,3*dir1,0,0) - 3*KRANC_GFOFFSET3D(u,-dir1,0,0) + 18*KRANC_GFOFFSET3D(u,dir1,0,0))*p1o12dx*dir1)
#else
#  define PDlopsided4th1(u) (PDlopsided4th1_impl(u,p1o12dx,cdj,cdk,dir1,dir2,dir3))
static CCTK_REAL PDlopsided4th1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o12dx, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDlopsided4th1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o12dx, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-10*KRANC_GFOFFSET3D(u,0,0,0) - 6*KRANC_GFOFFSET3D(u,2*dir1,0,0) + KRANC_GFOFFSET3D(u,3*dir1,0,0) - 3*KRANC_GFOFFSET3D(u,-dir1,0,0) + 18*KRANC_GFOFFSET3D(u,dir1,0,0))*p1o12dx*dir1;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDlopsided4th2(u) ((-10*KRANC_GFOFFSET3D(u,0,0,0) - 6*KRANC_GFOFFSET3D(u,0,2*dir2,0) + KRANC_GFOFFSET3D(u,0,3*dir2,0) - 3*KRANC_GFOFFSET3D(u,0,-dir2,0) + 18*KRANC_GFOFFSET3D(u,0,dir2,0))*p1o12dy*dir2)
#else
#  define PDlopsided4th2(u) (PDlopsided4th2_impl(u,p1o12dy,cdj,cdk,dir1,dir2,dir3))
static CCTK_REAL PDlopsided4th2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o12dy, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDlopsided4th2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o12dy, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-10*KRANC_GFOFFSET3D(u,0,0,0) - 6*KRANC_GFOFFSET3D(u,0,2*dir2,0) + KRANC_GFOFFSET3D(u,0,3*dir2,0) - 3*KRANC_GFOFFSET3D(u,0,-dir2,0) + 18*KRANC_GFOFFSET3D(u,0,dir2,0))*p1o12dy*dir2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDlopsided4th3(u) ((-10*KRANC_GFOFFSET3D(u,0,0,0) - 6*KRANC_GFOFFSET3D(u,0,0,2*dir3) + KRANC_GFOFFSET3D(u,0,0,3*dir3) - 3*KRANC_GFOFFSET3D(u,0,0,-dir3) + 18*KRANC_GFOFFSET3D(u,0,0,dir3))*p1o12dz*dir3)
#else
#  define PDlopsided4th3(u) (PDlopsided4th3_impl(u,p1o12dz,cdj,cdk,dir1,dir2,dir3))
static CCTK_REAL PDlopsided4th3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o12dz, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDlopsided4th3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o12dz, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDlopsided4th2_impl(u, p1o12dz, cdk, cdj, dir1, dir2, dir3);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDlopsided6th1(u) ((35*KRANC_GFOFFSET3D(u,0,0,0) - 2*KRANC_GFOFFSET3D(u,-2*dir1,0,0) + 30*KRANC_GFOFFSET3D(u,2*dir1,0,0) - 8*KRANC_GFOFFSET3D(u,3*dir1,0,0) + KRANC_GFOFFSET3D(u,4*dir1,0,0) + 24*KRANC_GFOFFSET3D(u,-dir1,0,0) - 80*KRANC_GFOFFSET3D(u,dir1,0,0))*pm1o60dx*dir1)
#else
#  define PDlopsided6th1(u) (PDlopsided6th1_impl(u,pm1o60dx,cdj,cdk,dir1,dir2,dir3))
static CCTK_REAL PDlopsided6th1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o60dx, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDlopsided6th1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o60dx, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (35*KRANC_GFOFFSET3D(u,0,0,0) - 2*KRANC_GFOFFSET3D(u,-2*dir1,0,0) + 30*KRANC_GFOFFSET3D(u,2*dir1,0,0) - 8*KRANC_GFOFFSET3D(u,3*dir1,0,0) + KRANC_GFOFFSET3D(u,4*dir1,0,0) + 24*KRANC_GFOFFSET3D(u,-dir1,0,0) - 80*KRANC_GFOFFSET3D(u,dir1,0,0))*pm1o60dx*dir1;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDlopsided6th2(u) ((35*KRANC_GFOFFSET3D(u,0,0,0) - 2*KRANC_GFOFFSET3D(u,0,-2*dir2,0) + 30*KRANC_GFOFFSET3D(u,0,2*dir2,0) - 8*KRANC_GFOFFSET3D(u,0,3*dir2,0) + KRANC_GFOFFSET3D(u,0,4*dir2,0) + 24*KRANC_GFOFFSET3D(u,0,-dir2,0) - 80*KRANC_GFOFFSET3D(u,0,dir2,0))*pm1o60dy*dir2)
#else
#  define PDlopsided6th2(u) (PDlopsided6th2_impl(u,pm1o60dy,cdj,cdk,dir1,dir2,dir3))
static CCTK_REAL PDlopsided6th2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o60dy, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDlopsided6th2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o60dy, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (35*KRANC_GFOFFSET3D(u,0,0,0) - 2*KRANC_GFOFFSET3D(u,0,-2*dir2,0) + 30*KRANC_GFOFFSET3D(u,0,2*dir2,0) - 8*KRANC_GFOFFSET3D(u,0,3*dir2,0) + KRANC_GFOFFSET3D(u,0,4*dir2,0) + 24*KRANC_GFOFFSET3D(u,0,-dir2,0) - 80*KRANC_GFOFFSET3D(u,0,dir2,0))*pm1o60dy*dir2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDlopsided6th3(u) ((35*KRANC_GFOFFSET3D(u,0,0,0) - 2*KRANC_GFOFFSET3D(u,0,0,-2*dir3) + 30*KRANC_GFOFFSET3D(u,0,0,2*dir3) - 8*KRANC_GFOFFSET3D(u,0,0,3*dir3) + KRANC_GFOFFSET3D(u,0,0,4*dir3) + 24*KRANC_GFOFFSET3D(u,0,0,-dir3) - 80*KRANC_GFOFFSET3D(u,0,0,dir3))*pm1o60dz*dir3)
#else
#  define PDlopsided6th3(u) (PDlopsided6th3_impl(u,pm1o60dz,cdj,cdk,dir1,dir2,dir3))
static CCTK_REAL PDlopsided6th3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o60dz, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDlopsided6th3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o60dz, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDlopsided6th2_impl(u, pm1o60dz, cdk, cdj, dir1, dir2, dir3);
}
#endif

