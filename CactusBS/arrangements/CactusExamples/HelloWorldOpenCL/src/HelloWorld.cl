// -*-C++-*-

// This is an OpenCL kernel, without the surrounding function
// declaration. Such a declaration will be added by OpenCLRunTime,
// passing the Cactus cctkGH in a way such that (most) standard Cactus
// grid information (the cctk_* variables) can be accessed. That is,
// the code below will be placed into a function.

if (get_global_id(0) == 0 && get_global_id(1) == 0 && get_global_id(2) == 0) {
  printf("Hello, World! This is iteration %d\n", cctk_iteration);
}
