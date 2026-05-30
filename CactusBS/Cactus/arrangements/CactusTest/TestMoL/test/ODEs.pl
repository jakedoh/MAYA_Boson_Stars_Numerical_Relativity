#! /usr/bin/env perl

use warnings;
use strict;

my %ODE;
# each entry of %ODE is a reference to a hash describing single ODE stepper.
# elements of the hash:
# method
# steps
# nscratch
# RHS
# RHSSlow - empty if no slow evolution
# gf - evolve a grid function
# ga - evolve a grid array
# RHS_zero - zero RHS before each step?
my @elems = ("method", "steps", "nscratch", "RHS", "RHSSlow", "gf", "ga", "RHS_zero");
$ODE{"RK4"} = {
  method=>"rk4", steps=>4, nscratch=>1, RHS=>"t**3", RHSSlow=>"", gf=>1, ga=>1, RHS_zero=>1
};
$ODE{"RK3"} = {
  method=>"rk3", steps=>3, nscratch=>1, RHS=>"t**3", RHSSlow=>"", gf=>1, ga=>1, RHS_zero=>1
};
$ODE{"RK2"} = {
  method=>"rk2", steps=>2, nscratch=>1, RHS=>"t", RHSSlow=>"", gf=>1, ga=>1, RHS_zero=>1
};
$ODE{"RK2-central"} = {
  method=>"rk2-central", steps=>2, nscratch=>1, RHS=>"t", RHSSlow=>"", gf=>1, ga=>1, RHS_zero=>1
};
$ODE{"ICN"} = {
  method=>"icn", steps=>3, nscratch=>0, RHS=>"t", RHSSlow=>"", gf=>1, ga=>1, RHS_zero=>1
};
$ODE{"ICN-avg"} = {
  method=>"icn-avg", steps=>3, nscratch=>0, RHS=>"1", RHSSlow=>"", gf=>1, ga=>1, RHS_zero=>1
};
$ODE{"RK45"} = {
  method=>"rk45", steps=>6, nscratch=>6, RHS=>"t**4", RHSSlow=>"", gf=>1, ga=>0, RHS_zero=>1
};
$ODE{"RK45CK"} = {
  method=>"rk45ck", steps=>6, nscratch=>6, RHS=>"t**4", RHSSlow=>"", gf=>1, ga=>0, RHS_zero=>1
};
$ODE{"RK65"} = {
  method=>"rk65", steps=>8, nscratch=>8, RHS=>"t**5", RHSSlow=>"", gf=>1, ga=>0, RHS_zero=>1
};
$ODE{"RK87"} = {
  method=>"rk87", steps=>13, nscratch=>13, RHS=>"t**7", RHSSlow=>"", gf=>1, ga=>0, RHS_zero=>1
};
$ODE{"RK4-RK2"} = {
  method=>"rk4-rk2", steps=>4, nscratch=>1, RHS=>"t**3", RHSSlow=>"t", gf=>1, ga=>0, RHS_zero=>1
};
$ODE{"RK2-MR-2-1"} = {
  method=>"rk2-mr-2:1", steps=>5, nscratch=>5, RHS=>"t", RHSSlow=>"t", gf=>1, ga=>0, RHS_zero=>0
};
$ODE{"RK4-MR-2-1"} = {
  method=>"rk4-mr-2:1", steps=>10, nscratch=>10, RHS=>"t**3", RHSSlow=>"t**3", gf=>1, ga=>0, RHS_zero=>0
};

foreach my $name (keys %ODE) {
  my ($ode_method, $MoL_Intermediate_Steps, $MoL_Num_Scratch_Levels,
      $RHSexpression, $RHSSlowexpression, $evolve_grid_function,
      $evolve_grid_array, $init_RHS_zero) = @{$ODE{$name}}{@elems};

  # construct a list of output variables
  # for multi rate methods include both slow and fast grid functions
  my @gfs = ('TestMoL::analytic_gf', 'TestMoL::diff_gf',
            'TestMoL::evolved_gf', 'TestMoL::constrained_gf',
            'TestMoL::sandr_gf');
  my @slowgfs = ('TestMoL::analyticslow_gf', 'TestMoL::diffslow_gf',
                 'TestMoL::evolvedslow_gf');
  push (@gfs, @slowgfs) if $RHSSlowexpression;
  my @gas = ('TestMoL::analytic_ga', 'TestMoL::diff_ga',
             'TestMoL::evolved_ga', 'TestMoL::constrained_ga',
             'TestMoL::sandr_ga');

  my $out_vars = "";
  $out_vars = join "\n        ",($out_vars, @gfs) if $evolve_grid_function;
  $out_vars = join "\n        ",($out_vars, @gas) if $evolve_grid_array;

  # norm of the difference to the analytic solution
  my $out_norms = "";
  if($evolve_grid_function) {
    $out_norms = join "\n        ",($out_norms, "TestMoL::diff_gf");
  }
  if($evolve_grid_function and $RHSSlowexpression) {
    $out_norms = join "\n        ",($out_norms, "TestMoL::diffslow_gf");
  }
  if($evolve_grid_array) {
    $out_norms = join "\n        ",($out_norms, "TestMoL::diff_ga");
  }

  # $RHSSlowexpression is unset unless a multi-rate method is used
  # however we need to provide some value for the parfile parser
  if (not $RHSSlowexpression) {
    $RHSSlowexpression = "1"; # provide a dummy expression
  }

  my $lines = <<EOF;
# !DESC verify convergence order of ODE steppers inside of MoL
# Do not modify this file by hand, instead modify ODEs.pl and regenerate this
# file.
ActiveThorns = "
        IOASCII
        IOBasic
        IOUtil
        LocalReduce
        MoL
        PUGH
        PUGHReduce
        PUGHSlab
        TestMoL
"

Cactus::terminate = "iteration"
Cactus::cctk_itlast = 3
PUGH::global_nx = 1
PUGH::global_ny = 1
PUGH::global_nz = 10

MoL::ode_method = "$ode_method"
MoL::MoL_Intermediate_Steps = $MoL_Intermediate_Steps
MoL::MoL_Num_Scratch_Levels = $MoL_Num_Scratch_Levels
MoL::init_RHS_zero = $init_RHS_zero

TestMoL::RHSexpression = "$RHSexpression"
TestMoL::RHSSlowexpression = "$RHSSlowexpression"
TestMoL::evolve_grid_function = $evolve_grid_function
TestMoL::evolve_grid_array = $evolve_grid_array

IO::out_dir			= \$parfile
IO::out_fileinfo                = "none"
IO::parfile_write               = "no"

IOASCII::out1D_every = 1
IOASCII::out1D_d = "no"
IOASCII::out1D_x = "yes"
IOASCII::out1D_y = "no"
IOASCII::out1D_z = "no"
IOASCII::out1D_vars  = "$out_vars"

IOBasic::outScalar_every = 1
IOBasic::outScalar_reductions = "norm_inf"
IOBasic::outScalar_vars  = "$out_norms"
EOF

  open(my $fh, ">", "$name.par");
  print $fh $lines;
  close $fh;
}
