#! /usr/bin/perl
use File::Copy;

die "Usage: ./MakeEEBwithAlts.perl <kranc exe location>\n", if ($#ARGV != 0); 
$kranc = $ARGV[0];

# create thorn
system "$kranc ExternalEtaBeta.m\n";

@file = ('ExternalEtaBeta/src/eta_beta_calc_etaBeta_exp.cc','ExternalEtaBeta/src/eta_beta_calc_etaBeta_poly.cc');

# loop over files with placeholders
for $file (@file)
{
   # create backup file
   system "mv $file $file.bak\n";
   
   # read in file, and get rid of placeholders with spherical surface
   # data
   open FILE, "$file.bak" or die $!;
   $outfile = '';
   while(<FILE>)
   {
       if ( /CCTK_REAL xx1 = 0;/ )
       {
           $lines = <<EOF;
        assert(surfindex1 >=0 && surfindex1 < nsurfaces);
        sf_valid[surfindex1]=1;
        sf_active[surfindex1]=1;
        CCTK_REAL xx1 = sf_centroid_x[surfindex1];
        CCTK_REAL yy1 = sf_centroid_y[surfindex1];
        CCTK_REAL zz1 = sf_centroid_z[surfindex1];
  
        assert(surfindex2 >=0 && surfindex2 < nsurfaces);
        sf_valid[surfindex2]=1;
        sf_active[surfindex2]=1;
        CCTK_REAL xx2 = sf_centroid_x[surfindex2];
        CCTK_REAL yy2 = sf_centroid_y[surfindex2];
        CCTK_REAL zz2 = sf_centroid_z[surfindex2];
EOF
        $outfile = $outfile . $lines;
       }
       elsif ( /CCTK_REAL yy1 = 0;/ ) {}
       elsif ( /CCTK_REAL zz1 = 0;/ ) {}
       elsif ( /CCTK_REAL xx2 = 0;/ ) {}
       elsif ( /CCTK_REAL yy2 = 0;/ ) {}
       elsif ( /CCTK_REAL zz2 = 0;/ ) {}
       else
       {
        $outfile = $outfile . $_;
       }
   }
   close FILE;

   # output new file
   open OUTFILE, ">$file";
   print OUTFILE $outfile;
   close OUTFILE;

   # remove backup file
   system "rm $file.bak";
}

# edit param.ccl to find spherical surfaces stuff
system "mv ExternalEtaBeta/param.ccl ExternalEtaBeta/param.ccl.bak";
open PAR, "ExternalEtaBeta/param.ccl.bak";
# append lines to the top param.ccl
$parout = "shares: SphericalSurface\n";
$parout = $parout . "USES CCTK_INT nsurfaces\n";
while(<PAR>)
{
  $parout = $parout . $_;
}
close PAR;

# output new param.ccl
open PAROUT, ">ExternalEtaBeta/param.ccl";
print PAROUT $parout;

# remove backup
system "rm ExternalEtaBeta/param.ccl.bak";
close PAROUT;

# put a copy of MakeEEBwithAlts.perl and ExternalEtaBeta.m into thorn
system "mkdir -p ExternalEtaBeta/m";
system "cp ExternalEtaBeta.m ExternalEtaBeta/m";
system "cp MakeEEBWithAlts.perl ExternalEtaBeta/m";
