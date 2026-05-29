#! /usr/bin/perl
use File::Copy;


@file = <./Kranc2BSSNChi/src/*>;

# loop over files
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
        my $replacestring = "CCTK_ATTRIBUTE_UNUSED"; 
        $_ =~ s/$replacestring//g;
        $outfile = $outfile . $_;
   }
   close FILE;

   # output new file
   open OUTFILE, ">$file";
   print OUTFILE $outfile;
   close OUTFILE;

   # remove backup file
   system "rm $file.bak";
}
