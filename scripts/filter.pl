#!/usr/bin/perl

=head1 NAME
 
 perl filter.pl  ------ the script to filter SNPs based on the reads counts of minor allele

=head1 Usage

  perl filter.pl  -d <threshold> <Reads count file>

       
=head1 Example

=cut


 use warnings;
 use strict;
 use Getopt::Long;
 use Cwd;
 my ($h,$d);
        GetOptions(
           "d=i" =>\$d,
           "help"=>\$h);
 die `pod2text $0` if ( $h);

      unless(defined($d)){
                     $d=3; 
                         }

      my @a;

    print "\#CHROM\tPOS\tH_REF\tH_ALT\tL_REF\tL_ALT\n";
      
 while(<>){
     chomp;
     @a=split/\t/,$_;
     unless($a[2] <= $d && $a[4] <= $d){
     unless($a[3] <= $d && $a[5]<= $d){
               print "$_\n";
                }
       }
       }


  
              
 
  
