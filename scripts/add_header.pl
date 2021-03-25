#!/usr/bin/perl



   use warnings;
   use strict;
  
  print "\#CHROM\tPOS\tREF\tALT\tP_HAP1\tP_HAP2\tP_NAME\tM_HAP1\tM_HAP2\tM_NAME\tH_HAP1\tH_HAP2\tH_NAME\tL_HAP1\tL_HAP2\tL_NAME\n";

  while(<>){
        chomp;
     print "$_\n";
           }  
