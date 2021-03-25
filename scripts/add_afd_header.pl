#!/usr/bin/perl


  use warnings;
  use strict;

  print "\#CHROM\tPOS\tH_REF\tH_ALT\tL_REF\tL_ALT\tH_REF_AF\tH_ALT_AF\tL_REF_AF\tL_ALT_AF\tG_VALUE\tGprimer_1M\tGprimer_0.75M\tGprimer_0.5M\tGprimer_0.25M\tCHROM\tPOS\tREF\tALT\tP_HAP1\tP_HAP2\tP_NAME\tM_HAP1\tM_HAP2\tM_NAME\tH_HAP1\tH_HAP2\tH_NAME\tL_HAP1\tL_HAP2\tL_NAME\n";
while(<>){
     chomp;
     print "$_\n";
        }
