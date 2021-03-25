#!/usr/bin/perl


   use warnings;
   use strict;

   my @a;
 while(<>){
    chomp;
  unless($_=~/#/){
    @a=split/\t/,$_;
    my $start=$a[1]-1;
   print "$a[0]\t$start\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]\t$a[8]\t$a[9]\n";
          }
            }
   
