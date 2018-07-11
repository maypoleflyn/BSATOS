#!/usr/bin/perl

  

 use warnings;
 use strict;

  my @a;
 while(<>){
    chomp;
   @a=split/\t/,$_;
  unless($a[2] <=3 && $a[4] <=3){
    unless($a[3]<=3 && $a[5]<=3){
               print "$_\n";
                }
       }
       }


  
              
 
  
