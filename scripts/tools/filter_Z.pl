#!/usr/bin/perl

 
    use warnings;
    use strict;
    my @a;
    my $g1;
    my $g2;
  while(<>){
     chomp;
     @a=split/\t/,$_;




     if($a[24]=~/(.*?):/){
           $g1=$1;
          }
     if($a[25]=~/(.*?):/){
           $g2=$1;
          }
      if( ($g1 eq "0/1") && ($g2 eq "0/0" || $g2 eq "1/1" )){

         if($_=~/dist=(\d+)\;dist=(\d+)/){
                   if($1 >2000  && $2 >2000){
                             next;
                              }
                         }       
                     unless($_=~/(\t)synonymous SNV/ || $_=~/intronic/ || $_=~/downstream/ ){
                                print "$_\n";
                           }
                       
                
                }                 }  
