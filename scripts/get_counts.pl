#!/usr/bin/perl


 
    use warnings;
    use strict;
    my @a;
    my $dp;
    my $ref;
    my $alt;
    my $all;
    my $snp_index;
 while(<>){
       chomp;
  
            
    unless($_=~/#/ ||$_=~/INDEL/ ){
      @a=split/\t/,$_;
          if($a[5]>=30 && $a[4]!~/,/){ 
         
        $_=~/DP=(\d+);/;
        $dp=$1; 
        $_=~/DP4=(\d+)\,(\d+)\,(\d+)\,(\d+)/;
        $ref=$1+$2;
        $alt=$3+$4;
        $all=$ref+$alt;
        $snp_index=$alt/$all;
        if($all>=30){
      
     
         print "$a[0]\t$a[1]\t$ref\t$alt\n";
          }}}
     }    
