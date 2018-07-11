#!/usr/bin/perl


  use warnings;
  use strict;
  
  open LOD,"<$ARGV[0]";
  open LOP,"<$ARGV[1]";
    my %hash;
    my @a;
    my @b;
  while(<LOD>){
       chomp;
     @a=split/\t/,$_;
    $hash{$a[0]}{$a[1]}=$_;
         }

   while(<LOP>){
          chomp;
      @b=split/\t/,$_;
    if(defined($hash{$b[0]}{$b[1]})){
           print "$_\t$hash{$b[0]}{$b[1]}\n";
          }else{
        print "$_\n";
             }
         }
        
