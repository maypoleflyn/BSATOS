#!/usr/bin/perl


  use warnings;
  use strict;
  
   open LOD,"<$ARGV[0]";
   open LOP,"<$ARGV[1]";
   open LOT,"<$ARGV[2]";
   
   my @a;
   my %hash; 
  while(<LOD>){
         chomp;
      @a=split/\t/,$_;
      $hash{$a[0]}++;
      }

    my @b;
    my %hash1;
   
  while(<LOP>){
         chomp;
        @b=split/\t/,$_;
        $hash1{$b[0]}++;
        }

   my @c;
   my %hash2;

 while(<LOT>){
         chomp;
       @c=split/\t/,$_;
       $hash2{$c[0]}++;
     }


 my @k=sort keys%hash;
 

  foreach my $index(@k){
           print "$index\t$hash{$index}\t$hash1{$index}\t$hash2{$index}\n";
          }



















  
