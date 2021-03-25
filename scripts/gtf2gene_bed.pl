#!/usr/bin/perl

 
  use warnings;
  use strict;
  my %hash;
  my @a;
  my $gn;
  my %chr;
  while(<>){
        chomp;
       @a=split/\t/,$_;
       if($_=~/gene_name "(\S+)";/){
          $gn=$1;
          $chr{$gn}=$a[0];
                                   }
                            
       if($a[2] eq "exon"){
             push@{$hash{$gn}},$a[3];
             push@{$hash{$gn}},$a[4];
                         }
           }

     
 my @gene=keys%hash;
 
    foreach my $index(@gene){
                 my @pos=sort {$a<=>$b} @{$hash{$index}};
                 print "$chr{$index}\t$pos[0]\t$pos[$#pos]\t$index\n";
                            }   


   
           
  
