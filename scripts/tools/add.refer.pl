#!/usr/bin/perl

  use warnings;
  use strict;
   my %hash;
   my @a;
  open LOD,"</home/shen/project/work/big_project/qtlseq/tools/refer.anno";
    
     while(<LOD>){
            chomp;
         @a=split/\t/,$_;
         $hash{$a[0]}=$a[1];
          }

 
  open LOP,"<$ARGV[0]";
     my @b;

   while(<LOP>){
        chomp;
     @b=split/\t/,$_;
   unless($b[6]=~/,/){
       if($b[6]=~/gene:(.*)/){
      if(defined($hash{$1})){
              print "$_\t$hash{$1}\n";
                  }else{
                   print "$_\n";
           }
       }
     }else{

         print "$_\n";
        }
 }
        
              
   
  
