#!/usr/bin/perl


   use warnings;
   use strict;

   my @a;
   my $name;
   my $count;
 while(<>){
     chomp;
  unless($_=~/\*/){ 
    if($_=~/BLOCK/){
       $count++;
       $name="block".$count;
                   }
     else{
        @a=split/\t/,$_;        
        print "$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[1]\t$a[2]\t$name\n";
        
                    }
       }
     } 
