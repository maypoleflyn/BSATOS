#!/usr/bin/perl


   use warnings;
   use strict;
 
 

  my @a;
  my %hash;
  my @c; 
 while(<>){
        chomp;
    @a=split/\t/,$_; 
    if($a[4] eq "0"){
          $hash{$a[6]}{"0"}++;
                   } 
     if($a[4] eq "1"){
          $hash{$a[6]}{"1"}++;
                     }
      push@c,$_; 
             }
  
 
     my @b;
  
     foreach (@c){ 
     @b=split/\t/,$_;
   # print "$hash{$b[6]}{"0"}\n";
     unless(defined($hash{$b[6]}{"0"})){
              $hash{$b[6]}{"0"}=0;
                                       }
      unless(defined($hash{$b[6]}{"1"})){
              $hash{$b[6]}{"1"}=0;
                                        }
     if($hash{$b[6]}{"0"} > $hash{$b[6]}{"1"}){
               
                         print "$b[0]\t$b[1]\t$b[2]\t$b[3]\t$b[4]\t$b[5]\t$b[6]\n";
                              }else{
                 
                         print "$b[0]\t$b[1]\t$b[2]\t$b[3]\t$b[5]\t$b[4]\t$b[6]\n";
                                }     
                              }
              
        
                 
               

              
      
                   
                        
 
     
  
  
