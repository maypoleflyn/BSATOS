#!/usr/bin/perl

   use warnings;
   use strict;
 
   my %hash1;
   my %hash2;
 
 open LOD,"<$ARGV[0]";
 open LOP,"<$ARGV[1]";
 open LOT,"<$ARGV[2]";
    my @a;
    my @b;
 
 while(<LOD>){
        chomp;
      @a=split/\t/,$_;
    $hash1{$a[0]}{$a[1]}=[$a[2],$a[3],$a[4]];
       }
 
        my $big;
        my $sm;
 while(<LOP>){
       chomp;
     @b=split/\t/,$_;
                   
         if($b[1] > $b[2]){
               $big=$b[1];
               $sm=$b[2];
              }
          if($b[1]< $b[2]){
               $big=$b[2];
               $sm=$b[1];
             }

      foreach my $index($sm .. $big){
                    

  
                    
    $hash2{$b[0]}{$index}=[$b[0],$b[1],$b[2],$b[3],$b[4],$b[5],$b[6]];
      # print "@{$hash2{$b[0]}{$index}}\n";  
     
          }
         }
       my @c;
       my @p;
       my @pp; 
  while(<LOT>){
       chomp;
     @c=split/\t/,$_; 
   if($c[11] >30){  
    if(defined($hash1{$c[0]}{$c[1]})){
    
     
          @p=@{$hash1{$c[0]}{$c[1]}};
      #   print "@p\n";
          }else{
             @p=("-","-","-");
              }
   if(defined($hash2{$c[0]}{$c[1]})){
     #    $pp="-"."\t"."-"."\t"."-"."\t"."-"."\t"."-"."\t"."-";
    
            @pp=@{$hash2{$c[0]}{$c[1]}};
             }else{
          @pp=("-","-","-","-","-","-");
           }
     print "$_\t@p\t@pp\n";
           }}



    
        
  




 
            
