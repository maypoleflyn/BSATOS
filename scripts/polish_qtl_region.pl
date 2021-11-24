#!/usr/bin/perl

=head1 NAME
 
 polish_qtl_region.pl   ------ polish qtl region and correct the QTL peaks

=head1 Usage

  perl polish_qtl_region.pl <afd.ad> <theshold 1M> <threhold 0.75M> <threshold 0.5M> <threshold 0.25M>

        
 

=head1 Example

=cut

   use warnings;
   use strict;
   use Getopt::Long;
   use Cwd;
    my ($h);
        GetOptions(
           "help"=>\$h);
 die `pod2text $0` if ( @ARGV==0 || $h);
    

      my @a;
      my %p;
      my %m;
      my $c;
      my %b1;
      my %k1;
      my %b2;
      my %k2;
      my %b3;
      my %k3;
      my %b4;
      my %k4;     
     open LOD,"<$ARGV[0]";
    
     while(<LOD>){
           chomp;
       if($_=~/#/){
             next;
                  } 
     
          @a=split/\t/,$_;
          $c=$a[7]-$a[9];
       
       if($a[10] >=$ARGV[1] ){
            $b1{$a[22]}{"T"}++;
            $k1{$a[25]}{"T"}++; 
                           }else{
             $b1{$a[22]}{"F"}++; 
             $k1{$a[25]}{"F"}++;
                                } 
       if($a[11] >=$ARGV[2]){
             $b2{$a[22]}{"T"}++;         
             $k2{$a[25]}{"T"}++;
                            }else{    
             $b2{$a[22]}{"F"}++;    
             $k2{$a[25]}{"F"}++;

                                 }
  
        if($a[12] >=$ARGV[3]){
             $b3{$a[22]}{"T"}++;    
             $k3{$a[25]}{"T"}++;
                            }else{
             $b3{$a[22]}{"F"}++; 
             $k3{$a[25]}{"F"}++;
                                 }
              
          if($a[13] >=$ARGV[4]){
             $b4{$a[22]}{"T"}++;    
             $k4{$a[25]}{"T"}++;
                            }else{
             $b4{$a[22]}{"F"}++; 
             $k4{$a[25]}{"F"}++;

                                 }

        
          if($a[21] eq "1"){
            if($c >0){           
            $p{$a[22]}{"H"}++;  
   
                     }else{
            $p{$a[22]}{"L"}++;
                          }
                    }
           if($a[21] eq "0"){
                 if($c <0){
                 $p{$a[22]}{"H"}++;
                          }else{
                 $p{$a[22]}{"L"}++;
                               }
                             }
            if($a[24] eq "1"){
                  if($c >0){
                    $m{$a[25]}{"H"}++;
                     }else{
                    $m{$a[25]}{"L"}++;
                          }
                    }
            if($a[24] eq "0"){
                  if($c <0){
                    $m{$a[25]}{"H"}++;
                     }else{
                    $m{$a[25]}{"L"}++;
                          }
                    }
                }
 
       close(LOD); 
              
    open LOP,"<$ARGV[0]";

    my @b;
    my $fp;
    my $lp;
    my $fm;
    my $mp;
    my $c2;
    my %remove;
    my ($n1,$n2,$n3,$n4,$n5,$n6,$n7,$n8);
   while(<LOP>){
         chomp;
    if($_=~/#/){
               next;
               }
 
       @b=split/\t/,$_;
   
      unless(defined($b1{$b[22]}{"T"})){
                  $b1{$b[22]}{"T"}=0;
                                       }        
      unless(defined($b1{$b[22]}{"F"})){
                  $b1{$b[22]}{"F"}=0;
                                       }
      unless(defined($b2{$b[22]}{"T"})){
                  $b2{$b[22]}{"T"}=0;
                                       }
      unless(defined($b2{$b[22]}{"F"})){
                  $b2{$b[22]}{"F"}=0;
                                       }
      unless(defined($b3{$b[22]}{"T"})){
                  $b3{$b[22]}{"T"}=0;
                                       }
      unless(defined($b3{$b[22]}{"F"})){
                  $b3{$b[22]}{"F"}=0;
                                       }
      unless(defined($b4{$b[22]}{"T"})){
                   $b4{$b[22]}{"T"}=0;
                                       }
      unless(defined($b4{$b[22]}{"F"})){
                  $b4{$b[22]}{"F"}=0;
                                       }

       unless(defined($k1{$b[25]}{"T"})){
                  $k1{$b[25]}{"T"}=0;
                                       }
      unless(defined($k1{$b[25]}{"F"})){
                  $k1{$b[25]}{"F"}=0;
                                       }
      unless(defined($k2{$b[25]}{"T"})){
                  $k2{$b[25]}{"T"}=0;
                                       }
      unless(defined($k2{$b[25]}{"F"})){
                  $k2{$b[25]}{"F"}=0;
                                       }
      unless(defined($k3{$b[25]}{"T"})){
                  $k3{$b[25]}{"T"}=0;
                                       }
      unless(defined($k3{$b[25]}{"F"})){
                  $k3{$b[25]}{"F"}=0;
                                       }
      unless(defined($k4{$b[25]}{"T"})){
                   $k4{$b[25]}{"T"}=0;
                                       }
      unless(defined($k4{$b[25]}{"F"})){
                  $k4{$b[25]}{"F"}=0;
                                       }  
       
   

      
 
      $n1=$b1{$b[22]}{"T"}/($b1{$b[22]}{"T"}+$b1{$b[22]}{"F"});
      $n2=$k1{$b[25]}{"T"}/($k1{$b[25]}{"T"}+$k1{$b[25]}{"F"}); 
      $n3=$b2{$b[22]}{"T"}/($b2{$b[22]}{"T"}+$b2{$b[22]}{"F"});
      $n4=$k2{$b[25]}{"T"}/($k2{$b[25]}{"T"}+$k2{$b[25]}{"F"});
      $n5=$b3{$b[22]}{"T"}/($b3{$b[22]}{"T"}+$b3{$b[22]}{"F"});
      $n6=$k3{$b[25]}{"T"}/($k3{$b[25]}{"T"}+$k3{$b[25]}{"F"});
      $n7=$b4{$b[22]}{"T"}/($b4{$b[22]}{"T"}+$b4{$b[22]}{"F"});
      $n8=$k4{$b[25]}{"T"}/($k4{$b[25]}{"T"}+$k4{$b[25]}{"F"});
  
     unless($n1 < 0.7 && $n2 < 0.7 && $n3 < 0.7 && $n4 < 0.7 && $n5 < 0.7 && $n6 < 0.7 && $n8 < 0.7){
               next;   
               
                                                                                          }
                                      
       $c2=$b[7]-$b[9];
       unless($p{$b[22]}{"H"}){
                   $p{$b[22]}{"H"}=0;
                              }
       unless($p{$b[22]}{"L"}){
                   $p{$b[22]}{"L"}=0;
                              }
       unless($m{$b[25]}{"H"}){
                   $m{$b[25]}{"H"}=0;
                              }
       unless($m{$b[25]}{"L"}){
                   $m{$b[25]}{"L"}=0;
                              }    
       $lp=$p{$b[22]}{"L"}+$p{$b[22]}{"H"};
       $mp=$m{$b[25]}{"L"}+$m{$b[25]}{"H"};
       if($lp >5){
       $fp=$p{$b[22]}{"H"}/($p{$b[22]}{"L"}+$p{$b[22]}{"H"});
                 }else{
                 $fp=0;
                      } 
 
       if($mp >5){    
       $fm=$m{$b[25]}{"H"}/($m{$b[25]}{"L"}+$m{$b[25]}{"H"});
                 }else{
                   $fm=0; 
                      }
                    
       if($b[21] eq "1"){
              if($fp >0.7){
                    if($c2<0){
                         # print "#$_\n";
                         $remove{$b[0]}{$b[1]}++;
                             }
                       } 
              if($fp <0.3 && $fp >0 ){
                    if($c2 > 0){
                        # print "#$_\n";
                         $remove{$b[0]}{$b[1]}++;
                             }
                       }
                         } 
         if($b[21] eq "0"){
                  if($fp >0.7){
                        if($c2 >0){
                          #  print "#$_\n";
                          $remove{$b[0]}{$b[1]}++;
                                  }
                            }
                  if($fp < 0.3 && $fp >0 ){
                        if( $c2 < 0){
                           # print "#$_\n";
                           $remove{$b[0]}{$b[1]}++;
                                  }
                           }
                        }

              if($b[24] eq "1"){
              if($fm >0.7){
                    if($c2<0){
                         # print "#$_\n";
                         $remove{$b[0]}{$b[1]}++;

                             }
                       }
              if($fm <0.3 && $fm >0 ){
                    if($c2 > 0){
                        # print "#$_\n";
                         $remove{$b[0]}{$b[1]}++;
                             }
                       }
                         }
  
            if($b[24] eq "0"){
                  if($fm >0.7){
                        if($c2 >0){
                          #  print "#$_\n";
                         $remove{$b[0]}{$b[1]}++;
                                  }
                            }
                  if($fm < 0.3 && $fm >0 ){
                        if( $c2 < 0){
                       #     print "#$_\n";
                         $remove{$b[0]}{$b[1]}++;
                                  }
                           }
                     }                          
                     }

       close(LOP);

       open LK,"<$ARGV[0]";
       my @c;
       while(<LK>){
              chomp;
                      
          if($_=~/#/){
                  next;
                     }   
             @c=split/\t/,$_;
             unless(defined($remove{$c[0]}{$c[1]})){
                  print "$c[0]\t$c[1]\t$c[2]\t$c[3]\t$c[4]\t$c[5]\n";
                                                    }
                 }
      
    
      
                            
                        
   
                     
                   
        
              
   
   
