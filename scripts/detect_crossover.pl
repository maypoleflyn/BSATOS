#!/usr/bin/perl

=head1 NAME
 
 detect_crossover.pl   ------ detect crossove between two parents  

=head1 Usage

  perl detect_crossover.pl [options] <overlapped> <haplotype file>

    
  --h:       help 

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
    my $fx;
    my %hash;
    my %hp;
    my $hn;
    my %ov;
    my $num=0;
    my $start;
    my $end;
    my $count=0;
    my %info;
    my @c=();
    my %ib;

  print "#SUB\tCHROM\tSTART\tEND\tHAP\tHAP_START\tHAP_END\n";

   while(<>){
          chomp;
     @a=split/\t/,$_;
    
     if($a[7] eq "-" || $a[10] eq "-"){
                next;
                                      } 
    if($a[7] eq $a[10]){
        # print "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[7]\t$a[10]\t+\n";
          $fx="+";
                       }else{

         # print "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[7]\t$a[10]\t-\n";
          $fx="-";    
                            }

    my $bl=$a[0]."_".$a[1]."_".$a[2];  
    $hp{$bl}++;
    
   if($hp{$bl} ==1){
          $num++;
          $hn="HAP".$num;
         
          $ov{$hn}=[$a[0],$a[1],$a[2]];    
           @c=();    
          my @end=();
          $count=1;
          %info=();             
                   }

 if($hp{$bl} >=1){   
                   
            push@c,$fx;   
       unless(defined($info{$count}{"s"})){ 
            $info{$count}{"s"}=$a[4]; 
       
                                          }
       unless(defined($info{$count}{"e"})){
             $info{$count}{"e"}= $info{$count}{"s"}; 
                                         }
           
         if($#c==1){           
           if($c[1] eq $c[0]){                   
            $info{$count}{"e"}=$a[4];
             if($ov{$hn}->[2] eq $a[4]){
 
               my $out1=$info{$count}{"s"};
               my $out2=$info{$count}{"e"};
               my $sub="SUB".$count;
            print "$sub\t$a[0]\t$out1\t$out2\t$hn\t$a[1]\t$a[2]\n";
                      
                                       }                       
                             }else{
                                     
               my $out1=$info{$count}{"s"};
               my $out2=$a[4];
            
                my $sub="SUB".$count;
            print "$sub\t$a[0]\t$out1\t$out2\t$hn\t$a[1]\t$a[2]\n"; 
                            
            $count++;
            $info{$count}{"s"}=$a[4]; 
            
                                }
            shift@c;
          

                   }
                  }
                   
        
  
         


                
         }       
     
          











