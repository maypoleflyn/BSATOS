#!/usr/bin/perl

=head1 NAME
 
  filter_snv_based_on_qtl.pl    ------ the script to filter SNVs based on the QTL type

=head1 Usage

  perl filter_snv_based_on_qtl.pl --i <QTL type> --d <promoter region> -q <snv quality>  --hap <enriched haplotype information> <annotated SNVs file from ANNOVAR>


        
=head1 Example

=cut

   use warnings;
   use strict;
 
   use Getopt::Long;
   use Cwd;
   my ($h,$i,$hap,$q,$d);
        GetOptions(
                  "hap=s" =>\$hap,
                  "i=s" =>\$i,
                  "q=i" =>\$q,
                  "d=i" =>\$d,            
                  "help"=>\$h);
 die `pod2text $0` if ( @ARGV==0 || $h || !defined($i));

      unless(defined($d)){
              $d=2000; 
                         }
      unless(defined($q)){
              $q=30;
                         }   
 





   open LKK,"<$hap";
   my %hap;
   my @p; 
   while(<LKK>){
           chomp;
           @p=split/\t/,$_;
           $hap{$p[0]}{$p[1]}=$_;
               }

    my @a;
    my $g1;
    my $g2;
  while(<>){
     chomp;
     @a=split/\t/,$_;
     if($a[18] >$d){
     if($a[22]=~/(.*?):/){
           $g1=$1;
                         }
     if($a[23]=~/(.*?):/){
           $g2=$1;
                         }
     if($i eq "P"){  
    
         if(($g1 eq "0/1") && ($g2 eq "0/0" || $g2 eq "1/1")){
             if($_=~/dist=(\d+)\;dist=(\d+)/){
                      if($1 >$d  && $2 >$d){
                               next;
                                }
                                }       
               unless($_=~/(\t)synonymous SNV/ || $_=~/intronic/ || $_=~/downstream/ ){

                  my $hp=$hap{$a[0]}{$a[1]}; 
                  unless(defined($hap{$a[0]}{$a[1]})){
                                 $hap{$a[0]}{$a[1]}="";
                                                     }
                  
               print "$a[0]\t$a[1]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]\t$a[8]\t$a[18]\t$a[22]\t$a[23]\t$hap{$a[0]}{$a[1]}\n";




                                }
                                }
                    }

      if($i eq "M"){

         if(($g2 eq "0/1") && ($g1 eq "0/0" || $g1 eq "1/1")){
             if($_=~/dist=(\d+)\;dist=(\d+)/){
                      if($1 >$d  && $2 >$d){
                               next;
                                }
                                }
                     unless($_=~/(\t)synonymous SNV/ || $_=~/intronic/ || $_=~/downstream/ ){

                   my $hp=$hap{$a[0]}{$a[1]};
                   unless(defined($hap{$a[0]}{$a[1]})){
                                 $hap{$a[0]}{$a[1]}="";
                                                     }
          
                     print "$a[0]\t$a[1]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]\t$a[8]\t$a[18]\t$a[22]\t$a[23]\t$hap{$a[0]}{$a[1]}\n";



                                }
                                }
                    } 

      if($i eq "H"){

         if(($g1 eq "0/1") && ($g2 eq "0/1")){
             if($_=~/dist=(\d+)\;dist=(\d+)/){
                      if($1 >$d  && $2 >$d){
                               next;
                                }
                                }
                     unless($_=~/(\t)synonymous SNV/ || $_=~/intronic/ || $_=~/downstream/ ){
                            
                     
                   my $hp=$hap{$a[0]}{$a[1]};
                  unless(defined($hap{$a[0]}{$a[1]})){
                                 $hap{$a[0]}{$a[1]}="";
                                                     }
           
                       print "$a[0]\t$a[1]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]\t$a[8]\t$a[18]\t$a[22]\t$a[23]\t$hap{$a[0]}{$a[1]}\n";



                                }
                                }
                    }
                  }  
                  }
     
