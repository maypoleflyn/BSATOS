#!/usr/bin/perl

=head1 NAME
 
  filter_sv_based_on_qtl.pl    ------ the script to filter SVs based on the QTL type

=head1 Usage

  perl filter_sv_based_on_qtl.pl --i <QTL type> -d <promoter region>  <annotated SNVs file from ANNOVAR>


        
=head1 Example

=cut

   use warnings;
   use strict;
 
   use Getopt::Long;
   use Cwd;
   my ($h,$i,$d);
        GetOptions(
                  "i=s" =>\$i,
                  "help"=>\$h,
                  "d=i" =>\$d);
 die `pod2text $0` if ( @ARGV==0 || $h || !defined($i));



    unless(defined($d)){
                  $d=2000;
                       }
       
    my @a;
    my $g1;
    my $g2;
  while(<>){
     chomp;
     @a=split/\t/,$_;
     if($a[19] eq "PASS"){
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

                
               print "$a[0]\t$a[1]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]\t$a[8]\t$a[18]\t$a[22]\t$a[23]\n";




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
          
                     print "$a[0]\t$a[1]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]\t$a[8]\t$a[18]\t$a[22]\t$a[23]\n";

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
                            
              
                       print "$a[0]\t$a[1]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]\t$a[8]\t$a[18]\t$a[22]\t$a[23]\n";



                                }
                                }
                    }
                  }  
                  }
     
