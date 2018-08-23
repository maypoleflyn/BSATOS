#!/usr/bin/perl


=head1 NAME
 
 perl get_counts.pl  ------ get reads counts data based on the SNP calling result. 

=head1 Usage

  perl get_counts.pl -d <dep> -q <quality> VCF file


=head1 Example

=cut

   use warnings;
   use strict;
   use Getopt::Long;
   use Cwd;
   my ($h,$d,$q);
        GetOptions(
           "d=i" =>\$d,
           "q=i" =>\$q,
           "help"=>\$h);
 die `pod2text $0` if ( $h);
       
       unless(defined($d)){
                      $d=30; 
                          }
       unless(defined($q)){
                      $q=30;
                          }
 
    my @a;
    my $dp;
    my $ref;
    my $alt;
    my $all;
    my $snp_index;
   

 while(<>){
       chomp;           
    unless($_=~/#/ ||$_=~/INDEL/ ){    
          @a=split/\t/,$_;
          if($a[5]>=$q && $a[4]!~/,/){ 
         
        $_=~/DP=(\d+);/;
        $dp=$1; 
        $_=~/DP4=(\d+)\,(\d+)\,(\d+)\,(\d+)/;
        $ref=$1+$2;
        $alt=$3+$4;
        $all=$ref+$alt;
        $snp_index=$alt/$all;
        if($all>=$d){
      
     
         print "$a[0]\t$a[1]\t$ref\t$alt\n";
          }}}
     }    
