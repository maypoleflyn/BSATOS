#!/usr/bin/perl 

=head1 NAME
 
 get_rds.pl  ------ the script to assign RDS to genes in the QTL region

=head1 Usage


  perl get_rds.pl <gene bed file> <chromosome> <start> <end> <peak position> <window size>


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
 


  my $chr=$ARGV[1];
  my $start=$ARGV[2];
  my $end=$ARGV[3];
  my $pos=$ARGV[4];
  my $w=$ARGV[5]; 
  my $rds;
   
  my @a; 

  open LOD,"<$ARGV[0]";

print "#CHROM\tSTART\tEND\tGENE\tPEAK\tRDS\n";

 while(<LOD>){
     chomp;
     @a=split/\t/,$_;

       if($_=~/#/){
	     next;
	          }

   if($a[0] eq $ARGV[1] && ($a[1] >= $ARGV[2]) && ($a[1] <=$ARGV[3])){

         $rds=(abs($a[1]-$pos)/($w/2))*100;
	  
       print "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$pos\t$rds\n";
          

       }
       }  
  
