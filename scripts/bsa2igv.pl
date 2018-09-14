#!/usr/bin/perl

=head1 NAME
 
   bsa2igv.pl  ------ the script to generate inputs files for IGV

=head1 Usage

   perl bsa2igv.pl --w <window size> -i <P/M/H/> 

   --h: help information
   --w: 1 or 0.5 or 0.25 or 0.125 [1]
   --i: maker type [P/M/H] 
        
=head1 Example

=cut

   use warnings;
   use strict;
 
   use Getopt::Long;
   use Cwd;
   my ($h,$w,$i);
        GetOptions(
                  "w=i" =>\$w,      
                  "i=s" =>\$i, 
                  "help"=>\$h);
 
 die `pod2text $0` if ( @ARGV==0 || $h || !defined($i));

  unless(defined($w)){
           $w=1;
                     }

    open IN,"<$ARGV[0]";
    my @a; 
    
    print "Chromosome\tStart\tEND\tFeature\t$i\n";

   while(<IN>){
          chomp;
    @a=split/\t/,$_;
   
    if($w==1){
        
          my $s=$a[1]-1;
          print "$a[0]\t$s\t$a[1]\t$i\t$a[11]\n";
            }
    if($w==0.5){
         
         my $s=$a[1]-1;
          print "$a[0]\t$s\t$a[1]\t$i\t$a[12]\n";
               }
    if($w==0.25){
         
           my $s=$a[1]-1;
          print "$a[0]\t$s\t$a[1]\t$i\t$a[13]\n";
                }
    if($w==0.125){
         
          my $s=$a[1]-1;
          print "$a[0]\t$s\t$a[1]\t$i\t$a[14]\n";
                 }
            }    

  
 
 


