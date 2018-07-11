#!/usr/bin/perl


=head1 NAME
  filter_vcf.pl ------ the script to get the genotype of 3-way 

=head1 Usage

  perl filter_get_genotype.pl --d <depth> --q <quality> --g <G1|G2|G3> --outprefix <prefix>

  --d: the theshold of coverage depth
  --q: the log-phred SNP/Indel quality
  --g: the genotype types; G1: the genotype of first parent (G1) is heterzygous; G2: the genotype of second parent (G2) is heterzyous; G3:both the parents is heterzygous; default: G1 & G2 & G3
  --outprefix <out prefix>
=head1 Example

=cut

   use warnings;
   use strict;
   use Getopt::Long;
   use Cwd;
   my ($q,$g,$d,$h,$o);
  GetOptions(
           "g=s"=>\$g,
           "d=i"=>\$d,
           "q=i"=>\$q,
           "outprefix=s"=>\$o, 
           "help"=>\$h);

 die `pod2text $0` if (@ARGV==0 || $h);
  
     unless(defined($d)){
                     $d=5;
                        }
     unless(defined($q)){
                     $q=30;
      
                   }
    
     
   
  
    my @a; 
    my $dep;
    while(<>){
        chomp;
       unless($_=~/#/){
          @a=split/\t/,$_; 
         $_=~/DP=(.*?);/;
         $dep=$1;   
          if($a[5] >= $q && $dep >=$d){
              if( $_!~/\.\/\./ &&  $_=~/0\/1/){
                    print "$_\n";              
                     } 
                      }
                      }
                    }

        
 


