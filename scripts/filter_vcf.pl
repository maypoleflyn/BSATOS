#!/usr/bin/perl


=head1 NAME

  filter_vcf.pl ------ the script to filter the VCF file

=head1 Usage

  perl filter_vcf.pl  --d <depth> --q <quality> <VCF file>

  --d: the theshold of coverage depth
  --q: the log-phred SNP/Indel quality

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
                     $d=10;
                        }
     unless(defined($q)){
                     $q=40;
      
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

        
 


