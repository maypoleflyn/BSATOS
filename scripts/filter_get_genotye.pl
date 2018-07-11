#!/usr/bin/perl


=head1 NAME
  filter_get_genotype.pl ------ the script to get the genotype of 3-way 

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
                     $d=10;
                        }
     unless(defined($q)){
                     $q=30;
      
                   }
     
    my $o1=$o."_G1";
    my $o2=$o."_G2";
    my $o3=$o."_G3";
     
    open LOA,">$o1";
    open LOB,">$o2";
    open LOC,">$o3";
    my @a; 
    my $dep;
    while(<>){
        chomp;
       unless($_=~/#/){
          @a=split/\t/,$_; 
         $_=~/DP=(.*?);/;
         $dep=$1;   
          if($a[5] >= $q && $dep >=$d){
              if($a[10]=~/0\/1/ && $a[9] =~/0\/0/){
                     print LOB "$a[0]\t$a[1]\t$a[3]\t$a[4]\n";
                      }
                 if($a[10]=~/0\/0/ && $a[9] =~/0\/1/){
                     print LOA "$a[0]\t$a[1]\t$a[3]\t$a[4]\n";
                      }
                if($a[10]=~/0\/1/ && $a[9] =~/0\/1/){
                     print LOC "$a[0]\t$a[1]\t$a[3]\t$a[4]\n";                  
                     } 
                      }
                      }
                    }

        
 


