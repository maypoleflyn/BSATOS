#!/usr/bin/perl


=head1 NAME
  filter_get_genotype.pl ------ the script to get the genotype of 3-way 

=head1 Usage

  perl filter_get_genotype.pl --d <depth> --q <quality> --g <G1|G2|G3> --outprefix <prefix>

  --d: the theshold of coverage depth [30]
  --q: the log-phred SNP/Indel quality [30]
  --g: the genotype types; G1: the genotype of first parent (G1) is heterzygous; G2: the genotype of second parent (G2) is heterzyous; G3:both the parents is heterzygous; default: G1 & G2 & G3
  --outprefix <out prefix>

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

 die `pod2text $0` if ($h);
  
     unless(defined($d)){
                     $d=10;
                        }
     unless(defined($q)){
                     $q=30;
      
                   }
     
    my $o1=$o.".p";
    my $o2=$o.".m";
    my $o3=$o.".pm";
     
    open LOA,">$o1";
    open LOB,">$o2";
    open LOC,">$o3";

   print LOA "\#CHROM\tPOS\tREF\tALT\n";
   print LOB "\#CHROM\tPOS\tREF\tALT\n";
   print LOC "\#CHROM\tPOS\tREF\tALT\n";


    my @a; 
    my $dep;
    my $indel_r=0;
    my $snp_r=0;
    my $indel_a=0;
    my $snp_a=0;    
    my $snv_m=0;
    my $snv_p=0;
    my $snv_pm=0;



    while(<>){
        chomp;
       unless($_=~/#/){
          @a=split/\t/,$_; 
     
      if($_=~/INDEL/){
       $indel_r++;
                     }else{
       $snp_r++;
                          }
                     
          $_=~/DP=(.*?);/;


         $dep=$1;   
       if($a[5] >= $q && $dep >=$d){

         if($_=~/INDEL/){
            $indel_a++;
                     }else{
            $snp_a++;
                          }
       

              if($a[10]=~/0\/1/ && $a[9] =~/0\/0/){

                    $snv_m++;
                     print LOB "$a[0]\t$a[1]\t$a[3]\t$a[4]\n";
                      }
                 if($a[10]=~/0\/0/ && $a[9] =~/0\/1/){


                     $snv_p++;   
                     print LOA "$a[0]\t$a[1]\t$a[3]\t$a[4]\n";
                      }
                if($a[10]=~/0\/1/ && $a[9] =~/0\/1/){

                     $snv_pm++;

                     print LOC "$a[0]\t$a[1]\t$a[3]\t$a[4]\n";                  
                     } 
                      }
                      }
                    }

        
 


 print "summary SNVs data\n\nSNPs(befor filtering):\t$snp_r\nINDELs(befor filtering):\t$indel_r\nSNPs(after filtering):\t$snp_a\nINDELs(after filtering):\t$indel_a\nP type SNVs:\t$snv_p\nM type SNVs:\t$snv_m\nPM type SNVs:\t$snv_pm\n";


    



