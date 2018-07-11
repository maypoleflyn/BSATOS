#!/usr/bin/perl

=head1 NAME
  preBSA.pl-----the pipeline used to prepare data for BSA

=head1 Usage

perl preBSA.pl 

-r1: paired1 reads
-r2: paired2 reads
-rg: reads group
-sp1: used to conduct QTL for parent1 (parent1 AB and parent2 AA)
-sp2: used to conduct QTL for parent2 (parent2 AB and parent1 AA)
-sp3: used to conduct QTL for parent2 (parent2 AB and parent1 AB)
-h: help





=head1 Example

=cut

   use warnings;
   use strict;
   use Getopt::Long;
   use Cwd;
    our ($r1,$r2,$rg,$sp1,$sp2,$sp3,$h,$cfg,$sp1_name,$sp2_name,$sp3_name,$ref,$pool,$sh);  
  GetOptions(
           "r1"=>\$r1,
           "r2"=>\$r2,
           "rg"=>\$rg,
           "sp1"=>\$sp1,
           "sp2"=>\$sp2,
           "cfg=s"=>\$cfg,   
           "help"=>\$h);

 die `pod2text $0` if (  $h);

   open LOD,"<$cfg";
   
   while(<LOD>){
           chomp;
    unless($_=~/#/){
       if($_=~/read1=(\S+)/){
            $r1=$1;
                   }
       if($_=~/read2=(\S+)/){
            $r2=$1;
             }
       if($_=~/readgroup=(\S+)/){
            $rg=$1;
             }
       if($_=~/parent1=(\S+)/){
            $sp1=$1;
           }
       if($_=~/pools_name=(\S+)/){
           $pool=$1;
          }  
  

        if($_=~/parent1_name=(\S+)/){
            $sp1_name=$1;
           }
        if($_=~/parent2_name=(\S+)/){
            $sp2_name=$1;
           }
         if($_=~/parent3_name=(\S+)/){
            $sp3_name=$1;
           }
           
       if($_=~/parent3=(\S+)/){
            $sp3=$1;
              }
         if($_=~/parent2=(\S+)/){
            $sp2=$1;
              }

          
       if($_=~/reference=(\S+)/){
                 $ref=$1;
           }
         
         }
        }   
     my $dir = getcwd;      

      $sh=$pool.".sh";

    open SH,">$sh";
 ######alignment using BWA                           
              
             my $bam=$pool.".bam";
             my $srt_bam=$pool.".srt.bam";
             my $rm_bam=$pool.".rm.bam";
             my $snp_res1=$pool."_".$sp1_name.".vcf";
             my $snp_res2=$pool."_".$sp2_name.".vcf";
             my $snp_res3=$pool."_".$sp3_name.".vcf";
             my $snp_res1_count=$pool."_".$sp1_name.".counts";
             my $snp_res2_count=$pool."_".$sp2_name.".counts";
             my $snp_res3_count=$pool."_".$sp3_name.".counts";
             my $get_count_pl=$dir."\/get_counts.pl";
     print SH "bwa mem -t 5 -R \"\@RG\\tID:$rg\\tSM:$rg\" $ref $r1 $r2 |samtools view -bS -> $bam
samtools sort -o $srt_bam $bam
samtools rmdup $srt_bam $rm_bam
samtools index $rm_bam
samtools mpileup -ugf $ref -q 30 -l $sp1 $rm_bam |bcftools call -m ->$snp_res1 
samtools mpileup -ugf $ref -q 30 -l $sp2 $rm_bam |bcftools call -m ->$snp_res2 
samtools mpileup -ugf $ref -q 30 -l $sp3 $rm_bam |bcftools call -m ->$snp_res3
perl $get_count_pl $snp_res1 > $snp_res1_count
perl $get_count_pl $snp_res2 > $snp_res2_count
perl $get_count_pl $snp_res3 > $snp_res3_count
gzip \*vcf 
mkdir $pool
mv \*bam \*counts \*vcf.gz *bam.bai $pool\n";                             
                 
        
