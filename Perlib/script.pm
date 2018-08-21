

package script;

BEGIN {
	$VERSION = '1.0';
}


use strict;
use Getopt::Long;
use base 'Exporter';
use File::Basename;
use FindBin qw($Bin);
use Cwd; 


our @EXPORT = qw(runscript);

my $G_USAGE = "
 
Usage: bsatos script [options]  

  Options: 
         --oprefix STR       the prefix of the result folder [script] 
         --r       FILE      the genome fasta file [fasta]
         --gtf     FILE      the GTF/GFF file of the genes
         --pf1     FILE      the paired-end1 fastq file of the pollen parent    | |    | --pb FILE  BAM file of the pollen parent
         --pf2     FILE      the paired-end2 fastq file of the pollen parent    | |    |
         --mf1     FILE      the paired-end1 fastq file of the maternal parent  | |    | --mb FILE  BAM file of the maternal parent 
         --mf2     FILE      the paired-end2 fastq file of the maternal parent  | | OR |     
         --hf1     FILE      the paired-end1 fastq file of the H pool reads     | |    | --hb FILE  BAM file of the H pool 
         --hf2     FILE      the paired-end2 fastq file of the H pool reads     | |    |
         --lf1     FILE      the paired-end1 fastq file of the L pool reads     | |    | --lb FILE  BAM file of the L pool
         --lf2:    FILE      the paired-end2 fastq file of the L pool reads     | |    |
 

1) prepar step options:  

Aligment Options:

           --t  INT        number of threads [1]

SNV Calling Options:

           --aq  INT       skip alignments with mapQ smaller than INT [30]

SV Calling Options:
          
           --vq INT        skip alignments with mapQ smaller than INT [30]

SNV filtering Options:

           --d  INT        skip SNVs with reads depth smaller than INT [10]   
           --sq INT        skip SNVs with phred-scaled quality smaller than [30]

2) prep step options:

Aligment Options:

           --t2   INT        number of threads [1]

SNV Calling Options:

           --aq2  INT        skip alignments with mapQ smaller than INT [30]

SNV filtering Options:

           --cov  INT        average sequencing coverage [30] 
           --sq2  INT        skip SNVs with phred-scaled quality smaller than [30]
           --pn   INT        reads counts of minor allele should greater than INT [3]

3) haplotype step options:

           --phase2   use samtools algorithm or HAPCUT2 algorithm to assembly haplotype [default:T; T: samtools; F:HAPCUT2]
          
 SNP genotyping Options:

           --dep  INT   skip SNPs with read depth smaller than INT [10]
           --aq3  INT   skip alignment with mapQ smaller than INT  [20]
           --vq3  INT   skip SNPs with phred-scaled quality smaller than INT [40]    

4) afd step options:

Statistics Options:

      --sd  STR    the statistic method: ED/g/si [g]    
      --w   INT    the sliding window size [1000000] 
      --th  STR    Threshold selection method [N] 
                                           1) N: Estimate parameters of the log-normal null-distribution. 
                                           2) Y: medium plus 3-times standard deviation   
      --m   STR      Smooth curve using multi-window size. Y: YES; N: NO [Y] 


5) polish step options:

Statistics Options:

      --sd  STR    the statistic method: ED/g/si [g]    
      --w   INT    the sliding window size [1000000] 
      --th  STR   Threshold selection method [N] 
                                           1) N: Estimate parameters of the log-normal null-distribution. 
                                           2) Y: medium plus 3-times standard deviation   
      --m   STR      Smooth curve using multi-window size. Y: YES; N: NO [Y] 


6) qtl_pick options:

      --q   INT        mininum phred-scaled quality score [30]         
      --pr  INT        promoter region [2000]
 


Example:

1) Use reads files to run bastos
 
  bastos script --o result --r genome.fasta --gtf gene.gtf --pf1 P_1.fastq.gz --pf2 P_2.fastq.gz --mf1 M_1.fastq.gz --mf2 M_2.fastq.gz --hf1 H_1.fastq.gz --hf2 H_2.fastq.gz --lf1 L_1.fastq.gz --lf2 L_2.fastq.gz > script.sh  

2) Use pre-aligned BAMs files to run bastos

 bastos script --o result --r genome.fasta --gtf gene.gtf --pb P.bam --mb M.bam --hb H.bam --lb L.bam > script.sh 
    
";
 
 
sub runscript {
        my $genome = undef;
        my $gtf = undef;
        my $pf1 = undef;
	my $pf2 = undef;
        my $mf1 = undef;
        my $mf2 = undef;
        my $hf1 = undef;
        my $hf2 = undef;
        my $lf1 = undef;
        my $lf2 = undef;  
        my $pb  = undef;
        my $mb  = undef;
        my $hb  = undef;
        my $lb  = undef;
        my $t   = undef;
        my $aq  = undef;
        my $vq  = undef;
        my $d   = undef;
        my $sq  = undef;
        my $t2  = undef;
        my $aq2 = undef;
        my $cov = undef;
        my $sq2 = undef;
        my $pn  = undef;
        my $phase2 = undef;
        my $dep    = undef;
        my $aq3    = undef;
        my $vq3    = undef;
        my $sd      = undef;
        my $w      = undef;
        my $th     = undef;        
        my $m      = undef; 
        my $q      = undef;
        my $pr     = undef;
         


     



        my $outputPrefix = undef;
	my $help = 0;
	
	GetOptions (
	
        "r=s" =>\$genome,
        "gtf=s" =>\$gtf,
        "pf1=s" =>\$pf1,
        "pf2=s" =>\$pf2,
        "mf1=s" =>\$mf1,
        "mf2=s" =>\$mf2,
        "hf1=s" =>\$hf1,
        "hf2=s" =>\$hf2,
        "lf1=s" =>\$lf1,
        "lf2=s" =>\$lf2, 
        "pb=s"  =>\$pb,
        "mb=s"  =>\$mb,
        "hb=s"  =>\$hb,
        "lb=s"  =>\$lb,
        "t=i" =>\$t,
        "aq=i" =>\$aq,
        "vq=i" =>\$vq,
        "d=i" =>\$d, 
        "sq=i" =>\$sq,
        "t2=i" =>\$t2,
        "aq2=i" =>\$aq2,
        "cov=i" =>\$cov,
        "sq2=i" =>\$sq2,
        "pn=i" =>\$pn,
        "phase2=s" =>\$phase2,
        "dep=i" =>\$dep,
        "aq3=i" =>\$aq3,
        "vq3=i" =>\$vq3,
        "sd=s" =>\$sd, 
        "w=s" =>\$w,
        "th=s" =>\$th,
        "m=s" =>\$m,
        "q=i" =>\$q,
        "pr=i" =>\$pr,  
        "oprefix=s"   => \$outputPrefix,
	"help=s" => \$help)
         
	or die("$G_USAGE");
	
	die "$G_USAGE" if ($help || !defined($genome));

	
    unless(defined($outputPrefix)){
                   $outputPrefix="script";
                                  }

     unless(defined($t)){
                    $t=1;
                        }
     unless(defined($aq)){
                    $aq=30;
                        }
     unless(defined($vq)){
                   $vq=30;
                         }
     unless(defined($d)){
                   $d=10;
                        }
     unless(defined($sq)){
                   $sq=30;
                         }
     unless(defined($t2)){
                   $t2=1;
                         }
     unless(defined($aq2)){
                   $aq2=30;
                          }
     unless(defined($cov)){
                   $cov=30;
                          }
     unless(defined($sq2)){
                   $sq2=30;
                         }
     unless(defined($pn)){
                   $pn=3;
                         }
     unless(defined($phase2)){
                    $phase2="T";
                             }
     unless(defined($dep)){
                    $dep=10;
                          }
     unless(defined($aq3)){
                    $aq3=20;
                          }
     unless(defined($vq3)){
                    $vq3=40;
                          }
     unless(defined($sd)){
                    $sd="g";
                         }
     unless(defined($w)){
                    $w=1000000;
                        }
     unless(defined($th)){
                    $th="N";
                         }
     unless(defined($m)){
                    $m="Y";
                        }                            
     unless(defined($q)){
                    $q=30;
                        }
     unless(defined($pr)){
                    $pr=2000;
                         }

 

         my $bsatos ="$FindBin::Bin/bsatos";        
         
         #my $cd = $outputPrefix.".cm.sh";
         my $work_dir=getcwd;
         my $dir = $work_dir."/".$outputPrefix."_dir";

       
         unless(-e $dir){
           print "mkdir $dir\n";
                        }   
       
           
      
      #    open (my $ofh, ">$script") or die $!;
     
      
        if(defined($pf1) && defined ($pf2) && defined($mf1) && defined($mf2)){
        print  "$bsatos prepar --t $t --aq $aq --vq $vq --d $d --sq $sq  --pf1 $pf1  --pf2 $pf2 --mf1 $mf1 --mf2 $mf2 --r $genome --gtf $gtf  --o prepar\n"  
                                                                             }else{
        print  "$bsatos prepar --t $t --aq $aq --vq $vq --d $d --sq $sq   --pb $pb --mb $mb --r $genome  --o prepar --gtf $gtf\n"; 
                                                                                  } 
        if(defined($hf1) && defined($hf2) && defined($lf1) && defined($lf2)){           
        print  "$bsatos prep --t2 $t2 --aq2 $aq2 --cov $cov --sq2 $sq2 --pn $pn --hf1 $hf1 --hf2 $hf2 --lf1 $lf1  --lf2 $lf2 --r $genome --g1 prepar_dir/P_M_G1 --g2 prepar_dir/P_M_G2 --g3 prepar_dir/P_M_G3  --o prep\n";  
                                                                            }else{ 
        print  "$bsatos prep --t2 $t2 --aq2 $aq2 --cov $cov --sq2 $sq2 --pn $pn --hb $hb --lb $lb --r $genome  --g1 prepar_dir/P_M_G1 --g2 prepar_dir/P_M_G2 --g3 prepar_dir/P_M_G3  --o prep\n";                                                                      
                                                   
                               }
           
          unless(defined($pb)){
                        $pb="$FindBin::Bin/prepar_dir/prepar_P_rm.bam";                      
                              }
          unless(defined($mb)){
                        $mb="$FindBin::Bin/prepar_dir/prepar_M_rm.bam";                    
                              } 
          unless(defined($hb)){
                        $hb="$FindBin::Bin/prep_dir/prep_H_rm.bam";
                              }
           unless(defined($lb)){
                        $lb="$FindBin::Bin/prep_dir/prep_L_rm.bam";
                              }

        print "$bsatos haplotype --phase2 $phase2 --dep $dep --aq3 $aq3 --vq3 $vq3  --pb $pb --mb $mb --hb $hb --lb $lb  --r $genome  --o haplotype --var prepar_dir/M_P.snp\n";
     
        print "$bsatos afd --sd $sd --w $w --th $th --m $m --o  afd  --g1 prep_dir/g1.res  --g2 prep_dir/g2.res  --g3 prep_dir/g3.res --h haplotype_dir/merged.block\n";

        print "$bsatos polish --sd $sd --w $w --th $th --m $m --o polish --g1 afd_dir/g1.res.ad --g2 afd_dir/g2.res.ad --g3 afd_dir/g3.res.ad -h haplotype_dir/merged.block\n";
  
        print "$bsatos qtl_pick --q $q --pr $pr --o qtl_pick --g1 polish_dir/g1.res.ad.ad --g2 polish_dir/g2.res.ad.ad --g3 polish_dir/g3.res.ad.ad --v prepar_dir/anno/snv.AT_multianno.txt --sv prepar_dir/anno/sv.AT_multianno.txt -gtf $gtf --h haplotype_dir/merged.block\n"; 
         
        print "mv prepar_dir prep_dir  haplotype_dir afd_dir polish_dir qtl_pick_dir *sh $dir\n";                                                                                                                                    
 #         &process_cmd("cat $script");
           exit(0);
 
              }


    sub process_cmd {
    my ($cmd) = @_;
   # print STDERR "CMD: $cmd\n";
    my $ret = system($cmd);
    if ($ret) {
        die "Error, cmd: $cmd died with ret $ret";
    }

    return;
                  }

 
  
