package all;

BEGIN {
	$VERSION = '1.0';
}


use strict;
use Getopt::Long;
use base 'Exporter';
use File::Basename;
use FindBin qw($Bin);
use Cwd; 


our @EXPORT = qw(runall);

my $G_USAGE = "
 
Usage: bsatos all [options]  

  Options: 
         --oprefix STR       the prefix of the result folder [all] 
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
         --log1    FILE      prepar module log file [prepar.log] 
         --log2    FILE      prep module log file [prep.log]
         --log3    FILE      haplotype log file [haplotype.log]
         --log4    FILE      afd module log [afd.log]
         --log5    FILE      polish module log file [polish.log]
         --log6    FILE      qtl_pick module log file [qtl_pick.log]
         --log     FILE      script module log [script.log]

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
      --fn  INT    batches for smoothing; the smaller the faster but more memory[20]

5) polish step options:

Statistics Options:

      --sd  STR    the statistic method: ED/g/si [g]    
      --w   INT    the sliding window size [1000000] 
      --fn  INT    batches for smoothing; the smaller the faster but more memory[20]

6) qtl_pick options:

      --q   INT        mininum phred-scaled quality score [30]         
      --pr  INT        promoter region [2000]
      


Example:

1) Use reads files to run bastos
 
  bastos all --o result --r genome.fasta --gtf gene.gtf --pf1 P_1.fastq.gz --pf2 P_2.fastq.gz --mf1 M_1.fastq.gz --mf2 M_2.fastq.gz --hf1 H_1.fastq.gz --hf2 H_2.fastq.gz --lf1 L_1.fastq.gz --lf2 L_2.fastq.gz 

2) Use pre-aligned BAMs files to run bastos

 bastos all --o result --r genome.fasta --gtf gene.gtf --pb P.bam --mb M.bam --hb H.bam --lb L.bam  
    
";
 
 
sub runall {
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
        my $fn  = undef;
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
        my $log1   = undef;
        my $log2   = undef;
        my $log3   = undef;
        my $log4   = undef;
        my $log5   = undef;
        my $log6   = undef; 
        my $log   = undef; 
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
        "fn=i" =>\$fn,
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
        "log1=s" => \$log1,
        "log2=s" =>\$log2,
        "log3=s" =>\$log3,
        "log4=s" =>\$log4,
        "log5=s" =>\$log5,
        "log6=s" =>\$log6,
        "log=s" =>\$log,
        "oprefix=s"   => \$outputPrefix,
	"help=s" => \$help)
         
	or die("$G_USAGE");
	
	die "$G_USAGE" if ($help || !defined($genome));

	
        unless(defined($outputPrefix)){
                   $outputPrefix="all";
                                      }


#################################################################################################

         my $work_dir=getcwd;
         my $dir = $work_dir."/".$outputPrefix."_dir";
         my $prepar_dir=$work_dir."/prepar_dir";
         my $prep_dir=$work_dir."/prep_dir";
         my $afd_dir=$work_dir."/AFD_dir";
         my $hap_dir=$work_dir."/haplotype_dir";
         my $qtl_dir=$work_dir."/qtl_pick_dir";      
         my $po_dir=$work_dir."/polish_dir";


        unless(defined($log)){                           
         $log=$work_dir."/".$outputPrefix.".log";
                             }

         my %hlog=();       

         if(-e $log){
           open LOD,"<$log"; 
           my $con=0;    

                   while(<LOD>){
                    chomp;
                    my @sp=split/\t/,$_;
                    if($sp[2] eq "\#done"){            
                    $hlog{$sp[0]}=1;
                                          }
                    if($sp[2] eq "\#die"){
                    $hlog{$sp[0]}=0;                           
                                         }
                    if($_=~/finish/){
                             $con++;       
                                    }                                                                     
                                } 
             if($con !=0 ){      
                print "
         
          
             Program: bsatos (Bulked segregant analysis tools for outbreeding species)
             Version: 1.0
             
               BSATOS has be finished!!!!!\n\n\n\n\n";
                 exit(0);
                             
                        }
                  }
  
#####################################################################################
          
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
                    $phase2="NO";
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



        unless(defined($fn)){
  
                       $fn=20;
                                     
                            }


############################################################################
############################################################################






    my $bsatos ="$FindBin::Bin/bsatos";               
    my $work_dir=getcwd;
    my $dir = $work_dir."/".$outputPrefix."_dir";

    open (my $lo,">>$log") or die $!;  
 
         unless(-e $dir){

         #print "mkdir $dir\n";

          print STDOUT "create igv result dir\n"; 
          print $lo "create igv result dir\nCMD1:\t";
          unless(defined($hlog{"CMD0:"})){
                         $hlog{"CMD0:"}=0;
                                         }
          my $mx1c = $hlog{"CMD0:"}; 
          my $cd1c = "mkdir $dir";
          &process_cmd($cd1c,$lo,$mx1c);

                        }          
       
           print STDOUT "

              \#\#\#\#\#Program: bsatos (Bulked segregant analysis tools for outbreeding species)
              \#\#\#\#\#Version: 1.0
              \#\#\#\#\#
              \#\#\#\#\#all         do prepar, prep, haplotype, afd, polish, qtl_pick, igv in turn\n\n\n\n\n\n";


           print $lo "

              \#\#\#\#\#Program: bsatos (Bulked segregant analysis tools for outbreeding species)
              \#\#\#\#\#Version: 1.0
              \#\#\#\#\#
              \#\#\#\#\#all         do prepar, prep, haplotype, afd, polish, qtl_pick, igv in turn\n\n\n\n\n\n";



########################################################################################################

       print STDOUT "\#\#\#\#\#The first module
              \#\#\#\#\#Program: bsatos (Bulked segregant analysis tools for outbreeding species)
              \#\#\#\#\#Version: 1.0
              \#\#\#\#\#
              \#\#\#\#\#prepar      prepare the parents data:\n";

         
       print $lo "\#\#\#\#\#The first module
              \#\#\#\#\#Program: bsatos (Bulked segregant analysis tools for outbreeding species)
              \#\#\#\#\#Version: 1.0
              \#\#\#\#\#
              \#\#\#\#\#prepar      prepare the parents data:\n";


#############################################################################################################

        
   
###########################################################################################################
         print STDOUT "stage1:\t";
         print $lo  "stage1:\t";


        if(defined($pf1) && defined ($pf2) && defined($mf1) && defined($mf2)){

        my $stage1 = "$bsatos prepar --t $t --aq $aq --vq $vq --d $d --sq $sq  --pf1 $pf1  --pf2 $pf2 --mf1 $mf1 --mf2 $mf2 --r $genome --gtf $gtf  --o prepar";
          
       unless(defined($hlog{"stage1"})){
                $hlog{"stage1"}=0;
                                         }
        my $m1=$hlog{"stage1"};
        
        &process_cmd($stage1,$lo,$m1);



                                                                             }else{
        my $stage1 = "$bsatos prepar --t $t --aq $aq --vq $vq --d $d --sq $sq   --pb $pb --mb $mb --r $genome  --o prepar --gtf $gtf";
          
        unless(defined($hlog{"stage1"})){
 
                $hlog{"stage1"}=0;
                                         }
               my $m1=$hlog{"stage1"};

              &process_cmd($stage1,$lo,$m1);

                            }

##################################################################################################################
 
 print STDOUT "\#\#\#\#\#The second module
              \#\#\#\#\#Program: bsatos (Bulked segregant analysis tools for outbreeding species)
              \#\#\#\#\#Version: 1.0
              \#\#\#\#\#
              \#\#\#\#\#prep      prepare the pool data:\n";

    print $lo "\#\#\#\#\#The second module
              \#\#\#\#\#Program: bsatos (Bulked segregant analysis tools for outbreeding species)
              \#\#\#\#\#Version: 1.0
              \#\#\#\#\#
              \#\#\#\#\#prep      prepare the pool data:\n";


       print STDOUT "stage2:\t";
       print $lo  "stage2:\t";




#####################################################################################################################

        if(defined($hf1) && defined($hf2) && defined($lf1) && defined($lf2)){


              my $stage2 = "$bsatos prep --t2 $t2 --aq2 $aq2 --cov $cov --sq2 $sq2 --pn $pn --hf1 $hf1 --hf2 $hf2 --lf1 $lf1  --lf2 $lf2 --r $genome --gp prepar_dir/P_M.m --gm prepar_dir/P_M.m --gpm prepar_dir/P_M.pm  --o prep";

              unless(defined($hlog{"stage2"})){

              $hlog{"stage2"}=0;
                                         }
              my $m2=$hlog{"stage2"};

              &process_cmd($stage2,$lo,$m2);

                     }else{

         my $stage2 = "$bsatos prep --t2 $t2 --aq2 $aq2 --cov $cov --sq2 $sq2 --pn $pn --hb $hb --lb $lb --r $genome  --gp prepar_dir/P_M.p --gm prepar_dir/P_M.m --gpm prepar_dir/P_M.pm  --o prep";
                       
             unless(defined($hlog{"stage2"})){

             $hlog{"stage2"}=0;
                                         }
             my $m2=$hlog{"stage2"};

              &process_cmd($stage2,$lo,$m2);
                                                                               
                         }

      
########################################################################################################################

print STDOUT "\#\#\#\#\#The third module
              \#\#\#\#\#Program: bsatos (Bulked segregant analysis tools for outbreeding species)
              \#\#\#\#\#Version: 1.0
              \#\#\#\#\#
              \#\#\#\#\#haplotype   construct haplotype block:\n";



print $lo    "\#\#\#\#\#The third module
              \#\#\#\#\#Program: bsatos (Bulked segregant analysis tools for outbreeding species)
              \#\#\#\#\#Version: 1.0
              \#\#\#\#\#
              \#\#\#\#\#haplotype   construct haplotype block:\n";


       print STDOUT "stage3:\t";
       print $lo  "stage3:\t";


########################################################################################################################

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


        my $stage3 = "$bsatos haplotype --hapcut2 $phase2 --dep $dep --aq3 $aq3 --vq3 $vq3  --pb $pb --mb $mb --hb $hb --lb $lb  --r $genome  --o haplotype --var prepar_dir/P_M.snv";

        unless(defined($hlog{"stage3"})){

             $hlog{"stage3"}=0;
                                         }
        my $m3=$hlog{"stage3"};

        &process_cmd($stage3,$lo,$m3);



###########################################################################################################################
 print STDOUT "\#\#\#\#\#The fourth module
              \#\#\#\#\#Program: bsatos (Bulked segregant analysis tools for outbreeding species)
              \#\#\#\#\#Version: 1.0
              \#\#\#\#\#
              \#\#\#\#\#afd         calculate and filter allele frequency difference between two extreme pools:\n";

print  $lo   "\#\#\#\#\#The fourth module
              \#\#\#\#\#Program: bsatos (Bulked segregant analysis tools for outbreeding species)
              \#\#\#\#\#Version: 1.0
              \#\#\#\#\#
              \#\#\#\#\#afd         calculate and filter allele frequency difference between two extreme pools:\n";


       print STDOUT "stage4:\t";
       print $lo  "stage4:\t";


###########################################################################################################################


      my $stage4 = "$bsatos afd --sd $sd --w $w --th $th --m $m --o  AFD  --pc prep_dir/P.counts  --mc prep_dir/M.counts  --pmc prep_dir/PM.counts --hap haplotype_dir/haplotype.block";

     unless(defined($hlog{"stage4"})){

             $hlog{"stage4"}=0;
                                         }
        my $m4=$hlog{"stage4"};

      &process_cmd($stage4,$lo,$m4);


###########################################################################################################################
          

 print STDOUT "\#\#\#\#\#The fifth module
              \#\#\#\#\#Program: bsatos (Bulked segregant analysis tools for outbreeding species)
              \#\#\#\#\#Version: 1.0
              \#\#\#\#\#
              \#\#\#\#\#polish      polish candidate QTLs region and remove nosiy makers based on haplotype information:\n";

 print $lo   "\#\#\#\#\#The fifth module
              \#\#\#\#\#Program: bsatos (Bulked segregant analysis tools for outbreeding species)
              \#\#\#\#\#Version: 1.0
              \#\#\#\#\#
              \#\#\#\#\#polish      polish candidate QTLs region and remove nosiy makers based on haplotype information:\n";

 print STDOUT "stage5:\t";
 print $lo  "stage5:\t";

 

##############################################################################################

        my $stage5 = "$bsatos polish --sd $sd --w $w --th $th --m $m --o polish --gs1 AFD_dir/P.AFD --gs2 AFD_dir/M.AFD --gs3 AFD_dir/PM.AFD -h haplotype_dir/haplotype.block";

        unless(defined($hlog{"stage5"})){

             $hlog{"stage5"}=0;
                                         }
        my $m5=$hlog{"stage5"};
 
       &process_cmd($stage5,$lo,$m5);

  
#################################################################################################################################
 print STDOUT "\#\#\#\#\#The sixth module
              \#\#\#\#\#Program: bsatos (Bulked segregant analysis tools for outbreeding species)
              \#\#\#\#\#Version: 1.0
              \#\#\#\#\#
              \#\#\#\#\#qtl_pick    judge and pick up QTLs from three types of peaks:\n";



 print $lo    "\#\#\#\#\#The sixth module
              \#\#\#\#\#Program: bsatos (Bulked segregant analysis tools for outbreeding species)
              \#\#\#\#\#Version: 1.0
              \#\#\#\#\#
              \#\#\#\#\#qtl_pick    judge and pick up QTLs from three types of peaks:\n";


 print STDOUT "stage6:\t";
 print $lo  "stage6:\t";


#########################################################################################################

        

###############################################################################################################################

        my $stage6 = "$bsatos qtl_pick --win $w --q $q --pr $pr --o qtl_pick --gp polish_dir/P.polished.afd --gm polish_dir/M.polished.afd --gpm polish_dir/PM.polished.afd --v prepar_dir/anno/snv.AT_multianno.txt --sv prepar_dir/anno/sv.AT_multianno.txt -gtf $gtf --h haplotype_dir/haplotype.block";

         
        unless(defined($hlog{"stage6"})){

             $hlog{"stage6"}=0;
                                        }
        my $m6=$hlog{"stage6"};

       &process_cmd($stage6,$lo,$m6);


###################################################################################################
##################################################################################################
 print $lo    "\#\#\#\#\#The last module
              \#\#\#\#\#Program: bsatos (Bulked segregant analysis tools for outbreeding species)
              \#\#\#\#\#Version: 1.0
              \#\#\#\#\#
              \#\#\#\#\#igv         generate files for Integrative Genomics Viewer\n";


 print STDOUT "stage7:\t";
 print $lo  "stage7:\t";



      my $stage7 = "$bsatos igv --o igv --r $genome --gtf $gtf --prepar $prepar_dir --hap $hap_dir --qtl $qtl_dir --po $po_dir";


        unless(defined($hlog{"stage7"})){

             $hlog{"stage7"}=0;
                                        }
        my $m7=$hlog{"stage7"};

       &process_cmd($stage7,$lo,$m7);




################################################################################################################################
#################################################################################################################################
       print STDOUT  "\#\#\#\#\#\#\#\#\#\#\#\#\#move all the data to result dir\n";
       print $lo    "\#\#\#\#\#\#\#\#\#\#\#\#\#move all the data to result dir\nstage8\t";

       my $stage8 = "mv  $work_dir/igv_dir $work_dir/prepar_dir $work_dir/prep_dir  $work_dir/haplotype_dir $work_dir/AFD_dir $work_dir/polish_dir $work_dir/qtl_pick_dir $dir";
  
       
        unless(defined($hlog{"stage8"})){

             $hlog{"stage8"}=0;
                                         }
        my $m8=$hlog{"stage8"};

       &process_cmd($stage8,$lo,$m8);
                                                                                                                                   
  
       print STDOUT  "\#\#\#\#\#\#\#\#\#\#\#\#\#BSATOS finished!!!!\n";
       print $lo    "\#\#\#\#\#\#\#\#\#\#\#\#\#BSATOS finished!!!\n";





       exit(0);

 
              }


########################################################################################

  sub process_cmd {
    my ($cmd,$lk,$hg) = @_;
    my $time="[".localtime()."]";
    if($hg == 1){
   

    print STDOUT "this command has be processed last time\n$cmd\t.........................................................\#done\t$time\n"; 
    print $lk "$cmd\t\#done\t$time\n";
                     }else{

    print STDOUT "$cmd\t........................................................."; 
    print $lk "$cmd\t";

    my $ret = system($cmd);
    if ($ret) {
           print $lk "\#die\t$time\n";  
           print STDOUT "\#die\t$time\n";       
           die "Error, cmd: $cmd died with ret $ret";
              }else{
        print $lk  "\#done\t$time\n";
        print STDOUT "\#done\t$time\n";
                   }
    return;
                  }                  
                  } 

##############################################################













##
##      BSATOS: Bulk Segregant Analysis Tools for Outbreeding Species
##
##      author: Fei Shen <352363546@qq.com>  
##
## Copyright (c) 2019 China Agricultural University 
##
## Permission is hereby granted, free of charge, to any person obtaining a copy
## of this software and associated documentation files (the "Software"), to deal
## in the Software without restriction, including without limitation the rights
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
## copies of the Software, and to permit persons to whom the Software is
## furnished to do so, subject to the following conditions:
## 
## The above copyright notice and this permission notice shall be included in all
## copies or substantial portions of the Software.
## 
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
## IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
## AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
## LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
## OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
## SOFTWARE.
