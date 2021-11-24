

package haplotype;

BEGIN {
	$VERSION = '1.0';
}


use strict;
use Getopt::Long;
use File::Basename;
use Cwd;
use base 'Exporter';
our @EXPORT = qw(runhaplotype);

my $G_USAGE = "

Usage: bastos haplotype 

Options:   
           --pb         FILE     pre-aligned bam file from pollen parent 
           --mb         FILE     pre-aligned bam file from maternal parent
           --hb         FILE     pre-aligned bam file from high extreme pool
           --lb         FILE     pre-aligned bam file from low extreme pool
           --r          FILE     reference genome [FASTA format]
           --o          STR      output dir name prefix 
           --log        FILE     haplotype log file [haplotype.log]
           --hapcut2    STF      use samtools algorithm or HAPCUT2 algorithm to assembly haplotype [NO]

SNP genotyping Options:

           --dep  INT   skip SNPs with read depth smaller than INT [10]
           --aq3  INT   skip alignment with mapQ smaller than INT  [20]
           --vq3  INT   skip SNPs with phred-scaled quality smaller than INT [40]    

Outputs:

  P_block            [FILE]      haplotype blocks of pollen parent
  M.block            [FILE]      haplotype blocks of maternal parent 
  haplotype.block    [FILE]      merged, corrected and patched haplotype blocks of pollen parent, maternal parent, High extreme pool and Low extreme pool
  P_haplotype.bed    [FILE]      BED format haplotype information of pollen parent
  M_haplotype.bed    [FILE]      BED format haplotype information of maternal parent
  overlapped.bed     [FILE]      BED format haplotype information classified from two parents and two pools
  sub_haplotype      [FILE]      haplotype sub-blocks information within haplotype block 


Example:

    bsatos haplotype --pb P.bam --mb M.bam --hb H.bam --lb L.bam --r genome.fasta --o haplotype  

";

sub runhaplotype {
	my $Pbam = undef;
 	my $Mbam = undef;
        my $Hbam = undef;
        my $Lbam = undef;
        my $genome = undef;
	my $outputPrefix = undef;
	my $help = 0;
	my $phase2 = undef;
        my $var = undef;
        my $var1 =undef;
        my $s = undef;
        my $dep = undef;
        my $aq = undef;
        my $vq = undef;
        my $log = undef;
	GetOptions (
	"pb=s" => \$Pbam,
	"mb=s" => \$Mbam,
        "hb=s" => \$Hbam,
        "lb=s" => \$Lbam,
        "hapcut2=s" =>\$phase2,
        "r=s" => \$genome,
        "var=s" => \$var, 
	"oprefix=s"   => \$outputPrefix,
        "s=s" => \$s,
        "dep=i" =>\$dep,
        "aq3=i" =>\$aq,
        "vq3=i" =>\$vq,
        "log=s" =>\$log,
	"help!" => \$help)
	or die("$G_USAGE");
	
	die "$G_USAGE" if ($help);
	
	die "$G_USAGE" if ((!defined ($outputPrefix)) || (!defined($genome)));


########################################################


        unless(defined($outputPrefix)){
                      $outputPrefix="haplotype";
                                     }
   
        unless(defined($phase2)){         
                        $phase2="NO";
                                }

        unless(defined($aq)){
                       $aq=20;
                             }
        unless(defined($vq)){
                       $vq=40;
                            }

       unless(defined($dep)){
                      $dep=10;
                            }
 


         my $samtools= "samtools";
         my $bcftools= "bcftools";
         my $bwa="bwa";
         my $filter="$FindBin::Bin/scripts/filter";
         my $gg="$FindBin::Bin/scripts/filter_get_genotype.pl";
         my $step1="extractHAIRS";
         my $step2="HAPCUT2";
         my $get_ref="$FindBin::Bin/scripts/get_refer.pl";
         my $get_block="$FindBin::Bin/scripts/get_block.pl";
         my $get_b="$FindBin::Bin/scripts/get_hap.pl";
         my $get_list="$FindBin::Bin/scripts/merge_haplotype.pl";
         my $sort="$FindBin::Bin/scripts/judge_and_sort_haplotype.pl";
         my $fill_homo="$FindBin::Bin/scripts/fill_homo.pl";
         my $cmp1="$FindBin::Bin/scripts/cmp_infer_block_step1.pl";
         my $cmp2="$FindBin::Bin/scripts/cmp_infer_block_step2.pl";
         my $sep= "$FindBin::Bin/scripts/classify_haplotype.pl";
         my $bedops = "bedops";
         my $filter_vcf="$FindBin::Bin/scripts/filter_vcf.pl";
         my $add_header="$FindBin::Bin/scripts/add_header.pl"; 
         my $work_dir=getcwd;
         my $dir = $work_dir."/".$outputPrefix."_dir";
         my $add_bl_header="$FindBin::Bin/scripts/add_bl_header.pl";            
         my $ch2bed="$FindBin::Bin/scripts/ch2bed.pl";     
         my $P_phase=$dir."/".basename($Pbam)."_phase";
         my $M_phase=$dir."/".basename($Mbam)."_phase";
         my $H_phase=$dir."/".basename($Hbam)."_phase";
         my $L_phase=$dir."/".basename($Lbam)."_phase";
         my $P_ref =$dir."/".basename($Pbam)."_ref";
         my $M_ref =$dir."/".basename($Mbam)."_ref";
         my $H_ref =$dir."/".basename($Hbam)."_ref";
         my $L_ref= $dir."/".basename($Lbam)."_ref";
         my $P_bl =$dir."/"."P_block";
         my $M_bl =$dir."/"."M_block";
         my $H_bl =$dir."/"."H_block";
         my $L_bl= $dir."/"."L_block";
         my $P_bl1 =$dir."/".basename($Pbam)."_block1";
         my $M_bl1 =$dir."/".basename($Mbam)."_block1";
         my $H_bl1 =$dir."/".basename($Hbam)."_block1";
         my $L_bl1= $dir."/".basename($Lbam)."_block1";
         my $merge =$dir."/"."merge.list";
         my $merge_bl=$dir."/"."merged.block";
         my $phase_res=$dir."/"."haplotype.block";
         my $P_out = $dir."/"."P.haplotype";
         my $P_srt = $dir."/"."P_haplotype.bed";
         my $M_out = $dir."/"."M.haplotype";
         my $M_srt = $dir."/"."M_haplotype.bed";
         my $overlap = $dir."/"."overlapped.bed";
         my $phase_sf =$dir."/"."phase_tmp";
      
        unless(defined($log)){
                           
         $log=$work_dir."/".$outputPrefix.".log";
                             
                             }


#######################################################

       my %hlog=();       
 
#########################

    if(-e $log){
        
           open LOD,"<$log"; 
           my $con=0;   

#####################################
                                 
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
                      
#########################################

            if($con !=0){                 
                 print "
             Program: bsatos (Bulked segregant analysis tools for outbreeding species)
             Version: 1.0
             The third module: haplotype
             haplotype   construct haplotype block\n\n\n\n 


             This module has be finished!!!!!\n\n\n";
                      
               exit(0);
                    
                        }
###########################################
                 }


###############################################################
   
         open (my $lo, ">>$log") or die $!;


         print STDOUT "

         Program: bsatos (Bulked segregant analysis tools for outbreeding species)
         Version: 1.0
         The third module: haplotype
         haplotype   construct haplotype block\n\n\n\n";

         print $lo "

         Program: bsatos (Bulked segregant analysis tools for outbreeding species)
         Version: 1.0
         The third module: haplotype
         haplotype   construct haplotype block\n\n\n\n";

      
######################################################

          unless(-e $dir){

          print STDOUT "create haplotype result dir\n";
          print $lo "create haplotype result dir\nCMD1:\n";
          unless(defined($hlog{"CMD1:"})){                
                         $hlog{"CMD1:"}=0;
                                         }
          my $mx1=$hlog{"CMD1:"};
          my $cd1 = "mkdir $dir";
          &process_cmd($cd1,$lo,$mx1);
  
                         }

########################################

  if($phase2 eq "NO" ){ 
      
##############################################
            
              unless(defined($var)){
                           $var=$dir."/"."var";

         print STDOUT "SNV calling based on BAM FILES\n";
         print $lo "SNV calling based on BAM FILES\nCMD2:\t";
         unless(defined($hlog{"CMD2:"})){                
                         $hlog{"CMD2:"}=0;
                                         }
         my $mx2=$hlog{"CMD2:"};
         my $cd2 = "$samtools mpileup -ugf $genome -q $aq $Pbam $Mbam $Hbam $Lbam |$bcftools call -mv -|$filter_vcf -d $dep -q $vq ->$var";
         &process_cmd($cd2,$lo,$mx2);
                                   }

#################################################
            
                    $var1=$dir."/"."var1";    

####################################

         print STDOUT "filtering SNVs based on reads depth and quality\n";
         print $lo "filtering SNVs based on reads depth and quality\nCMD3:\t";    
         unless(defined($hlog{"CMD3:"})){                
                         $hlog{"CMD3:"}=0;
                                         }
         my $mx3=$hlog{"CMD3:"};
         my $cd3= "$filter_vcf -d $dep -q $vq  $var > $var1";
         &process_cmd($cd3,$lo,$mx3);    
         
#######################################

         print STDOUT "samtools phase pollen parent\n";
         print $lo "samtools phase pollen parent\nCMD4:\t";
         unless(defined($hlog{"CMD4:"})){                
                         $hlog{"CMD4:"}=0;
                                        }
         my $mx4=$hlog{"CMD4:"};
         my $cd4 = "$samtools phase $Pbam >$P_phase";
         &process_cmd($cd4,$lo,$mx4);

#####################################

         print STDOUT "samtools phase maternal parent\n";
         print $lo "samtools phase maternal parent\nCMD5:\t";
         unless(defined($hlog{"CMD5:"})){                
                         $hlog{"CMD5:"}=0;
                                         }  
         my $mx5=$hlog{"CMD5:"};
         my $cd5 ="$samtools phase $Mbam >$M_phase";
         &process_cmd($cd5,$lo,$mx5);

#########################################

         print STDOUT "samtools phase High pool\n";
         print $lo "samtools phase High pool\nCMD6:\t";
         unless(defined($hlog{"CMD6:"})){                
                         $hlog{"CMD6:"}=0;
                                         } 
         my $mx6=$hlog{"CMD6:"};
         my $cd6 = "$samtools phase $Hbam >$H_phase";
         &process_cmd($cd6,$lo,$mx6);


############################

         print STDOUT "samtools phase Low pool\n";
         print $lo "samtools phase low pool\nCMD7:\t";
         unless(defined($hlog{"CMD7:"})){                
                         $hlog{"CMD7:"}=0;
                                      }
         my $mx7=$hlog{"CMD7:"};
         my $cd7 ="$samtools phase $Lbam >$L_phase";
         &process_cmd($cd7,$lo,$mx7);

#################################

         print STDOUT "get pollen parent reference genome loci\n";
         print $lo "get pollen parent reference genome loci\nCMD8:\t";
         my $cd8 = "$get_ref $P_phase $genome >$P_ref";
         unless(defined($hlog{"CMD8:"})){                
                         $hlog{"CMD8:"}=0;
                                        }
         my $mx8=$hlog{"CMD8:"};
         &process_cmd($cd8,$lo,$mx8);

##########################

         print STDOUT "get maternal parent reference genome loci\n";
         print $lo "get maternal parent reference genome loci\nCMD9:\t";
         unless(defined($hlog{"CMD9:"})){                
                         $hlog{"CMD9:"}=0;
                                         }  
         my $mx9=$hlog{"CMD9:"};
         my $cd9 = "$get_ref  $M_phase $genome >$M_ref";
         &process_cmd($cd9,$lo,$mx9);

#################################
       
         print STDOUT "get H pool reference genome loci\n";
         print $lo "get H pool reference genome loci\nCMD10:\t";
         unless(defined($hlog{"CMD10:"})){                
                         $hlog{"CMD10:"}=0;
                                         } 
         my $mx10=$hlog{"CMD10:"};
         my $cd10 = "$get_ref  $H_phase $genome >$H_ref";
         &process_cmd($cd10,$lo,$mx10);


###########################

         print STDOUT "get L pool reference genome loci\n";
         print $lo "get L pool reference genome loci\nCMD11:\t";
         my $cd11 = "$get_ref  $L_phase $genome >$L_ref";
         unless(defined($hlog{"CMD11:"})){                
                         $hlog{"CMD11:"}=0;
                                         }
         my $mx11=$hlog{"CMD11:"};
         &process_cmd($cd11,$lo,$mx11); 
 
###############################

         print STDOUT "construct pollen parent haplotype blocks\n";
         print $lo "construct pollen parent haplotype block\nCMD12:\t";             
         my $cd12 = "$get_block  $P_phase $P_ref | $sort -> $P_bl";
         unless(defined($hlog{"CMD12:"})){                
                         $hlog{"CMD12:"}=0;
                                         }
         my $mx12=$hlog{"CMD12:"};
         &process_cmd($cd12,$lo,$mx12); 


##########################

         print STDOUT "construct maternal parent haplotype blocks\n";
         print $lo "construct maternal parent haplotype blocks\nCMD13:\t";
         my $cd13 = "$get_block  $M_phase $M_ref | $sort -> $M_bl";
         unless(defined($hlog{"CMD13:"})){                
                         $hlog{"CMD13:"}=0;
                                         }
         my $mx13=$hlog{"CMD13:"};
         &process_cmd($cd13,$lo,$mx13);

############################
 
         print STDOUT "construct H pool haplotype blocks\n";
         print $lo "construct H pool haplotype blocks\nCMD14:\t";
         my $cd14 = "$get_block  $H_phase $H_ref | $sort -> $H_bl";
         unless(defined($hlog{"CMD14:"})){                
                         $hlog{"CMD14:"}=0;
                                         }  
         my $mx14=$hlog{"CMD14:"};
         &process_cmd($cd14,$lo,$mx14);


########################

         print STDOUT "construct L pool haplotype blocks\n";
         print $lo "contruct L pool haplotype blocks\nCMD15:\t";
         unless(defined($hlog{"CMD15:"})){                
                         $hlog{"CMD15:"}=0;
                                         }   
         my $mx15=$hlog{"CMD15:"};
         my $cd15 = "$get_block  $L_phase $L_ref | $sort -> $L_bl";
         &process_cmd($cd15,$lo,$mx15);

 #############################
 
                }else{


##########################    

          unless(defined($var)){

                           $var=$dir."/"."var";
         print STDOUT "get SNVs loci and filtering\n";
         print $lo "get SNVs loci and filtering\nCMD16:\t";
         unless(defined($hlog{"CMD16:"})){                
                         $hlog{"CMD16:"}=0;
                                         }
         my $mx16=$hlog{"CMD16:"};
         my $cd16= "$samtools mpileup -ugf $genome -q $aq $Pbam $Mbam $Hbam $Lbam |$bcftools call -mv - |$filter_vcf -d $dep -q $vq ->$var";
         &process_cmd($cd16,$lo,$mx16);
                
                               }

#####################################3

                             $var1=$dir."/"."var1";

#################################
         
         print STDOUT "filtering SNVs loci\n";
         print $lo "filtering SNVs loci\nCMD17:\t";
    
         unless(defined($hlog{"CMD17:"})){                
                         $hlog{"CMD17:"}=0;
                                         }   
         my $mx17=$hlog{"CMD17:"};
         my $cd17 = "$filter_vcf -d $dep -q $vq  $var > $var1";
         &process_cmd($cd17,$lo,$mx17);

######################          
         print STDOUT "extract pollen parent fragments\n";
         print $lo "extract pollen parent fragments\nCMD18:\t";
         unless(defined($hlog{"CMD18:"})){                
                         $hlog{"CMD18:"}=0;
                                         }     
         my $mx18=$hlog{"CMD18:"};
         my $cd18= "$step1 --bam $Pbam   --VCF $var1 --out $P_phase";
         &process_cmd($cd18,$lo,$mx18);

##########################3
         
         print STDOUT "extract maternal parent fragments\n"; 
         print $lo "extract maternal parent fragments\nCMD19:\t";  
         my $cd19 = "$step1 --bam $Mbam   --VCF $var1 --out $M_phase";
         unless(defined($hlog{"CMD19:"})){                
                         $hlog{"CMD19:"}=0;
                                         }  
         my $mx19=$hlog{"CMD19:"};
         &process_cmd($cd19,$lo,$mx19);        

##################### 
         print STDOUT "extract H pool fragments\n";
         print $lo "extract H pool fragmets\nCMD20:\t";
         my $cd20 = "$step1 --bam $Hbam   --VCF $var1 --out $H_phase";
         unless(defined($hlog{"CMD20:"})){                
                         $hlog{"CMD20:"}=0;
                                         }  
         my $mx20=$hlog{"CMD20:"};
         &process_cmd($cd20,$lo,$mx20);

#######################
        
         print STDOUT "extract L pool fragments\n";
         print $lo "extract L pool fragments\nCMD21:\t";
         my $cd21 = "$step1 --bam $Lbam   --VCF $var1 --out $L_phase";
         unless(defined($hlog{"CMD21:"})){                
                         $hlog{"CMD21:"}=0;
                                         }
           
         my $mx21=$hlog{"CMD21:"};
         &process_cmd($cd21,$lo,$mx21);
##########################33         
        
         print STDOUT "construct pollen parent haplotype blocks\n";
         print $lo "construct pollen parent haplotype blocks\nCMD22:\n";
         my $cd22 = "$step2 --fragments $P_phase  --VCF $var1 --output $P_bl1";
         unless(defined($hlog{"CMD22:"})){                
                         $hlog{"CMD22:"}=0;
                                         }
           
         my $mx22=$hlog{"CMD22:"};
         &process_cmd($cd22,$lo,$mx22);
##########################

         print STDOUT "construct maternal parent haplotype blocks\n";
         print $lo "construct maternal parent haplotype blocks\nCMD23:\n";
         my $cd23 = "$step2 --fragments $M_phase  --VCF $var1 --output $M_bl1";
         unless(defined($hlog{"CMD23:"})){                
                         $hlog{"CMD23:"}=0;
                                         }
         my $mx23=$hlog{"CMD23:"};
         &process_cmd($cd23,$lo,$mx23);
###################
             
         print STDOUT "construct H pool haplotype blocks\n";
         print $lo "construct H pool haplotype blocks\nCMD24:\t";
         my $cd24 = "$step2 --fragments $H_phase  --VCF $var1 --output $H_bl1";
         unless(defined($hlog{"CMD24:"})){                
                         $hlog{"CMD24:"}=0;
                                         } 
         my $mx24=$hlog{"CMD24:"};
         &process_cmd($cd24,$lo,$mx24);

#################
         print STDOUT "construct L pool haplotype blocks\n";
         print $lo "construct L pool haplotype blocks\nCMD25:\t";
         unless(defined($hlog{"CMD25:"})){                
                        $hlog{"CMD25:"}=0;
                                         }           
         my $mx25=$hlog{"CMD25:"};
         my $cd25 ="$step2 --fragments $L_phase  --VCF $var1 --output $L_bl1";
         &process_cmd($cd25,$lo,$mx25);

##################
         print STDOUT "sorting pollen parent haplotype blocks\n";
         print $lo "sorting pollen parent haplotype blocks\nCMD26:\t";  
         my $cd26 = "$get_b  $P_bl1 | $sort -> $P_bl";
         unless(defined($hlog{"CMD26:"})){                
                         $hlog{"CMD26:"}=0;
                                         }  
         my $mx26=$hlog{"CMD26:"};
         &process_cmd($cd26,$lo,$mx26);
################

         print STDOUT "sorting maternal parent haplotype blocks\n";
         print $lo "sorting maternal parent haplotype blocks\nCMD27:\t";
         unless(defined($hlog{"CMD27:"})){                
                         $hlog{"CMD27:"}=0;
                                         }  
         my $mx27=$hlog{"CMD27:"};         
         my $cd27 = "$get_b  $M_bl1 | $sort -> $M_bl";
         &process_cmd($cd27,$lo,$mx27); 

################ 
         print STDOUT "sorting H pool haplotype blocks\n";
         print $lo "sorting H pool haplotype blocks\nCMD28:\t";
         my $cd28 = "$get_b  $H_bl1 | $sort -> $H_bl";
         unless(defined($hlog{"CMD28:"})){                
                         $hlog{"CMD28:"}=0;
                                         }
         my $mx28=$hlog{"CMD28:"};
         &process_cmd($cd28,$lo,$mx28);


#################
 
         print STDOUT "sorting L pool haplotype blocks\n";
         print $lo "sorting L pool haplotype blocks\nCMD29:\t";
         my $cd29 = "$get_b  $L_bl1 | $sort -> $L_bl";
         unless(defined($hlog{"CMD29:"})){                
                         $hlog{"CMD29:"}=0;
                                         }           
         my $mx29=$hlog{"CMD29:"};
         &process_cmd($cd29,$lo,$mx29);

####################
#
                 }


######################                       
         print STDOUT "get all the haplotype loci\n";
         print $lo "get all the haplotype loci\nCMD30:\t";
         my $cd30 = "$get_list  --Pb $P_bl --Mb $M_bl --Hb $H_bl  --Lb $L_bl >$merge";
         unless(defined($hlog{"CMD30:"})){                
                         $hlog{"CMD30:"}=0;
                                         }
         my $mx30=$hlog{"CMD30:"};
         &process_cmd($cd30,$lo,$mx30);


####################  
         print STDOUT "merge all the haplptype blocks files\n";
         print $lo "merge all the haplotype blocks files\nCMD31:\t";
         unless(defined($hlog{"CMD31:"})){                
                         $hlog{"CMD31:"}=0;
                                         }  
         my $mx31=$hlog{"CMD31:"};
         my $cd31 = "$filter -k cn -A A -B E  -C E -D E -E E  $merge  $P_bl  $M_bl $H_bl $L_bl |cut -f  1-2,5-9,14-16,21-23,28-30 -> $merge_bl";
         &process_cmd($cd31,$lo,$mx31);
    
####################
         print STDOUT "stitch and compare the haplotype blocks\n";
         print $lo "stitch and compare the haplotype blocks\nCMD32:\t";
         unless(defined($hlog{"CMD32:"})){                
                         $hlog{"CMD32:"}=0;
                                         } 
         my $mx32=$hlog{"CMD32:"};
         my $cd32 = "$fill_homo $merge_bl $var1  |$cmp1 - |$cmp2 -|$cmp1 - |$cmp2 - |$cmp1 - |$add_header -> $phase_res";
         &process_cmd($cd32,$lo,$mx32);
    
#####################
         print STDOUT "separate the haplotype blocks from pollen and maternal parents\n";
         print $lo "separate the haplotype blocks from pollen and maternal parents\nCMD33:\t";        
         my $cd33 = "$sep $phase_res $P_out $M_out";
         unless(defined($hlog{"CMD33:"})){                
                         $hlog{"CMD33:"}=0;
                                         }         
         my $mx33=$hlog{"CMD33:"};
         &process_cmd($cd33,$lo,$mx33);

#######################
         print STDOUT "sorting and comparing the haplptype blocks from pollen and maternal parents\n";
         print $lo "sorting and comparing the haplotype blocks from pollen and maternal parents\nCMD34:\t";            
         my $cd34 = "sort -k 1,1 -k 2,2n $P_out > $P_srt\;sort -k 1,1 -k 2,2n $M_out  > $M_srt\;$bedops -i $P_srt $M_srt > $overlap"; 
         unless(defined($hlog{"CMD34:"})){                
                         $hlog{"CMD34:"}=0;
                                         }
           
         my $mx34=$hlog{"CMD34:"};
         &process_cmd($cd34,$lo,$mx34);   
########################################################
     
          my $tmp_1=$dir."/"."sub_haplotype";
         
        
          my $detect_cr="$FindBin::Bin/scripts/detect_crossover.pl";
         print STDOUT "detect crossovers\n";
         print $lo "detect crossovers\nCMD39:\t";
                   
         my $cd39 = "$ch2bed $phase_res > $phase_sf;bedtools intersect  -wao -a $overlap  -b $phase_sf  |perl $detect_cr ->$tmp_1"; 
         unless(defined($hlog{"CMD39:"})){
                         $hlog{"CMD39:"}=0;
                                          }
      

         my $mx39=$hlog{"CMD39:"};
         &process_cmd($cd39,$lo,$mx39);

          
         
################################################################

#######################
       if(-e $H_phase){
         print STDOUT "remove some files\n";
         print $lo "remove some files\nCMD35:\t";
         unless(defined($hlog{"CMD35:"})){         
                         $hlog{"CMD35:"}=0;
                                         }
         my $mx35=$hlog{"CMD35:"};               
         my $cd35 = "rm $var $var1 $P_out $M_out  $H_bl $L_bl  $H_phase $H_ref $L_phase $L_ref $M_phase $M_ref $merge $merge_bl $P_phase $P_ref $phase_sf";
         &process_cmd($cd35,$lo,$mx35);
        
                     }

         print STDOUT "haplotype module finished!!!\n";
         print $lo "haplotype module finished!!!\n";


         exit(0);

############################
                    
                     
                     }


#######################################
 

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
