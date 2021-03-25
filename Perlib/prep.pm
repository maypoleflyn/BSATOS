

package prep;

BEGIN {
	$VERSION = '1.0';
}


use strict;
use Getopt::Long;
use base 'Exporter';
use Cwd;
our @EXPORT = qw(runprep);

my $G_USAGE = "

Usage: bastos prep [options]
         
        
Options:  
          --hf1  FILE     paired-end1 fastq file from high extreme pool  [FASTQ format]
          --hf2  FILE     paired-end1 fastq file from high extreme pool  [FASTQ format]
          --lf1  FILE     paired-end1 fastq file from low extreme pool   [FASTQ format]
          --lf2  FILE     paired-end1 fastq file from low extreme pool   [FASTQ format]
          --hb   FILE     pre-aligned bam file from high extreme pool    [Not required, if --hf1 & --hf2] 
          --lb   FILE     pre-aligned bam file from low extreme pool     [Not required, if --lf1 & --lf2]
          --r    FILE     reference genome [FASTA format]
          --gp   FILE     P file from prepar step  
          --gm   FILE     M file from prepar step          
          --gpm  FILE     PM file from prepar step           
          --o    FILE     output dir name prefix [prep]
          --log  FILE     prep module log file [prep.log]


Aligment Options:

           --t2   INT        number of threads [1]

SNV Calling Options:

           --aq2  INT        skip alignments with mapQ smaller than INT [30]

SNV filtering Options:

           --cov  INT        average sequencing coverage [30] 
           --sq2  INT        skip SNVs with phred-scaled quality smaller than [30]
           --pn   INT        reads counts of minor allele should greater than INT [3]


Outputs:

 prefix_dir [DIR] 

     P.counts   [FILE]   read counts with different alleles from H & L pools in P  type loci 
     M.counts   [FILE]   read counts with different alleles from H & L pools in M  type loci 
     PM.counts  [FILE]   read counts with different alleles from H & L pools in PM type loci
     sum        [FILE]   summary of P.counts, M.counts and PM.counts file
    
 File format:


  Chromosome	Position     H_REF	H_ALT	L_REF	L_ALT		
   Chr01         11000       10           1      1        10
    .              .          .		  .      .        .
    .              .          .           .      .        .
    .              .          .           .      .        .


Chromosome: the chromosome of markers   
Position  : the positon in chromosome of markers
H_REF     : read counts with REF alleles in H pool
H_ALT     : read counts with ALT alleles in H pool 
L_REF     : read counts with REF alleles in L pool 
L_ALT     : read counts with ALT alleles in L pool

Example:

1)
    bsatos prep --hf1 H_1.fastq --hf2 H_2.fastq --lf1 L_1.fastq --lf2 L_2.fastq --r genome.fasta  --gp P_M.p --gm P_M.m --g3 P_M.pm --o prep

OR

2) 
    bastos prep --hb H.bam --lb L.bam --r genome.fasta --gp P_M.p --gm P_M.m --gpm P_M.pm --o prep 

";

sub runprep {
	my $Hp1 = undef;
 	my $Hp2 = undef;
        my $Lp1 = undef;
        my $Lp2 = undef;
        my $Hbam = undef;
        my $Lbam = undef;
        my $G1  = undef;
        my $G2  = undef;
        my $G3 = undef;
        my $genome = undef;
	my $outputPrefix = undef;
	my $help = 0;
	my $s = undef;
        my $t = undef;
        my $aq = undef;
        my $cov = undef;
        my $sq = undef;
        my $pn = undef; 
        my $log = undef;
 
	GetOptions (
	"hf1=s" => \$Hp1,
	"hf2=s" => \$Hp2,
        "lf1=s" => \$Lp1,
        "lf2=s" => \$Lp2,
        "hb=s" => \$Hbam,
        "lb=s" => \$Lbam,
        "r=s" => \$genome,
        "gp=s"  => \$G1,
        "gm=s"  => \$G2, 
        "gpm=s"  => \$G3,
	"o=s"   => \$outputPrefix,
        "s=s" => \$s,
        "t2=i" => \$t,
        "aq2=i" =>\$aq,
        "cov=i" =>\$cov,
        "sq2=i" =>\$sq,
        "pn=i" =>\$pn,
        "log=s" =>\$log,
	"help!" => \$help)
        
	or die("$G_USAGE");
	
	die "$G_USAGE" if ($help);
	
	die "$G_USAGE" if (!defined $outputPrefix);

##################################################################

      unless(defined($outputPrefix)){
                        $outputPrefix="prep";
                                    }

        unless(defined($t)){
                           $t=1;
                           }
        unless(defined($aq)){
                           $aq=30;
                            }
        unless(defined($sq)){
                           $sq=30;
                            }
        unless(defined($cov)){
                           $cov=30;
                             }
        unless(defined($pn)){
                           $pn=3;
                            }
      

#####################################################################




         my $samtools= "samtools";
         my $bcftools="bcftools";
         my $bwa="bwa";
         my $get_counts="$FindBin::Bin/scripts/get_counts.pl";
         my $get_res="$FindBin::Bin/scripts/get_counts.pl";
         my $filter_dep="$FindBin::Bin/scripts/filter.pl";

         my $work_dir=getcwd;
         my $dir = $work_dir."/".$outputPrefix."_dir";
         my $H_bam=$dir."/".$outputPrefix."_H.bam";
         my $H_srt_bam=$dir."/".$outputPrefix."_H_srt.bam";
         my $H_rm_bam=$dir."/".$outputPrefix."_H_rm.bam";
         my $H_counts_G1=$dir."/"."H_G1_counts";
         my $H_counts_G2=$dir."/"."H_G2_counts";
         my $H_counts_G3=$dir."/"."H_G3_counts";

         my $H_G1_snp = $dir."/"."H_G1_snp";
         my $H_G2_snp = $dir."/"."H_G2_snp";
         my $H_G3_snp = $dir."/"."H_G3_snp";
         my $g1 = $dir."/"."P.counts";
         my $g2 = $dir."/"."M.counts";
         my $g3 = $dir."/"."PM.counts";
   
       
         my $L_bam = $dir."/". $outputPrefix."_L.bam";
         my $L_srt_bam = $dir."/".$outputPrefix."_L_srt.bam";
         my $L_rm_bam = $dir."/".$outputPrefix."_L_rm.bam";
         my $L_counts_G1 = $dir."/"."L_G1_counts";
         my $L_counts_G2 = $dir."/"."L_G2_counts";
         my $L_counts_G3 = $dir."/"."L_G3_counts";

         my $L_G1_snp = $dir."/". "L_G1_snp";
         my $L_G2_snp = $dir."/"."L_G2_snp";
         my $L_G3_snp = $dir."/"."L_G3_snp";
         my $tmp1 = $dir."/"."tmp1";
         my $tmp2 = $dir."/"."tmp2";
         my $sum =  $dir."/"."summary";


        open SS,">>$sum";

###########################################

               unless(defined($log)){

         $log=$work_dir."/".$outputPrefix.".log";
                             }

###########################################

       my %hlog=();



       if(-e $log){
           open LOD,"<$log";
           my $con=0;

######################################################
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
   
##################################################################
         if($con !=0 ){
                 
             print  "

             Program: bsatos (Bulked segregant analysis tools for outbreeding species)
             Version: 1.0
             The first module: prep
             prep      prepare the pool data\n\n\n\n


             This module has be processed last time\n\n\n";
                
              exit(0);                                  
                        }

#########################################################################
                }



############################################################################################

         print STDOUT "

             Program: bsatos (Bulked segregant analysis tools for outbreeding species)
             Version: 1.0
             The second module: prep
             prep        prepare the pool data\n\n\n\n";

        
         
    
        open (my $lo, ">>$log") or die $!;


        print $lo "

             Program: bsatos (Bulked segregant analysis tools for outbreeding species)
             Version: 1.0
             The second module: prep
             prep        prepare the pool data\n\n\n\n";



          my $idx=$genome.".bwt";
          my $fai=$genome.".fai"; 
#########################################################

         unless( -e $dir){
         
                      
         print STDOUT "create prep result dir\n";
         print $lo "create prep result dir\nCMD1:\t";
         unless(defined($hlog{"CMD1:"})){                
                        $hlog{"CMD1:"}=0;
                                         }  
         my $mx1=$hlog{"CMD1:"};       
         my $cd1="mkdir $dir"; 
         &process_cmd($cd1,$lo,$mx1);          
              
                         }

###############################################

         unless( -e $fai){  
 
         print STDOUT "samtools index reference genome\n";
         print $lo "samtools index reference genome\nCMD2:\t";
         unless(defined($hlog{"CMD2:"})){                
                         $hlog{"CMD2:"}=0;
                                         }
         my $mx2=$hlog{"CMD2:"};

         my $cd2="$samtools faidx $genome";
         &process_cmd($cd2,$lo,$mx2);
     
                          }

##############################################################

         unless(defined($Hbam) || defined($Lbam)){  

##################################################################
             unless(-e $idx){

         print STDOUT "bwa index reference genome\n";
         print $lo "bwa index reference genome\nCMD2:\t";
         unless(defined($hlog{"CMD2b:"})){
                         $hlog{"CMD2b:"}=0;
                                         }         
         my $mx2b=$hlog{"CMD2b:"};

         my $cd2b="$bwa index $genome";
         &process_cmd($cd2b,$lo,$mx2b);

                          }

                               
####################################################################

         print STDOUT "align High pool reads to reference genome\n";
         print $lo "align High pool reads to reference genome\nCMD3:\t";     
         unless(defined($hlog{"CMD3:"})){                
                         $hlog{"CMD3:"}=0;
                                         } 
         my $mx3=$hlog{"CMD3:"};
         my $cd3="$bwa mem -t $t -R \"\@RG\\tID:H\\tSM:H\" $genome $Hp1 $Hp2 |$samtools view -bS -> $H_bam";
         &process_cmd($cd3,$lo,$mx3);
 

###################################################################################################

         print STDOUT "align Low pool reads to reference genome\n";
         print $lo "align Low pools reads to reference genome\nCMD4:\t";
         unless(defined($hlog{"CMD4:"})){                
                         $hlog{"CMD4:"}=0;
                                         }
         my $mx4=$hlog{"CMD4:"};

         my $cd4="$bwa mem -t $t -R \"\@RG\\tID:L\\tSM:L\" $genome $Lp1 $Lp2 |$samtools view -bS -> $L_bam";
         &process_cmd($cd4,$lo,$mx4);
    
#################################################################################


         print STDOUT "sort H pool BAM FILE\n";
         print $lo "sort H pool BAM FILE\nCMD5:\t";
         unless(defined($hlog{"CMD5:"})){                
                         $hlog{"CMD5:"}=0;
                                         }
           
         my $mx5=$hlog{"CMD5:"};
         my $cd5="$samtools sort -o $H_srt_bam $H_bam";
         &process_cmd($cd5,$lo,$mx5);
   
################################################################################


         print STDOUT "sort L pool BAM FILE\n";
         print $lo "sort L pool BAM FILE\nCMD6:\t";
         my $cd6= "$samtools sort -o $L_srt_bam $L_bam";
         unless(defined($hlog{"CMD6:"})){                
                         $hlog{"CMD6:"}=0;
                                         } 
         my $mx6=$hlog{"CMD6:"};
         &process_cmd($cd6,$lo,$mx6);
  


##############################################################################


         print STDOUT "remove duplicates reads of H BAM FILE\n";
         print $lo "remove duplicates reads of H BAM FILE\nCMD7:\t";
         unless(defined($hlog{"CMD7:"})){          
                         $hlog{"CMD7:"}=0;
                                         }
         my $mx7=$hlog{"CMD7:"};
         my $cd7="$samtools rmdup   $H_srt_bam $H_rm_bam";
         &process_cmd($cd7,$lo,$mx7);
    

########################################################################

         print STDOUT "remove duplicates reads of L BAM FILE\n";
         print $lo "remove duplicates reads of L BAM FILE\nCMD8:\t";
         my $cd8 = "$samtools rmdup   $L_srt_bam $L_rm_bam";
         unless(defined($hlog{"CMD8:"})){          
                         $hlog{"CMD8:"}=0;
                                         }
         my $mx8=$hlog{"CMD8:"};

         &process_cmd($cd8,$lo,$mx8);


########################################################################

         print STDOUT "samtools index L BAM FILE\n";
         print $lo "samtools index L BAM FILE\nCMD9:\t";
         unless(defined($hlog{"CMD9:"})){          
                         $hlog{"CMD9:"}=0;
                                         }
         my $mx9=$hlog{"CMD9:"};

         my $cd9 ="$samtools index  $L_rm_bam";
         &process_cmd($cd9,$lo,$mx9);

#########################################################################

         print STDOUT "samtools flagstat L BAM FILE\n";
         print $lo "samtools flagstat L BAM FILE\nCMD9b:\t";
         print SS "summary the alignment of low pool alignment\n";

         unless(defined($hlog{"CMD9b:"})){
                         $hlog{"CMD9b:"}=0;
                                         }
         my $mx9b=$hlog{"CMD9b:"};

         my $cd9b ="$samtools flagstat $L_rm_bam >>$sum";
         &process_cmd($cd9b,$lo,$mx9b);

####################################################################
######################################################################

         print STDOUT "samtools index H BAM FILE\n";
         print $lo "samtools index H BAM FILE\nCMD10:\t";
         unless(defined($hlog{"CMD10:"})){          
                         $hlog{"CMD10:"}=0;
                                         }
         my $mx10=$hlog{"CMD10:"};
         my $cd10="$samtools index $H_rm_bam";
         &process_cmd($cd10,$lo,$mx10);

######################################################################
######################################################################

         print STDOUT "samtools flagstat H BAM FILE\n";
         print $lo "samtools flagstat H BAM FILE\nCMD9b:\t";
         print SS "\nsummary the alignment of high pool alignment\n";

         unless(defined($hlog{"CMD10b:"})){
                         $hlog{"CMD10b:"}=0;
                                         }
         my $mx10b=$hlog{"CMD10b:"};

         my $cd10b ="$samtools flagstat $H_rm_bam >>$sum";
         &process_cmd($cd10b,$lo,$mx10b);

#################################################################
##################################################################
####################################################################################
       
      if( -e $H_bam){
         print STDOUT "remove some files\n";
         print $lo "remove some files\nCMD11:\t";
         unless(defined($hlog{"CMD11:"})){          
                         $hlog{"CMD11:"}=0;
                                         }
         my $mx11=$hlog{"CMD11:"};
         my $cd11="rm $H_bam $L_bam $H_srt_bam $L_srt_bam";
         &process_cmd($cd11,$lo,$mx11); 
                    }




################################################################################

         print STDOUT "P loci genotyping in High pool\n";
         print $lo "P loci genotyping in High pool\nCMD12:\t";   
         unless(defined($hlog{"CMD12:"})){          
                         $hlog{"CMD12:"}=0;
                                         }
         my $mx12=$hlog{"CMD12:"};
   
         my $cd12 ="$samtools mpileup -ugf $genome -q $aq -l $G1 $H_rm_bam |$bcftools call -m ->$H_G1_snp";
         &process_cmd($cd12,$lo,$mx12); 


################################################################################

         print STDOUT "M loci genotyping in High pool\n";
         print $lo "M loci genotyping in High pool\nCMD13:\t"; 
         unless(defined($hlog{"CMD13:"})){          
                         $hlog{"CMD13:"}=0;
                                         }
         my $mx13=$hlog{"CMD13:"};
         my $cd13= "$samtools mpileup -ugf $genome -q $aq -l $G2 $H_rm_bam |$bcftools call -m ->$H_G2_snp";
         &process_cmd($cd13,$lo,$mx13);
 

###############################################################################

         print STDOUT "PM loci genotyping in High pool\n";
         print $lo "PM loci genotyping in High pool\nCMD14:\t"; 
         unless(defined($hlog{"CMD14:"})){          
                         $hlog{"CMD14:"}=0;
                                         }
         my $mx14=$hlog{"CMD14:"};                    
         my $cd14= "$samtools mpileup -ugf $genome -q $aq -l $G3 $H_rm_bam |$bcftools call -m ->$H_G3_snp";
         &process_cmd($cd14,$lo,$mx14);
 

##############################################################################

         print STDOUT "P loci genotyping in Low pool\n";
         print $lo "P loci genotyping in Low pool\nCMD15:\t";
        unless(defined($hlog{"CMD15:"})){          
                         $hlog{"CMD15:"}=0;
                                         }
         my $mx15=$hlog{"CMD15:"};

         my $cd15="$samtools mpileup -ugf $genome -q $aq -l $G1 $L_rm_bam |$bcftools call -m ->$L_G1_snp";
         &process_cmd($cd15,$lo,$mx15);




      
#########################################################################

         print STDOUT "M loci genotyping in Low pool\n";
         print $lo "M loci genotyping in Low pool\nCMD16:\t";
         unless(defined($hlog{"CMD16:"})){          
                         $hlog{"CMD16:"}=0;
                                         }
         my $mx16=$hlog{"CMD16:"};

         my $cd16="$samtools mpileup -ugf $genome -q $aq -l $G2 $L_rm_bam |$bcftools call -m ->$L_G2_snp";
         &process_cmd($cd16,$lo,$mx16);


###############################################################
    
         print STDOUT "G3 loci genotyping in Low pool\n";
         print $lo "G3 loci genotyping in Low pool\nCMD17\t";
         unless(defined($hlog{"CMD17:"})){          
                         $hlog{"CMD17:"}=0;
                                         }
         my $mx17=$hlog{"CMD17:"};

         my $cd17= "$samtools mpileup -ugf $genome -q $aq -l $G3 $L_rm_bam |$bcftools call -m ->$L_G3_snp";
         &process_cmd($cd16,$lo,$mx17);



###########################################################

                                                  }else{ 
            
############################################################

         print STDOUT "P loci genotyping in High pool\n";
         print $lo "P loci genotyping in High pool\nCMD18:\t"; 
         
         unless(defined($hlog{"CMD18:"})){          
                         $hlog{"CMD18:"}=0;
                                         }
         my $mx18=$hlog{"CMD18:"};

         my $cd18 ="$samtools mpileup -ugf $genome -q $aq -l $G1 $Hbam |$bcftools call -m ->$H_G1_snp";
         &process_cmd($cd18,$lo,$mx18);

###########################################################


         print STDOUT "M loci genotyping in High pool\n";
         print $lo "M loci genotyping in High pool\nCMD19:\t";
         unless(defined($hlog{"CMD19:"})){          
                         $hlog{"CMD19:"}=0;
                                         }
         my $mx19=$hlog{"CMD19:"};
         my $cd19= "$samtools mpileup -ugf $genome -q $aq -l $G2 $Hbam |$bcftools call -m ->$H_G2_snp";
         &process_cmd($cd19,$lo,$mx19);


###########################################################

         print STDOUT "PM loci genotyping in High pool\n";
         print $lo "PM loci genotyping in High pool\nCMD20:\t";
         unless(defined($hlog{"CMD20:"})){          
                         $hlog{"CMD20:"}=0;
                                         }
         my $mx20=$hlog{"CMD20:"};

         my $cd20= "$samtools mpileup -ugf $genome -q $aq -l $G3 $Hbam |$bcftools call -m ->$H_G3_snp";
         &process_cmd($cd20,$lo,$mx20);
           
 
#########################################################


         print STDOUT "P loci genotyping in Low pool\n";
         print $lo "P loci genotyping in Low pool\nCMD21:\t";
         unless(defined($hlog{"CMD21:"})){          
                         $hlog{"CMD21:"}=0;
                                         }
         my $mx21=$hlog{"CMD21:"};

         my $cd21="$samtools mpileup -ugf $genome -q $aq -l $G1 $Lbam |$bcftools call -m ->$L_G1_snp";
         &process_cmd($cd21,$lo,$mx21);

###################################################3


         print STDOUT "M loci genotyping in Low pool\n";
         print $lo "M loci genotyping in Low pool\nCMD22:\t";
         unless(defined($hlog{"CMD22:"})){          
                         $hlog{"CMD22:"}=0;
                                         }
         my $mx22=$hlog{"CMD22:"};

         my $cd22="$samtools mpileup -ugf $genome -q $aq -l $G2 $Lbam |$bcftools call -m ->$L_G2_snp";
         &process_cmd($cd22,$lo,$mx22);

################################################

         print STDOUT "PM loci genotyping in Low pool\n";
         print $lo "PM loci genotyping in Low pool\nCMD23:\t";
         unless(defined($hlog{"CMD23:"})){          
                         $hlog{"CMD23:"}=0;
                                         }
         my $mx23=$hlog{"CMD23:"};

         my $cd23= "$samtools mpileup -ugf $genome -q $aq -l $G3 $Lbam |$bcftools call -m ->$L_G3_snp";
         &process_cmd($cd23,$lo,$mx23);
                                                       }

###########################################


         print STDOUT "filter and get P loci reads count data in H pool\n";
         print $lo "filter and get P loci reads count data in H pool\nCMD24:\t";
         unless(defined($hlog{"CMD24:"})){          
                         $hlog{"CMD24:"}=0;
                                         }
         my $mx24=$hlog{"CMD24:"};

         my $cd24= "$get_counts -d $cov -q $sq  $H_G1_snp > $H_counts_G1";
         &process_cmd($cd24,$lo,$mx24);
 
###############################################

         print STDOUT "filter and get M loci reads count data in H pool\n";  
         print $lo "filter and get M loci reads count data in H pool\nCMD25:\t";    
         unless(defined($hlog{"CMD25:"})){          
                         $hlog{"CMD25:"}=0;
                                         }
         my $mx25=$hlog{"CMD25:"};
         my $cd25= "$get_counts -d $cov -q $sq  $H_G2_snp > $H_counts_G2";
         &process_cmd($cd25,$lo,$mx25);

############################################

         print STDOUT "filter and get PM loci reads count data in H pool\n";         
         print $lo "filter and get PM loci reads count data in H pool\nCMD26:\t";
         unless(defined($hlog{"CMD26:"})){          
                         $hlog{"CMD26:"}=0;
                                         }
         my $mx26=$hlog{"CMD26:"};
         my $cd26= "$get_counts -d $cov -q $sq  $H_G3_snp > $H_counts_G3";
         &process_cmd($cd26,$lo,$mx26);
 
#############################################

         print STDOUT "filter and get P loci reads count data in L pool\n";
         print $lo "filter and get P loci reads count data in L pool\nCMD27:\t";
         unless(defined($hlog{"CMD27:"})){          
                         $hlog{"CMD27:"}=0;
                                         }
         my $mx27=$hlog{"CMD27:"};
         my $cd27= "$get_counts -d $cov -q $sq  $L_G1_snp > $L_counts_G1";
         &process_cmd($cd27,$lo,$mx27);

####################################################


         print STDOUT "filter and get M loci reads count data in L pool\n";
         print $lo "filter and get M loci reads count data in L pool\nCMD28:\t";
         my $cd28= "$get_counts -d $cov -q $sq  $L_G2_snp > $L_counts_G2";
         unless(defined($hlog{"CMD28:"})){          
                         $hlog{"CMD28:"}=0;
                                         }
         my $mx28=$hlog{"CMD28:"};

         &process_cmd($cd28,$lo,$mx28);
         
##########################################

         print STDOUT "filter and get PM loci reads count data in L pool\n";         
         print $lo "filter and get PM loci reads count data in L pool\nCMD29:\t";
         my $cd29= "$get_counts -d $cov -q $sq  $L_G3_snp > $L_counts_G3";
         unless(defined($hlog{"CMD29:"})){          
                         $hlog{"CMD29:"}=0;
                                         }
         my $mx29=$hlog{"CMD29:"};

         &process_cmd($cd29,$lo,$mx29);

#########################################

           
         print STDOUT "obtain P loci reads count data\n";
         print $lo "obtain P loci reads count data\nCMD30:\t";
         unless(defined($hlog{"CMD30:"})){          
                         $hlog{"CMD30:"}=0;
                                         }
         my $mx30=$hlog{"CMD30:"};

         my $cd30= "sort -k 1,1 -k 2,2n $H_counts_G1 > $tmp1;sort -k 1,1 -k 2,2n $L_counts_G1 >$tmp2;filter -k cn $tmp1 $tmp2 | cut -f 1-4,7-8 - |$filter_dep -d $pn -> $g1";
         &process_cmd($cd30,$lo,$mx30); 

############################################################################
   
     
         print STDOUT "summary P type loci data\n";
         print $lo "summary P loci reads count data\nCMD30b:\t";
         print SS "summary P type loci data\n";
         unless(defined($hlog{"CMD30b:"})){
                         $hlog{"CMD30b:"}=0;
                                         }
         my $mx30b=$hlog{"CMD30b:"};

         my $cd30b= "wc -l $g1 >>$sum";
         &process_cmd($cd30b,$lo,$mx30b);

######################################################################
#######################################################################

         print STDOUT "obtain M loci reads count data\n";
         print $lo "obtain M loci reads count data\nCMD31\t";
         unless(defined($hlog{"CMD31:"})){          
                         $hlog{"CMD31:"}=0;
                                         }
         my $mx31=$hlog{"CMD31:"};

         my $cd31= "sort -k 1,1 -k 2,2n $H_counts_G2 > $tmp1;sort -k 1,1 -k 2,2n $L_counts_G2 >$tmp2;filter -k cn $tmp1 $tmp2 | cut -f 1-4,7-8 - |$filter_dep -d $pn ->$g2";
         &process_cmd($cd31,$lo,$mx31); 
        
######################################################################
         print STDOUT "summary M type loci data\n";
         print $lo "summary M loci reads count data\nCMD31b:\t";
         print SS "summary M type loci data\n";
         unless(defined($hlog{"CMD31b:"})){
                         $hlog{"CMD31b:"}=0;
                                         }
         my $mx31b=$hlog{"CMD31b:"};

         my $cd31b= "wc -l $g2 >>$sum";
         &process_cmd($cd31b,$lo,$mx31b);


######################################################################
##########################################################################

         print STDOUT "obtain G3 loci reads count data\n";      
         print $lo "obtain G3 loci reads count data\nCMD32\t";
         unless(defined($hlog{"CMD32:"})){          
                         $hlog{"CMD32:"}=0;
                                         }
         my $mx32=$hlog{"CMD32:"};

         my $cd32= "sort -k 1,1 -k 2,2n $H_counts_G3 > $tmp1;sort -k 1,1 -k 2,2n $L_counts_G3 >$tmp2;filter -k cn $tmp1 $tmp2 | cut -f 1-4,7-8 - |$filter_dep -d $pn ->$g3";

         &process_cmd($cd32,$lo,$mx32);   

########################################################################
###############################################################
 
         print STDOUT "summary PM type loci data\n";
         print $lo "summary PM loci reads count data\nCMD32b:\t";
         print SS "summary PM type loci data\n";
         unless(defined($hlog{"CMD32b:"})){
                         $hlog{"CMD32b:"}=0;
                                         }
         my $mx32b=$hlog{"CMD32b:"};

         my $cd32b= "wc -l $g3 >>$sum";
         &process_cmd($cd32b,$lo,$mx32b);

#####################################
         
         if( -e $H_G1_snp){
         print STDOUT "remove some files\n";
         print $lo "remove some files\nCMD33\t";
         unless(defined($hlog{"CMD33:"})){          
                         $hlog{"CMD33:"}=0;
                                         }
         my $mx33=$hlog{"CMD33:"};

         my $cd33="rm $H_G1_snp $H_G2_snp $H_G3_snp $L_G1_snp $L_G2_snp $L_G3_snp $tmp1 $tmp2 $H_counts_G1 $H_counts_G2 $H_counts_G3 $L_counts_G1 $L_counts_G2 $L_counts_G3"; 
         &process_cmd($cd33,$lo,$mx33); 
                        }
######################################

         print STDOUT "prep module finished!!!\n";
         print $lo "prep module finished!!!\n";


         exit(0);
              
                         
   
 
              }


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
