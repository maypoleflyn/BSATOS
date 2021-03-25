

package prepar;

BEGIN {
	$VERSION = '1.0';
}


use strict;
use Getopt::Long;
use base 'Exporter';
use Cwd;
our @EXPORT = qw(runprepar);

my $G_USAGE = "

Usage:  bastos prepar [options] 

Options:  
          --pf1       FILE        paired-end1 fastq file from pollen parent  [FASTQ format]
          --pf2       FILE        paired-end2 fastq file from pollen parent  [FASTQ format]
          --mf1       FILE        paired-end1 fastq file from maternal parent [FASTQ format]
          --mf2       FILE        paired-end2 fastq file from maternal parent [FASTQ format]
          --pb        FILE        pre-aligned bam file from pollen parent  [Not required, if --Pp1 & --Pp2]
          --mb        FILE        pre-aligned bam file from maternal parent [Not required, if --Mp1 & --Mp2]  
          --r         FILE        reference genome [FASTA format]
          --o         FILE        output dir name prefix  [prepar] 
          --gtf/gff   FILE        gene GTF/GFF file [GTF2/GFF3 format]
          --log       FILE        prepar module log file [prepar.log]          
 


Aligment Options:

           --t  INT        number of threads [1]

SNV Calling Options:

           --aq  INT       skip alignments with mapQ smaller than INT [30]

SV Calling Options:
          
           --vq INT        skip alignments with mapQ smaller than INT [30]

SV filtering Options:

           --dp        INT        support reads number including paired-reads and split reads [10]            
           --mq        INT        average mapQ [20]
           --l         INT        the minimum SV length [1000]
           --precise   STR        only kept precise SVs [yes]






SNV filtering Options:

           --d  INT        skip SNVs with reads depth smaller than INT [10]   
           --sq INT        skip SNVs with phred-scaled quality smaller than [30]
                       

Outputs:

   prefix_dir  [DIR]

   P_M.p  [FILE]  the genotype of the markers are homozygous in maternal parent but are heterozygous in pollen parent       
   P_M.m  [FILE]  the genotype of the markers are homozygous in pollen parent but is heterozygous in maternal parent
   P_M.pm [FILE]  the genotype of the markers are both heterozygous 
   M_P.snv [FILE] the SNVs vcf file of two parents   
   sv.vcf [FILE]  the SVs vcf file of two  parents

   prefix_dir/anno [DIR]

   AT_refGene.txt [FILE]        GenePred file
   AT_refGeneMrna.fa [FILE]     transcript FASTA file
   snv.AT_multianno.txt [FILE]  SNVs multianno file 
   snv.AT_multianno.vcf [FILE]  SNVs annotated VCF file
   snv.avinput [FILE]           SNVs input file
   sv.AT_multianno.txt [FILE]   SVs multianno file 
   sv.AT_multianno.vcf [FILE]   SVs annotated VCF file
   sv.avinput [FILE]            SVs input file
      

         
Example:
 
 1) bsatos prepar --pf1 P_1.fastq --pf2 P_2.fastq --mf1 M_1.fastq --mf2 M_fastq --r genome.fasta --gtf gene.gtf --o prepar
 OR
 2) bsatos prepar --pb P.bam --mb M.bam --r genome.fasta --gtf gene.gtf --o prepar 

";

sub runprepar {
	my $Pp1 = undef;
 	my $Pp2 = undef;
        my $Mp1 = undef;
        my $Mp2 = undef;
        my $genome = undef;
        my $Pbam = undef;
        my $Mbam = undef;   
	my $outputPrefix = undef;
	my $help = 0;
	my $gtf = undef;
        my $s = undef;
        my $t = undef;
        my $aq = undef;
        my $vq = undef;
        my $dep = undef;   
        my $sq = undef;
        my $gff = undef;
        my $dp  = undef;
        my $l   = undef;
        my $mq  = undef;
        my $precise = undef;
        my $log = undef;
	GetOptions (
	"pf1=s" => \$Pp1,
	"pf2=s" => \$Pp2,
        "mf1=s" => \$Mp1,
        "mf2=s" => \$Mp2,
        "pb=s" => \$Pbam,
        "mb=s" => \$Mbam,
        "r=s" => \$genome,
	"o=s"   => \$outputPrefix,
        "gtf=s" => \$gtf,
        "gff=s" => \$gff,
        "s=s" =>\$s,
        "t=s" =>\$t,
        "aq=i" =>\$aq,
        "vq=i" =>\$vq,
        "d=i" =>\$dep, 
        "sq=i" =>\$sq,
        "dp=i" =>\$dp,
        "precise=s" =>\$precise,
        "mq=i" =>\$mq,
        "l=i" =>\$l,
        "log=s" =>\$log,
	"help!" => \$help)
	or die("$G_USAGE");
	
	die "$G_USAGE" if ($help);
	
	die "$G_USAGE" if (!defined $outputPrefix);


####################################################################

    unless(defined($outputPrefix)){
                   $outputPrefix="prepar";
                                  }

        unless(defined($dp)){
              $dp=10;
                           }
       unless(defined($mq)){
                $mq=20;
                           }
       unless(defined($precise)){
                 $precise="yes";
                                }
       unless(defined($l)){
                 $l=1000;
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
        unless(defined($dep)){
                        $dep=10;
                             } 
        unless(defined($sq)){
                        $sq=30;
                            }  
 
#####################################################################
           
         my $samtools= "samtools";
         my $bwa="bwa";
         my $bcftools="bcftools";
         my $delly ="delly";
         my $gg="$FindBin::Bin/scripts/filter_get_genotye.pl";
         my $gtfpred="gtfToGenePred";
         my $gffread="gffread";      
         my $filter_sv="$FindBin::Bin/scripts/filter_raw_sv.pl"; 
         my $work_dir= getcwd; 
         my $dir=$work_dir."/".$outputPrefix."_dir";


        unless(defined($log)){
                           
         $log=$work_dir."/".$outputPrefix.".log";
                             
                             }

##########################################################

       my %hlog=();       
 
       if(-e $log){
        
           open LOD,"<$log"; 
           my $con=0;  

################################################
                                  
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
                      
########################################################
            if($con !=0 ){                 
                 print "
             Program: bsatos (Bulked segregant analysis tools for outbreeding species)
             Version: 1.0
             The first module: prepar
             prepar      prepare the parents data\n\n\n\n


             This module has be finished!!!!!\n\n\n";

                   exit(0);
                             
                        }
##########################################################     
              }
                    


##############################################################################
   
                                                     
         my $anno_dir=$dir."/anno";

         my $P_bam=$dir."/".$outputPrefix."_P.bam";
         my $P_srt_bam=$dir."/".$outputPrefix."_P_srt.bam";
         my $P_rm_bam=$dir."/".$outputPrefix."_P_rm.bam";
               
         my $M_bam = $dir."/".$outputPrefix."_M.bam";
         my $M_srt_bam = $dir."/".$outputPrefix."_M_srt.bam";
         my $M_rm_bam = $dir."/".$outputPrefix."_M_rm.bam";
         my $M_P_snp=$dir."/"."P_M.snv";    
         my $M_P_del_bcf=$dir."/"."P_M.del.bcf";
         my $M_P_del_bcf_csi=$dir."/"."P_M.del.bcf.csi"; 
         my $M_P_del_vcf=$dir."/"."P_M.del.vcf";
         my $M_P_ins_bcf=$dir."/"."P_M.ins.bcf";
         my $M_P_ins_bcf_csi=$dir."/"."P_M.ins.bcf.csi";
         my $M_P_ins_vcf=$dir."/"."P_M.ins.vcf";
         my $sv=$dir."/"."sv.vcf";       
         my $sum=$dir."/"."summary";         
         
       
         open SS,">>$sum";
       
#####################################################################################       
     
       print STDOUT "
             Program: bsatos (Bulked segregant analysis tools for outbreeding species)
             Version: 1.0
             The first module: prepar
             prepar      prepare the parents data\n\n\n\n";
                




        open (my $lo,">>$log") or die $!;
        

        print $lo "
             Program: bsatos (Bulked segregant analysis tools for outbreeding species)
             Version: 1.0
             The first module: prepar
             prepar      prepare the parents data\n";
            
        my $idx=$genome.".bwt";
        my $fai=$genome.".fai";
#############################################################

        if(defined($gff)){

          $gtf=$dir."/"."gtf";  

          print STDOUT "converting GFF to GTF file:\n";
          print $lo "converting GFF to GTF file:\nCMD1:\t";
          unless(defined($hlog{"CMD1:"})){                
                         $hlog{"CMD1:"}=0;
                                         }
           
          my $mx1=$hlog{"CMD1:"};
           
          my $cd1 ="$gffread $gff -o $gtf -T";
          &process_cmd($cd1,$lo,$mx1);
                         }

############################################################

#############################################
          unless( -f $fai){
 
          print STDOUT "samtools index reference genome\n";
          print $lo "samtools index genome\nCMD3:\n";
          unless(defined($hlog{"CMD3:"})){
                         $hlog{"CMD3:"}=0;
                                         }
          my $mx3=$hlog{"CMD3:"};
          my $cd3="$samtools faidx $genome";
          &process_cmd($cd3,$lo,$mx3);

###################################################                
                          }
           

#######################################

        unless( -e $dir){


          print STDOUT "creat result dir\n";
          print $lo "creat result dir\nCMD4:\t"; 
           unless(defined($hlog{"CMD4:"})){
                         $hlog{"CMD4:"}=0;
                                         }
          my $mx4=$hlog{"CMD4:"};
          my $cd4="mkdir $dir";
          &process_cmd($cd4,$lo,$mx4);                    
                         }


##############################################
 
         my $P_M=$dir."/"."P_M";    
  unless(defined($Pbam) || defined($Mbam)){     

######################################################

          unless( -f $idx){

          print STDOUT "bwa index reference genome\n";
          print $lo "bwa index reference genome\nCMD2:\t";
          unless(defined($hlog{"CMD2:"})){
                         $hlog{"CMD2:"}=0;
                                         }
          my $mx2=$hlog{"CMD2:"};
          my $cd2="$bwa index $genome";
          &process_cmd($cd2,$lo,$mx2);
                          }

################################################################################################

         print STDOUT "align pollen parent reads to reference genome and generate BAM file\n"; 
         print $lo "align pollen parent reads to reference genome and generate BAM file\nCMD5:\t";
          unless(defined($hlog{"CMD5:"})){
                         $hlog{"CMD5:"}=0;
                                         }
          my $mx5=$hlog{"CMD5:"};
         my $cd5="$bwa mem -t $t -R \"\@RG\\tID:P\\tSM:P\" $genome $Pp1 $Pp2 |$samtools view -bS -> $P_bam";
         &process_cmd($cd5,$lo,$mx5);
      
#####################################################################################################
   
         print STDOUT "align maternal parent reads to reference genome and generate BAM file\n";  
         print $lo "align maternal parent reads to reference genome and generate BAM file\nCMD6:\t";
         unless(defined($hlog{"CMD6:"})){
                         $hlog{"CMD6:"}=0;
                                         }
          my $mx6=$hlog{"CMD6:"};
         my $cd6="$bwa mem -t $t -R \"\@RG\\tID:M\\tSM:M\" $genome $Mp1 $Mp2 |$samtools view -bS -> $M_bam";
         &process_cmd($cd6,$lo,$mx6);
        
############################################################################

         print STDOUT "sort pollen parent BAM file\n";
         print $lo "sort pollen parent BAM file\nCMD7:\t";   
         unless(defined($hlog{"CMD7:"})){
                         $hlog{"CMD7:"}=0;
                                         }
         my $mx7=$hlog{"CMD7:"};                
         my $cd7="$samtools sort -o $P_srt_bam $P_bam"; 
         &process_cmd($cd7,$lo,$mx7);
         

#######################################################################

         print STDOUT "sort maternal parent BAM file\n"; 
         print $lo "sort maternal parent BAM file\nCMD8:\t";
         unless(defined($hlog{"CMD8:"})){
                         $hlog{"CMD8:"}=0;
                                         }
         
         my $mx8=$hlog{"CMD8:"};
         my $cd8="$samtools sort -o $M_srt_bam $M_bam";
         &process_cmd($cd8,$lo,$mx8);


#######################################################################

         print $STDOUT "remove duplicates pollen parent BAM file\n";
         print $lo "remove duplicates pollen parent BAM file\nCMD9:\t";
         unless(defined($hlog{"CMD9:"})){
                         $hlog{"CMD9:"}=0;
                                         }
         my $mx9=$hlog{"CMD9:"};
         my $cd9="$samtools rmdup   $P_srt_bam $P_rm_bam";
         &process_cmd($cd9,$lo,$mx9);
       

#####################################################################

         print STDOUT "remove duplicates maternal parent BAM file\n";
         print $lo "remove duplicates maternal parent BAM file\nCMD10:\t";
         unless(defined($hlog{"CMD10:"})){
                         $hlog{"CMD10:"}=0;
                                         }
         my $mx10=$hlog{"CMD10:"};
         my $cd10="$samtools rmdup   $M_srt_bam $M_rm_bam";
         &process_cmd($cd10,$lo,$mx10);
         

#################################################################

         print STDOUT "samtools index pollen parent rmdup BAM file\n"; 
         print $lo "samtools index pollen parent rmdup BAM file\nCMD11:\t";
         unless(defined($hlog{"CMD11:"})){
                         $hlog{"CMD11:"}=0;
                                         }
         my $mx11=$hlog{"CMD11:"};
         my $cd11="$samtools index  $P_rm_bam";
         &process_cmd($cd11,$lo,$mx11);

###################################################################
        print STDOUT "summary the alignment result\n";
        print $lo "summary the alignment result\nCMDs1:\t";
        print SS "summary the alignment result of pollen parent:\n\n";
        unless(defined($hlog{"CMDs1:"})){
                         $hlog{"CMDs1:"}=0;
                                        }
       my $name_s1=$hlog{"CMDs1:"};
       my $cd_s1="$samtools flagstat $P_rm_bam >>$sum";
       &process_cmd($cd_s1,$lo,$name_s1); 
##########################################################################

         print STDOUT "samtools index maternal parent rmdup BAM file\n";   
         print $lo "samtools index maternal parent rmdup BAM file\nCMD12:\t"; 
         unless(defined($hlog{"CMD12:"})){
                         $hlog{"CMD12:"}=0;
                                         }
         
         my $mx12=$hlog{"CMD12:"};
         my $cd12="$samtools index $M_rm_bam";
         &process_cmd($cd12,$lo,$mx12);  

############################################################################
############################################################################
        print STDOUT "summary the alignment result\n";
        print $lo "summary the alignment result\nCMDs2:\t";
        print SS "\nsummary the alignment result of maternal parent:\n\n";
        unless(defined($hlog{"CMDs2:"})){
                         $hlog{"CMDs2:"}=0;
                                        }
       my $name_s2=$hlog{"CMDs2:"};
       my $cd_s2="$samtools flagstat $M_rm_bam >>$sum";
       &process_cmd($cd_s2,$lo,$name_s2);

#############################################################

         print STDOUT "parent SNV calling\n";
         print $lo "parent SNV calling\nCMD13:\t";
         unless(defined($hlog{"CMD13:"})){
                         $hlog{"CMD13:"}=0;
                                         }
         my $mx13=$hlog{"CMD13:"};
         my $cd13 ="$samtools mpileup -ugf $genome -q $aq $P_rm_bam $M_rm_bam |$bcftools call -mv ->$M_P_snp";
         &process_cmd($cd13,$lo,$mx13);

################################################################
###############################################################


         print STDOUT  "DEL SV calling\n";
         print $lo "DEL SV calling\nCMD14:\t";
         unless(defined($hlog{"CMD14:"})){
                         $hlog{"CMD14:"}=0;
                                         }
         my $mx14=$hlog{"CMD14:"};
         my $cd14="$delly call -q $vq -g $genome -t DEL -o $M_P_del_bcf $P_rm_bam $M_rm_bam";
         &process_cmd($cd14,$lo,$mx14);
       

###############################################################

         print STDOUT "convert DEL BCF file to DEL VCF file\n";
         print $lo "convert DEL BCF file to DEL VCF file\nCMD15:\t";
         unless(defined($hlog{"CMD15:"})){
                         $hlog{"CMD15:"}=0;
                                         }
         
         my $mx15=$hlog{"CMD15:"};
        #my $cd15 = "$bcftools view $M_P_del_bcf > $M_P_del_vcf";
         my $cd15="$bcftools view $M_P_del_bcf |$filter_sv --d $dp --mq 20 --l 1000 --precise $precise --header yes -> $M_P_del_vcf";
         &process_cmd($cd15,$lo,$mx15);     



########################################################
   
         print STDOUT "INS SV calling\n";
         print $lo "INS SV calling\nCMD16:\t";
         unless(defined($hlog{"CMD16:"})){
                         $hlog{"CMD16:"}=0;
                                         }
         my $mx16=$hlog{"CMD16:"};
         my $cd16="$delly call -q $vq -g $genome -t INS -o $M_P_ins_bcf $P_rm_bam $M_rm_bam";
        # my $cd16= "$bcftools view $M_P_ins_bcf | $filter_sv --d $dp --mq $mq --l $l --precise  $precise --header no -> $M_P_ins_vcf";
         &process($cd16,$lo,$mx16);
         
#####################################################################


         print STDOUT "convert INS BCF file to VCF file\n"; 
         print $lo "convert INS BCF file to VCF file\nCMD17:\t";
         unless(defined($hlog{"CMD17:"})){
                         $hlog{"CMD17:"}=0;
                                         }         
         my $mx17=$hlog{"CMD17:"};
         #my $cd17="$bcftools view $M_P_ins_bcf > $M_P_ins_vcf";
         my $cd17= "$bcftools view $M_P_ins_bcf | $filter_sv --d $dp --mq $mq --l $l --precise  $precise --header no -> $M_P_ins_vcf";
         &process_cmd($cd17,$lo,$mx17);
       

###########################################################

         print STDOUT "convert INS BCF file to VCF file\n"; 
         print $lo "merge and filter SV\nCMD18:\t";      
         unless(defined($hlog{"CMD18:"})){
                         $hlog{"CMD18:"}=0;
                                         }
         my $mx18=$hlog{"CMD18:"};
         my $cd18="cat $M_P_del_vcf $M_P_ins_vcf |grep -v \"#\" >$sv";
         &process_cmd($cd18,$lo,$mx18); 
         

#######################################

   if(-e $M_P_del_bcf){
         print STDOUT "remove files\n";
         print $lo "remove files\nCMD19:\t";
         unless(defined($hlog{"CMD19:"})){
                         $hlog{"CMD19:"}=0;
                                         }
         my $mx19=$hlog{"CMD19:"};
         my $cd19="rm $M_P_del_bcf $M_P_ins_bcf $M_P_del_bcf_csi $M_P_ins_bcf_csi";
         &process_cmd($cd19,$lo,$mx19);
                      }
########################################################                     

                    }else{
        
#############################################################

        print STDOUT "samtools index pollen parent BAM file\n"; 
        print $lo "samtools index pollen parent BAM file\nCMD20:\t";
        unless(defined($hlog{"CMD20:"})){
                         $hlog{"CMD20:"}=0;
                                         }
        my $mx20=$hlog{"CMD20:"};
        my $cd20="$samtools index $Pbam";
        &process_cmd($cd20,$lo,$mx20);

#############################################################

        print STDOUT "summary the alignment result\n";
        print $lo "summary the alignment result\nCMDs1:\t";
        print SS "summary the alignment result of pollen parent\n";
        unless(defined($hlog{"CMDs1:"})){
                         $hlog{"CMDs1:"}=0;
                                        }
       my $name_s1=$hlog{"CMDs1:"};
       my $cd_s1="$samtools flagstat $Pbam >>$sum";
       &process_cmd($cd_s1,$lo,$name_s1);

############################################################
############################################################
        print STDOUT "samtools index maternal parent BAM file\n";
        print $lo "samtools index maternal parent BAM file\nCMD21:\t";
        unless(defined($hlog{"CMD21:"})){
                         $hlog{"CMD21:"}=0;
                                         }
        my $mx21=$hlog{"CMD21:"};
        my $cd21="$samtools index $Mbam";
        &process_cmd($cd21,$lo,$mx21);

############################################################

        print STDOUT "summary the alignment result\n";
        print $lo "summary the alignment result\nCMDs2:\t";
        print SS "summary the alignment result of maternal parent\n";
        unless(defined($hlog{"CMDs2:"})){
                         $hlog{"CMDs2:"}=0;
                                        }
       my $name_s2=$hlog{"CMDs2:"};
       my $cd_s2="$samtools flagstat $Mbam >>$sum";
       &process_cmd($cd_s2,$lo,$name_s2);

############################################################
############################################################

        
        print STDOUT "SNV calling\n";  
        print $lo "SNV calling\nCMD22:\t";
        unless(defined($hlog{"CMD22:"})){
                         $hlog{"CMD22:"}=0;
                                         }
        my $mx22=$hlog{"CMD22:"};
        my $cd22="$samtools mpileup -ugf $genome -q $aq $Pbam $Mbam |$bcftools call -mv ->$M_P_snp";
        &process_cmd($cd22,$lo,$mx22);
  
########################################################

        print STDOUT "DEL SV calling\n";      
        print $lo "DEL SV calling\nCMD23:\t";
        unless(defined($hlog{"CMD23:"})){
                         $hlog{"CMD23:"}=0;
                                         }
         my $mx23=$hlog{"CMD23:"};

        my $cd23="$delly call -q $vq -g $genome -t DEL -o $M_P_del_bcf $Pbam $Mbam";
        &process_cmd($cd23,$lo,$mx23);


###################################################

        print STDOUT "convert DEL BCF to VCF\n";
        print $lo "convert DEL BCF to VCF\nCMD24:\t";
        unless(defined($hlog{"CMD24:"})){
                         $hlog{"CMD24:"}=0;
                                         }
        my $mx24=$hlog{"CMD24:"};
        my $cd24="$bcftools view $M_P_del_bcf |$filter_sv --d $dp --mq 20 --l 1000 --precise $precise --header yes -> $M_P_del_vcf";
        &process_cmd($cd24,$lo,$mx24);

#####################################################


        print STDOUT "INS SV calling\n";
        print $lo "INS SV calling\nCMD25:\t";
        unless(defined($hlog{"CMD25:"})){
                         $hlog{"CMD25:"}=0;
                                         }
         my $mx25=$hlog{"CMD25:"};

        my $cd25= "$delly call -q $vq -g $genome -t INS -o $M_P_ins_bcf $Pbam $Mbam";
        &process_cmd($cd25,$lo,$mx25);


######################################################

        print STDOUT "convert INS BCF to VCF\n";
        print $lo "convert INS BCF to VCF\nCMD26:\t";
        unless(defined($hlog{"CMD26:"})){
                         $hlog{"CMD26:"}=0;
                                         }
         my $mx26=$hlog{"CMD26:"};

        my $cd26= "$bcftools view $M_P_ins_bcf | $filter_sv --d $dp --mq $mq --l $l --precise  $precise --header no -> $M_P_ins_vcf";
        &process_cmd($cd26,$lo,$mx26);

##############################################
        print STDOUT "merge and filter SV file\n";
        print $lo "merge and filter SV file\nCMD27:\t";
        my $cd27="cat $M_P_del_vcf $M_P_ins_vcf |grep -v \"#\" >$sv";
        unless(defined($hlog{"CMD27:"})){
                         $hlog{"CMD27:"}=0;
                                         }
         my $mx27=$hlog{"CMD27:"};

        &process_cmd($cd27,$lo,$mx27);



##################################################

       if(-e $M_P_del_bcf){
        print STDOUT "remove SV files\n"; 
        print $lo "remove SV files\nCMD28:\t";
        unless(defined($hlog{"CMD28:"})){
                         $hlog{"CMD28:"}=0;
                                         }
        my $mx28=$hlog{"CMD28:"};
        my $cd28 = "rm $M_P_del_bcf $M_P_ins_bcf $M_P_del_bcf_csi $M_P_ins_bcf_csi $M_P_del_vcf $M_P_ins_vcf";
        &process_cmd($cd28,$lo);
                          }    

 ########################################################
                          }

######################################################

        print STDOUT "filter SNV based on reads depth and phred-scaled quality\n"; 
        print $lo "filter SNV based on reads depth and phred-scaled quality\nCMD29:\t";
        my $cd29="$gg -d $dep -q $sq -o $P_M $M_P_snp >>$sum";
        unless(defined($hlog{"CMD29:"})){
                         $hlog{"CMD29:"}=0;
                                         }
         my $mx29=$hlog{"CMD29:"};

        &process_cmd($cd29,$lo,$mx29);  

###############################################################
#######################################################

        print STDOUT "create SNV and SV annotation dir\n";
        print $lo "create SNV and SV annotation dir\nCMD30:\t";
        unless(defined($hlog{"CMD30:"})){
                         $hlog{"CMD30:"}=0;
                                         }
         my $mx30=$hlog{"CMD30:"};
        my $cd30= "mkdir $anno_dir";
        my $mx30="CMD30:";
    unless(-e $anno_dir){
        &process_cmd($cd30,$lo,$mx30);
                        }
        my $out_snv=$anno_dir."/snv";
        my $out_sv= $anno_dir."/sv";

##############################################
        
        print STDOUT "create genome annotation database step1\n";
        print $lo "create genome annotation database step1\nCMD31:\t";    
        unless(defined($hlog{"CMD31:"})){
                         $hlog{"CMD31:"}=0;
                                         }
         my $mx31=$hlog{"CMD31:"};
        my $cd31 = "$gtfpred -genePredExt $gtf $anno_dir/AT_refGene.txt";
        &process_cmd($cd31,$lo,$mx31);
      
################################################


        print STDOUT "creat genome annotation database step2\n";
        print $lo "creat genome annotation database step2\nCMD32:\t";   
        unless(defined($hlog{"CMD32:"})){
                         $hlog{"CMD32:"}=0;
                                         }
         my $mx32=$hlog{"CMD32:"};
        my $cd32= "retrieve_seq_from_fasta.pl --format refGene --seqfile $genome $anno_dir/AT_refGene.txt --outfile  $anno_dir/AT_refGeneMrna.fa";
        &process_cmd($cd32,$lo,$mx32);

################################################### 


       print STDOUT "annotate SNV file\n"; 
       print $lo "annotate SNV file\nCMD33:\t";  
       unless(defined($hlog{"CMD33:"})){
                         $hlog{"CMD33:"}=0;
                                         }
         my $mx33=$hlog{"CMD33:"};

       my $cd33 = "table_annovar.pl $M_P_snp $anno_dir -buildver AT -out $out_snv -remove   -protocol refGene -operation g -nastring . -vcfinput";
       &process_cmd($cd33,$lo,$mx33);   


############################################
              
       print STDOUT "annotate SV file\n";
       print $lo "annotate SV file\nCMD34:\t";
       unless(defined($hlog{"CMD34:"})){
                         $hlog{"CMD34:"}=0;
                                         }
       my $mx34=$hlog{"CMD34:"};
       my $cd34 ="table_annovar.pl $sv $anno_dir -buildver AT -out $out_sv -remove   -protocol refGene -operation g -nastring . -vcfinput";
       &process_cmd($cd34,$lo,$mx34);

#############################################################

      




###############################################

       print STDOUT "prepar module finished!!!";               
       print $lo "prepar module finished!!!";

###################################

                          }




####################################################################

 sub process_cmd {
    my ($cmd,$lk,$hg) = @_;
    my $time="[".localtime()."]";
#######################
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

########################################

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
