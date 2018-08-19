

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

Options:  --pf1       FILE        paired-end1 fastq file from pollen parent  [FASTQ format]
          --pf2       FILE        paired-end2 fastq file from pollen parent  [FASTQ format]
          --mf1       FILE        paired-end1 fastq file from maternal parent [FASTQ format]
          --mf2       FILE        paired-end2 fastq file from maternal parent [FASTQ format]
          --pb        FILE        pre-aligned bam file from pollen parent  [Not required, if --Pp1 & --Pp2]
          --mb        FILE        pre-aligned bam file from maternal parent [Not required, if --Mp1 & --Mp2]  
          --r         FILE        reference genome [FASTA format]
          --o         FILE        output dir name prefix  [prepar] 
          --gtf/gff   FILE        gene GTF/GFF file [GTF2/GFF3 format]
          



Aligment Options:

           --t  INT        number of threads [1]

SNV Calling Options:

           --aq  INT       skip alignments with mapQ smaller than INT [30]

SV Calling Options:
          
           --vq INT        skip alignments with mapQ smaller than INT [30]

SNV filtering Options:

           --d  INT        skip SNVs with reads depth smaller than INT [10]   
           --sq INT        skip SNVs with phred-scaled quality smaller than [30]
                       

Outputs:

   prefix_dir  [DIR]

   P_M_G1 [FILE]  the genotype of the markers are homozygous in maternal parent but are heterozygous in pollen parent       
   P_M_G2 [FILE]  the genotype of the markers are homozygous in pollen parent but is heterozygous in maternal parent
   P_M_G3 [FILE]  the genotype of the markers are both heterozygous 
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
 2) bastos prepar --pb P.bam --mb M.bam --r genome.fasta --gtf gene.gtf --o prepar 

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
	"help!" => \$help)
	or die("$G_USAGE");
	
	die "$G_USAGE" if ($help);
	
	die "$G_USAGE" if (!defined $outputPrefix);


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
 
           
         my $samtools= "samtools";
         my $bwa="bwa";
         my $bcftools="bcftools";
         my $delly ="delly";
         my $gg="$FindBin::Bin/scripts/filter_get_genotye.pl";
         my $gtfpred="gtfToGenePred";
         my $gffread="gffread";      
         
         my $work_dir= getcwd; 
         my $dir=$work_dir."/prepar_dir";
         my $anno_dir=$dir."/anno";

         my $P_bam=$dir."/".$outputPrefix."_P.bam";
         my $P_srt_bam=$dir."/".$outputPrefix."_P_srt.bam";
         my $P_rm_bam=$dir."/".$outputPrefix."_P_rm.bam";
               
         my $M_bam = $dir."/".$outputPrefix."_M.bam";
         my $M_srt_bam = $dir."/".$outputPrefix."_M_srt.bam";
         my $M_rm_bam = $dir."/".$outputPrefix."_M_rm.bam";
         my $M_P_snp=$dir."/"."M_P.snv";    
         my $M_P_del_bcf=$dir."/"."M_P.del.bcf";
         my $M_P_del_bcf_csi=$dir."/"."M_P.del.bcf.csi"; 
         my $M_P_del_vcf=$dir."/"."M_P.del.vcf";
         my $M_P_ins_bcf=$dir."/"."M_P.ins.bcf";
         my $M_P_ins_bcf_csi=$dir."/"."M_P.ins.bcf.csi";
         my $M_P_ins_vcf=$dir."/"."M_P.ins.vcf";
         my $sv=$dir."/"."sv.vcf";       
         
         
       

        my $script = $outputPrefix.".sh";
       
    
        open (my $ofh, ">$script") or die $!;
                            
          my $idx=$genome.".bwt";

        if(defined($gff)){

          $gtf=$dir."/"."gtf";           
          print $ofh "$gffread $gff -o $gtf -T\n";   
                         }

          unless( -f $idx){  
         print $ofh "$bwa index $genome\n$samtools faidx $genome\n";
                          }
           
              unless( -e $dir){
         print $ofh "mkdir $dir\n";
                          }
 
         my $P_M=$dir."/"."P_M";    
     unless(defined($Pbam) || defined($Mbam)){               
         print $ofh "$bwa mem -t $t -R \"\@RG\\tID:P\\tSM:P\" $genome $Pp1 $Pp2 |$samtools view -bS -> $P_bam\n";
         print $ofh "$bwa mem -t $t -R \"\@RG\\tID:M\\tSM:M\" $genome $Mp1 $Mp2 |$samtools view -bS -> $M_bam\n";
                              
          
         print $ofh "$samtools sort -o $P_srt_bam $P_bam\n";
         print $ofh "$samtools sort -o $M_srt_bam $M_bam\n";
         print $ofh "$samtools rmdup   $P_srt_bam $P_rm_bam\n";
         print $ofh "$samtools rmdup   $M_srt_bam $M_rm_bam\n";
         print $ofh "$samtools index   $P_rm_bam\n";
         print $ofh "$samtools index   $M_rm_bam\n";      
         print $ofh "$samtools mpileup -ugf $genome -q $aq $P_rm_bam $M_rm_bam |$bcftools call -mv ->$M_P_snp\n";
         print $ofh "$delly call -q $vq -g $genome -t DEL -o $M_P_del_bcf $P_rm_bam $M_rm_bam\n$bcftools view $M_P_del_bcf > $M_P_del_vcf\n";
         print $ofh "$delly call -q $vq -g $genome -t INS -o $M_P_ins_bcf $P_rm_bam $M_rm_bam\n$bcftools view $M_P_ins_bcf > $M_P_ins_vcf\n";
         print $ofh "cat $M_P_del_vcf $M_P_ins_vcf |grep -v \"#\" >$sv";
         print $ofh "rm $M_P_del_bcf $M_P_ins_bcf $M_P_del_bcf_csi $M_P_ins_bcf_csi\n"
                                             }else{
         print $ofh "$samtools index   $Pbam\n";
         print $ofh "$samtools index   $Mbam\n";
         print $ofh "$samtools mpileup -ugf $genome -q $aq $Pbam $Mbam |$bcftools call -mv ->$M_P_snp\n";
         print $ofh "$delly call -q $vq -g $genome -t DEL -o $M_P_del_bcf $Pbam $Mbam\n$bcftools view $M_P_del_bcf > $M_P_del_vcf\n";
         print $ofh "$delly call -q $vq -g $genome -t INS -o $M_P_ins_bcf $Pbam $Mbam\n$bcftools view $M_P_ins_bcf > $M_P_ins_vcf\n";
         print $ofh "cat $M_P_del_vcf $M_P_ins_vcf |grep -v \"#\" >$sv\n";
         print $ofh "rm $M_P_del_bcf $M_P_ins_bcf $M_P_del_bcf_csi $M_P_ins_bcf_csi\n"
                                                 }
         print $ofh "$gg -d $dep -q $sq -o $P_M $M_P_snp\n";  
         print $ofh "mkdir $anno_dir\n";  
         my $out_snv=$anno_dir."/snv";
         my $out_sv= $anno_dir."/sv";
                  
         print $ofh "$gtfpred -genePredExt $gtf $anno_dir/AT_refGene.txt\n";
         print $ofh "retrieve_seq_from_fasta.pl --format refGene --seqfile $genome $anno_dir/AT_refGene.txt --outfile  $anno_dir/AT_refGeneMrna.fa\n"; 
         print $ofh "table_annovar.pl $M_P_snp $anno_dir -buildver AT -out $out_snv -remove   -protocol refGene -operation g -nastring . -vcfinput\n";          
         print $ofh "table_annovar.pl $sv $anno_dir -buildver AT -out $out_sv -remove   -protocol refGene -operation g -nastring . -vcfinput\n";             
         close $ofh;
     
      unless(defined($s)){
         &process_cmd("sh $script");
                        }else{
         &process_cmd("cat $script");                      
                            }

         exit(0);
 
              }


    sub process_cmd {
    my ($cmd) = @_;
    print STDERR "CMD: $cmd\n";
    my $ret = system($cmd);
    if ($ret) {
        die "Error, cmd: $cmd died with ret $ret";
    }

    return;
                  }
       

	
