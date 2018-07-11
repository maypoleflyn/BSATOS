

package prepar;

BEGIN {
	$VERSION = '1.0';
}


use strict;
use Getopt::Long;
use base 'Exporter';
our @EXPORT = qw(runprepar);

my $G_USAGE = "

Usage:  bastos prepar [options] 

Options:  --Pp1          paired-end1 fastq file from pollen parent  [FASTQ format]
          --Pp2          paired-end2 fastq file from pollen parent  [FASTQ format]
          --Mp1          paired-end1 fastq file from maternal parent [FASTQ format]
          --Mp2          paired-end2 fastq file from maternal parent [FASTQ format]
          --Pbam         pre-aligned bam file from pollen parent  [Not required, if --Pp1 & --Pp2]
          --Mbam         pre-aligned bam file from maternal parent [Not required, if --Mp1 & --Mp2]  
          --genome       reference genome [FASTA format]
          --oprefix      output dir name prefix   
          --gtf          gene gtf file [GTF format]Example:

Example:
 
 1) bsatos prepar --Pp1 P_1.fastq --Pp2 P_2.fastq --Mp1 M_1.fastq --Mp2 M_fastq --genome genome.fasta --gtf gene.gtf --prefix prepar
 OR
 2) bastos prepar --Pbam P.bam --Mbam M.bam --genome genome.fasta --gtf gene.gtf --prefix prepar 

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
	GetOptions (
	"Pp1=s" => \$Pp1,
	"Pp2=s" => \$Pp2,
        "Mp1=s" => \$Mp1,
        "Mp2=s" => \$Mp2,
        "Pbam=s" => \$Pbam,
        "Mbam=s" => \$Mbam,
        "genome=s" => \$genome,
	"oprefix=s"   => \$outputPrefix,
        "gtf=s" => \$gtf,
        "s=s" =>\$s, 
	"help!" => \$help)
	or die("$G_USAGE");
	
	die "$G_USAGE" if ($help);
	
	die "$G_USAGE" if (!defined $outputPrefix);
           
         my $samtools= "samtools";
         my $bwa="bwa";
         my $bcftools="bcftools";
         my $delly ="delly";
         my $gg="$FindBin::Bin/scripts/filter_get_genotye.pl";
         my $gtfpred="gtfToGenePred";
         
       
     
         my $dir="prepar_dir";
         my $anno_dir=$dir."/anno";

         my $P_bam=$dir."/".$outputPrefix."_P.bam";
         my $P_srt_bam=$dir."/".$outputPrefix."_P_srt.bam";
         my $P_rm_bam=$dir."/".$outputPrefix."_P_rm.bam";
               
         my $M_bam = $dir."/".$outputPrefix."_M.bam";
         my $M_srt_bam = $dir."/".$outputPrefix."_M_srt.bam";
         my $M_rm_bam = $dir."/".$outputPrefix."_M_rm.bam";
         my $M_P_snp=$dir."/"."M_P.snp";    
         my $M_P_del_bcf=$dir."/"."M_P.del.bcf";
         my $M_P_del_vcf=$dir."/"."M_P.del.vcf";
         my $M_P_ins_bcf=$dir."/"."M_P.ins.bcf";
         my $M_P_ins_vcf=$dir."/"."M_P.ins.vcf";
         my $sv=$dir."/"."sv.vcf";       
         
         
       

        my $script = $outputPrefix.".sh";
       
    
        open (my $ofh, ">$script") or die $!;
                            
          my $idx=$genome.".bwt";
          unless( -f $idx){  
         print $ofh "$bwa index $genome\n$samtools faidx $genome\n";
                          }
           
              unless( -e $dir){
         print $ofh "mkdir $dir\n";
                          }
 
         my $P_M=$dir."/"."P_M";    
     unless(defined($Pbam) || defined($Mbam)){               
         print $ofh "$bwa mem -t 5 -R \"\@RG\\tID:P\\tSM:P\" $genome $Pp1 $Pp2 |$samtools view -bS -> $P_bam\n";
         print $ofh "$bwa mem -t 5 -R \"\@RG\\tID:M\\tSM:M\" $genome $Mp1 $Mp2 |$samtools view -bS -> $M_bam\n";
                              
          
         print $ofh "$samtools sort -o $P_srt_bam $P_bam\n";
         print $ofh "$samtools sort -o $M_srt_bam $M_bam\n";
         print $ofh "$samtools rmdup   $P_srt_bam $P_rm_bam\n";
         print $ofh "$samtools rmdup   $M_srt_bam $M_rm_bam\n";
         print $ofh "$samtools index   $P_rm_bam\n";
         print $ofh "$samtools index   $M_rm_bam\n";      
         print $ofh "$samtools mpileup -ugf $genome -q 30 $P_rm_bam $M_rm_bam |$bcftools call -mv ->$M_P_snp\n";
         print $ofh "$delly call -g $genome -t DEL -o $M_P_del_bcf $P_rm_bam $M_rm_bam\n$bcftools view $M_P_del_bcf > $M_P_del_vcf\n";
         print $ofh "$delly call -g $genome -t INS -o $M_P_ins_bcf $P_rm_bam $M_rm_bam\n$bcftools view $M_P_ins_bcf > $M_P_ins_vcf\n";
         print $ofh "cat $M_P_del_vcf $M_P_ins_vcf |grep -v \"#\" >$sv";
                                             }else{
         print $ofh "$samtools index   $Pbam\n";
         print $ofh "$samtools index   $Mbam\n";
         print $ofh "$samtools mpileup -ugf $genome -q 30 $Pbam $Mbam |$bcftools call -mv ->$M_P_snp\n";
         print $ofh "$delly call -g $genome -t DEL -o $M_P_del_bcf $Pbam $Mbam\n$bcftools view $M_P_del_bcf > $M_P_del_vcf\n";
         print $ofh "$delly call -g $genome -t INS -o $M_P_ins_bcf $Pbam $Mbam\n$bcftools view $M_P_ins_bcf > $M_P_ins_vcf\n";
         print $ofh "cat $M_P_del_vcf $M_P_ins_vcf |grep -v \"#\" >$sv\n";
                                                 }
         print $ofh "$gg -o $P_M $M_P_snp\n";  
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
       

	
