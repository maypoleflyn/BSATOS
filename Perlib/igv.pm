

package igv;

BEGIN {
	$VERSION = '1.0';
}


use strict;
use Getopt::Long;
use base 'Exporter';
our @EXPORT = qw(runigv);

my $G_USAGE = "

Usage: bastos igv [options]
         
        
Options:  --hf1  FILE     paired-end1 fastq file from high extreme pool  [FASTQ format]
          --hf2  FILE     paired-end1 fastq file from high extreme pool  [FASTQ format]
          --lf1  FILE     paired-end1 fastq file from low extreme pool   [FASTQ format]
          --lf2  FILE     paired-end1 fastq file from low extreme pool   [FASTQ format]
          --hb   FILE     pre-aligned bam file from high extreme pool    [Not required, if --Hp1 & --Hp2] 
          --lb   FILE     pre-aligned bam file from low extreme pool     [Not required, if --Lp1 & --Lp2]
          --r    FILE     reference genome [FASTA format]
          --g1   FILE     G1 file from prepar step  
          --g2   FILE     G2 file from prepar step          
          --g3   FILE     G3 file from prepar step           
          --o    FILE     output dir name prefix [prep]



Aligment Options:

           --t  INT        number of threads [1]

SNV Calling Options:

           --aq  INT        skip alignments with mapQ smaller than INT [30]

SNV filtering Options:

           --cov INT        average sequencing coverage [30] 
           --sq  INT        skip SNVs with phred-scaled quality smaller than [30]
           --pn  INT        reads counts of minor allele should greater than INT [3]


Example:

1)
    bsatos prep --hf1 H_1.fastq --hf2 H_2.fastq --lf1 L_1.fastq --lf2 L_2.fastq --genome genome.fasta  --g1 P_M_G1 --g2 P_M_G2 --g3 P_M_G3 --prefix prep

OR

2) 
    bastos prep --hb H.bam --lb L.bam --genome genome.fasta --g1 P_M_G1 --g2 P_M_G2 --g3 P_M_G3 --prefix prep 

";

sub runigv {
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

	GetOptions (
	"hf1=s" => \$Hp1,
	"hf2=s" => \$Hp2,
        "lf1=s" => \$Lp1,
        "lf2=s" => \$Lp2,
        "hb=s" => \$Hbam,
        "lb=s" => \$Lbam,
        "r=s" => \$genome,
        "g1=s"  => \$G1,
        "g2=s"  => \$G2, 
        "g3=s"  => \$G3,
	"o=s"   => \$outputPrefix,
        "s=s" => \$s,
        "t=i" => \$t,
        "aq=i" =>\$aq,
        "cov=i" =>\$cov,
        "sq=i" =>\$sq,
        "pn=i" =>\$pn,
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
        unless(defined($sq)){
                           $sq=30;
                            }
        unless(defined($cov)){
                           $cov=30;
                             }
        unless(defined($pn)){
                           $pn=3;
                            }
      
         my $samtools= "samtools";
         my $bcftools="bcftools";
         my $bwa="bwa";
         my $get_counts="$FindBin::Bin/scripts/get_counts.pl";
         my $get_res="$FindBin::Bin/scripts/get_counts.pl";
         my $filter_dep="$FindBin::Bin/scripts/filter.pl";
         my $dir = "prep_dir";
         my $H_bam=$dir."/".$outputPrefix."_H.bam";
         my $H_srt_bam=$dir."/".$outputPrefix."_H_srt.bam";
         my $H_rm_bam=$dir."/".$outputPrefix."_H_rm.bam";
         my $H_counts_G1=$dir."/"."H_G1_counts";
         my $H_counts_G2=$dir."/"."H_G2_counts";
         my $H_counts_G3=$dir."/"."H_G3_counts";

         my $H_G1_snp = $dir."/"."H_G1_snp";
         my $H_G2_snp = $dir."/"."H_G2_snp";
         my $H_G3_snp = $dir."/"."H_G3_snp";
         my $g1 = $dir."/"."g1.res";
         my $g2 = $dir."/"."g2.res";
         my $g3 = $dir."/"."g3.res";
   
       
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
        
         my $script = $outputPrefix.".sh";
    
        open (my $ofh, ">$script") or die $!;
          my $idx=$genome.".bwt";
         unless( -e $dir){
         print $ofh "mkdir $dir\n";             
                         }
         unless( -e $idx){  
         print $ofh "$bwa index $genome\n$samtools faidx $genome\n";
                          }
         unless(defined($Hbam) || defined($Lbam)){       
         print $ofh "$bwa mem -t $t -R \"\@RG\\tID:H\\tSM:H\" $genome $Hp1 $Hp2 |$samtools view -bS -> $H_bam\n";
         print $ofh "$bwa mem -t $t -R \"\@RG\\tID:L\\tSM:L\" $genome $Lp1 $Lp2 |$samtools view -bS -> $L_bam\n";
         print $ofh "$samtools sort -o $H_srt_bam $H_bam\n";
         print $ofh "$samtools sort -o $L_srt_bam $L_bam\n";
         print $ofh "$samtools rmdup   $H_srt_bam $H_rm_bam\n";
         print $ofh "$samtools rmdup   $L_srt_bam $L_rm_bam\n";
         print $ofh "$samtools index   $L_rm_bam\n";
         print $ofh "$samtools index   $H_rm_bam\n";         
         print $ofh "$samtools mpileup -ugf $genome -q $aq -l $G1 $H_rm_bam |$bcftools call -m ->$H_G1_snp\n";      
         print $ofh "$samtools mpileup -ugf $genome -q $aq -l $G2 $H_rm_bam |$bcftools call -m ->$H_G2_snp\n";
         print $ofh "$samtools mpileup -ugf $genome -q $aq -l $G3 $H_rm_bam |$bcftools call -m ->$H_G3_snp\n";
         print $ofh "$samtools mpileup -ugf $genome -q $aq -l $G1 $L_rm_bam |$bcftools call -m ->$L_G1_snp\n";     
         print $ofh "$samtools mpileup -ugf $genome -q $aq -l $G2 $L_rm_bam |$bcftools call -m ->$L_G2_snp\n";
         print $ofh "$samtools mpileup -ugf $genome -q $aq -l $G3 $L_rm_bam |$bcftools call -m ->$L_G3_snp\n";
                                                  }else{ 
         print $ofh "$samtools mpileup -ugf $genome -q $aq -l $G1 $Hbam |$bcftools call -m ->$H_G1_snp\n";
         print $ofh "$samtools mpileup -ugf $genome -q $aq -l $G2 $Hbam |$bcftools call -m ->$H_G2_snp\n";
         print $ofh "$samtools mpileup -ugf $genome -q $aq -l $G3 $Hbam |$bcftools call -m ->$H_G3_snp\n";
         print $ofh "$samtools mpileup -ugf $genome -q $aq -l $G1 $Lbam |$bcftools call -m ->$L_G1_snp\n";
         print $ofh "$samtools mpileup -ugf $genome -q $aq -l $G2 $Lbam |$bcftools call -m ->$L_G2_snp\n";
         print $ofh "$samtools mpileup -ugf $genome -q $aq -l $G3 $Lbam |$bcftools call -m ->$L_G3_snp\n";
                                                       } 
         print $ofh "$get_counts -d $cov -q $sq  $H_G1_snp > $H_counts_G1\n"; 
         print $ofh "$get_counts -d $cov -q $sq  $H_G2_snp > $H_counts_G2\n";
         print $ofh "$get_counts -d $cov -q $sq  $H_G3_snp > $H_counts_G3\n";
         print $ofh "$get_counts -d $cov -q $sq  $L_G1_snp > $L_counts_G1\n";
         print $ofh "$get_counts -d $cov -q $sq  $L_G2_snp > $L_counts_G2\n";
         print $ofh "$get_counts -d $cov -q $sq $L_G3_snp > $L_counts_G3\n";


      print $ofh "sort -k 1,1 -k 2,2n $H_counts_G1 > $tmp1;sort -k 1,1 -k 2,2n $L_counts_G1 >$tmp2;filter -k cn $tmp1 $tmp2 | cut -f 1-4,7-8 - |$filter_dep -d $pn -> $g1\n";
      print $ofh "sort -k 1,1 -k 2,2n $H_counts_G2 > $tmp1;sort -k 1,1 -k 2,2n $L_counts_G2 >$tmp2;filter -k cn $tmp1 $tmp2 | cut -f 1-4,7-8 - |$filter_dep -d $pn ->$g2\n";   

     print $ofh "sort -k 1,1 -k 2,2n $H_counts_G3 > $tmp1;sort -k 1,1 -k 2,2n $L_counts_G3 >$tmp2;filter -k cn $tmp1 $tmp2 | cut -f 1-4,7-8 - |$filter_dep -d $pn ->$g3\n";


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

	
