

package haplotype;

BEGIN {
	$VERSION = '1.0';
}


use strict;
use Getopt::Long;
use File::Basename;
use base 'Exporter';
our @EXPORT = qw(runhaplotype);

my $G_USAGE = "

Usage: bastos haplotype 

   Option: --Pbam     pre-aligned bam file from pollen parent 
           --Mbam     pre-aligned bam file from maternal parent
           --Hbam     pre-aligned bam file from high extreme pool
           --Lbam     pre-aligned bam file from high extreme pool
           --genome   reference genome [FASTA format]
           --phase2   use samtools algorithm or HAPCUT2 algorithm to assembly haplotype [default:T; T: samtools; F:HAPCUT2]   
           --oprefix  output dir name prefix 

Example:

    bsatos haplotype --Pbam P.bam --Mbam M.bam --Hbam H.bam --Lbam L.bam --genome genome.fasta --oprefix haplotype  

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
	GetOptions (
	"Pbam=s" => \$Pbam,
	"Mbam=s" => \$Mbam,
        "Hbam=s" => \$Hbam,
        "Lbam=s" => \$Lbam,
        "phase2=s" =>\$phase2,
        "genome=s" => \$genome,
        "var=s" => \$var, 
	"oprefix=s"   => \$outputPrefix,
        "s=s" => \$s,
	"help!" => \$help)
	or die("$G_USAGE");
	
	die "$G_USAGE" if ($help);
	
	die "$G_USAGE" if ((!defined ($outputPrefix)) || (!defined($genome)));
   
         my $samtools= "samtools";
         my $bcftools= "bcftools";
         my $bwa="bwa";
         my $filter="filter";
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
         my $filter_vcf="$FindBin::Bin//scripts/filter_vcf.pl";
         my $dir = "haplotype_dir";
             
        
         my $P_phase=$dir."/".basename($Pbam)."_phase";
         my $M_phase=$dir."/".basename($Mbam)."_phase";
         my $H_phase=$dir."/".basename($Hbam)."_phase";
         my $L_phase=$dir."/".basename($Lbam)."_phase";
         my $P_ref =$dir."/".basename($Pbam)."_ref";
         my $M_ref =$dir."/".basename($Mbam)."_ref";
         my $H_ref =$dir."/".basename($Hbam)."_ref";
         my $L_ref= $dir."/".basename($Lbam)."_ref";
         my $P_bl =$dir."/".basename($Pbam)."_block";
         my $M_bl =$dir."/".basename($Mbam)."_block";
         my $H_bl =$dir."/".basename($Hbam)."_block";
         my $L_bl= $dir."/".basename($Lbam)."_block";
         my $P_bl1 =$dir."/".basename($Pbam)."_block1";
         my $M_bl1 =$dir."/".basename($Mbam)."_block1";
         my $H_bl1 =$dir."/".basename($Hbam)."_block1";
         my $L_bl1= $dir."/".basename($Lbam)."_block1";
         my $merge =$dir."/"."merge.list";
         my $merge_bl=$dir."/"."merged.block";
         my $phase_res=$dir."/"."phase_res";
         my $P_out = $dir."/"."phase_P";
         my $P_srt = $dir."/"."phase_P.bed";
         my $M_out = $dir."/"."phase_M";
         my $M_srt = $dir."/"."phase_M.bed";
         my $overlap = $dir."/"."overlapped.bed";
       
         my $script = $outputPrefix.".sh";
   
         open (my $ofh, ">$script") or die $!;

          unless(-e $dir){
           print $ofh "mkdir $dir\n";
                        }
        unless(defined($phase2)){ 
                  
              unless(defined($var)){
                           $var=$dir."/"."var";
         print $ofh "$samtools mpileup -ugf $genome $Pbam $Mbam $Hbam $Lbam |$bcftools call -mv -|$filter_vcf ->$var\n";
                               }
            
                      
         print $ofh "$samtools phase $Pbam >$P_phase\n";
         print $ofh "$samtools phase $Mbam >$M_phase\n";
         print $ofh "$samtools phase $Hbam >$H_phase\n";
         print $ofh "$samtools phase $Lbam >$L_phase\n";
         print $ofh "$get_ref  $P_phase $genome >$P_ref\n";
         print $ofh "$get_ref  $M_phase $genome >$M_ref\n";
         print $ofh "$get_ref  $H_phase $genome >$H_ref\n";
         print $ofh "$get_ref  $L_phase $genome >$L_ref\n";
         print $ofh "$get_block  $P_phase $P_ref | $sort -> $P_bl\n";
         print $ofh "$get_block  $M_phase $M_ref | $sort -> $M_bl\n";
         print $ofh "$get_block  $H_phase $H_ref | $sort -> $H_bl\n";
         print $ofh "$get_block  $L_phase $L_ref | $sort -> $L_bl\n";
                }else{
           
          unless(defined($var)){
                           $var=$dir."/"."var";
         print $ofh "$samtools mpileup -ugf $genome $Pbam $Mbam $Hbam $Lbam |$bcftools call -mv  |$filter_vcf ->$var\n";
                
                               }

                             $var1=$dir."/"."var1";
         print $ofh "$filter_vcf $var > $var1\n";
         print $ofh "$step1 --bam $Pbam   --VCF $var1 --out $P_phase\n";
         print $ofh "$step1 --bam $Mbam   --VCF $var1 --out $M_phase\n";
         print $ofh "$step1 --bam $Hbam   --VCF $var1 --out $H_phase\n";
         print $ofh "$step1 --bam $Lbam   --VCF $var1 --out $L_phase\n";
         print $ofh "$step2 --fragments $P_phase  --VCF $var1 --output $P_bl1\n";
         print $ofh "$step2 --fragments $M_phase  --VCF $var1 --output $M_bl1\n";
         print $ofh "$step2 --fragments $H_phase  --VCF $var1 --output $H_bl1\n";
         print $ofh "$step2 --fragments $L_phase  --VCF $var1 --output $L_bl1\n";
         print $ofh "$get_b  $P_bl1 | $sort -> $P_bl\n";
         print $ofh "$get_b  $M_bl1 | $sort -> $M_bl\n";
         print $ofh "$get_b  $H_bl1 | $sort -> $H_bl\n";
         print $ofh "$get_b  $L_bl1 | $sort -> $L_bl\n";
            
                    }
                       
         print $ofh "$get_list  --Pb $P_bl --Mb $M_bl --Hb $H_bl  --Lb $L_bl >$merge\n";
         print $ofh "$filter -k cn -A A -B E  -C E -D E -E E  $merge  $P_bl  $M_bl $H_bl $L_bl |cut -f  1-2,5-9,14-16,21-23,28-30 -> $merge_bl\n"; 
         print $ofh "$fill_homo $merge_bl $var1  |$cmp1 - |$cmp2 -|$cmp1 - |$cmp2 - |$cmp1 -> $phase_res\n";  
         print $ofh "$sep $phase_res $P_out $M_out\n";
         print $ofh  "sort -k 1,1 -k 2,2n $P_out > $P_srt\nsort -k 1,1 -k 2,2n $M_out  > $M_srt\n$bedops -i $P_srt $M_srt > $overlap\n";    
          
              
        

  
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

   sub local_phase {
               

                
                   }	

  sub cmp_phase{
      
      
               }
 
  sub merge_hp{
   
              }
