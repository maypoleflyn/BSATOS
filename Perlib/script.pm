

package script;

BEGIN {
	$VERSION = '1.0';
}


use strict;
use Getopt::Long;
use base 'Exporter';
use File::Basename;
use FindBin qw($Bin);

our @EXPORT = qw(runscript);

my $G_USAGE = "
 
Usage: bsatos script [options]  

   
 Options: 
         --oprefix: the prefix of the result folder
         --genome:  the genome fasta file
         --gtf:     the gtf file of the genes
         --pf1:     the paired-end1 fastq file of the pollen parent    | |    | --pb  :  BAM file of the pollen parent
         --pf2:     the paired-end2 fastq file of the pollen parent    | |    |
         --mf1:     the paired-end1 fastq file of the maternal parent  | |    | --mb  :  BAM file of the maternal parent 
         --mf2:     the paired-end2 fastq file of the maternal parent  | | OR |     
         --hf1:     the paired-end1 fastq file of the H pool reads     | |    | --hb  :  BAM file of the H pool 
         --hf2:     the paired-end2 fastq file of the H pool reads     | |    |
         --lf1:     the paired-end1 fastq file of the L pool reads     | |    | --lb  :  BAM file of the L pool
         --lf2:     the paired-end2 fastq file of the L pool reads     | |    |
 

Example:

1) Use reads files to run bastos
 
  bastos script --oprefix result --genome genome.fasta --gtf gene.gtf --pf1 P_1.fastq.gz --pf2 P_2.fastq.gz --mf1 M_1.fastq.gz --mf2 M_2.fastq.gz --hf1 H_1.fastq.gz --pf2 H_2.fastq.gz --lf1 L_1.fastq.gz --lf2 L_2.fastq.gz > script.sh  

2) Use pre-aligned BAMs files to run bastos

 bastos script --oprefix result --genome genome.fasta --gtf gene.gtf --pb P.bam --mb M.bam --hb H.bam --lh L.bam > script.sh 
    
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
        
        my $outputPrefix = undef;
	my $help = 0;
	
	GetOptions (
	
        "genome=s" =>\$genome,
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
	"oprefix=s"   => \$outputPrefix,
	"help=s" => \$help)
         
	or die("$G_USAGE");
	
	die "$G_USAGE" if ($help || !defined($genome));
	
             
         my $bsatos ="$FindBin::Bin/bsatos.pl";        
         
         #my $cd = $outputPrefix.".cm.sh";
         my $dir = $outputPrefix."_dir";
         unless(-e $dir){
           print "mkdir $dir\n";
                        }   
       
           
      
      #    open (my $ofh, ">$script") or die $!;
     
      
        if(defined($pf1) && defined ($pf2) && defined($mf1) && defined($mf2)){
        print  "perl $bsatos prepar  --Pp1 $pf1  --Pp2 $pf2 --Mp1 $mf1 --Mp2 $mf2 --genome $genome --gtf $gtf  --oprefix prepar\n"  
                                                                             }else{
        print  "perl $bsatos prepar --Pbam $pb --Mbam $mb --genome $genome  --oprefix prepar --gtf $gtf\n"; 
                                                                                  } 
        if(defined($hf1) && defined($hf2) && defined($lf1) && defined($lf2)){           
        print  "perl $bsatos prep --Hp1 $hf1 --Hp2 $hf2 --Lp1 $lf1  --Lp2 $lf2 --genome $genome --g1 prepar_dir/P_M_G1 --g2 prepar_dir/P_M_G2 --g3 prepar_dir/P_M_G3  --oprefix prep\n";  
                                                                            }else{ 
        print  "perl $bsatos prep  --Hbam $hb --Lbam $lb --genome $genome   --g1 prepar_dir/P_M_G1 --g2 prepar_dir/P_M_G2 --g3 prepar_dir/P_M_G3  --oprefix prep\n";                                                                      
                                                   
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

        print "perl $bsatos haplotype --Pbam $pb --Mbam $mb --Hbam $hb --Lbam $lb  --genome $genome  --oprefix haplotype --var prepar_dir/M_P.snp\n";     
        print "perl $bsatos afd  --oprefix afd  -g1 prep_dir/g1.res  -g2 prep_dir/g2.res  -g3 prep_dir/g3.res  -phase haplotype_dir/merged.block\n";
        print "perl $bsatos polish  --oprefix polish  -g1 afd_dir/g1.res.ad -g2 afd_dir/g2.res.ad -g3  afd_dir/g3.res.ad -phase  haplotype_dir/merged.block\n";  
        print "perl $bsatos qtl_pick  --oprefix qtl_pick -g1 polish_dir/g1.res.ad.ad -g2 polish_dir/g2.res.ad.ad  -g3 polish_dir/g3.res.ad.ad --v prepar_dir/anno/snv.AT_multianno.txt --sv prepar_dir/anno/sv.AT_multianno.txt -gtf $gtf --hap haplotype_dir/merged.block\n"; 
         
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

 
  
