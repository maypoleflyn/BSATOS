

package qtl_pick;

BEGIN {
	$VERSION = '1.0';
}


use strict;
use Getopt::Long;
use base 'Exporter';
use File::Basename;
our @EXPORT = qw(runqtl_pick);

my $G_USAGE = "

Usage: bastos qtl_pick [options]   
     
     Options: --oprefix          output dir name prefix
              --g1               G1 file from polish step 
              --g2               G2 file from polish step
              --g3               G3 file from polish step
              --v                annotated SNVs file from prepar step
              --sv               annotated SVs file from prepar step
              --gtf              gene.gtf file                             
              --hap              haplotye file file from haplotype step 

Example:

  bsatos qtl_pick --oprefix qtl_pick --g1 g1.res.ad.ad --g2 g2.res.ad.ad --g3 g3.res.ad.ad -v snv.AT_multianno.txt --sv sv.AT_multianno.txt --gtf gene.gtf --hap merged.block

";
 
 
sub runqtl_pick {
        my $g1 = undef;
        my $g2 = undef;
        my $g3 = undef;
	my $v = undef;
        my $sv = undef;
        my $gtf = undef;
	my $outputPrefix = undef;
	my $help = 0;
	my $s = undef;
        my $hap = undef;
	GetOptions (
	"hap=s"=>\$hap,
        "g1=s" =>\$g1,
        "g2=s" =>\$g2,
        "g3=s" =>\$g3,
        "gtf=s" => \$gtf, 
        "v=s" => \$v, 
        "s=s" => \$s,
        "sv=s" => \$sv,
	"oprefix=s"   => \$outputPrefix,
	"help!" => \$help)
	or die("Error in command line arguments\n$G_USAGE");
	
	die "$G_USAGE" if ($help);
	
	die "Please specify the output prefix\n$G_USAGE" if (!defined $outputPrefix);
             
       
            
                 
         my $script = $outputPrefix.".sh";
         my $dir = "$FindBin::Bin/qtl_pick_dir";
         my $gn1 = basename($g1);
         my $gn2 = basename($g2);
         my $gn3 = basename($g3); 
         
         my $sum1=$dir."/".$gn1.".cal.out";
         my $sum2=$dir."/".$gn2.".cal.out";
         my $sum3=$dir."/".$gn3.".cal.out";
         my $ad1=$dir."/".$gn1.".ad"; 
         my $ad2=$dir."/".$gn2.".ad";
         my $ad3=$dir."/".$gn3.".ad";       

         my $filter = "filter";        
         my $plot = "$FindBin::Bin/R/plot.R";
         my $peak = "$FindBin::Bin/R/peak.R";
         my $thres = "$FindBin::Bin/R/thres.R";
         my $main_peak="$FindBin::Bin/scripts/main_peak.pl";     
         my $assign="$FindBin::Bin/scripts/assign_enriched_haplotype.pl"; 
         my @t1 =split/\n/,`Rscript $thres $g1`;
         my @t2 =split/\n/,`Rscript $thres $g2`;
         my @t3 =split/\n/,`Rscript $thres $g3`;
         my $gtf2bed="$FindBin::Bin/scripts/gtf2gene_bed.pl";
         my $getbed="$FindBin::Bin/scripts/get.bed.pl"; 
         my $filter_pa="$FindBin::Bin/scripts/filter_snv_based_on_qtl.pl"; 
         my $filter_sv="$FindBin::Bin/scripts/filter_sv_based_on_qtl.pl"; 
         my $subplot = "$FindBin::Bin/R/sub_plot.R";  

      open (my $ofh, ">$script") or die $!;
                    
           my @chr=&get_chr($g1);
           my @chrom;
           my $pl;     

           unless(-e $dir){
           print $ofh "mkdir $dir\n";
                          }
           print $ofh "perl $assign $g1 $t1[0] $t1[1] $t1[2] $t1[3]  $hap > $dir/g1_hap\n";
           print $ofh "perl $assign $g2 $t2[0] $t2[1] $t2[2] $t2[3]  $hap > $dir/g2_hap\n"; 
           print $ofh "perl $assign $g3 $t3[0] $t3[1] $t3[2] $t3[3]  $hap > $dir/g3_hap\n";         
           print $ofh "perl $gtf2bed $gtf |sort -k 1,1 -k 2,2n -  > $dir/gene.bed\n";
          
                  
           foreach my $cl(@chr){  
                  $pl=$cl.".pdf";
                  push@chrom,$pl;     
           print $ofh "Rscript $plot $g1 $g2 $g3 $cl $t1[0] $t2[0] $t3[0]\n";      
           print $ofh "perl $main_peak --th1 $t1[0] --th2 $t2[0] --th3 $t3[0] --g1 $g1 --g2 $g2 --g3 $g3 >> $dir/qtl\n";       
                               }

           open LJ,"<$dir\/qtl";
           my @lk;
           my $out_name;         
           while(<LJ>){
                  chomp;
               @lk=split/\t/,$_; 
               my $out_gene=$dir."/".$lk[6].".gene";
               my $out_snv=$dir."/".$lk[6].".snv";
               my $out_sv=$dir."/".$lk[6].".sv";
               my $out_hap=$dir."/".$lk[6].".hap";
  
               
            print $ofh "perl $getbed $dir/gene.bed $lk[1] $lk[2] $lk[3] >$out_gene\n";        
             
            if($lk[0] eq "P"){  
            print $ofh "perl $getbed $v $lk[1] $lk[2] $lk[3] |perl $filter_pa --i $lk[0] --hap $dir/g1_hap - > $out_snv\n";
            print $ofh "perl $getbed $dir/g1_hap $lk[1] $lk[2] $lk[3] > $out_hap\n"; 
            print $ofh "perl $getbed $sv $lk[1] $lk[2] $lk[3] |perl $filter_sv --i $lk[0] - > $out_sv\n"; 
            print $ofh "Rscript $subplot $g1 $lk[1] $lk[2] $lk[3] $lk[6]\n";  
                             }
            if($lk[0] eq "M"){
            print $ofh "perl $getbed $v $lk[1] $lk[2] $lk[3] |perl $filter_pa --i $lk[0] --hap $dir/g2_hap - > $out_snv\n";
            print $ofh "perl $getbed $dir/g2_hap $lk[1] $lk[2] $lk[3] > $out_hap\n"; 
            print $ofh "perl $getbed $sv $lk[1] $lk[2] $lk[3] |perl $filter_sv --i $lk[0] - > $out_sv\n";
            print $ofh "Rscript $subplot $g2 $lk[1] $lk[2] $lk[3] $lk[6]\n";
                             }
            if($lk[0] eq "H"){
            print $ofh "perl $getbed $v $lk[1] $lk[2] $lk[3] |perl $filter_pa --i $lk[0] --hap $dir/g3_hap - > $out_snv\n";
            print $ofh "perl $getbed $dir/g3_hap $lk[1] $lk[2] $lk[3] > $out_hap\n";
            print $ofh "perl $getbed $sv $lk[1] $lk[2] $lk[3] |perl $filter_sv --i $lk[0] - > $out_sv\n";
            print $ofh "Rscript $subplot $g3 $lk[1] $lk[2] $lk[3] $lk[6]\n";
                             }
             

                       }
             
        
           print $ofh "mv @chrom *pdf  $dir\n";
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

 
  
 sub  get_chr {
        my @a=@_;
        my %h;
        my @e;
        my @k;                    
        open LOD,"<$a[0]";
       
        while(<LOD>){
                  chomp;
              @e=split/\t/,$_;
              $h{$e[0]}++;
                    }
             @k=keys%h;
          return(@k);
                   
             }
