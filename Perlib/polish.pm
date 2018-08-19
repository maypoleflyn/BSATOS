

package polish;

BEGIN {
	$VERSION = '1.0';
}


use strict;
use Getopt::Long;
use base 'Exporter';
use File::Basename;
use Cwd;
our @EXPORT = qw(runpolish);

my $G_USAGE = "

Usage: bastos polish 
 
    Options: 
             --o    STR       output dir name prefix [polish]
             --gs1  FILE      smoothed curve base on G1 type loci across genome with haplotype information from afd module 
             --gs2  FILE      smoothed curve base on G2 type loci across genome with haplotype information from afd module
             --gs3  FILE      smoothed curve base on G3 type loci across genome with haplotype information from afd module
             --h    FILE      merged haplotype block file from haplotype module 


Statistics Options:

      --sd  STR   the statistic method: ED/g/si [g]    
      --w   INT   the sliding window size [1000000] 
      --th  STR   Threshold selection method [N] 
                                           1) N: Estimate parameters of the log-normal null-distribution. 
                                           2) Y: medium plus 3-times standard deviation   
      --m STR      Smooth curve using multi-window size. Y: YES; N: NO [Y] 



Outputs:

     g1.res.ad.polish [FILE]  read counts with different alleles in G1 type loci after polish 
     g2.res.ad.polish [FILE]  read counts with different alleles in G2 type loci after polish
     g3.res.ad.polish [FILE]  read counts with different alleles in G3 type loci after polish

     *_g1.res.ad_afd  [FILE]  smoothed curve with different window in each chromosome based on polished read counts of G1 type loci
     *_g2.res.ad_afd  [FILE]  smoothed curve with different window in each chromosome based on polished read counts of G2 type loci
     *_g3.res.ad_afd  [FILE]  smoothed curve with different window in each chromosome based on polished read counts of G3 type loci

     g1.res.ad.cal.out [FILE] smoothed curve across genome based on polished read counts of G1 type loci 
     g2.res.ad.cal.out [FILE] smoothed curve across genome based on polished read counts of G2 type loci
     g3.res.ad.cal.out [FILE] smoothed curve across genome based on polished read counts of G3 type loci
      
     g1.res.ad.ad  [FILE] smoothed curve with haplotype information based on polished read counts of G1 type loci      
     g2.res.ad.ad  [FILE] smoothed curve with haplotype information based on polished read counts of G2 type loci        
     g3.res.ad.ad  [FILE] smoothed curve with haplotype information based on polished read counts of G3 type loci 


Example:

bsatos polish --sd g --gs1 g1.res.ad --gs2 g2.res.ad --gs3 g3.res.ad --h merged.block --o polish
                 
";
  
sub runpolish {
        my $g1 = undef;
        my $g2 = undef;
        my $g3 = undef;
        my $phase = undef;
	my $d = undef;
	my $outputPrefix = undef;
	my $help = 0;
	my $s = undef;
        my $t = undef;
        my $w = undef;
        my $m = undef;
        my $md = undef;
	GetOptions (
	"sd=s" =>\$d,
        "gs1=s" =>\$g1,
        "gs2=s" =>\$g2,
        "gs3=s" =>\$g3,
        "h=s" => \$phase, 
	"o=s" => \$outputPrefix,
        "s=s" => \$s,
        "w=i" =>\$w,
        "m=s" =>\$m,
        "th=s" =>\$t, 
	"help!" => \$help)
	or die("$G_USAGE");	
	die "$G_USAGE" if ($help);
	die "$G_USAGE" if (!defined $outputPrefix);
             
      unless(defined($w)){
           $w=1000000;
                           }
        unless(defined($m)){
           $m= "Y";   
                           }
        if($m eq "Y"){
               $md =1;
                     }else{
               $md =0;
                          }
   

          my $cal;  
          my $polish="$FindBin::Bin/scripts/polish_qtl_region.pl";
          my $thes="$FindBin::Bin/R/thres.R"; 
     
        if(!defined($d) || $d eq "g" ){         
          $cal = "$FindBin::Bin/R/g.R";  
        
                                         }
       if($d eq "ed"){  
            $cal = "$FindBin::Bin/R/ed.R";
                     }
       if($d eq "si"){
           $cal = "$FindBin::Bin/R/si.R";            
                     } 
                 
         my $script = $outputPrefix.".sh";
         my $work_dir=getcwd;
         my $dir = $work_dir."/polish_dir";
         my $gn1 = basename($g1);
         my $gn2 = basename($g2);
         my $gn3 = basename($g3); 
         my @t1 =split/\n/,`Rscript $thes $g1`;
         my @t2 =split/\n/,`Rscript $thes $g2`;
         my @t3 =split/\n/,`Rscript $thes $g3`;
        
 
         my $sum1=$dir."/".$gn1.".cal.out";
         my $sum2=$dir."/".$gn2.".cal.out";
         my $sum3=$dir."/".$gn3.".cal.out";
         my $ad1=$dir."/".$gn1.".ad"; 
         my $ad2=$dir."/".$gn2.".ad";
         my $ad3=$dir."/".$gn3.".ad";       
         my $filter="filter";        
         my $polish1=$dir."/".$gn1.".polish";
         my $polish2=$dir."/".$gn2.".polish";       
         my $polish3=$dir."/".$gn3.".polish";   
           
              
            
      open (my $ofh, ">$script") or die $!;
                    
                my @chr=&get_chr($g1);
       

        unless(-e $dir){        
        print $ofh "mkdir $dir\n";  
                      }

        print $ofh "$polish $g1 $t1[0] $t1[1] $t1[2] $t1[3] >$polish1\n";
        print $ofh "$polish $g2 $t2[0] $t2[1] $t2[2] $t2[3] >$polish2\n";
        print $ofh "$polish $g3 $t3[0] $t3[1] $t3[2] $t3[3] >$polish3\n";
            
           foreach my $cl(@chr){                  
           print $ofh "Rscript $cal $polish1 $cl $gn1 $w $md\n"; 
           print $ofh "Rscript $cal $polish2 $cl $gn2 $w $md\n";
           print $ofh "Rscript $cal $polish3 $cl $gn3 $w $md\n";                              
                               }
           print $ofh "mkdir $dir\n";
           print $ofh "cat \*$gn1\_afd >>$sum1\n";
           print $ofh "cat \*$gn2\_afd >>$sum2\n";
           print $ofh "cat \*$gn3\_afd >>$sum3\n";
           print $ofh "$filter -k cn -A A -B E $sum1 $phase > $ad1\n";
           print $ofh "$filter -k cn -A A -B E $sum2 $phase > $ad2\n";
           print $ofh "$filter -k cn -A A -B E $sum3 $phase > $ad3\n";     
                 
 
           print $ofh "mv \*afd $dir\n";  
            
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
