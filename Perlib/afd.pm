

package afd;

BEGIN {
	$VERSION = '1.0';
}


use strict;
use Getopt::Long;
use base 'Exporter';
use File::Basename;
use Cwd;
our @EXPORT = qw(runafd);

my $G_USAGE = "

Usage:  bastos afd [options]


Options: 
           --g1 FILE        read counts with different alleles from H & L pools in G1 type loci from prep module [required]
           --g2 FILE        read counts with different alleles from H & L pools in G2 type loci from prep module [required]
           --g3 FILE        read counts with different alleles from H & L pools in G3 type loci from prep module [required]
           --h  FILE        merged, corrented and patched haplotype block file of two parents and two pools from haplotype module [required]
           --o  STR         output dir name prefix [afd]

Statistics Options:

      --sd STR    the statistic method: ED/g/si [g]    
      --w  INT    the sliding window size [1000000] 
      --th  STR    Threshold selection method [N] 
                                           1) N: Estimate parameters of the log-normal null-distribution. 
                                           2) Y: medium plus 3-times standard deviation   
      --m STR      Smooth curve using multi-window size. Y: YES; N: NO [Y] 

Outputs:

    *_g1.res_afd   [FILE]  G value (ED/SI) based G1 type loci and smoothed curve with different window in each chromosome 
    *_g2.res_afd   [FILE]  G value (ED/SI) based G2 type loci and smoothed curve with different window in each chromosome
    *_g3.res_afd   [FILE]  G value (ED/SI) based G3 type loci and smoothed curve with different window in each chromosome

    g1.res.cal.out [FILE] G value (ED/SI) based G1 type loci and smoothed curve with different window across genome
    g1.res.ad      [FILE] G value (ED/SI) based G1 type loci and smoothed curve with different window across genome with haplotype information 
  
    g2.res.cal.out [FILE] G value (ED/SI) based G2 type loci and smoothed curve with different window across genome
    g2.res.ad      [FILE] G value (ED/SI) based G2 type loci and smoothed curve with different window across genome with haplotype information 

    g3.res.cal.out [FILE] G value (ED/SI) based G3 type loci and smoothed curve with different window across genome
    g3.res.ad      [FILE] G value (ED/SI) based G3 type loci and smoothed curve with different window across genome with haplotype information 
     

Example:

   bsatos afd --sd g --g1 g1.res --g2 g2.res --g3 g3.res --h merged.block --o afd --w 1000000
 
";
 
 
sub runafd {
        my $g1 = undef;
        my $g2 = undef;
        my $g3 = undef;
        my $phase = undef;
	my $d = undef;
	my $outputPrefix = undef;
	my $help = 0;
        my $s = undef;	
        my $md = undef;
        my $w = undef;
        my $t = undef;
        my $m = undef;
	GetOptions (
	"sd=s" =>\$d,
        "g1=s" =>\$g1,
        "g2=s" =>\$g2,
        "g3=s" =>\$g3,
        "h=s" => \$phase, 
	"o=s"   => \$outputPrefix,
        "s=s" => \$s, 
        "w=i" =>\$w,
        "m=s" =>\$m,
        "th=s" =>\$t,
	"help!" => \$help)
	or die("$G_USAGE");
	
	die "$G_USAGE" if ($help);
	
	die "$G_USAGE" if (!defined $outputPrefix);
             
          my $cal;  
            
        if(!defined($d) || $d eq "g" ){         
          $cal = "$FindBin::Bin/R/g.R";  
        
                                         }
       if($d eq "ed"){  
            $cal = "$FindBin::Bin/R/ed.R";
                     }
       if($d eq "si"){
           $cal = "$FindBin::Bin/R/si.R";            
                     } 
 
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
                   
         my $script = $outputPrefix.".sh";
         my $work_dir=getcwd;
         my $dir = $work_dir."/afd_dir";
         my $gn1 = basename($g1);
         my $gn2 = basename($g2);
         my $gn3 = basename($g3); 
         
         my $sum1=$dir."/".$gn1.".cal.out";
         my $sum2=$dir."/".$gn2.".cal.out";
         my $sum3=$dir."/".$gn3.".cal.out";
         my $ad1=$dir."/".$gn1.".ad"; 
         my $ad2=$dir."/".$gn2.".ad";
         my $ad3=$dir."/".$gn3.".ad";       
         my $filter="filter";        
         
            
      open (my $ofh, ">$script") or die $!;
                    
             my @chr=&get_chr($g1);
    
           foreach my $cl(@chr){                  
           print $ofh "Rscript $cal $g1 $cl $gn1 $w $md \n"; 
           print $ofh "Rscript $cal $g2 $cl $gn2 $w $md\n";
           print $ofh "Rscript $cal $g3 $cl $gn3 $w $md\n";                              
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
