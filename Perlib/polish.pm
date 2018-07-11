

package polish;

BEGIN {
	$VERSION = '1.0';
}


use strict;
use Getopt::Long;
use base 'Exporter';
use File::Basename;
our @EXPORT = qw(runpolish);

my $G_USAGE = "

Usage: bastos polish 
 
    Options: --d         The statistic in the pipeline: ED/G value/SNP index; [ed/g/si]; default:g
             --oprefix   output dir name prefix
             --g1        G1 file from afd step
             --g2        G2 file from afd step
             --g3        G3 file from afd step 
             --phase     phase block file from phase module 

Example:

bsatos polish --d g --g1 g1.res.ad.ad --g2 g2.res.ad.ad --g3 g3.res.ad.ad --phase merged.block --prefix polish
                 
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
	GetOptions (
	"d=s" =>\$d,
        "g1=s" =>\$g1,
        "g2=s" =>\$g2,
        "g3=s" =>\$g3,
        "phase=s" => \$phase, 
	"oprefix=s"   => \$outputPrefix,
        "s=s" => \$s,
	"help!" => \$help)
	or die("$G_USAGE");
	
	die "$G_USAGE" if ($help);
	
	die "$G_USAGE" if (!defined $outputPrefix);
             
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
         my $dir = "$FindBin::Bin/polish_dir";
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
           print $ofh "Rscript $cal $polish1 $cl $gn1\n"; 
           print $ofh "Rscript $cal $polish2 $cl $gn2\n";
           print $ofh "Rscript $cal $polish3 $cl $gn3\n";                              
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
