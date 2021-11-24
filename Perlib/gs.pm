

package gs;

BEGIN {
	$VERSION = '1.0';
}


use strict;
use Getopt::Long;
use base 'Exporter';
use Cwd;
our @EXPORT = qw(rungs);

my $G_USAGE = "

Usage: bastos gs [options]
         
        
Options:  
          --gen  FILE     Genotype file
          --phe  FILE     Phenotype file
          --rou  INT      The rounds of calculate [1000]
          --pre FLOAT    The trainning set in all the populations [0.6]
          --o    STR      outputPrefix           
Outputs:

 gs_dir [DR]



";

sub rungs {
	my $gen = undef;
 	my $phe = undef;
        my $rou = undef;
        my $pre = undef;
	my $outputPrefix = undef;
	my $help = 0;
        my $log = undef;
 
	GetOptions (
	"gen=s" => \$gen,
	"phe=s" => \$phe,
        "pre=f" => \$pre,
        "rou=i" => \$rou,
	"o=s"   => \$outputPrefix,
	"help!" => \$help)
        
	or die("$G_USAGE");
	
	die "$G_USAGE" if ($help);
	
	die "$G_USAGE" if (!defined $outputPrefix);

##################################################################

      unless(defined($outputPrefix)){
                        $outputPrefix="gs";
                                    }

        unless(defined($rou)){
                           $rou=1000;
                           }
        unless(defined($pre)){
                           $pre=0.6;
                            }
        
      
#####################################################################

         my $gs2valid="$FindBin::Bin/scripts/gs2valid.pl"; 
         my $pr="$FindBin::Bin/scripts/pearson.R";
         my $work_dir=getcwd;
         my $dir = $work_dir."/".$outputPrefix."_dir"; 
         my $acc= $dir."/accurancy";
      unless( -e $dir){
         system("mkdir $dir");
                      }

         system("cp $pr $work_dir");
 
         foreach my $index (1 .. $rou){            
                        system("perl $gs2valid -p $pre $gen $phe > $acc ");   
                                      }
              system("mv EPV_file marker_effect.txt pearson.R $dir");        

                  }
