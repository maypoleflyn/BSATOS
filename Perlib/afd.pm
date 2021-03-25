

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
           --pc    FILE        read counts with different alleles from H & L pools in G1 type loci from prep module [required]
           --mc    FILE        read counts with different alleles from H & L pools in G2 type loci from prep module [required]
           --pmc   FILE        read counts with different alleles from H & L pools in G3 type loci from prep module [required]
           --hap   FILE        merged, corrented and patched haplotype block file of two parents and two pools from haplotype module [required]
           --o     STR         output dir name prefix [AFD]
           --log   FILE        afd module log [afd.log]     


Statistics Options:

      --sd  STR    the statistic method: ED/g/si [g]    
      --w   INT    the sliding window size [1000000] 
      --fn  INT    batches for smoothing ;ther smaller the faster, but more memory [20]  


Outputs:

    P.AFD       [FILE] G value (ED/SI) based P type loci and smoothed curve with different window across genome with haplotype information
    M.AFD       [FILE] G value (ED/SI) based M type loci and smoothed curve with different window across genome with haplotype information
    PM.AFD      [FILE] G value (ED/SI) based PM type loci and smoothed curve with different window across genome with haplotype information 

Example:

   bsatos afd --sd g --pc P.counts --mc M.counts --pmc PM.counts --hap haplotype.block --o AFD --w 1000000
 
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
        my $log = undef;
        my $fn = undef;

	GetOptions (
	"sd=s" =>\$d,
        "pc=s" =>\$g1,
        "mc=s" =>\$g2,
        "pmc=s" =>\$g3,
        "hap=s" => \$phase, 
	"o=s"   => \$outputPrefix,
        "s=s" => \$s, 
        "w=i" =>\$w,
        "m=s" =>\$m,
        "th=s" =>\$t,
        "fn=i" =>\$fn,
        "log=s" =>\$log,
	"help!" => \$help)
	or die("$G_USAGE");
	
	die "$G_USAGE" if ($help);
	
	die "$G_USAGE" if (!defined $outputPrefix);
             
        my $cal;  

#############################################################


       my $add_header= "$FindBin::Bin/scripts/add_afd_header.pl";

             
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

       unless(defined($fn)){
       
             $fn=20;
                           }     

################################################################                   
        
      unless(defined($outputPrefix)){
              $outputPrefix="AFD";        
                                    }
         my $work_dir=getcwd;
         my $dir = $work_dir."/".$outputPrefix."_dir";
        
     
###########################################################
        unless(defined($log)){
                           
         $log=$work_dir."/".$outputPrefix.".log";
                             
                             }

#############################################
       my %hlog=();       

       if(-e $log){
           open LOD,"<$log"; 
           my $con=0;    

########################################

           while(<LOD>){
                    chomp;
                    my @sp=split/\t/,$_;
                    if($sp[2] eq "\#done"){            
                    $hlog{$sp[0]}=1;
                                          }
                    if($sp[2] eq "\#die"){
                    $hlog{$sp[0]}=0;                           
                                         }
                    if($_=~/finish/){
                             $con++;       
                                    }
                                                                     
                      }  
#####################################          
            if($con !=0 ){                 
                print "
             Program: bsatos (Bulked segregant analysis tools for outbreeding species)
             Version: 1.0
             The third module: afd
             afd         calculate and filter allele frequency difference between two extreme pools\n\n\n\n


             This module has be finished!!!!!\n\n\n";
                  
                 exit(0);
                             
                        }
####################################

                  }

######################################################################################
            

        print STDOUT "

             Program: bsatos (Bulked segregant analysis tools for outbreeding species)
             Version: 1.0
             The fourth module: afd
             afd         calculate and filter allele frequency difference between two extreme pools\n\n\n\n";
                



        open (my $lo,">>$log") or die $!;
        
        print $lo "

             Program: bsatos (Bulked segregant analysis tools for outbreeding species)
             Version: 1.0
             The fourth module: afd
             afd         calculate and filter allele frequency difference between two extreme pools\n\n\n\n";


#################################################################################################################

       unless( -e $dir){

       
          print STDOUT "creat result dir\n";
          print $lo "creat result dir\nCMD0:\t"; 
          my $cd1="mkdir $dir"; 
          unless(defined($hlog{"CMD0"})){
                         $hlog{"CMD0"}=0;
                                         }
           my $mfd=$hlog{"CMD0"};

          &process_cmd($cd1,$lo,$mfd);    

                      }
              
###################################################

         my $gn1 = basename($g1);
         my $gn2 = basename($g2);
         my $gn3 = basename($g3); 
         
         my $sum1=$dir."/".$gn1.".cal.out";
         my $sum2=$dir."/".$gn2.".cal.out";
         my $sum3=$dir."/".$gn3.".cal.out";
        # my $ad1=$dir."/".$gn1.".ad"; 
        # my $ad2=$dir."/".$gn2.".ad";
        # my $ad3=$dir."/".$gn3.".ad";       
        
         my $ad1=$dir."/"."P.AFD";
         my $ad2=$dir."/"."M.AFD";
         my $ad3=$dir."/"."PM.AFD";




         my $filter="filter";        
     
########################################################
               
         my @chr=&get_chr($g1);
         my $idx=0;
 
        foreach my $cl(@chr){         
           $idx++;
           my $cdn1 = "CMD".$idx.":";

########################################################
              
           print STDOUT "calculate P allele frequency difference between two extreme pools\n";
           print $lo "calculate P allele frequency difference between two extreme pools\t$cdn1\n";
           unless(defined($hlog{$cdn1})){
                         $hlog{$cdn1}=0;
                                         }
           my $mx1=$hlog{$cdn1}; 

           my $g1c = "Rscript $cal $g1 $cl $gn1 $w $md $fn"; 
           my $g1c2 = "Rscript $cal $g1 $cl $gn1 $w 0 $fn";
           &process_cmd2($g1c,$lo,$mx1,$g1c2);
          
###########################################################

           $idx++;
           my $cdn2 = "CMD".$idx.":";
           print STDOUT "calculate M allele frequency difference between two extreme pools\n";
           print $lo "calculate M allele frequency difference between two extreme pools\t$cdn2\n";
           unless(defined($hlog{$cdn2})){
                         $hlog{$cdn2}=0;
                                         }
           my $mx2=$hlog{$cdn2};
           my $g2c = "Rscript $cal $g2 $cl $gn2 $w $md $fn";
           my $g2c2 = "Rscript $cal $g2 $cl $gn2 $w 0 $fn";
           &process_cmd2($g2c,$lo,$mx2,$g2c2); 
   
############################################################

           $idx++;
           my $cdn3 = "CMD".$idx.":"; 
           print STDOUT "calculate PM allele frequency difference between two extreme pools\n";
           print $lo "calculate PM allele frequency difference between two extreme pools\t$cdn3\n";
           unless(defined($hlog{$cdn3})){
                         $hlog{$cdn3}=0;
                                         }
           my $mx3=$hlog{$cdn3};
           my $g3c = "Rscript $cal $g3 $cl $gn3 $w $md $fn"; 
           my $g3c2 = "Rscript $cal $g3 $cl $gn3 $w 0 $fn"; 
           &process_cmd2($g3c,$lo,$mx3,$g3c2);  
 
##########################################################

                               }

###########################################
          
           $idx++;        

           my $cdn4 = "CMD".$idx.":";   
           print STDOUT "merge G1 AFD file of different chromosomes\n";
           print $lo "merge G1 AFD file of differenr chromosomes\n$cdn4\t";
           unless(defined($hlog{$cdn4})){
                         $hlog{$cdn4}=0;
                                         }
           my $mx4=$hlog{$cdn4};
           my $m1c= "cat $work_dir/\*_$gn1\_afd >>$sum1";
           &process_cmd($m1c,$lo,$mx4);
           
####################################################

           $idx++;

           my $cdn5 = "CMD".$idx.":";
           print STDOUT "merge G2 AFD file of different chromosomes\n";
           print $lo "merge G2 AFD file of differenr chromosomes\n$cdn5\t";
           unless(defined($hlog{$cdn5})){
                         $hlog{$cdn5}=0;
                                         }
           my $mx5=$hlog{$cdn5};
           my $m2c= "cat $work_dir/\*_$gn2\_afd >>$sum2";
           &process_cmd($m2c,$lo,$mx5);

###############################################

           $idx++;
           my $cdn6 = "CMD".$idx.":";
           print STDOUT "merge G3 AFD file of different chromosomes\n";
           print $lo "merge G3 AFD file of differenr chromosomes\n$cdn6\t";
           unless(defined($hlog{$cdn6})){
                         $hlog{$cdn6}=0;
                                        }
           my $mx6=$hlog{$cdn6};
           my $m3c= "cat $work_dir/\*_$gn3\_afd >>$sum3";
           &process_cmd($m3c,$lo,$mx6);

#################################################

           $idx++;

           my $cdn7 = "CMD".$idx.":";
           print STDOUT "adding haplotype block information to G1 AFD file\n";
           print $lo "adding haplotype block information to G1 AFD file\n$cdn7\t";
           unless(defined($hlog{$cdn7})){
                         $hlog{$cdn7}=0;
                                        }
           my $mx7=$hlog{$cdn7};
           my $adc1 = "$filter -k cn -A A -B E $sum1 $phase  |$add_header -> $ad1";
           &process_cmd($adc1,$lo,$mx7);

####################################################

           $idx++;
           my $cdn8 = "CMD".$idx.":";
           print STDOUT "adding haplotype block information to G2 AFD file\n";
           print $lo "adding haplotype block information to G2 AFD file\n$cdn8\t";
           unless(defined($hlog{$cdn8})){
                         $hlog{$cdn8}=0;
                                        }
           my $mx8=$hlog{$cdn8};
           my $adc2 ="$filter -k cn -A A -B E $sum2 $phase |$add_header -> $ad2";
           &process_cmd($adc2,$lo,$mx8);

##################################################

           $idx++;
           my $cdn9 = "CMD".$idx.":";
           print STDOUT "adding haplotype block information to G3 AFD file\n";
           print $lo "adding haplotype block information to G3 AFD file\n$cdn9\t";
           unless(defined($hlog{$cdn9})){
                         $hlog{$cdn9}=0;
                                        }
           my $mx9=$hlog{$cdn9};
           my $adc3 = "$filter -k cn -A A -B E $sum3 $phase  |$add_header -> $ad3";
           &process_cmd($adc3,$lo,$mx9);     

#################################################   

           $idx++;
           my $cdn10 ="CMD".$idx.":";
           print STDOUT "move all the files to the result dir\n";
           print $lo "move all the files to the result dir\n$cdn10\t";
           unless(defined($hlog{$cdn10})){
                         $hlog{$cdn10}=0;
                                        }
           my $mx10=$hlog{$cdn10};
           my $mc = "rm $work_dir/\*counts_afd\nrm $sum1 $sum2 $sum3";
           &process_cmd($mc,$lo,$mx10);  


          print STDOUT "AFD module finished!!!\n";
          print $lo "AFD module finished!!!\n";



##################################################            
       
 
#############################
              }

########################################################################################

 sub process_cmd {
    my ($cmd,$lk,$hg) = @_;
    my $time="[".localtime()."]";
    if($hg == 1){
   

    print STDOUT "this command has be processed last time\n$cmd\t.........................................................\#done\t$time\n"; 
    print $lk "$cmd\t\#done\t$time\n";
                     }else{

    print STDOUT "$cmd\t........................................................."; 
    print $lk "$cmd\t";

    my $ret = system($cmd);
    if ($ret) {
           print $lk "\#die\t$time\n";  
           print STDOUT "\#die\t$time\n";       
           die "Error, cmd: $cmd died with ret $ret";
              }else{
        print $lk  "\#done\t$time\n";
        print STDOUT "\#done\t$time\n";
                   }

    return;
                  }                  
                  } 


##############################################################################

 sub  get_chr {
        my @a=@_;
        my %h;
        my @e;
        my @k;                    
        open LOD,"<$a[0]";
       
        while(<LOD>){
                  chomp;
           if($_=~/#/){
                    next;
                      }
              @e=split/\t/,$_;
              $h{$e[0]}++;
                    }
             @k=keys%h;
          return(@k);
                   
             }

##################################################################################


 sub process_cmd2 {
    my ($cmd,$lk,$hg,$cmd2) = @_;
    my $time="[".localtime()."]";
    if($hg == 1){


    print STDOUT "this command has be processed last time\n$cmd\t.........................................................\#done\t$time\n";
    print $lk "$cmd\t\#done\t$time\n";
                     }else{

    print STDOUT "$cmd\t.........................................................";
    print $lk "$cmd\t";

    my $ret = system($cmd);

    if ($ret) {
         my $ret2 = system($cmd2);
        print $lk  "\#fail and try other method\t$time\n";
        print STDOUT "\#fail and try other method\t$time\n";
        print STDOUT "$cmd2\t.........................................................";
        print $lk "$cmd2\t";
      
         if ($ret2) {
           
           print $lk "\#die\t$time\n";
           print STDOUT "\#die\t$time\n";
           die "Error, cmd: $cmd2 died with ret $ret";
              }else{
        print $lk  "\#done\t$time\n";
        print STDOUT "\#done\t$time\n";
          
           
              }
             }else{
        print $lk  "\#done\t$time\n";
        print STDOUT "\#done\t$time\n";
                   }

      return;
                  }
                  }



##
##      BSATOS: Bulk Segregant Analysis Tools for Outbreeding Species
##
##      author: Fei Shen <352363546@qq.com>  
##
## Copyright (c) 2019 China Agricultural University 
##
## Permission is hereby granted, free of charge, to any person obtaining a copy
## of this software and associated documentation files (the "Software"), to deal
## in the Software without restriction, including without limitation the rights
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
## copies of the Software, and to permit persons to whom the Software is
## furnished to do so, subject to the following conditions:
## 
## The above copyright notice and this permission notice shall be included in all
## copies or substantial portions of the Software.
## 
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
## IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
## AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
## LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
## OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
## SOFTWARE.
