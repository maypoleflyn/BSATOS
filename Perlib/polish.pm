

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
             --gs1  FILE      smoothed curve base on P type loci across genome with haplotype information from afd module 
             --gs2  FILE      smoothed curve base on M type loci across genome with haplotype information from afd module
             --gs3  FILE      smoothed curve base on PM type loci across genome with haplotype information from afd module
             --h    FILE      merged haplotype block file from haplotype module 
             --log  FILE      polish module log file [polish.log]
             --fdr  INT       FDR threshold for the polishing [0.01]

Statistics Options:

      --sd  STR   the statistic method: ED/g/si [g]    
      --w   INT   the sliding window size [1000000] 
      --fn  INT   batches for smoothing; the smaller the faster but more memory needed; [20]  
 

Outputs:

      
     P.polished.afd   [FILE] smoothed curve with haplotype information based on polished read counts of P type loci      
     M.polished.afd   [FILE] smoothed curve with haplotype information based on polished read counts of M type loci        
     PM.polished.afd  [FILE] smoothed curve with haplotype information based on polished read counts of PM type loci 

     p.igv        [FILE] P IGV format file for Integrative Genomics Viewer    
     m.igv        [FILE] M IGV format file for Integrative Genomics Viewer
     pm.igv       [FILE] PM IGV format file for Integrative Genomics Viewer


Example:

bsatos polish --sd g --gs1 P.AFD --gs2 M.AFD --gs3 PM.AFD --h haplotype.block --o polish
                 
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
        my $fn = undef;
        my $fdr = undef;
		my $log = undef;
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
             "fdr=f" =>\$fdr, 
             "fn=i" =>\$fn,
	     "log=s" =>\$log, 
	     "help!" => \$help)
	or die("$G_USAGE");	
	die "$G_USAGE" if ($help);
	die "$G_USAGE" if (!defined $outputPrefix);

##############################################################
	
       unless(defined($outputPrefix)){
                    $outputPrefix="polish";
                                     }
       unless(defined($fdr)){
                     $fdr="0.01";
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


############################################

         my $work_dir=getcwd;
         my $dir = $work_dir."/".$outputPrefix."_dir";


        unless(defined($log)){
                           
         $log=$work_dir."/".$outputPrefix.".log";
                             
                             }
##############################################
      my %hlog=();       
 
       if(-e $log){
        
           open LOD,"<$log"; 
           my $con=0;  
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

         if($con !=0 ){                 
                 print "
             Program: bsatos (Bulked segregant analysis tools for outbreeding species)
             Version: 1.0
             The fourth module: polish
             polish      polish candidate QTLs region and remove nosiy makers based on haplotype information\n\n\n\n


             This module has be finished!!!!!\n\n\n";
                      
              exit(0);       
                       }
                       }
  						  
################################################################
          my $cal;  
          my $polish="$FindBin::Bin/scripts/polish_qtl_region.pl";
          my $bsa2igv="$FindBin::Bin/scripts/bsa2igv.pl";
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
                 
        
 
         

         my $gn1 = basename($g1);
         my $gn2 = basename($g2);
         my $gn3 = basename($g3); 
    
       
        my $g1_th = $work_dir."/".$outputPrefix."_dir"."/"."g1.thres";
        my $g2_th = $work_dir."/".$outputPrefix."_dir"."/"."g2.thres";
        my $g3_th = $work_dir."/".$outputPrefix."_dir"."/"."g3.thres";

        
 
         my $sum1=$dir."/".$gn1.".cal.out";
         my $sum2=$dir."/".$gn2.".cal.out";
         my $sum3=$dir."/".$gn3.".cal.out";
         my $ad1=$dir."/"."P.polished.afd"; 
         my $ad2=$dir."/"."M.polished.afd";
         my $ad3=$dir."/"."PM.polished.afd";       
         my $filter="filter";        
         my $polish1=$dir."/".$gn1.".polish";
         my $polish2=$dir."/".$gn2.".polish";       
         my $polish3=$dir."/".$gn3.".polish";   
         my $igv1=$dir."/p.igv";
         my $igv2=$dir."/m.igv";
         my $igv3=$dir."/pm.igv";

         my $add_header= "$FindBin::Bin/scripts/add_afd_header.pl";
              
############################################################
            
      open (my $lo, ">$log") or die $!;

 
      print  STDOUT "
     
             Program: bsatos (Bulked segregant analysis tools for outbreeding species)
             Version: 1.0
             The fifth module: polish
             polish      polish candidate QTLs region and remove nosiy makers based on haplotype information\n\n\n\n";
                

     
        
      print $lo "
             
             Program: bsatos (Bulked segregant analysis tools for outbreeding species)
             Version: 1.0
             The fifth module: polish
             polish      polish candidate QTLs region and remove nosiy makers based on haplotype information\n\n\n\n";

			 
			 
#####################################################################################			 
            
     my @chr=&get_chr($g1);
 
########################################################
 
   unless(-e $dir){        
 
        print STDOUT "create result dir\n";
        print $lo "create result dir\nCMD0:\t";
        unless(defined($hlog{"CMD0:"})){
                         $hlog{"CMD0:"}=0;
                                         }
        my $mx0=$hlog{"CMD0:"};
        my $cd1 =  "mkdir $dir";
        &process_cmd($cd1,$lo,$mx0);  

                  }
#############################################################

        print STDOUT "calculate threshold g1\n";
        print $lo "calculate g1 threshold\nCMD B1:\t";
        my $cdb1 = "Rscript $thes $g1 $fdr $g1_th";
                unless(defined($hlog{"CMD B1:"})){
                         $hlog{"CMD B1:"}=0;
                                       }
         my $mxb1=$hlog{"CMD B1:"};
        &process_cmd($cdb1,$lo,$mxb1);

################################################################
#
        print STDOUT "calculate threshold g2\n";
        print $lo "calculate g2 threshold\nCMD B2:\t";
        my $cdb2 = "Rscript $thes $g2 $fdr $g2_th";
                unless(defined($hlog{"CMD B2:"})){
                         $hlog{"CMD B2:"}=0;
                                       }
         my $mxb2=$hlog{"CMD B2:"};
        &process_cmd($cdb2,$lo,$mxb2);

#####################################################################

        print STDOUT "calculate threshold g3\n";
        print $lo "calculate g1 threshold\nCMD B3:\t";
        my $cdb3 = "Rscript $thes $g3 $fdr $g3_th";
                unless(defined($hlog{"CMD B3:"})){
                         $hlog{"CMD B1:"}=0;
                                       }
         my $mxb3=$hlog{"CMD B3:"};
        &process_cmd($cdb3,$lo,$mxb3);

##################################################################################
        my @t1 =split/\n/,`cat $g1_th`;
        my @t2 =split/\n/,`cat $g2_th`;
        my @t3 =split/\n/,`cat $g3_th`;


####################################################################
            
        print STDOUT "polish G1 loci data\n";
        print $lo "polish G1 loci data\nCMD1:\t";         
        my $cd2 = "$polish $g1 $t1[0] $t1[1] $t1[2] $t1[3] >$polish1";
		unless(defined($hlog{"CMD1:"})){
                         $hlog{"CMD1:"}=0;
                                       }
         my $mx1=$hlog{"CMD1:"};	
        &process_cmd($cd2,$lo,$mx1);
      
#######################################################################	  
	  
	  
        print STDOUT "polish G2 loci data\n";
        print $lo "polish G2 loci data\nCMD2:\t";
        my $cd3 = "$polish $g2 $t2[0] $t2[1] $t2[2] $t2[3] >$polish2";
		unless(defined($hlog{"CMD2:"})){
                         $hlog{"CMD2:"}=0;
                                       }
        my $mx2=$hlog{"CMD2:"};	
        &process_cmd($cd3,$lo,$mx2);
      
	 
############################################################################
 
        print STDOUT "polish G3 loci data\n";
        print $lo "polish G3 loci data\nCMD3:\t";
        my $cd4= "$polish $g3 $t3[0] $t3[1] $t3[2] $t3[3] >$polish3";
		
		unless(defined($hlog{"CMD3:"})){
                         $hlog{"CMD3:"}=0;
                                       }
        my $mx3=$hlog{"CMD3:"};	
		
        &process_cmd($cd4,$lo,$mx3);
         

#############################################################################		 
        
		my $idx =3;
		
     foreach my $cl(@chr){   

###################################
	 
           $idx++;
           my $cdn1="CMD".$idx.":";		   
           print STDOUT "smooth G1 loci in $cl\n";
           print $lo "smooth the G1 loci in $cl\n$cdn1\t";
           unless(defined($hlog{$cdn1})){
                         $hlog{$cdn1}=0;
                                       }
           my $mxc1=$hlog{$cdn1};	
           my $g1c = "Rscript $cal $polish1 $cl $gn1 $w $md $fn";
           my $g1c2 = "Rscript $cal $polish1 $cl $gn1 $w 0 $fn";
           &process_cmd2($g1c,$lo,$mxc1,$g1c2); 
		   
#################################################
		   
           $idx++;
           my $cdn2= "CMD".$idx.":";
           print STDOUT "smooth G2 loci in $cl\n";
           print $lo "smooth G2 loci in $cl\n$cdn2\t"; 
           unless(defined($hlog{$cdn2})){
                         $hlog{$cdn2}=0;
                                       }
           my $mxc2=$hlog{$cdn2};
           my $g2c = "Rscript $cal $polish2 $cl $gn2 $w $md $fn";
           my $g2c2 = "Rscript $cal $polish2 $cl $gn2 $w 0 $fn";
           &process_cmd2($g2c,$lo,$mxc2,$g2c2);
  
#######################################################  
  
           $idx++;
           my $cdn3= "CMD".$idx.":";
           print STDOUT "smooth G3 loci in $cl\n";
           print $lo "smooth G3 loci in $cl\n$cdn3\t";
           unless(defined($hlog{$cdn3})){
                         $hlog{$cdn3}=0;
                                        }
           my $mxc3=$hlog{$cdn3};
           my $g3c = "Rscript $cal $polish3 $cl $gn3 $w $md $fn";
           my $g3c2 = "Rscript $cal $polish3 $cl $gn3 $w 0 $fn";
           &process_cmd2($g3c,$lo,$mxc3,$g3c2);
                            
                         }

           $idx++;
           my $cdn4= "CMD".$idx.":";
           print STDOUT "merge smoothed G1 loci\n";
           print $lo "merge smoothed G1 loci\n$cdn4\t";   
           unless(defined($hlog{$cdn4})){
                         $hlog{$cdn4}=0;
                                        }
           my $mxc4=$hlog{$cdn4};
           my $m1c = "cat $work_dir/\*_$gn1\_afd >>$sum1";
           &process_cmd($m1c,$lo,$mxc4);


           $idx++;
           my $cdn5= "CMD".$idx.":";
           print STDOUT "merge smoothed G2 loci\n";
           print $lo "merge smoothed G2 loci\n$cdn5\t";
           unless(defined($hlog{$cdn5})){
                         $hlog{$cdn5}=0;
                                        }
           my $mxc5=$hlog{$cdn5};
           my $m2c = "cat $work_dir/\*_$gn2\_afd >>$sum2";
           &process_cmd($m2c,$lo,$mxc5);

           $idx++;
           my $cdn6= "CMD".$idx.":";
           print STDOUT "merge smoothed G3 loci\n";
           print $lo "merge smoothed G3 loci\n$cdn6\t";
           unless(defined($hlog{$cdn6})){
                         $hlog{$cdn6}=0;
                                        }
           my $mxc6=$hlog{$cdn6};
           my $m3c = "cat $work_dir/\*_$gn3\_afd >>$sum3";
           &process_cmd($m3c,$lo,$mxc6);

           $idx++;
           my $cdn7= "CMD".$idx.":";
           print STDOUT "adding haplotype block information to smoothed G1 loci\n";
           print $lo "adding haplotype block information to smoothed G1 loci\n$cdn7\t";
           unless(defined($hlog{$cdn7})){
                         $hlog{$cdn7}=0;
                                        }
           my $mxc7=$hlog{$cdn7};
           my $adc1 = "$filter -k cn -A A -B E $sum1 $phase  |$add_header -> $ad1";
           &process_cmd($adc1,$lo,$mxc7);
           
           $idx++;
           my $cdn8= "CMD".$idx.":";
           print STDOUT "adding haplotype information to smoothed G2 loci\n";
           print $lo "adding haplotype information to smoothed G2 loci\n$cdn8\t";
           my $adc2 = "$filter -k cn -A A -B E $sum2 $phase |$add_header -> $ad2";
           unless(defined($hlog{$cdn8})){
                         $hlog{$cdn8}=0;
                                        }
           my $mxc8=$hlog{$cdn8};
           &process_cmd($adc2,$lo,$mxc8);
        
           $idx++; 
           my $cdn9= "CMD".$idx.":";
           print STDOUT "adding haplotype information to smoothed G3 loci\n";
           print $lo "adding haplotype information to smooothed G3 loci\n$cdn9\t";  
           unless(defined($hlog{$cdn9})){
                         $hlog{$cdn9}=0;
                                        }
           my $mxc9=$hlog{$cdn9};
           my $adc3 = "$filter -k cn -A A -B E $sum3 $phase |$add_header -> $ad3";
           &process_cmd($adc3,$lo);
        
           $idx++;
           my $cdn10= "CMD".$idx.":";
           print STDOUT "convert G1 BSA FILE to IGV format FILE\n";
           print $lo "convert G1 BSA FILE to IGV format FILE\n$cdn10\t";
           unless(defined($hlog{$cdn10})){
                         $hlog{$cdn10}=0;
                                        }
           my $mxc10=$hlog{$cdn10};
           my $igvc1 = "$bsa2igv -i P $ad1 > $igv1";
           &process_cmd($igvc1,$lo,$mxc10);
            
           $idx++;
           my $cdn11= "CMD".$idx.":";
           print STDOUT "convert G2 BSA FILE to IGV format FILE\n";
           print $lo "convert G2 BSA FILE to IGV format FILE\nCMD$cdn11\t";   
           unless(defined($hlog{$cdn11})){
                         $hlog{$cdn11}=0;
                                        }
           my $mxc11=$hlog{$cdn11};  
           my $igvc2 = "$bsa2igv -i M $ad2 > $igv2";
           &process_cmd($igvc2,$lo,$mxc11);

    
           $idx++;
           my $cdn12= "CMD".$idx.":";
           print STDOUT "convert G3 BSA FILE to IGV format FILE\n";
           print $lo "convert G3 BSA FILE to IGV format FILE\n$cdn12\t";
           unless(defined($hlog{$cdn12})){
                         $hlog{$cdn12}=0;
                                        }
           my $mxc12=$hlog{$cdn12};
           my $igvc3 = "$bsa2igv -i PM $ad3 > $igv3";
           &process_cmd($igvc3,$lo,$mxc12);
                   

           $idx++;  
           my $cdn13= "CMD".$idx.":";
           print STDOUT "move all the files to result dir\n";
           print $lo "move all the files to result dir\nCMD$cdn13\t";
           unless(defined($hlog{$cdn13})){
                         $hlog{$cdn13}=0;
                                        }
           my $mxc13=$hlog{$cdn13};
           my $mc = "rm $work_dir/\*AFD_afd  $sum1 $sum2 $sum3 $polish1 $polish2 $polish3 $g1_th $g2_th $g3_th";
           &process_cmd($mc,$lo,$mxc13);  
      
           print STDOUT "polish module finished!!!\n";
           print $lo "polish module finished!!!\n";          
            



              }


  
####################################################################

sub process_cmd {
    my ($cmd,$lk,$hg) = @_;
    my $time="[".localtime()."]";
#######################
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
  
 ############################################################
  
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
#################################################################
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
