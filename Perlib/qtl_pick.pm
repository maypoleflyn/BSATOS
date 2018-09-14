

package qtl_pick;

BEGIN {
	$VERSION = '1.0';
}


use strict;
use Getopt::Long;
use base 'Exporter';
use File::Basename;
use Cwd;
our @EXPORT = qw(runqtl_pick);

my $G_USAGE = "

Usage: bastos qtl_pick [options]   
     
     
     Options: --o        STR        output dir name prefix
              --gp1      FILE       smoothed curve base on G1 type loci across genome with haplotype information from polish module 
              --gp2      FILE       smoothed curve base on G2 type loci across genome with haplotype information from polish module
              --gp3      FILE       smoothed curve base on G3 type loci across genome with haplotype information from polish module
              --v        FILE       annotated SNVs file from prepar step
              --sv       FILE       annotated SVs file from prepar step
              --gtf      FILE       gene.gtf file                             
              --h        FILE       haplotye file file from haplotype step 
              --q        INT        mininum phred-scaled quality score [30]         
              --pr       INT        promoter region [2000]
              --per      FLOAT      minimum percentage of blocks in QTLs region [0.7]
              --bl       INT        minimum enriched block length [2] 
              --log      FILE       qtl_pick module log file [qtl_pick.log]
              --fdr      FLOAT      FDR threshold for the qtl [0.01]    
Outputs:

       qtl     [FILE]         detected QTLs list file 
       *.pdf   [FILE]         G value/ED/SI profiles across each chromosome (*:chromosome) 
       g1_hap  [FILE]         haplotype information in G1 type loci 
       g2_hap  [FILE]         haplotype information in G2 type loci
       g2_hap  [FILE]         haplotype information in G3 type loci
       gene.bed [FILE]        gene bed file 
       *.gene   [FILE]        gene list located in the QTL regions (*: QTL accession)
       *.hap  [FILE]          haplotype information located in the QTL regions (*:QTL accession)
       *.snv [FILE]           screened SNVs based on genetic rules located in the QTL regions (*: QTL accession)  
       *.sv [FIE]             screened SNVs based on genetic rules located in the QTL regions （*： QTL accession)  
    


Example:

  bsatos qtl_pick --o qtl_pick --gp1 g1.res.ad.ad --gp2 g2.res.ad.ad --gp3 g3.res.ad.ad -v snv.AT_multianno.txt --sv sv.AT_multianno.txt --gtf gene.gtf --h merged.block

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
        my $pr = undef;
        my $qu = undef;
        my $per = undef;
        my $bl = undef;      
        my $log = undef;
        my $fdr = undef;
	GetOptions (
	"h=s"=>\$hap,
        "gp1=s" =>\$g1,
        "gp2=s" =>\$g2,
        "gp3=s" =>\$g3,
        "gtf=s" => \$gtf, 
        "v=s" => \$v, 
        "s=s" => \$s,
        "pr=i" =>\$pr,
        "q=i" =>\$qu,
        "sv=s" => \$sv,
        "fdr=f" =>\$fdr,
	"o=s"   => \$outputPrefix,
        "per=f" =>\$per,
        "bl=i" =>\$bl,
        "log=s" =>\$log,
	"help!" => \$help)
	or die("$G_USAGE");
	
	die "$G_USAGE" if ($help);
	
	die "$G_USAGE" if (!defined $outputPrefix);
             

#################################################

        unless(defined($outputPrefix)){
                     $outputPrefix="qtl_pick";
                                      }

 
         unless(defined($fdr)){
 
                     $fdr=0.01;
                             
                             }







        unless(defined($pr)){
                      $pr=2000;
                            }
        unless(defined($qu)){
                      $qu=30;
                           }
        unless(defined($per)){
                       $per=0.7;
                             } 

        unless(defined($bl)){
                       $bl=2;
                            }
         
                                  
         my $work_dir=getcwd;
         my $dir = $work_dir."/".$outputPrefix."_dir";
         

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
        # my @t1 =split/\n/,`Rscript $thres $g1`;
        # my @t2 =split/\n/,`Rscript $thres $g2`;
        # my @t3 =split/\n/,`Rscript $thres $g3`;
         my $gtf2bed="$FindBin::Bin/scripts/gtf2gene_bed.pl";
         my $getbed="$FindBin::Bin/scripts/get.bed.pl"; 
         my $filter_pa="$FindBin::Bin/scripts/filter_snv_based_on_qtl.pl"; 
         my $g_igv="$FindBin::Bin/scripts/generate_vcf_and_maf.pl";
         my $filter_sv="$FindBin::Bin/scripts/filter_sv_based_on_qtl.pl"; 

         my $subplot = "$FindBin::Bin/R/sub_plot.R";  

############################################################################

        my $g1_th = $work_dir."/".$outputPrefix."_dir"."/"."g1.thres";
        my $g2_th = $work_dir."/".$outputPrefix."_dir"."/"."g2.thres";
        my $g3_th = $work_dir."/".$outputPrefix."_dir"."/"."g3.thres";


##########################################################################


         unless(defined($log)){
              $log=$work_dir."/".$outputPrefix.".log";
                             }

############################################################################

             my %hlog=();
   

           if(-e $log){
 
               open LOD,"<$log";
               my $con=0;
 ##############################################
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

################################################    
                  if($con !=0 ){
                 
                 print "

             Program: bsatos (Bulked segregant analysis tools for outbreeding species)
             Version: 1.0
             The sixth module: qtl_pick 
             qtl_pick    judge and pick up QTLs from three types of peaks\n\n\n\n


             This module has be processed last time\n\n\n";
                       
                   exit(0);           
                            }
 #####################################
                        }
           
 #############################################

         print STDOUT "

             Program: bsatos (Bulked segregant analysis tools for outbreeding species)
             Version: 1.0
             The sixth module: qtl_pick 
             qtl_pick    judge and pick up QTLs from three types of peaks\n\n\n\n";

        
         
    
        open (my $lo, ">>$log") or die $!;


        print $lo "

             Program: bsatos (Bulked segregant analysis tools for outbreeding species)
             Version: 1.0
             The sixth module: qtl_pick 
             qtl_pick    judge and pick up QTLs from three types of peaks\n\n\n\n";

###########################################################################

                    
           my @chr=&get_chr($g1);
           my @chrom;
           my $pl;     

           unless(-e $dir){
    
           print STDOUT "create result dir\n";
           print $lo "create result dir\nCMD0:\t";
           unless(defined($hlog{"CMD0:"})){                
                        $hlog{"CMD0:"}=0;
                                          }  
           my $mx0=$hlog{"CMD0:"};       
           my $cd0 = "mkdir $dir";
           &process_cmd($cd0,$lo,$mx0);
         
                          }
          
########################################################################
        print STDOUT "calculate threshold g1\n";
        print $lo "calculate g1 threshold\nCMD B1:\t";
        my $cdb1 = "Rscript $thres $g1 $fdr $g1_th";
                unless(defined($hlog{"CMD B1:"})){
                         $hlog{"CMD B1:"}=0;
                                       }
         my $mxb1=$hlog{"CMD B1:"};
        &process_cmd($cdb1,$lo,$mxb1);
###########################################################################
        print STDOUT "calculate threshold g2\n";
        print $lo "calculate g2 threshold\nCMD B2:\t";
        my $cdb2 = "Rscript $thres $g2 $fdr $g2_th";
                unless(defined($hlog{"CMD B2:"})){
                         $hlog{"CMD B2:"}=0;
                                       }
         my $mxb2=$hlog{"CMD B2:"};
        &process_cmd($cdb2,$lo,$mxb2);
##########################################################################
        print STDOUT "calculate threshold g3\n";
        print $lo "calculate g1 threshold\nCMD B3:\t";
        my $cdb3 = "Rscript $thres $g3 $fdr $g3_th";
                unless(defined($hlog{"CMD B3:"})){
                         $hlog{"CMD B1:"}=0;
                                       }
         my $mxb3=$hlog{"CMD B3:"};
        &process_cmd($cdb3,$lo,$mxb3);
###################################################################33
        my @t1 =split/\n/,`cat $g1_th`;
        my @t2 =split/\n/,`cat $g2_th`;
        my @t3 =split/\n/,`cat $g3_th`;
################################################################3

############################################################################

      
           print STDOUT "assign enriched haplotype block in G1 loci\n";
           print $lo "assign enriched haplotype block in G1 loci\nCMD1:\t"; 
           my $cd1="perl $assign --per $per --block $bl $g1 $t1[0] $t1[1] $t1[2] $t1[3]  $hap > $dir/g1_hap";
           unless(defined($hlog{"CMD1:"})){                
                        $hlog{"CMD1:"}=0;
                                         }  
           my $mx1=$hlog{"CMD1:"};       

           &process_cmd($cd1,$lo,$mx1);



#######################################################################################


 
           print STDOUT "assign enriched haplotype block in G2 loci\n";
           print $lo "assign enriched haplotype block in G2 loci\nCMD2:\t";
           my $cd2 = "perl $assign --per $per --block $bl $g2 $t2[0] $t2[1] $t2[2] $t2[3]  $hap > $dir/g2_hap";
           unless(defined($hlog{"CMD2:"})){                
                        $hlog{"CMD2:"}=0;
                                         }  
           my $mx2=$hlog{"CMD2:"};       

           &process_cmd($cd2,$lo,$mx2);



#####################################################################################

           print STDOUT "assign enriched haplotype block in G3 loci\n";
           print $lo "assign enriched haplotype block in G3 loci\nCMD3:\t";
           my $cd3 = "perl $assign --per $per --block $bl  $g3 $t3[0] $t3[1] $t3[2] $t3[3]  $hap > $dir/g3_hap"; 
           unless(defined($hlog{"CMD3:"})){                
                        $hlog{"CMD3:"}=0;
                                         }  
           my $mx3=$hlog{"CMD3:"};       

           &process_cmd($cd3,$lo,$mx3);
         

#################################################################################################

           print STDOUT "convert gene GTF FILE to BED FILE\n";
           print $lo "convert gene GTF FILE to BED FILE\nCMD4:\t";
           my $cd4 = "perl $gtf2bed $gtf |sort -k 1,1 -k 2,2n -  > $dir/gene.bed";
           unless(defined($hlog{"CMD4:"})){                
                        $hlog{"CMD4:"}=0;
                                          }  
           my $mx4=$hlog{"CMD4:"};       

           &process_cmd($cd4,$lo,$mx4);


######################################################################################### 
        
           my $idx=4;        
      foreach my $cl(@chr){  
                  $pl=$dir."/".$cl.".pdf";
                  push@chrom,$pl;    
  
#######################################################################
           $idx++;
           my $i1="CMD".$idx.":";

           print STDOUT "plot AFD profiles in $cl\n";
           print $lo "plot AFD profiles in $cl\n$i1\t";
           unless(defined($hlog{$i1})){                
                        $hlog{$i1}=0;
                                          }  
           my $mx5=$hlog{$i1};       
           my $cdc = "Rscript $plot $g1 $g2 $g3 $cl $t1[0] $t2[0] $t3[0]";
           &process_cmd($cdc,$lo,$mx5);

####################################################################
             
           $idx++;
           my $i2="CMD".$idx.":";
           print STDOUT "obtain QTL peak region in $cl\n";
           print $lo "obtain QTL peak region in $cl\n$i2\t";
           unless(defined($hlog{$i2})){             
                        $hlog{$i2}=0;
                                       }  
           my $mx6=$hlog{$i2};                        
           my $cdd = "perl $main_peak --th1 $t1[0] --th2 $t2[0] --th3 $t3[0] --g1 $g1 --g2 $g2 --g3 $g3 >> $dir/qtl";
           &process_cmd($cdd,$lo,$mx6);
           
####################################################################################
                               }


           open LJ,"<$dir\/qtl";
           my @lk;
           my $out_name;
           my $qi=0; 
                   
        while(<LJ>){
                  chomp;
               @lk=split/\t/,$_; 
               my $out_gene=$dir."/".$lk[6].".gene";
               my $out_snv=$dir."/".$lk[6].".snv";
               my $out_sv=$dir."/".$lk[6].".sv";
               my $out_hap=$dir."/".$lk[6].".hap";
               $qi++;
               $idx++;

            my $i3="CMD".$idx.":";

#####################################################

            print STDOUT "get QTL$qi region genes\n";
            print $lo "get QTL$qi region genes\n$i3\t:";
            unless(defined($hlog{$i3})){             
                        $hlog{$i3}=0;
                                       }  
            my $mx7=$hlog{$i3};

            my $co1 = "perl $getbed $dir/gene.bed $lk[1] $lk[2] $lk[3] >$out_gene";
            &process_cmd($co1,$lo,$mx7);        
###################################################################
             
            if($lk[0] eq "P"){ 
           
            $idx++;
            my $i4 = "CMD".$idx.":";
            print STDOUT "obtain and filter SNVs in QTL$qi region based on genetic laws\n";
            print $lo "obtain and filter SNVs in QTL$qi region based on genetic laws\n$i4\t";         
            my $gsnv = "perl $getbed $v $lk[1] $lk[2] $lk[3] |perl $filter_pa --d $pr --q $qu --i $lk[0] --hap $dir/g1_hap - > $out_snv";
            unless(defined($hlog{$i4})){             
                        $hlog{$i4}=0;
                                          }  
            my $mx8=$hlog{$i4};

            &process_cmd($gsnv,$lo,$mx8);


##############################################################
           
            $idx++;
            my $i5 ="CMD".$idx.":";

            print STDOUT "convert filtered SNVs file to IGV format\n";
            print $lo "convert filtered SNVs file to IGV format\n$i5\t";
            unless(defined($hlog{$i5})){             
                        $hlog{$i5}=0;
                                          }  
            my $mx9=$hlog{$i5};
            my $gigv = "perl $g_igv -var snv $out_snv";
            &process_cmd($gigv,$lo,$mx9);
  


##############################################################################

            $idx++;
            my $i6 ="CMD".$idx.":";
            print STDOUT "obtain enriched haplotype information in QTL$qi region\n";
            print $lo "obtain enriched haplotype information in QTL$qi region\n$i6\t";
            unless(defined($hlog{$i6})){            
                        $hlog{$i6}=0;
                                          }
            my $mx10=$hlog{$i6};
            my $ocd = "perl $getbed $dir/g1_hap $lk[1] $lk[2] $lk[3] > $out_hap";
            &process_cmd($ocd,$lo,$mx10);
 


##########################################################

            $idx++;
            my $i7 ="CMD".$idx.":";
            print STDOUT "obtain SVs in QTL$qi region\n";
            print $lo "obtain SVs in QTL$qi region\n$i7\t";
            unless(defined($hlog{$i7})){
                     $hlog{$i7}=0;
                                       }
            my $mx11=$hlog{$i7};
            my $svc = "perl $getbed $sv $lk[1] $lk[2] $lk[3] |perl $filter_sv --d $pr  --i $lk[0] - > $out_sv";
            &process_cmd($svc,$lo,$mx11);
 


################################################################
            $idx++;
            my $i8 ="CMD".$idx.":";
            print STDOUT "convert SVs in QTL$qi region to IGV format\n";
            print $lo "convert SVs in QTL$qi region to IGV format\n$i8\t";
            unless(defined($hlog{$i8})){
                     $hlog{$i8}=0;
                                       }
            my $mx12=$hlog{$i8};
            my $igv_sv = "perl $g_igv -var sv $out_sv";
            &process_cmd($igv_sv,$lo,$mx12);
  

############################################################################

            $idx++;
            my $i9 ="CMD".$idx.":";
            print STDOUT "plot profiles\n";
            print $lo "plot profiles\n$i9\t";
            unless(defined($hlog{$i9})){
                     $hlog{$i9}=0;
                                       }
            my $mx13=$hlog{$i9};
            my $cmdd1 = "Rscript $subplot $g1 $lk[1] $lk[2] $lk[3] $lk[6]";
            &process_cmd($cmdd1,$lo,$mx13);
  
#####################################################################
                        
               }

###########################################################


            if($lk[0] eq "M"){

#############################################################
            $idx++;
            my $i10="CMD".$idx.":";
            print STDOUT "obtain and filter SNVs in QTL$qi region based on genetic laws\n";
            print $lo "obtain and filter SNVs in QTL$qi region based on genetic laws\n$i10\t";
            my $gsnv = "perl $getbed $v $lk[1] $lk[2] $lk[3] |perl $filter_pa --d $pr --q $qu --i $lk[0] --hap $dir/g2_hap - > $out_snv";
            unless(defined($hlog{$i10})){
                     $hlog{$i10}=0;
                                       }
            my $mx14=$hlog{$i10};
          
            &process_cmd($gsnv,$lo,$mx14);

#############################################################

            $idx++;
            my $i11="CMD".$idx.":";
            print STDOUT "convert filtered SNVs file to IGV format\n";
            print $lo "convert filtered SNVs file to IGV format\n$i11\t";
            unless(defined($hlog{$i11})){
                     $hlog{$i11}=0;
                                       }
            my $mx15=$hlog{$i11};

            my $gigv = "perl $g_igv -var snv $out_snv";
            &process_cmd($gigv,$lo,$mx15);

###################################################################


            $idx++;
            my $i12 = "CMD".$idx.":";
            print STDOUT "obtain enriched haplotype information in QTL$qi region\n";
            print $lo "obtain enriched haplotype information in QTL$qi region\n$i12\t";
            unless(defined($hlog{$i12})){
                     $hlog{$i12}=0;
                                       }
            my $mx16=$hlog{$i12};
            my $ocd = "perl $getbed $dir/g2_hap $lk[1] $lk[2] $lk[3] > $out_hap";
            &process_cmd($ocd,$lo,$mx16);
            
###################################################################

            $idx++;
            my $i13="CMD".$idx.":";
            print STDOUT "obtain SVs in QTL$qi region\n";
            print $lo "obtain SVs in QTL$qi region\n$i13\t";
            unless(defined($hlog{$i13})){
                     $hlog{$i13}=0;
                                       }
            my $mx17=$hlog{$i13};
            my $svc = "perl $getbed $sv $lk[1] $lk[2] $lk[3] |perl $filter_sv --d $pr  --i $lk[0] - > $out_sv";
            &process_cmd($svc,$lo,$mx17);



########################################################################
            $idx++;
            my $i14="CMD".$idx.":";
            print STDOUT "convert SVs in QTL$qi region to IGV format\n";
            print $lo "convert SVs in QTL$qi region to IGV format\n$14\t";
            unless(defined($hlog{$i14})){
                     $hlog{$i14}=0;
                                       }
            my $mx18=$hlog{$i14};
            my $igv_sv = "perl $g_igv -var sv $out_sv";
            &process_cmd($igv_sv,$lo,$mx18);

#########################################################
            $idx++;
            my $i15 = "CMD".$idx.":";
            print STDOUT "plot profiles\n";
            print $lo "plot profiles\n$i15\t";
            unless(defined($hlog{$i15})){
                     $hlog{$i15}=0;
                                       }
            my $mx19=$hlog{$i15};

            my $cmdd1 = "Rscript $subplot $g1 $lk[1] $lk[2] $lk[3] $lk[6]";
            &process_cmd($cmdd1,$lo,$mx19);
            

################################################
                             }

            if($lk[0] eq "H"){

##################################################

            $idx++;
            my $i16 = "CMD".$idx.":";
            print STDOUT "obtain and filter SNVs in QTL$qi region based on genetic laws\n";
            print $lo "obtain and filter SNVs in QTL$qi region based on genetic laws\n$i16\t";
            unless(defined($hlog{$i16})){
                     $hlog{$i16}=0;
                                       }
            my $mx20=$hlog{$i16};
            my $gsnv = "perl $getbed $v $lk[1] $lk[2] $lk[3] |perl $filter_pa --d $pr --q $qu --i $lk[0] --hap $dir/g3_hap - > $out_snv";
            &process_cmd($gsnv,$lo,$mx20);

##################################################

            $idx++;
            my $i17 = "CMD".$idx.":";
            print STDOUT "convert filtered SNVs file to IGV format\n";
            print $lo "convert filtered SNVs file to IGV format\n$i17\t";
            unless(defined($hlog{$i17})){
                     $hlog{$i17}=0;
                                       }
            my $mx21=$hlog{$i17};
            my $gigv = "perl $g_igv -var snv $out_snv";
            &process_cmd($gigv,$lo,$mx21);

#################################################

            $idx++;
            my $i18 = "CMD".$idx.":";
            print STDOUT "obtain enriched haplotype information in QTL$qi region\n";
            print $lo "obtain enriched haplotype information in QTL$qi region\n$i18\t";
            my $ocd = "perl $getbed $dir/g3_hap $lk[1] $lk[2] $lk[3] > $out_hap";
            unless(defined($hlog{$i18})){
                     $hlog{$i18}=0;
                                       }
            my $mx22=$hlog{$i18};

            &process_cmd($ocd,$lo,$mx22);


####################################################################
            $idx++;
            my $i19 = "CMD".$idx.":";            
            print STDOUT "obtain SVs in QTL$qi region\n";
            print $lo "obtain SVs in QTL$qi region\n$i19\t";
            my $svc = "perl $getbed $sv $lk[1] $lk[2] $lk[3] |perl $filter_sv --d $pr  --i $lk[0] - > $out_sv";
            unless(defined($hlog{$i19})){
                     $hlog{$i19}=0;
                                       }
            my $mx23=$hlog{$i19};
            &process_cmd($svc,$lo,$mx23);

################################################################

            $idx++;
            my $i20 = "CMD".$idx.":"; 
            print STDOUT "convert SVs in QTL$qi region to IGV format\n";
            print $lo "convert SVs in QTL$qi region to IGV format\n$i20\t";
            my $igv_sv = "perl $g_igv -var sv $out_sv";
            unless(defined($hlog{$i20})){
                     $hlog{$i20}=0;
                                       }
            my $mx24=$hlog{$i20};
            &process_cmd($igv_sv,$lo,$mx24);

######################################################

            $idx++;
            my $i21 = "CMD".$idx.":";
            print STDOUT "plot profiles\n";
            print $lo "plot profiles\n$i21\t";
            my $cmdd1 = "Rscript $subplot $g1 $lk[1] $lk[2] $lk[3] $lk[6]";
            unless(defined($hlog{$i21})){
                     $hlog{$i21}=0;
                                       }
            my $mx25=$hlog{$i21};
            &process_cmd($cmdd1,$lo,$mx25);
 


##############################################################
                             }
             
                       }

           $idx++;  
           my $i22 = "CMD".$idx.":";
           print STDOUT "move all the files to result dir\n";
           print $lo "move all the files to result dir\nCMD$i22\t";
           unless(defined($hlog{$i22})){
                     $hlog{$i22}=0;
                                      }
           my $mx26=$hlog{$i22};
           my $mvc = "mv *pdf $dir";
            &process_cmd($mvc,$lo,$mx26);


         print STDOUT "qtl_pick module finished!!!\n";
         print $lo "qtl_pick_module finished!!!\n";
         
              }


######################################################################################

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
  
######################################################################################

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
########################################################################################### 
