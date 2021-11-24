

package igv;

BEGIN {
	$VERSION = '1.0';
}


use strict;
use Getopt::Long;
use base 'Exporter';
use Cwd;
our @EXPORT = qw(runigv);

my $G_USAGE = "

Usage: bastos igv [options]
         
        
Options:  
          --r            FILE     reference genome [FASTA format]
          --gtf/gff      FILE     gene GTF/GFF file  
          --prepar       DIR      result DIR from prepar module [prepar_dir]
          --hap          DIR      result DIR from haplotype module [haplotype_dir]
          --qtl          DIR      result DIR from qtl_pick module [qtl_pick_dir]
          --po           DIR      result DIR from polish module   [polish_dir]
          --o            STR      prefix name of the outputs [igv] 
          --log          FILE     igv log file [igv.log]



";


##################################################################


sub runigv {
        my $genome = undef;
	my $outputPrefix = undef;
	my $help = 0;
        my $prepar = undef;
        my $prep = undef;
        my $hap = undef;
        my $gtf = undef;
        my $gff = undef;
        my $s = undef;
        my $qtl = undef;
        my $po = undef;
        my $log = undef;      
 
 	GetOptions (
        "r=s" => \$genome,
        "gtf=s"  => \$gtf,
        "gff=s"  => \$gff, 
        "prepar=s" => \$prepar,
        "prep=s" =>\$prep,
	"o=s"   => \$outputPrefix,
        "hap=s" => \$hap,
        "s=s" =>\$s,
        "qtl=s" =>\$qtl,
        "po=s" =>\$po,
        "log=s" =>\$log,
	"help!" => \$help)
        
	or die("$G_USAGE");
	die "$G_USAGE" if ($help);
	die "$G_USAGE" if (!defined $outputPrefix);

#############################################################
#############################################################

      my $work_dir=getcwd;
     


       unless(defined($outputPrefix)){
                    $outputPrefix="igv";
                                     } 

       my $dir = $work_dir."/".$outputPrefix."_dir";
         
                     
       unless(defined($log)){
                    $log=$work_dir."/".$outputPrefix.".log";
                            }
       
       unless(defined($prepar)){
                     $prepar = $work_dir."/prepar_dir";
                                 }
       unless(defined($hap)){
                     $hap = $work_dir."/haplotype_dir";
                              }
       unless(defined($qtl)){
                     $qtl =$work_dir."/qtl_pick_dir";
                              }
      
##############################################################
##############################################################

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
                    if($_=~/finished/){
                             $con++;       
                                     }                                                                     
                      } 

##############################################################
##############################################################

       if($con !=0 ){                 
                 print "
             Program: bsatos (Bulked segregant analysis tools for outbreeding species)
             Version: 1.0
             The last module: igv
             igv         generate files for Integrative Genomics Viewer\n\n\n\n


             This module has be finished!!!!!\n\n\n";
                   exit(0);                            
                    }
##############################################################
              }   
##############################################################

             print STDOUT "
             Program: bsatos (Bulked segregant analysis tools for outbreeding species)
             Version: 1.0
             The last module: igv
             igv         generate files for Integrative Genomics Viewer\n\n\n\n";

             open (my $lo,">>$log") or die $!;
        
             print $lo "
             Program: bsatos (Bulked segregant analysis tools for outbreeding species)
             Version: 1.0
             The last module: igv
             igv         generate files for Integrative Genomics Viewer\n\n\n\n";

##############################################################
##############################################################
 
       unless(-e $dir){

  
          print STDOUT "create igv result dir\n"; 
          print $lo "create igv result dir\nCMD1:\t";
          unless(defined($hlog{"CMD1:"})){
                         $hlog{"CMD1:"}=0;
                                         }
          my $mx1 = $hlog{"CMD1:"}; 
          my $cd1 = "mkdir $dir";
          &process_cmd($cd1,$lo,$mx1);
                  }
###############################################################           
###############################################################
 
         my $ref=$dir."/reference.fasta";
         my $gff1=$dir."/gene.gff";
         my $gtf1=$dir."/gene.gtf";
         my $snv_vcf=$dir."/snv.vcf";
         my $snv_maf=$dir."/snv.maf";    
         my $sv_vcf=$dir."/sv.vcf"; 
         my $sv_maf=$dir."/sv.maf";
         my $haplotype_bed=$dir."/haplotye.bed";
         my $H_enriched_bed=$dir."/H_enriched.bed";
         my $L_enriched_bed=$dir."/L_enriched.bed";
         my $H_enriched_vcf=$dir."/H_enriched.vcf";
         my $L_enriched_vcf=$dir."/L_enriched.vcf";
         my $P_profile=$dir."/P.igv"; 
         my $M_profile=$dir."/M.igv";
         my $H_profile=$dir."/PM.igv";

###############################################################   
###############################################################
       
          print STDOUT "create igv genome file\n";
          print $lo "create igv genome file\nCMD2:\t";
          unless(defined($hlog{"CMD2:"})){
                         $hlog{"CMD2:"}=0;
                                         }            
          my $cd2 ="cp $genome $ref";
          my $mx2 = $hlog{"CMD2:"};
          &process_cmd($cd2,$lo,$mx2);

################################################################
################################################################

         
        if(defined($gtf)){ 
        
          print STDOUT "create igv genome gtf\n";
          print $lo "create igv genome gtf\nCMD3:\t";
          unless(defined($hlog{"CMD3:"})){
                         $hlog{"CMD3:"}=0;
                                         }
          my $cd3 = "cp $gtf $gtf1";
          my $mx3 = $hlog{"CMD3:"};
          &process_cmd($cd3,$lo,$mx3);

                         }

################################################################           
################################################################



        if(defined($gff)){ 

          print STDOUT "create igv genome gff\n";
          print $lo "create igv genome gff\nCMD4:\t";
          unless(defined($hlog{"CMD4:"})){
                         $hlog{"CMD4:"}=0;
                                         }
          my $cd4 = "cp $gff $gff1";
          my $mx4 = $hlog{"CMD4:"};
          &process_cmd($cd4,$lo,$mx4);

                         }

################################################################
################################################################

          print STDOUT "copy pollen parent profiles\n";
          print $lo "copy pollen parent profiles\nCMD4:\t";
          unless(defined($hlog{"CMD5:"})){
                         $hlog{"CMD5:"}=0;
                                         }
          my $cd5 = "cp $po\/p.igv $P_profile\n"; 
          my $mx5 = $hlog{"CMD5:"};
          &process_cmd($cd5,$lo,$mx5);

#################################################################
#################################################################

          print STDOUT "copy maternal parent profiles\n";
          print $lo "copy maternal parent profiles\nCMD6:\t";
          unless(defined($hlog{"CMD6:"})){
                         $hlog{"CMD6:"}=0;
                                         }
          my $cd6 = "cp $po\/m.igv $M_profile\n";    
          my $mx6 = $hlog{"CMD6:"};
          &process_cmd($cd6,$lo,$mx6);

#################################################################
################################################################# 
 
          print STDOUT "copy pollen and maternal parent profiles\n";
          print $lo "copy pollen and maternal parent profiles\nCMD7:\t";
          unless(defined($hlog{"CMD7:"})){
                         $hlog{"CMD7:"}=0;
                                         }
          my $cd7 = "cp $po\/pm.igv $H_profile\n";
          my $mx7 = $hlog{"CMD6:"};
          &process_cmd($cd7,$lo,$mx7);

##################################################################
##################################################################              

          print STDOUT "copy and merge IGV files form qtl_pick_dir\n";
          print $lo "copy and merge IGV files from qtl_pick_dir\nCMD8:\t";
          unless(defined($hlog{"CMD8:"})){
                         $hlog{"CMD8:"}=0;
                                         }
          my $cd8 = "cp $qtl\/*igv* $dir\n";
          my $mx8 = $hlog{"CMD8:"};
          &process_cmd($cd8,$lo,$mx8);
  
####################################################################
####################################################################
             
          print STDOUT "copy and merge IGV files form qtl_pick_dir\n";
          print $lo "copy and merge IGV files from qtl_pick_dir\nCMD9:\t";
          unless(defined($hlog{"CMD9:"})){
                         $hlog{"CMD9:"}=0;
                                         }
          my $cd9 = "cat $qtl\/*snv.igv.vcf >>$snv_vcf\n";
          my $mx9 = $hlog{"CMD9:"};
          &process_cmd($cd9,$lo,$mx9);

#######################################################################


          print STDOUT "copy and merge IGV files form qtl_pick_dir\n";
          print $lo "copy and merge IGV files from qtl_pick_dir\nCMD10:\t";
          unless(defined($hlog{"CMD10:"})){
                         $hlog{"CMD10:"}=0;
                                         }
          my $cd10 = "cat $qtl\/*sv.igv.vcf >>$sv_vcf\n";
          my $mx10 = $hlog{"CMD10:"};
          &process_cmd($cd10,$lo,$mx10);

#########################################################################

          print STDOUT "copy and merge IGV files form qtl_pick_dir\n";
          print $lo "copy and merge IGV files from qtl_pick_dir\nCMD11:\t";
          unless(defined($hlog{"CMD11:"})){
                         $hlog{"CMD11:"}=0;
                                         }
          my $cd11 = "cat $qtl\/*sv.igv.mut >>$sv_maf\n";
          my $mx11 = $hlog{"CMD11:"};
          &process_cmd($cd11,$lo,$mx11);

##########################################################################

          print STDOUT "copy and merge IGV files form qtl_pick_dir\n";
          print $lo "copy and merge IGV files from qtl_pick_dir\nCMD12:\t";
          unless(defined($hlog{"CMD12:"})){
                         $hlog{"CMD12:"}=0;
                                         }
          my $cd12 = "cat $qtl\/*snv.igv.mut >>$snv_maf\n";
          my $mx12 = $hlog{"CMD12:"};
          &process_cmd($cd12,$lo,$mx12);

############################################################


          print STDOUT "IGV module finished!!!!\n";
          print $lo "IGV module finished!!!!\n";




################################################################################




##################################################################
################################################################## 

         exit(0);
 
              }


##################################################################  	


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

####################################################

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
