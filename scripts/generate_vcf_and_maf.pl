#!/usr/bin/perl

=head1 NAME
 
   generate_vcf_and_maf.pl    ------ the script to generate VCF and MAF files based on the variation and annotation file 

=head1 Usage

  perl generate_vcf_and_maf.pl -var <variation type>  -header <yes>

   -h: help information
   -var: variation type snv or sv [snv]
   -header: add header information (yes/no) [yes]
        
=head1 Example

=cut

   use warnings;
   use strict;
 
   use Getopt::Long;
   use Cwd;
   my ($h,$hap,$header,$var);
        GetOptions(
                 
                  "var=s" =>\$var,    
                
                  "header=s" =>\$header,      
                  "help"=>\$h);
 
 die `pod2text $0` if ( @ARGV==0 || $h);

  
                      
      unless(defined($var)){
               $var="snv";
                           }
  
      unless(defined($header)){
               $header="yes";
                              }        



  open IN,"<$ARGV[0]";
  my $vcf_out=$ARGV[0].".igv.vcf"; 
  my $mut_out=$ARGV[0].".igv.mut";
  open VO,">$vcf_out";
  open VV,">$mut_out";

   my @k;
   my $c=0;
    
  while(<IN>){
          chomp; 
         $c++;
    if($c==1 && $header eq "yes" ){
    if($var eq "snv"){
         print VO "\#\#fileformat=VCFv4.2\n\#\#fileDate=20100501\n\#\#\#reference\n\#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP_SNV\tM_SNV\n";            
         print VV "chr\tstart\tend\tsample\ttype\n";
           
                     }
    if($var eq "sv"){
         print VO "\#\#fileformat=VCFv4.2\n\#\#fileDate=20100501\n\#\#\#reference\n\#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP_SV\tM_SV\n";      
         print VV "chr\tstart\tend\tsample\ttype\n";          
                    }
             }else{          
      if($var eq "snv"){  
                    @k=split/\t/,$_;   
                    print VO "$k[13]\t$k[14]\t$k[15]\t$k[16]\t$k[17]\t$k[18]\t$k[19]\t$k[20]\t$k[21]\t$k[22]\t$k[23]\n"; 
                    print VV "$k[0]\t$k[1]\t$k[2]\tparents\t$k[8]\n";        
                       }
      if($var eq "sv"){
                    @k=split/\t/,$_;
                    print VO "$k[13]\t$k[14]\t$k[15]\t$k[16]\t$k[17]\t$k[18]\t$k[19]\t$k[20]\t$k[21]\t$k[22]\t$k[23]\n"; 
                     print VV "$k[0]\t$k[1]\t$k[2]\tparents\t$k[8]\n";   
                       }          
                   }
            } 
  

    










     
