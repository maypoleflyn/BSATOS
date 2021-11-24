#!/usr/bin/perl

=head1 NAME
 
    igv_header.pl ------ the script to add IGV header to VCF and MAF files based on the variation and annotation file 

=head1 Usage

  perl igv_header.pl -file <file type> --var <variation type>

   -h: help information
   -file: variation file type [vcf or mut]
   -var: variation type [snv/sv]      

=head1 Example

=cut

   use warnings;
   use strict;
 
   use Getopt::Long;
   use Cwd;
   my ($h,$file,$var);
        GetOptions(  
                  "file=s" =>\$file,   
                  "var=s" =>\$var,
                  "help"=>\$h);
 
 die `pod2text $0` if (@ARGV==0 || $h || (!defined($file)) || (!defined($var)));
                      
    
  open IN,"<$ARGV[0]";
   
    my @k;
    my $c=0;
    
  while(<IN>){
       chomp; 
  if($_=~/#/){
          next;
             } 
       $c++;
  if($c==1){
    if($file eq "vcf"){                  
            if($var eq "snv"){                 
       print "\#\#fileformat=VCFv4.2\n\#\#fileDate=20100501\n\#\#\#reference\n\#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP_SNV\tM_SNV\n";      
              }else{
       print "\#\#fileformat=VCFv4.2\n\#\#fileDate=20100501\n\#\#\#reference\n\#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP_SV\tM_SV\n";                     }
                       }
     if($file eq "mut"){
               print "chr\tstart\tend\tsample\ttype\n";
                       }
             }else{
      unless($_=~/#/ || $_=~/type/){
          print "$_\n";
                    }
                  }
             }                         
 











