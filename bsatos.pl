#!/usr/bin/perl

use strict;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/Perlib";
use prep;
use prepar;
use afd;
use phase;

my $G_USAGE = "
$0 <command> -h
 
 -- prepar  prepare the parents data   
 -- prep    prepare the input data
 -- phase   construct haplotype block 
 -- afd     calculate and filter allele frequency difference between two extreme pools [ED/g value/SNP index] and smooth the afd/ed/si/g value
 -- qtl_pick judge and pick up QTLs from three types of peaks  
 -- narrow_peak narrow down the peaks based on the multiple G value peaks and assign RDS to regions
 -- multi   select genes based on multi-omics data 
 -- qh      estimate QTL heritability based on AFD
 -- script   write a bash-script for all the data processing

bsatos (bulked segregant analysis tools for outbreeding species) Version: 1.0.1; Author: Fei Shen (flynn123\@cau.edu.cn);

";

my $command = undef;
if (!defined $ARGV[0] || substr($ARGV[0],0,1) eq '-' || substr($ARGV[0],0,2) eq '--') {
        die("Please specify the command.\n",$G_USAGE);
}
$command = shift @ARGV; $command = lc $command;

# auto-flush for both stderr and stdout
 select(STDERR);
 $| = 1;
 select(STDOUT);
 $| = 1;

if('prepar' eq $command){
         runprepar();
                       }
if ('prep' eq $command) {
        runprep();
} elsif ('phase' eq $command) {
         runphase();
} elsif ('qtl_pick' eq $command){
                                                       
} elsif ('afd' eq $command) {
         runafd();
} elsif ('ed' eq $command) {
        runed();
} elsif ('si' eq $command) {
        runsi();
} elsif ('smooth' eq $command){
        runsmooth();
} elsif ('narrow_peak narrow' eq $command){
        runnarrow_peak narrow();
} elsif ( 'multi' eq $command){
        runmulti();
} elsif  ('qh' eq $command) {
        runqh();
} else {
  
    print $G_USAGE;
       }
 
      
  exit 0;
       
















                                                           
