#!/usr/bin/perl

=head1 NAME

  get_ref.pl ---- get reference genotype from mass of phase files

=head1 Usage

  perl get_ref.pl <phase file> <genome file>

  
=head1 Example

=cut

   use warnings;
   use strict;
   use Getopt::Long;
   use Cwd; 
   
   my $h; 
   my %hash;
  GetOptions( 
           "help"=>\$h);

 die `pod2text $0` if (@ARGV==0 || $h);

   my @a;
  open LOD, "<$ARGV[0]";
  open LOP,"<$ARGV[1]";
 
 while(<LOD>){
        chomp;
    @a=split/\t/,$_;
   if($_=~/^M/){
    $hash{$a[1]}{$a[3]}++;
       }}

      my $chr;
      my @seq;
      my $pos;
  while(<LOP>){
          chomp;
        if($_=~/>(\S+)/){
              $chr=$1;
              $pos=0; 
                        }else{
        @seq=split//,$_;

    foreach my $index(@seq){
                      $pos++;
                  if(defined($hash{$chr}{$pos})){
                                  print "$chr\t$pos\t$index\n";
                                              }
                         }
           }                      
         }

  
  
 



