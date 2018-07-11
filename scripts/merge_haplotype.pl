#!/usr/bin/perl


=head1 NAME
  merge_haplotype.pl ------ the script to merge haplotype of the files

=head1 Usage

  perl merge_haplotype.pl --outprefix <prefix> --Pb <> --Mb <> --Hb <> --Lb <>

  --outprefix <out prefix>
  --Pb block file from pollen parent
  --Mb block file from maternal parent
  --Hb block file from high pool 
  --Lb block file from low  pool




=head1 Example

=cut

   use warnings;
   use strict;
   use Getopt::Long;
   use Cwd;
   my ($h,$o,$Pb,$Mb,$Hb,$Lb);
    GetOptions(
           "Pb=s" =>\$Pb,
           "Mb=s" =>\$Mb,
           "Hb=s" =>\$Hb,
           "Lb=s" =>\$Lb,
           "outprefix=s"=>\$o, 
           "help"=>\$h);

 die `pod2text $0` if ($h);
   
    my %hash;
    
    open LOD,"<$Pb";
    open LOP,"<$Mb";
    open LOT,"<$Hb";
    open LOG,"<$Lb";
    my @a;

  while(<LOD>){
           chomp;
       @a=split/\t/,$_;
       $hash{$a[0]}{$a[1]}++;
              }

  while(<LOP>){
           chomp;
       @a=split/\t/,$_;
       $hash{$a[0]}{$a[1]}++;
              }

   while(<LOT>){
           chomp;
       @a=split/\t/,$_;
       $hash{$a[0]}{$a[1]}++;
              }

    while(<LOG>){
           chomp;
       @a=split/\t/,$_;
       $hash{$a[0]}{$a[1]}++;
              }

  
  my @chr=sort %hash;
  
  foreach my $index(@chr){
             my @pos= sort {$a<=>$b} keys%{$hash{$index}};
            foreach my $index2 (@pos){
                    print "$index\t$index2\n";
                                     }
                         }



             
                        








     



