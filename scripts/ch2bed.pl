#!/usr/bin/perl

=head1 NAME
 
  ch2bed.pl ------ change to BED file

=head1 Usage

  perl ch2bed.pl    


=head1 Example

=cut

   use warnings;
   use strict;
   use Getopt::Long;
   use Cwd;
  my ($q,$g,$d,$h,$o);
  GetOptions(
        
           "help"=>\$h);

 die `pod2text $0` if (@ARGV==0 || $h);


my @a;  
open LOD,"<$ARGV[0]";

while(<LOD>){
     chomp;
     @a=split/\t/,$_;
 if($_=~/#/){
      next;
            }

  my $s=$a[1]-1;
  print "$a[0]\t$s\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]\t$a[8]\t$a[9]\t$a[10]\t$a[11]\t$a[12]\t$a[13]\t$a[14]\t$a[15]\n";
 
           } 
     


