#!/usr/bin/perl

=head1 NAME
 
  get_block.pl ------ get the phase blocks from the phase file 

=head1 Usage

  perl get_block.pl <result from get_ref.pl> <phase file>   


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
    my $count=0;
    my $name;

open LOD,"<$ARGV[0]";
open LOP,"<$ARGV[1]";
 my %hash2;
 my @c;

while(<LOP>){
     chomp;
   @c=split/\t/,$_;
  $hash2{$c[0]}{$c[1]}=$c[2];
   }


   my $ref;
   my $alt1;
   my $alt2;
while(<LOD>){
      chomp;
     if($_=~/PS/){
            $count++;
        $name="block".$count;
        }
     if($_=~/^M/){
     #          print "$name\t$_\n";
    @a=split/\t/,$_;

         $ref=$hash2{$a[1]}{$a[3]}; 
         if($ref eq $a[4]){
              $alt1="0";
              $alt2="1";
  print "$a[1]\t$a[3]\t$ref\t$a[5]\t$alt1\t$alt2\t$name\n";
                         }
         if($ref eq $a[5]){
              $alt1="1";
              $alt2="0";
                          
                          
    print "$a[1]\t$a[3]\t$ref\t$a[4]\t$alt1\t$alt2\t$name\n";
              }
         }
     }








