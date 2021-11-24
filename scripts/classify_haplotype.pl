#!/usr/bin/perl


=head1 NAME
 
 classify_haplotype.pl   ------ compare and classify  phase blocks from the merged file 

=head1 Usage

  perl classify_haplotype.pl <file>

        
 

=head1 Example

=cut

   use warnings;
   use strict;
   use Getopt::Long;
   use Cwd;
    my ($h);
        GetOptions(
           "help"=>\$h);
 die `pod2text $0` if ( @ARGV==0 || $h);
   
     my @a;
     my %hash1; 
     my %hash2; 
     my %hash3;
     my %hash4;
    open LOD,"<$ARGV[0]";
    while(<LOD>){
      
       if($_=~/#/){
              next;
                 }

          chomp;
       @a=split/\t/,$_;
     unless($a[6] eq "-"){
       $hash1{$a[6]}=$a[0];
                         }
     unless($a[9] eq "-"){
       $hash3{$a[9]}=$a[0];
                         }
     unless($a[6] eq "-"){
       push@{$hash2{$a[6]}},$a[1];
                          }
     unless($a[9] eq "-"){
       push@{$hash4{$a[9]}},$a[1];
            }
           }


  my @k1=keys %hash1;
  my @k2=keys %hash3;
  
   my $outp=$ARGV[1];
   my $outm=$ARGV[2];
     open LOL,">$outp";
     open LOK,">$outm";

     foreach my $index (@k1){
                  my @d=sort {$a<=>$b} @{$hash2{$index}};
                if($#d >5){  
                  my $start =shift@d;
                  my $end =pop@d; 
              print LOL "$hash1{$index}\t$start\t$end\t$index\n";
                             }
                            }
      
      foreach my $index1 (@k2){
                  my @e=sort {$a<=>$b} @{$hash4{$index1}};
                if($#e >5){
                  my $start1 =shift@e;
                  my $end1 =pop@e;
              print LOK "$hash3{$index1}\t$start1\t$end1\t$index1\n";
                             }
                            }
   





