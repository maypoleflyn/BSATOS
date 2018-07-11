#!/usr/bin/perl

=head1 NAME
 
 perl cmp_infer_block.pl  ------ compare and infer  phase blocks from the merged file 

=head1 Usage

  perl cmp_infer_block.pl  

        
 

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

 
###first step infer missing 

        my @a;
        my $c=0;
        my @b;
        my %hash;
        my @z;
        my @x;         
        while(<>){
              chomp;
            push@a,$_;
            $c=$#a;
          if($c==1){
                my %pb=();
                my $p_s1="";
                my $p_s2="";
              
                my $m_s1="";
                my $m_s2="";
              
                my $h_s1="";
                my $h_s2="";
             
                my $l_s1="";
                my $l_s2="";
          
                @z=split/\t/,$a[0];
                @x=split/\t/,$a[1];

            foreach my $index(@a){
               @b=split/\t/,$index;
               $hash{$b[0]}{$b[1]}=[@b]; 
               $p_s1 .=$b[4];
               $p_s2 .=$b[5];
               $m_s1 .=$b[7];
               $m_s2 .=$b[8]; 
               $h_s1 .=$b[10];
               $h_s2 .=$b[11];
               $l_s1 .=$b[13];
               $l_s2 .=$b[14];                  
                               }
       
            if(($z[12] eq $x[12]) && ($z[12] ne "-") && ($z[15] eq $x[15]) && ($z[15] ne "-")){ 
                            if(($h_s1 ne $l_s1) && ($h_s1 ne $l_s2)){ 
                                                        shift@a;                                                        
                                                        next;
                                                                   }                  
                                                 }


                
             if(($z[12] eq $x[12]) && ($z[12] ne "-")){
                     if($z[4] eq $z[10] && (($x[4] eq "-") || ($x[6] ne $z[6]))){
                               $x[4]=$x[10];
                               $x[5]=$x[11];
                               $x[6]=$z[6];
                                                       }
                     if($z[4] eq $z[11] && (($x[4] eq "-") || ($x[6] ne $z[6]))){
                               $x[4]=$x[11];
                               $x[5]=$x[10];
                               $x[6]=$z[6];
                                                       }
                     if($z[7] eq $z[10] && (($x[7] eq "-") || ($x[9] ne $z[9]))){
                               $x[7]=$x[10];
                               $x[8]=$x[11];
                               $x[9]=$z[9];
                               
                                                       }                                                
                     if($z[7] eq $z[11] && (($x[7] eq "-") || ($x[9] ne $z[9]))) {
                               $x[7]=$x[11];
                               $x[8]=$x[10];
                               $x[9]=$z[9];
                                                       } 
                                                                
                                                   } 
                                       
                 if(($z[15] eq $x[15]) && ($z[15] ne "-")){
                     if($z[4] eq $z[13] && (($x[4] eq "-") ||($x[6] ne $z[6]))){
                               $x[4]=$x[13];
                               $x[5]=$x[14];
                               $x[6]=$z[6];
                                                       }
                     if($z[4] eq $z[14] && (($x[4] eq "-") || ($x[6] ne $z[6]))){
                               $x[4]=$x[14];
                               $x[5]=$x[13];
                               $x[6]=$z[6];
                                                       }
                     if($z[7] eq $z[13] && (($x[7] eq "-") || ($x[9] ne $z[9]))){
                               $x[7]=$x[13];
                               $x[8]=$x[14];
                               $x[9]=$z[9];

                                                       }
                     if($z[7] eq $z[14] && (($x[7] eq "-") || ($x[9] ne $z[9]))) {
                               $x[7]=$x[14];
                               $x[8]=$x[13];
                               $x[9]=$z[9];
                                                   }
                    
                                                 }                                                      
                      $hash{$x[0]}{$x[1]}=[@x];
                        
                         my $op=join("\t",@x);
                         $a[1]=$op;
                          shift@a; 
                                              
                          }
                          }
    my @chr=keys%hash;
                
           foreach my $index(@chr){
                     my @pos=sort {$a<=>$b} keys%{$hash{$index}};
                       foreach my $i(@pos){
                             print "@{$hash{$index}{$i}}\n";
                                }
                                          }  

                                               
