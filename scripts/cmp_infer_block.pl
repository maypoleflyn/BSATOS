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
 die `pod2text $0` if ( $ARGV[0] || $h);

 
###first step infer missing 

        my @a;
        my $c=0;
        my @b;
        my %hash;         
        while(<>){
              chomp;
            push@a,$_;
            $c=$#a;
          if($c==2){
                my %pb=();
                my $p_s1="";
                my $p_s2="";
                my %mb=();
                my $m_s1="";
                my $m_s2="";
                my %hb=();
                my $h_s1="";
                my $h_s2="";
                my %lb=();
                my $l_s1="";
                my $l_s2="";
                my %bl;
            foreach my $index(@a){
               @b=split/\t/,$index;
               $hash{$a[0]}{$a[1]}=[@b];
               $pb{$hash{$a[0]}{$a[1]}->[6]}++;
               $p_s1 .=$b[4];
               $p_s2 .=$b[5];
               $mb{$hash{$a[0]}{$a[1]}->[9]}++;
               $m_s1 .=$b[7];
               $m_s2 .=$b[8]; 
               $hb{$hash{$a[0]}{$a[1]}->[12]}++;
               $h_s1 .=$b[10];
               $h_s2 .=$b[11];
               $lb{$hash{$a[0]}{$a[1]}->[15]}++;
               $l_s1 .=$b[13];
               $l_s2 .=$b[14];                  
                                 }
             
          
                 my $l_p= $#{keys%pb};
                 my $l_m= $#{keys%mb}; 
                 my $l_h= $#{keys%hb};
                 my $l_l= $#{keys%lb}; 
                 my @seq=($p_s1,$p_s2,$m_s1,$m_s2,$h_s1,$h_s2,$l_s1,$l_s2);
                unless($p_s1=~/-/){
                if($l_p==0){  
                $bl{$p_s1}++;
                           }
                               }
             unless($p_s2=~/-/){
                if($l_p==0){
                $bl{$p_s2}++;
                           }
                               }
             unless($m_s1=~/-/){
                if($l_m ==0){
                $bl{$m_s1}++;
                            }
                               }
             unless($m_s2=~/-/){
                if($l_m ==0){
                $bl{$m_s2}++;
                            }
                               }
             unless($h_s1=~/-/){
                if($l_h ==0){
                $bl{$h_s1}++;
                            }
                               }
             unless($h_s2=~/-/){
               if($l_h ==0){   
                $bl{$h_s2}++;
                          }
                               }
             unless($l_s1=~/-/){
                 if($l_l==0){
                $bl{$l_s1}++;
                            }
                                }
              unless($l_s2=~/-/){
                if($l_l==0){
                $bl{$l_s2}++;
                           }
                               }
           
              if($l_p >0 && $l_m >0 && $l_h >0 && $l_l >0){
                            shift@a;  
                            next; 
                                      }else{
                               my @values=sort {$b<=>$a} values %bl;
                               my @key=keys %bl;
                               my @max;    
                               foreach my $index3(@key){          
                                      if($bl{$index3} == $values[0]){
                                          push@max,,$index3;                
                                            }
                                                       }
                              foreach my $index4(@seq){
                                                    
                                                     }   
