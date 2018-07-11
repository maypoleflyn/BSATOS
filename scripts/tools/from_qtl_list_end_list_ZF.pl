#!/usr/bin/perl


   use warnings;
   use strict;
   open LOD,"<$ARGV[0]";
     my %hash;
     my @a;
  while(<LOD>){
         chomp;
      @a=split/\t/,$_;
      $hash{$a[0]}=[$a[1],$a[2],$a[3]];
        }
 
  my @key=keys%hash;
 
  foreach my $index(@key){
            my $dir = $index."_qtl";
            my $out1 =$dir."/".$index."_Gvalue";
            my $out2 =$dir."/".$index."_variation";
            my $out3=$dir."/".$index."_RE";
            my $out4=$dir."/".$index."_domain";
            my $out5=$dir."/".$index."_variation_func";
            my $out6=$dir."/".$index."_variation_func.anno";
            my $out7=$dir."/".$index."_Gvalue.anno";
            my $out8=$dir."/".$index."_varition_filtered";
            
            my $out9=$dir."/".$index."_sv_filtered";
            my $out10=$dir."/".$index."_sv_out"; 
       system("mkdir $dir;perl /home/shen/project/work/big_project/qtlseq/tools/get.bed.pl sum_all $hash{$index}->[0] $hash{$index}->[1] $hash{$index}->[2] $index > $out1;perl /home/shen/project/work/big_project/qtlseq/tools/get.bed.pl  /home/shen/project/work/big_project/project3/new/G_vs_J_genome_variations/SNV/domestic.vcf_out.apple_multianno.txt $hash{$index}->[0] $hash{$index}->[1] $hash{$index}->[2] $index > $out2");
       system("perl /home/shen/project/work/big_project/qtlseq/tools/get.bed.pl  /home/shen/project/work/big_project/project3/new/function_domain/domain $hash{$index}->[0] $hash{$index}->[1] $hash{$index}->[2] $index > $out4; perl /home/shen/project/work/big_project/qtlseq/tools/get.bed.pl  /home/shen/project/work/big_project/project3/new/function_domain/RE.bed $hash{$index}->[0] $hash{$index}->[1] $hash{$index}->[2] $index > $out3");  
      system("perl /home/shen/project/work/big_project/qtlseq/tools/add_re_and_domain.pl $out4 $out3 $out2 > $out5");
      system("perl /home/shen/project/work/big_project/qtlseq/tools/add.refer.pl $out5 >$out6" );
      system("perl /home/shen/project/work/big_project/qtlseq/tools/Gvalue_add_anno.pl $out6 $out1 >$out7");  
         
     system("perl /home/shen/project/work/big_project/qtlseq/tools/get.bed.pl /home/shen/project/work/big_project/project3/new/G_vs_J_genome_variations/SV/sv.filtered $hash{$index}->[0] $hash{$index}->[1] $hash{$index}->[2] $index > $out10");  
 

      if($index=~/H/){
                 system("perl /home/shen/project/work/big_project/qtlseq/tools/filter_ZF.pl $out6 > $out8"); 
                 system("perl /home/shen/project/work/big_project/qtlseq/tools/filter_ZF.pl $out10 > $out9");
                  }
        if($index=~/Z/){
                 system("perl /home/shen/project/work/big_project/qtlseq/tools/filter_Z.pl $out6 > $out8");
                
                 system("perl /home/shen/project/work/big_project/qtlseq/tools/filter_Z.pl $out10 > $out9");

     
                  }
if($index=~/F/){
                 system("perl /home/shen/project/work/big_project/qtlseq/tools/filter_F.pl $out6 > $out8"); 

              system("perl /home/shen/project/work/big_project/qtlseq/tools/filter_F.pl $out10 > $out9");
    
                  }
                            
        system("cat $out7 >> qtl_Gvalue_sum");
        system("cat $out8 >> qtl_variation_sum");
        system("cat $out9 >> qtl_sv_sum");   
 




        }    
                  
      
          
   
