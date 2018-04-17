# BSATOS
Bulked segregation analysis tools for outbreeding species
Welcome to the BSATOS wiki! bulked-segregant analyis tools for outbreeding species(BSATOS) is developed by Dr.Shen Fei,from institute for horticultue plants,China Argicultural University.BSATOS is designed for NGS-based Bulked-segregant analysis for outbreeding species,includings fruits trees,such as apple,cirtus.

To improve the gene mapping efficiency of next generation sequencing-based segregant analysis (BSA) in outbreeding species and realize the rapid candidate gene mining based on multi-omics data, bulked segregant analysis tools for outbreeding species (BSATOS) was developed. Different from the classic two-way pseudo-testcross (PT) strategy, Local phased 3-way BSA mapping strategy was introduced to use three type of segregate makers unbiased. After comparing different statistical methods to evaluate the allele frequency difference, G value method was selected(Paul M. Magwene 2011). Multi-round G value screening and multi-omics data based gene mining were also integrated.****





Program: bsatos (bulked segregant analysis tools for outbreeding species)
Version: 1.0.1

Usage:  bsatos <command> [options]

Commands:
  -- prep    prepare the input data
  -- phase   construct haplotype block 
  -- afd     calculate and filter allele frequency difference between two extreme pools
  -- ed      calculate ED based on the AFD  
  -- si      calculate |SNP index| based on the AFD
  -- smooth  smooth the afd/ed/si/g value 
  -- narrow_peak narrow down the peaks based on the multiple G value peaks and assign RDS to regions
  -- multi select genes based on multi-omics data 




