

load(file="Chr01.Rdata")
chr1<-data;
load(file="Chr02.Rdata")
chr2<-data;
load(file="Chr03.Rdata")
chr3<-data;
load(file="Chr04.Rdata")
chr4<-data;
load(file="Chr05.Rdata")
chr5<-data;
load(file="Chr06.Rdata")
chr6<-data;
load(file="Chr07.Rdata")
chr7<-data;
load(file="Chr08.Rdata")
chr8<-data;
load(file="Chr09.Rdata")
chr9<-data;
load(file="Chr10.Rdata")
chr10<-data;
load(file="Chr11.Rdata")
chr11<-data;
load(file="Chr12.Rdata")
chr12<-data;
load(file="Chr13.Rdata")
chr13<-data;
load(file="Chr14.Rdata")
chr14<-data;
load(file="Chr15.Rdata")
chr15<-data;
load(file="Chr16.Rdata")
chr16<-data;
load(file="Chr17.Rdata")
chr17<-data;
load(file="Chr00.Rdata")
chr0<-data;
data<-rbind(chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr0)

save(data,file="sum.Rdata")
write.table(data,file="sum_all",append = FALSE, quote =F, sep ="\t",row.names = F,col.names = F);
