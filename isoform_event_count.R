library(stringr)
library(ggplot2)

speList = c('Human','PteAle','EonSpe','Pacbio')
speDir = c('hg38','pteAle','eonSpe','pacbio')
gffList = c( 'hg38/Homo_sapiens.GRCh38.84.chr.gtf','pteAle/GCF_000325575.1_ASM32557v1_genomic.gtf','./eonSpe/maker.gff3.MAKER.rename.gtf','./eonSpe/gffcmp.annotated.gtf')

allCount = data.frame(Species=character(),SplicingEvents=character(), Count=integer())
i=1
j=1
for (i in 1:length(speList)){
    ioeFille = list.files(speDir[i], pattern="*.ioe.count")
    for (j in 1:length(ioeFille)){
        count =  read.table(paste(speDir[i],ioeFille[j],sep="/"), stringsAsFactors=F,header=F) [,1]
        variantType=gsub("_","", sub("gtf", "", str_match(ioeFille[j], "gtf_.*?_"),"gtf"))
        tmpCount = data.frame(Species=speList[i], SplicingEvents=variantType, Count=count)
        allCount = rbind(allCount,tmpCount)
    }
}

pdf("isoform_event_count_barplot.pdf", width=7,height=5,bg="white")
ggplot(data=allCount, aes(x=SplicingEvents, y=Count, fill=Species)) + geom_bar(stat="identity", position=position_dodge())
dev.off()
