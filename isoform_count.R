library(stringr)
library(ggplot2)

speList = c('Human','PteAle','EonSpe','Pacbio')
gfflist = c( 'hg38/Homo_sapiens.GRCh38.84.chr.gff3','pteAle/GCF_000325575.1_ASM32557v1_genomic.gff','./eonSpe/maker.gff3.MAKER.rename','./pacbio/gffcmp.annotated.gtf')

allCount = data.frame(Species=character(),IsoformCount=numeric(), Count=integer())
i=1
for (i in 1:length(speList)){
    gff = read.table(gfflist[i],header=F,sep="\t", stringsAsFactors=F,fill=T)
    gff = gff[gff[,3] == "mRNA" | gff[,3] == 'transcript' | gff[,3] == 'pseudogene'|gff[,3] == 'lincRNA'|gff[,3]=='processed_transcript' |gff[,3]=='processed_pseudogene',]
    parents = sapply(gff[,9], function(x) str_match(x, "Parent=.*?;|gene=.*?;|gene_id.*?;"))
    print(length(parents))
    count = table(table(parents))
    count = c(count[1:15], sum(count[16:end(count)[1]]))
    names(count)[16]= "16+"
    tmpCount = data.frame(Species=rep(speList[i],length(count)), IsoformCount=names(count), Count=count)
    allCount = rbind(allCount,tmpCount)
}

theTable <- within(allCount, IsoformCount <- factor(IsoformCount, levels=names(sort(table(IsoformCount), decreasing=TRUE))))

positions = names(count)
pdf("isoform_count_barplot.pdf", width=7,height=5,bg="white")
ggplot(data=theTable, aes(x=IsoformCount, y=Count, fill=Species)) + geom_bar(stat="identity", position=position_dodge()) + scale_x_discrete(limits = positions)
dev.off()
