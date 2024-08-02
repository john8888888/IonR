#function to convert bed file from WG or ampliseq.com to a bed file with gene symbol
add.gene.id.to.bed<-function(bedfile, out=paste0(bedfile, ".withGene.bed")){
  #bed file is the .Designed.bed file from ampliseq.com
bedfile='~/Downloads/IAD193572_361_results/IAD193572_361_Designed.bed'
header=t(as.vector(read.delim(bedfile, nrows=1, sep=' ', header=F, stringsAsFactors = F)))
genome=strsplit(header[grep(header,pattern='hg')], split="=")[[1]][2]
if(!genome %in% c('hg19', 'hg38')) stop('The design was not human genome!!!')
targets<-read.delim(file=bedfile, header=F, skip=1, sep='\t', stringsAsFactors = F)
targets$chromsome_name<-substr(targets[,1], rep(4, nrow(targets)), nchar(targets[,1]))
temp<-apply(targets[1:2,c('chromsome_name', 'V2', 'V3')], 1, function(x)as.list(x, names=NULL))
library("biomaRt")
if(genome=='hg19'){
  hostURL="http://grch37.ensembl.org"
  ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host=hostURL)
  attributes<-c('hgnc_symbol','chromosome_name', 'start_position', 'end_position')
  filters = c("chromosome_name", "start", "end")
  temp1<-getBM(attributes=attributes, filters =filters, values =targets[1:2,c('chromsome_name', 'V2', 'V3')] , mart = ensembl)
  
} else {
  hostURL="http://grch37.ensembl.org"
}
listMarts()
#listMarts(host="dec2013.archive.ensembl.org") #make sure we use the same ensembl version
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl") #for v74 hg19
#ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl) #latest v95
temp<-listFilters(ensembl)
filters<-'hgnc_symbol'
attributes<-c('hgnc_symbol','chromosome_name', 'start_position', 'end_position', 'strand',  'entrezgene')
attributes_l<-c('hgnc_symbol','cds_length')
annot<-getBM(attributes=attributes, filters =filters, values = rownames(Cdat1), mart = ensembl)
annot<-annot[order(annot[,2],annot[,3]),] #sort
annot1<-annot[!duplicated(annot[,1]),]
annot.f<-merge(x=data.frame(hgnc_symbol=rownames(Cdat1)), y=annot1, all.x=T)
rownames(annot.f)<-as.vector(annot.f[,1])
annot.f<-annot.f[match(rownames(Cdat1), rownames(annot.f)),]
}
