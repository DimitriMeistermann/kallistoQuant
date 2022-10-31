args <- commandArgs(trailingOnly = TRUE)
options(stringsAsFactors = FALSE)

require(data.table)

#############
# FUNCTIONS #
#############
fastRead <- function(fileName, sep = '\t',row.names = 1,as.matrix=FALSE,stringsAsFactors=FALSE,...){
	require(data.table)
	dat <- as.data.frame(data.table::fread(fileName,stringsAsFactors=stringsAsFactors, sep = sep,...))
	if(!is.null(row.names)){
	  rownames(dat) <- dat[,row.names]
	  dat <- dat[,-row.names,drop=FALSE]
	}
	if(as.matrix) dat<-as.matrix(dat)
	return(dat)
}

fastWrite <- function(x, fileName = "default.tsv", headRow="Name",row.names=TRUE,col.names=TRUE, dec=".",sep="\t",...) {
	require(data.table)
	if(row.names){
		x=cbind(rownames(x),x)
		colnames(x)[1]<-headRow
	}
	fwrite(x=x,file=fileName,sep="\t", row.names = FALSE, col.names = col.names, quote = FALSE, dec=dec,...)
}


########
# MAIN #
########
# Construit la matrice de comptage des samples (en colonne) pour chaque gène (en ligne)

if (length(args)==0) stop("The fastq path argument is missing .n", call.=FALSE)

out_path <- args[1]
samples<- list.files(paste0(out_path,"/KALLISTO"))
countsList<-list()
tpmList<-list()
for(sample in samples){
	counts<-fastRead(paste0(out_path,"/KALLISTO/",sample,"/abundance.tsv"))
	countsList[[sample]]<-counts[,"est_counts",drop=FALSE]
	tpmList[[sample]]<-counts[,"tpm",drop=FALSE]
}
countsTable<-do.call(cbind,countsList)
colnames(countsTable)<-samples

tpmTable<-do.call(cbind,tpmList)
colnames(tpmTable)<-samples


fastWrite(round(countsTable),paste0(out_path,"/results/rawCountsTable.tsv"))
fastWrite(tpmTable,paste0(out_path,"/results/TranscriptPerMillion.tsv"))