if(!"tripsAndDipR" %in% installed.packages()){
  if(!"devtools" %in% installed.packages()){
    install.packages("devtools")
  }  
  library(devtools)
  install_github("delomast/tripsAndDipR",force=TRUE)
}
suppressMessages(require(tripsAndDipR))
#suppressMessages(require(parallel))


args <- commandArgs()
minP <- as.numeric(args[6])
maxP <- as.numeric(args[7])
LLRthreshold <- as.numeric(args[8])
MRCPL <- as.numeric(args[9])
MIN <- as.numeric(args[10])
OTthreshold <- as.numeric(args[11])
wd <- args[12]
library_name <- args[13]



## wd to be suppplied as argument
# wd<-getwd()
genosFiles <- dir(path = wd,pattern = "\\.genos$")
geno_fastq<-as.data.frame(matrix(nrow=length(genosFiles),ncol=0))
geno_fastq$geno<-genosFiles
wd2<-paste(wd,"/",sep="")
genosFiles <- paste(wd2,genosFiles,sep="")

#set geno classes, to be supplied as argument or take default
# minP<-4
# maxP<-6
ploidies<-seq(minP,maxP)


##automatically set #loci; read num.lines in first geno file
firstgeno<-read.delim(genosFiles[1],header=FALSE)
nlines<-nrow(firstgeno)
nucs<-lapply(lapply(as.list(as.character(firstgeno[2:nlines,])), function(x) strsplit(x,",")), function(x) unlist(strsplit(unlist(x)[c(2,3)],"="))[c(1,3)])
names(nucs)<-lapply(lapply(as.list(as.character(firstgeno[2:nlines,])), function(x) strsplit(x,",")), function(x) unlist(x)[1])

refCounts <- matrix(nrow = 0, ncol = nlines-1)
altCounts <- matrix(nrow = 0, ncol = nlines-1)

#f=1
fastq<-NULL
headers<-NULL
for(f in genosFiles){
	rReads <- c() # ref
	aReads <- c() # alt
	mNames <-c() # locus names
	gFile <- file(f, "r")
	line <- readLines(gFile, n = 1)
	headers<-c(headers,line)
	fastq<-c(fastq,strsplit(line, ",")[[1]][1])
	sName <- gsub("\\.fastq$", "", strsplit(line, ",")[[1]][1]) # get sample name
	line <- readLines(gFile, n = 1) # read first marker line
	while(length(line) > 0){
		sep <- strsplit(line, ",")[[1]]
		mNames <- c(mNames, sep[1])
		rReads <- c(rReads, as.numeric(gsub("^[ACTG-]=", "", sep[2])))
		aReads <- c(aReads, as.numeric(gsub("^[ACTG-]=", "", sep[3])))
		line <- readLines(gFile, n = 1)
	}
	close(gFile)

	# save data, use names to make sure all in same order
	names(rReads) <- mNames
	names(aReads) <- mNames
	refCounts <- rbind(refCounts, rReads)
	altCounts <- rbind(altCounts, aReads)
	rownames(refCounts)[nrow(refCounts)] <- sName
	rownames(altCounts)[nrow(altCounts)] <- sName

}
geno_fastq$fastq<-fastq
geno_fastq$headers<-headers
print(paste(c("R script counted reads from ",nrow(refCounts)," files"),collapse=""))

# fit models and calculate LLRs
print("R script estimating ploidy: PLEASE WAIT")
fp <- funkyPloid(refCounts, altCounts, ploidy = ploidies, maxIter = 150, maxDiff = .01, model = "BB_noise")
# est_ploid<-function(i){
#   funkyPloid(t(as.matrix(refCounts[i,])), t(as.matrix(altCounts[i,])), ploidy = ploidies, maxIter = 150, maxDiff = .01, model = "BB_noise")
# }
# time<-system.time(save<-mclapply(1:nrow(refCounts),est_ploid))
# print(time)
# fp<-data.frame(matrix(unlist(save), nrow=length(save), byrow=T),stringsAsFactors=FALSE)
# colnames(fp)<-colnames(save[[1]])
# fp$Ind<-rownames(refCounts)
# fp[,2:ncol(fp)]<-apply(fp[,2:ncol(fp)],2,function(x) as.numeric(x))
print("R script finished esimating ploidy; proceeding to genotyping")

##get LLR & OT read count thresholds as arguments; determine ploidy from LLR excluding those below LLR and OT read count threshold
# LLRthreshold=3
# OTthreshold=25000
LLRs <- split(fp[3:(2+length(ploidies))], seq(nrow(fp[3:(2+length(ploidies))])))
names(LLRs)<-fp$Ind
fp$minAltLLR<-unlist(lapply(LLRs,function(x) ifelse(as.logical(unique(x %in% "NA")),NA,x[which(x==min(x[!x==0]))])))
fp$ploidy<-unlist(lapply(LLRs,function(x) ifelse(as.logical(unique(x %in% "NA")),NA,ploidies[which(x==min(x))])))
LLRs<-NULL
fp$OTreads<-rowSums(refCounts)+rowSums(altCounts)
fp[fp$minAltLLR<LLRthreshold|fp$OTreads<OTthreshold,"ploidy"]<-NA

#create matrix of read ratios; exclude those below minimum read count per locus (MRCPL) threshold; get RCPL from arguments
# MRCPL<-20
ratios<-refCounts/altCounts
counts<-refCounts+altCounts
ratios[ratios=="Inf"]<-100
ratios[counts<MRCPL]<--9

##calculate read ratio limits based on 95% confidence intervals for each ploidy, with homozygotes adjusted to approx. maintain maximum 1in20 ratio
ploidy.ratio.limits<-list()
for(p in 1:length(ploidies)){
  div<-ploidies[p]
  sd<-0.05*(4/div)
  het.quants<-NULL
  for(h in 1:(div-1)){
    het.quants<-c(het.quants,qnorm(0.025, mean = h*(1/div), sd = sd),qnorm(0.975, mean = h*(1/div), sd = sd))
  }
  quants95<-c(0,qnorm(0.975, mean = 0, sd = sd),het.quants,qnorm(0.025, mean = 1, sd = sd))
  ratio.limits<-(quants95)/(1-quants95)
  ratio.limits[1]<-ratio.limits[1]/(1+(1/(div-3)))
  ratio.limits[length(ratio.limits)]<-ratio.limits[length(ratio.limits)]*(1+(1/(div-3)))
  ploidy.ratio.limits[[div]]<-ratio.limits
}
#ploidy.ratio.limits

##create genotype placeholders to correspond to read ratio limits
ploidy.genotype.classes<-list()
for(p in 1:length(ploidies)){
  n=0
  genotype.classes<-NULL
  for(g in 1:(1+ploidies[p])) {
    genotype.classes<-c(genotype.classes,c(paste(c(rep("M",times=(n)),rep("N",times=(ploidies[p]-n))),collapse=""),paste(rep("0",times=(ploidies[p])),collapse="")))
    n=n+1
  }
  genotype.classes<-genotype.classes[1:length(ploidy.ratio.limits[[ploidies[p]]])]
  ploidy.genotype.classes[[ploidies[p]]]<-genotype.classes
}
#ploidy.genotype.classes

genotypes<-as.data.frame(matrix(nrow=0,ncol=ncol(refCounts)))
for(i in 1:nrow(ratios)){
  i.ploid<-fp[fp$Ind==rownames(ratios)[i],"ploidy"]
  if(is.na(i.ploid)){
    genotypes<-rbind(genotypes,rep(paste(rep("0",times=(min(ploidies))),collapse=""),times=ncol(refCounts)),stringsAsFactors = FALSE)
  } else if(!is.na(i.ploid)){
    geno<-NULL
    for(l in 1:length(ratios[i,])) { 
      if (ratios[i,l]<0){
        geno<-c(geno,paste(rep("0",times=(i.ploid)),collapse=""))
      } else if (as.numeric(ratios[i,l])>=0) {
        keep<-ploidy.genotype.classes[[i.ploid]][(as.numeric(ratios[i,l])+1E-20)>ploidy.ratio.limits[[i.ploid]]]
        class(keep)
        geno<-c(geno,keep[length(keep)])
      }
    }
    genotypes<-rbind(genotypes,geno,stringsAsFactors = FALSE)
  }
}
colnames(genotypes)<-colnames(refCounts)
rownames(genotypes)<-rownames(refCounts)

for(l in 1:ncol(genotypes)){
  genotypes.list<-lapply(as.list(genotypes[,l]),function(x) unlist(strsplit(x,"")))
  gsub("M",nucs[[l]][1],genotypes.list[[1]])
  genotypes.list<-lapply(genotypes.list,function(x) gsub("M",nucs[[l]][1],x))
  genotypes.list<-lapply(genotypes.list,function(x) gsub("N",nucs[[l]][2],x))
  genotypes[,l]<-unlist(lapply(genotypes.list,function(x) paste(x,collapse="")))
}

rownames(fp)<-fp$Ind
genotypes<-merge(fp,genotypes,by="row.names")
genotypes<-genotypes[,!colnames(genotypes) %in% c("Row.names")]
loci<-grep("Atr",colnames(genotypes),value=TRUE)
missing.genos<-NULL
for(N in 1:maxP){
  missing.genos<-c(missing.genos,paste(rep("0",times=N),collapse=""))
}
genotypes$X.GT<-((length(loci)-apply(genotypes[,loci], 1, function(x) sum(x %in% missing.genos)))/length(loci))*100;
rownames(genotypes)<-rownames(fp)

geno_fastq$fastq_name<-unlist(lapply(as.list(geno_fastq$headers),function(x) gsub("\\.fastq$", "",unlist(strsplit(x,","))[1]) ))
geno_fastq$raw<-unlist(lapply(as.list(geno_fastq$headers),function(x) unlist(strsplit(unlist(strsplit(x,","))[2],":"))[2] ))
geno_fastq$OTheader<-unlist(lapply(as.list(geno_fastq$headers),function(x) unlist(strsplit(unlist(strsplit(x,","))[3],":"))[2] ))
geno_fastq$OTpercent<-unlist(lapply(as.list(geno_fastq$headers),function(x) unlist(strsplit(unlist(strsplit(x,","))[4],":"))[2] ))
geno_fastq$IFI<-unlist(lapply(as.list(geno_fastq$headers),function(x) unlist(strsplit(unlist(strsplit(x,","))[5],":"))[2] ))
geno_fastq$genos_name<-unlist(lapply(as.list(geno_fastq$headers),function(x) paste(na.omit(unlist(strsplit(gsub("\\.fastq$", "",unlist(strsplit(x,","))[1]),"_"))[4:100]),collapse="_") ))
rownames(geno_fastq)<-rownames(genotypes)
csv_dump<-merge(geno_fastq,genotypes,by="row.names")

#header for Lxxxx_GenosNF.csv: Sample	Raw Reads	On-Target Reads	%On-Target	%GT --> Loci
#header for Lxxxx_ProgGenos90.csv: Sample	Raw Reads	On-Target Reads	%On-Target	%GT	IFI --> Loci
#add: ploidy, minaltLLR

csv_GenosNF<-csv_dump[,c("fastq_name","raw","OTheader","OTpercent","X.GT","ploidy","minAltLLR",loci)]
header_csv_GenosNF<-c("Sample","Raw Reads","On-Target Reads","%On-Target","%GT","Ploidy Result","Ploidy Confidence",names(nucs))
csv_GenosNF<-rbind(header_csv_GenosNF,csv_GenosNF)
#library_name<-"Lxxxx"
csv_GenosNF_name<-paste(c(library_name,"_ploidy_GenosNF.csv"),collapse="")
write.table(csv_GenosNF,file=csv_GenosNF_name,col.names=FALSE,row.names=FALSE,quote=FALSE,sep=",")

#MIN<-90
csv_GenosMIN<-csv_dump[,c("genos_name","raw","OTheader","OTpercent","X.GT","IFI","ploidy","minAltLLR",loci)]
for(r in 1:nrow(csv_GenosMIN)){
  if(csv_GenosMIN[r,"X.GT"]<MIN&(!is.na(csv_GenosMIN[r,"ploidy"]))){csv_GenosMIN[r,loci]<-missing.genos[as.numeric(csv_GenosMIN[r,"ploidy"])] }
}
header_csv_GenosMIN<-c("Sample","Raw Reads","On-Target Reads","%On-Target","%GT","IFI","Ploidy Result","Ploidy Confidence",names(nucs))
csv_GenosMIN<-rbind(header_csv_GenosMIN,csv_GenosMIN)
csv_GenosMIN_name<-paste(c(library_name,"_ploidy_ProgGenos",MIN,".csv"),collapse="")
write.table(csv_GenosMIN,file=csv_GenosMIN_name,col.names=FALSE,row.names=FALSE,quote=FALSE,sep=",")

if(nrow(genotypes)>0&nrow(csv_dump)>0){
  print("R script completed successfully")
}
