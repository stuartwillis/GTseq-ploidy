
setwd("~/Desktop/");


#######----------------------------------------------------------------------------------------------------------------
##ploidy estimation

if(!"tripsAndDipR" %in% installed.packages()){
  if(!"devtools" %in% installed.packages()){
    install.packages("devtools")
  }  
  library(devtools)
  install_github("delomast/tripsAndDipR",force=TRUE)
}
library(tripsAndDipR)
library(parallel)

# pull out read counts from sample sturgeon data, wd to be suppplied as argument
wd<-getwd()
genosFiles <- dir(path = wd,pattern = "\\.genos$")
tail(genosFiles)

length(genosFiles)

geno_fastq<-as.data.frame(matrix(nrow=length(genosFiles),ncol=0))
geno_fastq$geno<-genosFiles
wd2<-paste(wd,"/",sep="")
genosFiles <- paste(wd2,genosFiles,sep="")

#set geno classes, to be supplied as argument or take default
minP<-4
maxP<-6
ploidies<-seq(minP,maxP)


##automatically set #loci; read num.lines in first geno file
firstgeno<-read.delim(genosFiles[1],header=FALSE)
nlines<-nrow(firstgeno)
head(firstgeno)
class(firstgeno)
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
	# get sample name
	line <- readLines(gFile, n = 1)
	headers<-c(headers,line)
	fastq<-c(fastq,strsplit(line, ",")[[1]][1])
	sName <- gsub("\\.fastq$", "", strsplit(line, ",")[[1]][1])
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
head(geno_fastq)
head(refCounts)
nrow(refCounts)
sum(colnames(refCounts) %in% names(nucs))

print(paste(c("R script counted reads from ",nrow(refCounts)," files"),collapse=""))

# fit models and calculate LLRs
print("R script estimating ploidy: PLEASE WAIT")
#fp <- funkyPloid(refCounts, altCounts, ploidy = ploidies, maxIter = 150, maxDiff = .01, model = "BB_noise")#non-parallel
est_ploid<-function(i){
  funkyPloid(t(as.matrix(refCounts[i,])), t(as.matrix(altCounts[i,])), ploidy = ploidies, maxIter = 150, maxDiff = .01, model = "BB_noise")
}
time<-system.time(save<-mclapply(1:nrow(refCounts),est_ploid))
#print(time)
fp<-data.frame(matrix(unlist(save), nrow=length(save), byrow=T),stringsAsFactors=FALSE)
colnames(fp)<-colnames(save[[1]])
fp$Ind<-rownames(refCounts)
fp[,2:ncol(fp)]<-apply(fp[,2:ncol(fp)],2,function(x) as.numeric(x))
print(paste("R script finished esimating ploidy in ",floor(as.numeric(time[3]))," seconds; proceeding to genotyping",sep=""))
#head(fp)
nrow(fp)

##determine ploidy from LLR (excluding those below LLR and OT read count threshold, if applicable)
# LLRthreshold=3
# OTthreshold=25000
LLRs <- split(fp[3:(2+length(ploidies))], seq(nrow(fp[3:(2+length(ploidies))])))
names(LLRs)<-fp$Ind
#as.logical(unique(LLRs[[2453]] %in% "NA"))
fp$minAltLLR<-unlist(lapply(LLRs,function(x) ifelse(as.logical(unique(x %in% "NA")),NA,x[which(x==min(x[!x==0]))])))
fp$ploidy<-unlist(lapply(LLRs,function(x) ifelse(as.logical(unique(x %in% "NA")),NA,ploidies[which(x==min(x))])))
LLRs<-NULL
fp$OTreads<-rowSums(refCounts)+rowSums(altCounts)
#fp[fp$minAltLLR<LLRthreshold|fp$OTreads<OTthreshold,"ploidy"]<-NA
head(fp)
#write.table(fp,file="ploidy-test-results.txt",col.names=T,row.names = T)




#######----------------------------------------------------------------------------------------------------------------
##genotyping

#create matrix of read ratios; exclude those below minimum read count per locus (MRCPL) threshold
MRCPL<-20
ratios<-refCounts/altCounts
counts<-refCounts+altCounts
ratios[ratios=="Inf"]<-100
ratios[counts<MRCPL]<--9

##calculate read ratio limits based on 95% confidence intervals for each ploidy, with homozygotes adjusted to approx. maintain maximum 1in20 ratio
#p=1
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
ploidy.ratio.limits

##create genotype placeholders to correspond to read ratio limits; change for Polygene/Colony2/Adegent or ProgenyCSV output
p=1
ploidy.genotype.classes<-list()
for(p in 1:length(ploidies)){
  n=0
  genotype.classes<-NULL
  for(g in 1:(1+ploidies[p])) {
    genotype.classes<-c(genotype.classes,c(paste(c(rep("1",times=(n)),rep("2",times=(ploidies[p]-n))),collapse=","),""))
    n=n+1
  }
  genotype.classes<-genotype.classes[1:length(ploidy.ratio.limits[[ploidies[p]]])]
  ploidy.genotype.classes[[ploidies[p]]]<-genotype.classes
}
ploidy.genotype.classes

genotypes<-as.data.frame(matrix(nrow=0,ncol=ncol(refCounts)))
for(i in 1:nrow(ratios)){
  i.ploid<-fp[fp$Ind==rownames(ratios)[i],"ploidy"]
  if(is.na(i.ploid)){
    genotypes<-rbind(genotypes,rep(paste(rep("0",times=(min(ploidies))),collapse=""),times=ncol(refCounts)),stringsAsFactors = FALSE)
  } else if(!is.na(i.ploid)){
    geno<-NULL
    for(l in 1:length(ratios[i,])) { 
      if (ratios[i,l]<0){
        geno<-c(geno,"")
      } else if (as.numeric(ratios[i,l])>=0) {
        keep<-ploidy.genotype.classes[[i.ploid]][(as.numeric(ratios[i,l])+1E-20)>ploidy.ratio.limits[[i.ploid]]]
        geno<-c(geno,keep[length(keep)])
      }
    }
    genotypes<-rbind(genotypes,geno,stringsAsFactors = FALSE)
  }
}
colnames(genotypes)<-colnames(refCounts)
rownames(genotypes)<-rownames(refCounts)
head(genotypes[,1:10])

genos_files<-genosFiles
length(genos_files)
coverage<-NULL
OTcoverage<-NULL
fastqs<-NULL
heterozygosity<-NULL
for (f in 1:length(genos_files)){
  genos<-read.delim(genos_files[f],header=FALSE,sep=" ")
  #head(genos)
  fastq<-unlist(strsplit(as.character(genos[1,1]),","))[1]
  reads<-unlist(strsplit(unlist(strsplit(as.character(genos[1,1]),","))[2],":"))[2]
  counts<-lapply(genos[2:nrow(genos),1],function(x) unlist(strsplit(unlist(strsplit(as.character(x),","))[2:3],"="))[c(2,4)])
  ontarget<-unlist(strsplit(unlist(strsplit(as.character(genos[1,2]),","))[1],":"))[2]
  ##what is a heterozygote, i.e. not a homozygote, if the lower-count allele is still >0
  ##depends on ploidy and expected read ratios, i.e. hexaploid has ABBBBB ratio of 1/6 (0.16667)
  ##so..use half of this as threshold (1/6)/2, or (1/maxP)/2 (completely aside from gentoyping; only for determining heterozygosity)
  lapply(counts,function(x) (min(as.numeric(x))/max(as.numeric(x))))
  hetero<-sum(unlist(lapply(counts,function(x) ifelse(min(as.numeric(x))/max(as.numeric(x)+0.0000001)<((1/maxP)/2),0,1))))/length(counts)
  #nrow(stats[stats$V2 %in% depth$Pos,])
  OTcoverage<-c(OTcoverage,ontarget)
  coverage<-c(coverage,reads)
  fastqs<-c(fastqs,fastq)
  heterozygosity<-c(heterozygosity,hetero)
}
covs<-as.data.frame(cbind(fastqs,coverage,OTcoverage,heterozygosity))
nrow(covs)
rownames(covs)<-rownames(refCounts)
rownames(fp)<-rownames(refCounts)
covs<-merge(covs,fp,by="row.names")
covs<-covs[!colnames(covs) %in% "Row.names"]
covs$heterozygosity<-as.numeric(as.character(covs$heterozygosity))
covs$OTcoverage<-as.numeric(as.character(covs$OTcoverage))
rownames(covs)<-covs$Ind
genotypes<-merge(covs,genotypes,by="row.names")
genotypes<-genotypes[,!colnames(genotypes) %in% c("Row.names")]
loci<-grep("Atr",colnames(genotypes),value=TRUE)
genotypes$X.GT<-((length(loci)-apply(genotypes[,loci], 1, function(x) sum(x=="")))/length(loci))*100;

geno_fastq$fastq_name<-unlist(lapply(as.list(geno_fastq$headers),function(x) gsub("\\.fastq$", "",unlist(strsplit(x,","))[1]) ))
geno_fastq$raw<-unlist(lapply(as.list(geno_fastq$headers),function(x) unlist(strsplit(unlist(strsplit(x,","))[2],":"))[2] ))
geno_fastq$OTheader<-unlist(lapply(as.list(geno_fastq$headers),function(x) unlist(strsplit(unlist(strsplit(x,","))[3],":"))[2] ))
geno_fastq$OTpercent<-unlist(lapply(as.list(geno_fastq$headers),function(x) unlist(strsplit(unlist(strsplit(x,","))[4],":"))[2] ))
geno_fastq$IFI<-unlist(lapply(as.list(geno_fastq$headers),function(x) unlist(strsplit(unlist(strsplit(x,","))[5],":"))[2] ))
geno_fastq$genos_name<-unlist(lapply(as.list(geno_fastq$headers),function(x) paste(na.omit(unlist(strsplit(gsub("\\.fastq$", "",unlist(strsplit(x,","))[1]),"_"))[4:100]),collapse="_") ))
rownames(genotypes)<-genotypes$Ind
rownames(geno_fastq)<-rownames(genotypes)
csv_dump<-merge(geno_fastq,genotypes,by="row.names")
csv_dump<-csv_dump[,!names(csv_dump)=="Row.names"]


####create CSV file for import into Progeny (or parentage) 
csv_dumpProgeny<-csv_dump

##convert missing data to ploidy-accurate '0's
missing.genos<-NULL
for(N in 1:maxP){
  missing.genos<-c(missing.genos,paste(rep("0",times=N),collapse=""))
}

i=1
for (i in 1:nrow(csv_dumpProgeny)){
  csv_dumpProgeny[i,loci]<-gsub(",","",gsub("^$",missing.genos[as.numeric(as.character(csv_dumpProgeny[i,"ploidy"]))],csv_dumpProgeny[i,loci]))
}
#head(csv_dumpProgeny[,8:25])

#convert REF/ALT placeholders (1,2) to nucleotide genotypes (A,C,T,G)
l=1
for(l in 1:length(loci)){
  genotypes.list<-lapply(as.list(csv_dumpProgeny[,loci[l]]),function(x) unlist(strsplit(x,"")))
  #gsub("M",nucs[[l]][1],genotypes.list[[1]])
  genotypes.list<-lapply(genotypes.list,function(x) gsub("1",nucs[[l]][1],x))
  genotypes.list<-lapply(genotypes.list,function(x) gsub("2",nucs[[l]][2],x))
  csv_dumpProgeny[,loci[l]]<-unlist(lapply(genotypes.list,function(x) paste(x,collapse="")))
}
head(csv_dumpProgeny[,8:25])

csv_GenosNF<-csv_dumpProgeny[,c("genos_name","raw","OTheader","OTpercent","X.GT","ploidy","minAltLLR",loci)]
header_csv_GenosNF<-c("Sample","Raw Reads","On-Target Reads","%On-Target","%GT","Ploidy Result","Ploidy Confidence",loci)
csv_GenosNF<-rbind(header_csv_GenosNF,csv_GenosNF)
head(csv_GenosNF[,1:10])
#write.table(csv_GenosNF,file= "ploidy-genotypes-test-results.csv" ,col.names=FALSE,row.names=FALSE,quote=FALSE,sep=",")
csv_dumpProgeny<-NULL

##--->Optionally: skip to parentage at the end of this script






##output for Polygene

##Polygene format:
#(tab-delimited)
#ID Pop Ploidy Loc1 Loc2
#Ind1  pop1  4  2,3,4  1,2,4
#Ind2  pop1  4  4      2,3,4
#Ind3  pop2  4  1,4    1,2,3,4

loci<-grep("Atr",colnames(csv_dump),value=TRUE)
csv_keep<-csv_dump[!(is.na(csv_dump$ploidy))&csv_dump$X.GT>0,]
csv_keep$pop<-paste("pop",csv_keep$ploidy,sep="")
csv_polygene<-csv_keep[,c("genos_name","pop","ploidy",loci)]
csv_polygene<-apply(csv_polygene,2,function(x) as.character(x))
header_csv_polygene<-c("ID","Pop","Ploidy",loci)
csv_polygene<-rbind(header_csv_polygene,csv_polygene)
#write.table(csv_polygene,file="polygene.txt",col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")






####make files for Colony, pseudo-dominant coding of bi-allelic SNPs
csv_keep<-csv_dump[!(is.na(csv_dump$ploidy))&csv_dump$X.GT>0,]
loci<-grep("Atr",colnames(csv_keep),value=TRUE)
csv_colony<-csv_keep[,c("genos_name",loci)]
ncol(csv_colony)
csv_colony<-as.data.frame(apply(csv_colony,2,function(x) as.character(x)),stringsAsFactors = FALSE)
csv_colony[csv_colony==""]<-0
head(csv_colony)

colony_genos<-as.data.frame(matrix(nrow=nrow(csv_colony),ncol=0))
colony_genos$sample<-csv_colony$genos_name
#l=1
for(l in 1:length(loci)){
  #  ifelse(sum(1%in%unlist(strsplit(as.character(csv_GenosNF[1,loci[l]]),",")))==1,1,2)
  #  ifelse(sum(2%in%unlist(strsplit(as.character(csv_GenosNF[1,loci[l]]),",")))==1,1,2)
  #  lapply(csv_GenosNF[,loci[l]],function(x) unlist(strsplit(x,",")) )
  a1<-unlist(lapply(csv_colony[,loci[l]],function(x) ifelse(sum(1%in%unlist(strsplit(as.character( x ),",")))==1,1,2) ))
  a2<-unlist(lapply(csv_colony[,loci[l]],function(x) ifelse(sum(2%in%unlist(strsplit(as.character( x ),",")))==1,1,2) ))
  a1[csv_colony[,loci[l]]=="0"]<-"0"
  a2[csv_colony[,loci[l]]=="0"]<-"0"
  allele1<-paste(loci[l],"a",sep="_")
  allele2<-paste(loci[l],"b",sep="_")
  colony_genos[,allele1]<-a1
  colony_genos[,allele2]<-a2
}
head(colony_genos)

colnames(colony_genos)<-gsub("-","_",colnames(colony_genos))
dom_loci<-colnames(colony_genos)[!colnames(colony_genos)%in%c("sample","pop")]
length(dom_loci)


##NOTE: you will need to separate genotypes of adults/offspring before importing to colony
##optionally, import meta-information as "colony_genos$pop" indicating male/female/adult or offspring
#meta<-read.delim("colony_samples_meta_lifestage_or_sex.csv",header=TRUE,sep=",")
#colony_genos$pop<-meta[,1]#NOTE: assumes samples are in the same order!!


###sample names must be max 20 characters; truncate from the left (saving text on the right)
if(!"stringr" %in% installed.packages()){
  install.packages("stringr")
}  
require(stringr)
colony_genos$sample<-unlist(lapply(colony_genos$sample,function(x) str_trunc(x, 20, "left",ellipsis = "")))

##uncomment for lifestage or sex-specific exports
# #output offspring
# offspring_genos<-colony_genos[colony_genos$pop=="offspring",c("sample",dom_loci)]
# nrow(offspring_genos)
# write.table(offspring_genos,file="offspring_colony2.txt",col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
# #output dams
# mom_genos<-colony_genos[colony_genos$pop=="dam",c("sample",dom_loci)]
# head(mom_genos)
# write.table(mom_genos,file="dams_colony2.txt",col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
# #output sires
# pop_genos<-colony_genos[colony_genos$pop=="sire",c("sample",dom_loci)]
# head(pop_genos)
# write.table(pop_genos,file="sires_colony2.txt",col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
# #output all adults together
# mompop_genos<-colony_genos[colony_genos$pop%in%c("sire","dam"),c("sample",dom_loci)]
# nrow(mompop_genos)
# write.table(mompop_genos,file="adults_colony2.txt",col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")

#output all individuals together
export_genos<-colony_genos[,c("sample",dom_loci)]
nrow(export_genos)
#write.table(export_genos,file="samples_colony2_rate12loci.txt",col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")

##output marker meta information (required by Colony)
marker_meta<-as.data.frame(matrix(nrow=0,ncol=length(dom_loci)))
marker_meta[1,]<-dom_loci#marker names
marker_meta[2,]<-rep(1,times=length(dom_loci))#indicate dominant coding
marker_meta[3,]<-rep("0.0000",times=length(dom_loci))#estimated rate of null alleles
marker_meta[4,]<-rep("0.05",times=length(dom_loci))#estimated rate of genotyping error
#write.table(marker_meta,file="marker_meta_colony2_error05.txt",col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")











#######----------------------------------------------------------------------------------------------------------------
##estimate genotyping error using "named" duplicates genotyped at >90% (X.GT)

samples<-unlist(lapply(as.list(csv_dump$headers),function(x) paste(na.omit(unlist(strsplit(gsub("\\.fastq$", "",unlist(strsplit(x,","))[1]),"_"))[4:100]),collapse="_") ))
plates<-unlist(lapply(as.list(csv_dump$headers),function(x) paste(na.omit(unlist(strsplit(gsub("\\.fastq$", "",unlist(strsplit(x,","))[1]),"_"))[3]),collapse="_") ))
known_ploidy_index<-grep("P",plates,invert=TRUE)
samples[known_ploidy_index]<-unlist(lapply(as.list(geno_fastq$headers[known_ploidy_index]),function(x) paste(na.omit(unlist(strsplit(gsub("\\.fastq$", "",unlist(strsplit(x,","))[1]),"_"))[5:100]),collapse="_") ))
plates[known_ploidy_index]<-unlist(lapply(as.list(geno_fastq$headers[known_ploidy_index]),function(x) paste(na.omit(unlist(strsplit(gsub("\\.fastq$", "",unlist(strsplit(x,","))[1]),"_"))[4]),collapse="_") ))
sample_table<-as.data.frame(table(samples))
nrow(sample_table[sample_table$Freq>1,])
nrow(sample_table[sample_table$Freq==2,])

csv_dump$samples<-samples
csv_dump$plates<-plates
dgenos<-sample_table[sample_table$Freq>1,"samples"]
length(dgenos)
g=1
dup_genos<-as.data.frame(matrix(nrow=0,ncol=12))
for (g in 1:length(dgenos)){
  df<-csv_dump[csv_dump$samples %in% dgenos[g],]
  df2<-head(df[order(df$X.GT,decreasing=TRUE),],n=2)
  if(min(df2$X.GT)>=90){
    samp1=df2$fastq[1];samp2=df2$fastq[2]
    errors<-sum( as.integer(df2[df2$fastq==samp1,loci]==df2[df2$fastq==samp2,loci]) [(!df2[df2$fastq==samp1,loci]=="")&(!df2[df2$fastq==samp2,loci]=="")] )/length( as.integer(df2[df2$fastq==samp1,loci]==df2[df2$fastq==samp2,loci]) [(!df2[df2$fastq==samp1,loci]=="")&(!df2[df2$fastq==samp2,loci]=="")] )
    dup_genos<-rbind(dup_genos, cbind(as.character(dgenos[g]),samp1,df2[df2$fastq==samp1,"plates"],df2[df2$fastq==samp1,"X.GT"],df2[df2$fastq==samp1,"OTreads"],samp2,df2[df2$fastq==samp2,"plates"],df2[df2$fastq==samp2,"X.GT"],df2[df2$fastq==samp2,"OTreads"],errors,min(df2$X.GT),min(df2$OTreads)) )
  }
}
colnames(dup_genos)<-c("sample","fastq1","plate1","X.GT1","OTreads1","fastq2","plate2","X.GT2","OTreads2","X.GTmatch","minX.GT","minOTreads")
dup_genos$X.GTmatch<-as.numeric(as.character(dup_genos$X.GTmatch))
dup_genos$minX.GT<-as.numeric(as.character(dup_genos$minX.GT))
dup_genos$minOTreads<-as.numeric(as.character(dup_genos$minOTreads))

mean(dup_genos$X.GTmatch)





#######----------------------------------------------------------------------------------------------------------------
##filter samples to keep best duplicate among "named" duplicates
genos<-rownames(refCounts)
samples<-unlist(lapply(as.list(geno_fastq$headers),function(x) paste(na.omit(unlist(strsplit(gsub("\\.fastq$", "",unlist(strsplit(x,","))[1]),"_"))[4:100]),collapse="_") ))
plates<-unlist(lapply(as.list(geno_fastq$headers),function(x) paste(na.omit(unlist(strsplit(gsub("\\.fastq$", "",unlist(strsplit(x,","))[1]),"_"))[3]),collapse="_") ))
known_ploidy_index<-grep("P",plates,invert=TRUE)
samples[known_ploidy_index]<-unlist(lapply(as.list(geno_fastq$headers[known_ploidy_index]),function(x) paste(na.omit(unlist(strsplit(gsub("\\.fastq$", "",unlist(strsplit(x,","))[1]),"_"))[5:100]),collapse="_") ))
plates[known_ploidy_index]<-unlist(lapply(as.list(geno_fastq$headers[known_ploidy_index]),function(x) paste(na.omit(unlist(strsplit(gsub("\\.fastq$", "",unlist(strsplit(x,","))[1]),"_"))[4]),collapse="_") ))
reads<-rowSums(refCounts)+rowSums(altCounts)
df<-as.data.frame(cbind(genos,samples,plates,reads))
#rownames(df)
ugenos<-unique(samples)
nrow(refCounts)-length(ugenos)
new_genos<-as.data.frame(matrix(nrow=0,ncol=ncol(df)))
for (g in 1:length(ugenos)){
  df2<-df[df$samples %in% ugenos[g],]
  new_genos<-rbind(new_genos,head(df2[order(df2$reads,decreasing=TRUE),],n=1))
}
nrow(new_genos)

csv_dump_new_genos<-csv_dump[csv_dump$Row.names %in% as.character(new_genos$genos),]
nrow(csv_dump_new_genos)

##output numbers of sequenced reads for all and 90+% genotyped
min(csv_dump_new_genos$OTreads)
mean(csv_dump_new_genos$OTreads)
max(csv_dump_new_genos$OTreads)
min(as.numeric(csv_dump_new_genos$OTpercent))
mean(as.numeric(csv_dump_new_genos$OTpercent))
max(as.numeric(csv_dump_new_genos$OTpercent))
mean(csv_dump_new_genos[csv_dump_new_genos$X.GT>=90,"OTreads"])
median(csv_dump_new_genos[csv_dump_new_genos$X.GT>=90,"OTreads"])
min(csv_dump_new_genos[csv_dump_new_genos$X.GT>=90,"OTreads"])
quantile(csv_dump_new_genos[csv_dump_new_genos$X.GT>=90,"OTreads"],probs=c(0.01,0.05,0.1))


#######----------------------------------------------------------------------------------------------------------------
##plot relationship between reads, genotyping, heterozygosity, and ploidy confidence
library(ggplot2)
p1<-ggplot(csv_dump_new_genos,aes(x=OTcoverage/1E6,y=minAltLLR,color=X.GT))+geom_point()+xlab("Million Reads (On Target)")+ylab("Minimum Alternate LLR")+
  scale_color_gradientn(limits=c(0,100),breaks=c(seq(0,100,by=25)),colours = rainbow(5))+
  theme_classic()+theme(legend.title = element_text(size=16, face="bold"),
                        legend.position=c(0.94,0.89),
                        legend.text = element_text(size=16, face="bold"),
                        axis.text=element_text(size=24,color="black"),
                        axis.title=element_text(size=28,face="bold"),
                        axis.title.y = element_text(margin = margin(t = 0, r = 6, b = 0, l = 0)),
                        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))+
                        labs(color = "% Geno")
p2<-ggplot(csv_dump_new_genos,aes(x=OTcoverage/1E6,y=minAltLLR,color=X.GT))+geom_point()+
  scale_color_gradientn(colours = rainbow(5))+
  scale_x_continuous(breaks=c(seq(0,0.1,by=0.05)),limits=c(0,0.1))+
  scale_y_continuous(breaks=c(seq(0,20,by=10)),limits=c(0,20))+
  theme(axis.text=element_text(size=14,color="black"),
                               panel.background = element_rect(fill="white"),
                               panel.border = element_blank(),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               axis.line = element_line(size = 0.5, linetype = "solid",colour = "black"),
                               legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())#+
p1+annotation_custom(ggplotGrob(p2), xmin = 1.1, xmax = 1.85, ymin = -10, ymax = 50)

p1<-ggplot(csv_dump_new_genos,aes(x=OTcoverage/1E6,y=minAltLLR,color=heterozygosity))+geom_point()+xlab("Million Reads (On Target)")+ylab("Minimum Alternate LLR")+
  scale_color_gradientn(colours = rainbow(5))+
#  scale_color_gradientn(limits=c(0,1),breaks=c(seq(0,1,by=0.25)),colours = rainbow(5))+
  theme_classic()+theme(legend.title = element_text(size=16, face="bold"),
                        legend.position=c(0.92,0.89),
                        legend.text = element_text(size=16, face="bold"),
                        axis.text=element_text(size=24,color="black"),
                        axis.title=element_text(size=28,face="bold"),
                        axis.title.y = element_text(margin = margin(t = 0, r = 6, b = 0, l = 0)),
                        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))+
                        labs(color = "Heterozyg.")
p2<-ggplot(csv_dump_new_genos,aes(x=OTcoverage/1E6,y=minAltLLR,color=heterozygosity))+geom_point()+
  scale_color_gradientn(colours = rainbow(5))+
  scale_x_continuous(breaks=c(seq(0,0.1,by=0.05)),limits=c(0,0.1))+
  scale_y_continuous(breaks=c(seq(0,20,by=10)),limits=c(0,20))+
  theme(axis.text=element_text(size=14,color="black"),
        panel.background = element_rect(fill="white"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid",colour = "black"),
        legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())#+
p1+annotation_custom(ggplotGrob(p2), xmin = 1.1, xmax = 1.85, ymin = -10, ymax = 50)

ggplot(csv_dump_new_genos,aes(x=OTcoverage/1E6,y=X.GT,color=heterozygosity))+geom_point()+
  xlab("Million Reads (On Target)")+ylab("Percent Genotyping Success")+
  scale_x_continuous(breaks=c(seq(0,1.75,by=0.25)),limits=c(0,1.75),labels=c("0.0","",0.5,"","1.0","",1.5,""))+
  scale_color_gradientn(name="Heterozygosity",colours = rainbow(5))+
  theme_classic()+theme(legend.title = element_text(size=16, face="bold"),
                        legend.position=c(0.7,0.5),
                        legend.text = element_text(size=16, face="bold"),
                        axis.text=element_text(size=24,color="black"),
                        axis.title=element_text(size=28,face="bold"),
                        axis.title.y = element_text(margin = margin(t = 0, r = 6, b = 0, l = 0)),
                        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))+
  geom_vline(xintercept=0.1,color="black",linetype="dotted",size=1)

cor.test(csv_dump_new_genos$OTreads,csv_dump_new_genos$minAltLLR,method="pearson")
cor.test(csv_dump_new_genos$OTreads,csv_dump_new_genos$X.GT,method="pearson")
cor.test(csv_dump_new_genos$OTreads,csv_dump_new_genos$heterozygosity,method="pearson")
cor.test(csv_dump_new_genos$X.GT,csv_dump_new_genos$minAltLLR,method="pearson")
cor.test(csv_dump_new_genos$heterozygosity,csv_dump_new_genos$minAltLLR,method="pearson")




#######----------------------------------------------------------------------------------------------------------------
####allele plots with 4N individuals

names(csv_dump_new_genos)[!names(csv_dump_new_genos)%in%loci]
keep4Ninds<-csv_dump_new_genos[csv_dump_new_genos$OTcoverage>50000&csv_dump_new_genos$minAltLLR>25&csv_dump_new_genos$ploidy==4,"fastq_name"]
length(keep4Ninds)
sum(rownames(refCounts) %in% keep4Ninds)

refCounts4Ninds<-refCounts[rownames(refCounts) %in% as.character(keep4Ninds),]
ncol(refCounts4Ninds)
nrow(refCounts4Ninds)
altCounts4Ninds<-altCounts[rownames(altCounts) %in% keep4Ninds,]
ncol(altCounts4Ninds)
nrow(altCounts4Ninds)

#n=1
library(ggplot2)
for (n in 1:length(loci)){
  df<-as.data.frame(cbind(refCounts4Ninds[,loci[n]],altCounts4Ninds[,loci[n]]))
  colnames(df)<-c("ref","alt")
  p1<-ggplot(data=df,aes(x=ref))+geom_point(data=df,aes(y=alt))+
    ylim(0,max(c(df$ref,df$alt)))+xlim(0,max(c(df$ref,df$alt)))+
    xlab("Number REF allele reads")+ylab("Number ALT allele reads")+
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill="white"),
          axis.line = element_line(size = 0.5, linetype = "solid",colour = "black"),
          axis.text.y=element_text(size=22, color="black"),
          axis.text.x=element_text(size=22, color="black"),
          axis.title=element_text(size=28,face="bold"),
          axis.title.y = element_text(margin = margin(t = 0, r = 6, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
  plot1name<-paste(fpLoc$Ind[n],".pdf",sep="")
  pdf(plot1name,width = 10, height = 10) 
  plot(p1, size = 1, fast = FALSE, labels = TRUE,cex=0.5)
  dev.off()
  
  div<-8#set boundaries for octosomic to compare to observed (putatively tetrasomic) allele ratios
  sd<-0.05*(4/div)
  dfOCTO<-df
  dfOCTO$x<-df$ref
  dfOCTO$A1HOMtarget<-dfOCTO$x*0
  dfOCTO$HET1target<-dfOCTO$x*((1/div)/(1-(1/div)))
  dfOCTO$HET2target<-dfOCTO$x*((2/div)/(1-(2/div)))
  dfOCTO$HET3target<-dfOCTO$x*((3/div)/(1-(3/div)))
  dfOCTO$HET4target<-dfOCTO$x*((4/div)/(1-(4/div)))
  dfOCTO$HET5target<-dfOCTO$x*((5/div)/(1-(5/div)))
  dfOCTO$HET6target<-dfOCTO$x*((6/div)/(1-(6/div)))
  dfOCTO$HET7target<-dfOCTO$x*((7/div)/(1-(7/div)))
  #  head(dfOCTO)
  
  octo.95quants<-c(qnorm(0.975, mean = 0, sd = sd),
                   qnorm(0.025, mean = 1*(1/div), sd = sd),
                   qnorm(0.975, mean = 1*(1/div), sd = sd),
                   qnorm(0.025, mean = 2*(1/div), sd = sd),
                   qnorm(0.975, mean = 2*(1/div), sd = sd),
                   qnorm(0.025, mean = 3*(1/div), sd = sd),
                   qnorm(0.975, mean = 3*(1/div), sd = sd),
                   qnorm(0.025, mean = 4*(1/div), sd = sd),
                   qnorm(0.975, mean = 4*(1/div), sd = sd),
                   qnorm(0.025, mean = 5*(1/div), sd = sd),
                   qnorm(0.975, mean = 5*(1/div), sd = sd),
                   qnorm(0.025, mean = 6*(1/div), sd = sd),
                   qnorm(0.975, mean = 6*(1/div), sd = sd),
                   qnorm(0.025, mean = 7*(1/div), sd = sd),
                   qnorm(0.975, mean = 7*(1/div), sd = sd),
                   qnorm(0.025, mean = 1, sd = sd))
  
  dfOCTO$A2HOMb<-octo.95quants[1]/(1-octo.95quants[1])*dfOCTO$x
  dfOCTO$A2HOMbALT<-(octo.95quants[1]/(1-octo.95quants[1])*dfOCTO$x)/(1+(1/(div-3)))
  dfOCTO$HET7Lb<-octo.95quants[2]/(1-octo.95quants[2])*dfOCTO$x
  dfOCTO$HET7Ub<-octo.95quants[3]/(1-octo.95quants[3])*dfOCTO$x
  dfOCTO$HET6Lb<-octo.95quants[4]/(1-octo.95quants[4])*dfOCTO$x
  dfOCTO$HET6Ub<-octo.95quants[5]/(1-octo.95quants[5])*dfOCTO$x
  dfOCTO$HET5Lb<-octo.95quants[6]/(1-octo.95quants[6])*dfOCTO$x
  dfOCTO$HET5Ub<-octo.95quants[7]/(1-octo.95quants[7])*dfOCTO$x
  dfOCTO$HET4Lb<-octo.95quants[8]/(1-octo.95quants[8])*dfOCTO$x
  dfOCTO$HET4Ub<-1/(dfOCTO$HET4Lb/dfOCTO$x)*dfOCTO$x
  dfOCTO$HET3Lb<-1/(dfOCTO$HET5Ub/dfOCTO$x)*dfOCTO$x
  dfOCTO$HET3Ub<-1/(dfOCTO$HET5Lb/dfOCTO$x)*dfOCTO$x
  dfOCTO$HET2Lb<-1/(dfOCTO$HET6Ub/dfOCTO$x)*dfOCTO$x
  dfOCTO$HET2Ub<-1/(dfOCTO$HET6Lb/dfOCTO$x)*dfOCTO$x
  dfOCTO$HET1Lb<-1/(dfOCTO$HET7Ub/dfOCTO$x)*dfOCTO$x
  dfOCTO$HET1Ub<-1/(dfOCTO$HET7Lb/dfOCTO$x)*dfOCTO$x
  dfOCTO$A1HOMb<-1/(dfOCTO$A2HOMb/dfOCTO$x)*dfOCTO$x
  dfOCTO$A1HOMbALT<-1/(dfOCTO$A2HOMbALT/dfOCTO$x)*dfOCTO$x
  dfOCTO[1,]<-0
  #head(dfOCTO)
  
  p2<-ggplot(data=dfOCTO,aes(x=x))+geom_point(data=dfOCTO,aes(y=alt))+
    ylim(0,max(c(dfOCTO$x,dfOCTO$alt)))+xlim(0,max(c(dfOCTO$x,dfOCTO$alt)))+
    xlab("Number REF allele reads")+ylab("Number ALT allele reads")+
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill="white"),
          axis.line = element_line(size = 0.5, linetype = "solid",colour = "black"),
          axis.text.y=element_text(size=22, color="black"),
          axis.text.x=element_text(size=22, color="black"),
          axis.title=element_text(size=28,face="bold"),
          axis.title.y = element_text(margin = margin(t = 0, r = 6, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))+
    geom_line(data=dfOCTO,aes(y=A1HOMtarget),color="black",linetype="dotted",size=1)+
    geom_line(data=dfOCTO,aes(y=HET1target),color="black",linetype="dotted",size=1)+
    geom_line(data=dfOCTO,aes(y=HET2target),color="black",linetype="dotted",size=1)+
    geom_line(data=dfOCTO,aes(y=HET3target),color="black",linetype="dotted",size=1)+  
    geom_line(data=dfOCTO,aes(y=HET4target),color="black",linetype="dotted",size=1)+  
    geom_line(data=dfOCTO,aes(y=HET5target),color="black",linetype="dotted",size=1)+  
    geom_line(data=dfOCTO,aes(y=HET6target),color="black",linetype="dotted",size=1)+  
    geom_line(data=dfOCTO,aes(y=HET7target),color="black",linetype="dotted",size=1)+  
    geom_vline(xintercept=0,color="black",linetype="dotted",size=1)+
    #geom_line(data=dfOCTO,aes(y=A2HOMb),color="blue",linetype="dashed")+
    geom_line(data=dfOCTO,aes(y=A2HOMbALT),color="blue",linetype="dotted")+
    geom_line(data=dfOCTO,aes(y=A1HOMbALT),color="blue",linetype="dotted")+
    #geom_line(data=dfOCTO,aes(y=A1HOMb),color="blue",linetype="dashed")+
    geom_line(data=dfOCTO,aes(y=HET1Ub),color="red",linetype="dashed")+
    geom_line(data=dfOCTO,aes(y=HET1Lb),color="red",linetype="dashed")+
    geom_line(data=dfOCTO,aes(y=HET2Lb),color="blue",linetype="dashed")+
    geom_line(data=dfOCTO,aes(y=HET2Ub),color="blue",linetype="dashed")+
    geom_line(data=dfOCTO,aes(y=HET3Lb),color="red",linetype="dashed")+
    geom_line(data=dfOCTO,aes(y=HET3Ub),color="red",linetype="dashed")+
    geom_line(data=dfOCTO,aes(y=HET4Lb),color="blue",linetype="dashed")+
    geom_line(data=dfOCTO,aes(y=HET4Ub),color="blue",linetype="dashed")+
    geom_line(data=dfOCTO,aes(y=HET5Lb),color="red",linetype="dashed")+
    geom_line(data=dfOCTO,aes(y=HET5Ub),color="red",linetype="dashed")+
    geom_line(data=dfOCTO,aes(y=HET6Lb),color="blue",linetype="dashed")+
    geom_line(data=dfOCTO,aes(y=HET6Ub),color="blue",linetype="dashed")+
    geom_line(data=dfOCTO,aes(y=HET7Lb),color="red",linetype="dashed")+
    geom_line(data=dfOCTO,aes(y=HET7Ub),color="red",linetype="dashed")
  
  plot2name<-paste(fpLoc$Ind[n],"-limits.pdf",sep="")
  pdf(plot2name,width = 10, height = 10) 
  plot(p2, size = 1, fast = FALSE, labels = TRUE,cex=0.5)
  dev.off()
}



#######----------------------------------------------------------------------------------------------------------------
#####get sample meta information, create combined data frames

##get site/locale information
progeny<-read.csv("Archived Sturgeon Samples in Hagerman-FDR.csv",header=TRUE)
progeny$Individual.name<-as.character(progeny$Individual.name)
progeny$newStage<-rep("UNK",times=nrow(progeny))
progeny[(!is.na(progeny$Length..Fork..mm.))&progeny$Length..Fork..mm.<=300,"newStage"]<-"YOY"
progeny[(!is.na(progeny$Length..Fork..mm.))&progeny$Length..Fork..mm.>300&progeny$Length..Fork..mm.<1360,"newStage"]<-"JUV"
progeny[(!is.na(progeny$Length..Fork..mm.))&progeny$Length..Fork..mm.>=1360,"newStage"]<-"Adult"
#fix some names
progeny[grep("SIRE",progeny$Individual.name,value=FALSE),"Individual.name"]<-gsub("SIRE ", "SIRE-",progeny[grep("SIRE",progeny$Individual.name,value=FALSE),"Individual.name"])
progeny[grep("DAM",progeny$Individual.name,value=FALSE),"Individual.name"]<-gsub("DAM ", "DAM-",progeny[grep("DAM",progeny$Individual.name,value=FALSE),"Individual.name"])

rownames(new_genos)<-new_genos$samples
rownames(progeny)<-progeny$Individual.name
#create data frame of samples filtered for duplicates, with site information
new_genos_sites<-merge(new_genos,progeny[progeny$Individual.name %in% new_genos$samples,],by="row.names")

#write.table(new_genos_sites,file="sequenced_Atr325_samples_2019.txt",row.names=FALSE,col.names=TRUE,sep="\t")

#classified sites by reach and region by hand; re-import
#new_genos_sites<-read.csv("sequenced_Atr325_samples_2019.csv",header=TRUE)

##import additional meta-data
progeny<-read.csv("WST_progeny_export_April2020.csv",header=TRUE)
names(progeny)
progeny$Individual.name<-as.character(progeny$Individual.name)

#fix some names
progeny[grep("SIRE",progeny$Individual.name,value=FALSE),"Individual.name"]<-gsub("SIRE ", "SIRE-",progeny[grep("SIRE",progeny$Individual.name,value=FALSE),"Individual.name"])
progeny[grep("DAM",progeny$Individual.name,value=FALSE),"Individual.name"]<-gsub("DAM ", "DAM-",progeny[grep("DAM",progeny$Individual.name,value=FALSE),"Individual.name"])
progeny[grep("DonBP",progeny$Individual.name,value=FALSE),"Individual.name"]<-gsub("DonBP", "DONBP",progeny[grep("DonBP",progeny$Individual.name,value=FALSE),"Individual.name"])

row.names(progeny)<-progeny$Individual.name
row.names(new_genos_sites)<-new_genos_sites$samples
new_genos_sites<-new_genos_sites[!colnames(new_genos_sites)%in%"Row.names"]
##create dataframe with samples filtered for duplicates, sites, and meta-data
new_genos_sites_meta<-merge(new_genos_sites,progeny[progeny$Individual.name %in% new_genos_sites$samples,],by="row.names")
new_genos_sites_meta<-new_genos_sites_meta[!colnames(new_genos_sites_meta) %in% "Row.names"]

csv_dump<-csv_dump[!colnames(csv_dump)%in%"Row.names"]
row.names(new_genos_sites_meta)<-new_genos_sites_meta$genos
row.names(csv_dump)<-csv_dump$fastq_name
##create data frame of samples filtered for duplicates, sites, meta-data, and genotypes 
new_genos_sites_meta_seq<-merge(new_genos_sites_meta,csv_dump,by="row.names")
nrow(new_genos_sites_meta_seq)

#write.table(new_genos_sites_meta_seq,file="all325_new_genos_sites_meta_seq.txt",col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")



#######----------------------------------------------------------------------------------------------------------------
##subset and output for Polygene with true population information

##Polygene format:
#(tab-delimited)
#ID Pop Ploidy Loc1 Loc2
#Ind1  pop1  4  2,3,4  1,2,4
#Ind2  pop1  4  4      2,3,4
#Ind3  pop2  4  1,4    1,2,3,4

# new_genos_sites_meta_seq<-read.delim(file="all325_new_genos_sites_meta_seq.txt",header=TRUE)
loci<-grep("Atr",colnames(new_genos_sites_meta_seq),value=TRUE)
csv_keep<-new_genos_sites_meta_seq[!(is.na(new_genos_sites_meta_seq$ploidy)|new_genos_sites_meta_seq$X.GT<80|new_genos_sites_meta_seq$minAltLLR<10),]
csv_keep$pop<-paste("pop",csv_keep$Reach,csv_keep$ploidy,sep="-")
csv_GenosNF<-csv_keep[,c("samples","pop","ploidy",loci)]
csv_GenosNF<-apply(csv_GenosNF,2,function(x) as.character(x))
header_csv_GenosNF<-c("ID","Pop","Ploidy",names(nucs))
csv_GenosNF<-rbind(header_csv_GenosNF,csv_GenosNF)
#write.table(csv_GenosNF,file="all325_80p.polygeneB.txt",col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")

###--->calculate relatedness, within population, for filtering, in Polygene


#######----------------------------------------------------------------------------------------------------------------
##import ML relatedness from Polygene; filter by relatedness; write another Polygene file to calculate population statistics

pop_matrix_all<- read.table("all325_80p.polygeneB-o_relatednessB.txt",header=FALSE,row.names=1,sep="\t")
pops<-grep("pop",rownames(pop_matrix_all),value=TRUE)
pops_lines<-which(rownames(pop_matrix_all)%in%pops)
pops_lines<-c(pops_lines,(nrow(pop_matrix_all)+1))

#get list of hatchery ("stocked") samples from meta-data
names(new_genos_sites_meta_seq)
table(new_genos_sites_meta_seq$Stocked)
stocked_samples<-new_genos_sites_meta_seq[new_genos_sites_meta_seq$Stocked=="yes","samples"]
stocked_samples<-stocked_samples[as.character(stocked_samples) %in% rownames(pop_matrix_all)[!rownames(pop_matrix_all)%in%pops]]

pop_matrices<-list()
##store population matrices in list indexed by population name
for(p in 1:length(pops)){
  pop_matrix<-as.data.frame(pop_matrix_all[pops_lines[p]:(pops_lines[p+1]-1),1:((pops_lines[p+1]-1)-pops_lines[p])],stringsAsFactors=FALSE)
  colnames(pop_matrix)<-as.character(unlist(pop_matrix[1,]))
  pop_matrix<-tail(pop_matrix,((pops_lines[p+1]-1)-pops_lines[p]))
  pop_matrices[[pops[p]]]<-pop_matrix
}

##loop through population matrices and find individuals related above "cutoff"; save in "keep" or "toss" lists indexed by population name; save all exclusions in "exclude"
dup_cut<-0.8#unknown duplicates
cutoff<-0.2
exclude<-NULL
keep<-list()
toss<-list()
duplicates<-list()
for(p in 1:length(pops)){
  pop_matrix<-pop_matrices[[pops[p]]]
  exclude.pop<-rownames(pop_matrix)[rownames(pop_matrix)%in%stocked_samples]
  exclude<-unique(c(exclude,exclude.pop))
  #pop_matrix[!(rownames(pop_matrix)%in%stocked_samples),!(colnames(pop_matrix)%in%stocked_samples)]
  pop_matrix2<-as.data.frame(pop_matrix[!(rownames(pop_matrix)%in%stocked_samples),!(colnames(pop_matrix)%in%stocked_samples)],stringsAsFactors=FALSE)
  colnames(pop_matrix2)<-colnames(pop_matrix)[!(colnames(pop_matrix)%in%stocked_samples)]
  rownames(pop_matrix2)<-rownames(pop_matrix)[!(rownames(pop_matrix)%in%stocked_samples)]
  pop_matrix<-pop_matrix2
  if(nrow(pop_matrix)>1){
    pop_matrix[upper.tri(pop_matrix,diag=TRUE)]<-NA#only keep one, not symmetric, comparison; also remove self on diagonal
    pop_matrix<-apply(pop_matrix[,1:ncol(pop_matrix)],2,function(x) as.numeric(as.character(x)))
    rownames(pop_matrix)<-colnames(pop_matrix)
    relatedness_pop_stack<-as.data.frame(matrix(ncol=3,nrow=0))
    for(c in 1:ncol(pop_matrix)){
      relatedness_pop_stack<-rbind(relatedness_pop_stack,cbind(unlist(pop_matrix[,c]),as.character(rownames(pop_matrix)),rep(colnames(pop_matrix)[c],times=nrow(pop_matrix))))
    }
    nrow(relatedness_pop_stack)
    #relatedness_pop_stack<-cbind(stack(pop_matrix),rep(colnames(pop_matrix),times=ncol(pop_matrix)))
    colnames(relatedness_pop_stack)<-c("value","row","col")
    relatedness_pop_stack$value<-as.numeric(as.character(relatedness_pop_stack$value))
    relatedness_pop_stack<-relatedness_pop_stack[!is.na(relatedness_pop_stack$value),]
    nrow(relatedness_pop_stack)
    #head(relatedness_pop_stack)
    #ggplot(relatedness_pop_stack,aes(x=value))+geom_density()
    
    samples<-unique(c(as.character(relatedness_pop_stack$col),as.character(relatedness_pop_stack$row)))
#    s=1
    for(s in 1:length(samples)){
      #identify duplicates
      df<-relatedness_pop_stack[relatedness_pop_stack$col%in%samples[s]&relatedness_pop_stack$value>=dup_cut,]
      if(nrow(df)>0){
        duplicates[[pops[p]]]<-rbind(duplicates[[pops[p]]],df)
        exclude.pop<-unique(c(exclude.pop,samples[s]))#should only exclude one duplicate b/c lower diagonal of matrix only
      }
      
      # ##if stocked sample; now already above
      # if(samples[s]%in%stocked_samples){
      #   exclude.pop<-unique(c(exclude.pop,samples[s]))
      # }
      
      ##add sample to exclude list if matrix contains samples above cutoff compared to sample
      if(!samples[s]%in%exclude.pop){
        df<-relatedness_pop_stack[relatedness_pop_stack$col%in%samples[s]&relatedness_pop_stack$value>=cutoff,]
        if(nrow(df)>0){
          exclude.pop<-unique(c(exclude.pop,samples[s]))
        }
      }
    }
    #length(exclude.pop)
    exclude<-unique(c(exclude,exclude.pop))
    keep[[pops[p]]]<-samples[!samples%in%exclude.pop]
    toss[[pops[p]]]<-exclude.pop
  } else if(nrow(pop_matrix)==1){
    exclude<-unique(c(exclude,exclude.pop))
    keep[[pops[p]]]<-as.character(rownames(pop_matrix)[1])
    toss[[pops[p]]]<-exclude.pop
  } else if(nrow(pop_matrix)==0){
    exclude<-unique(c(exclude,exclude.pop))
#    keep[[pops[p]]]<-NULL
    toss[[pops[p]]]<-exclude.pop
  }
}
length(exclude)
exclude0.2<-exclude
#write.table(exclude0.2,file="all325_80p-exclude0.2.txt",col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")


#keep only samples genotyped at >80% and ploidy confidence of 10
csv_keep<-new_genos_sites_meta_seq[!(is.na(new_genos_sites_meta_seq$ploidy)|new_genos_sites_meta_seq$X.GT<80|new_genos_sites_meta_seq$minAltLLR<10),]

#exclude0.2<-read.table(file="all325_80p-exclude0.2.txt",row.names=1,sep="\t")
nrow(exclude0.2)
exclude0.2<-exclude0.2$x
sum(csv_keep$samples%in%exclude0.2)
#keep only 4N/8N individuals and those related below 0.2
csv_GenosNF<-csv_keep[(!csv_keep$samples%in%exclude0.2)&csv_keep$ploidy<5,c("samples","pop","ploidy",loci)]
csv_GenosNF<-apply(csv_GenosNF,2,function(x) as.character(x))
header_csv_GenosNF<-c("ID","Pop","Ploidy",names(nucs))
csv_GenosNF<-rbind(header_csv_GenosNF,csv_GenosNF)
##write Polygene file
write.table(csv_GenosNF,file="all325_filtered.polygene.txt",col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")



#######----------------------------------------------------------------------------------------------------------------
#####adegenet

#keep only samples genotyped at >80% and ploidy confidence of 10
#csv_keep<-new_genos_sites_meta_seq[!(is.na(new_genos_sites_meta_seq$ploidy)|new_genos_sites_meta_seq$X.GT<80|new_genos_sites_meta_seq$minAltLLR<10),]
loci<-grep("Atr",names(csv_keep),value=T)
csv_GenosNF_genepop<-csv_keep[(!csv_keep$samples%in%exclude0.2)&csv_keep$ploidy<5,c("samples",loci)]
csv_GenosNF_genepop[,loci]<-apply(csv_GenosNF_genepop[,loci],2,function(x)as.character(x))
csv_GenosNF_genepop[csv_GenosNF_genepop==""]<-"0,0,0,0"
pops_ass<-droplevels(csv_keep[(!csv_keep$samples%in%exclude0.2)&csv_keep$ploidy<5,"Reach"])
table(pops_ass)

library(adegenet)
rownames(csv_GenosNF_genepop)<-csv_GenosNF_genepop$samples
csv_GenosNF_genepop<-csv_GenosNF_genepop[,loci]
loci<-gsub("\\.","_",loci)
colnames(csv_GenosNF_genepop)<-loci
mydata<-df2genind(csv_GenosNF_genepop,ploidy=4,sep=",", NA.char = "0")
mydata@pop<-as.factor(pops_ass)
mydata

#PCA
X <- scaleGen(mydata, NA.method="mean")
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
# pca1
# pca1$li
write.table(pca1$li,file="pca_table.txt",append=FALSE,quote=FALSE,sep="\t",col.names=T,row.names=T)
write.table(pca1$co,file="pca_table_loci.txt",append=FALSE,quote=FALSE,sep="\t",col.names=T,row.names=T)

barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))
pca1$eig[1]/sum(pca1$eig)
pca1$eig[2]/sum(pca1$eig)
col <- funky(nPop(mydata)*2)[sample(nPop(mydata)*2,nPop(mydata))]
s.class(pca1$li, pop(mydata), col=transp(col,1), axesell=FALSE, cpoint=0, grid=FALSE, label=levels(mydata@pop))
s.class(pca1$li, pop(mydata), xax = 1, yax = 2, col=transp(col,.8), axesell=FALSE, cstar=0, cpoint=1, grid=FALSE, label=NULL)
#s.label(pca1$li)

df<-cbind(pca1$li,mydata@pop)
library(ggplot2)
ggplot(df,aes(x=Axis1,y=Axis2,color=mydata@pop))+geom_point()+ylim(min(df$Axis1),max(df$Axis1))

##Heterozygosity
sum_mydata<-summary(mydata)
write.table(sum_mydata$Hobs,file="Hobs_table.txt",append=FALSE,quote=FALSE,sep="\t",col.names=T,row.names=T)
write.table(sum_mydata$Hexp,file="Hexp_table.txt",append=FALSE,quote=FALSE,sep="\t",col.names=T,row.names=T)

#inbreeding
inb_mydata<-inbreeding(mydata,pop=pop(mydata),res.type="estimate",N=500,M=N*100)
ggplot(as.data.frame(inb_mydata),aes(x=inb_mydata))+geom_density()
df_inb<-as.data.frame(cbind(as.numeric(as.character(inb_mydata)),as.character(pop(mydata)),rownames(mydata@tab)))
names(df_inb)<-c("inbreeding","pop","sample")
df_inb$inbreeding<-as.numeric(as.character(df_inb$inbreeding))
ggplot(df_inb,aes(x=pop,y=as.numeric(inbreeding)))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90))

#Fst-analog
library(hierfstat)
matFst <- pairwise.fst(mydata)#Nei 1973
as.vector(as.matrix(matFst))


#jackknife/bootstrap Fst-analog
prop<-0.5
Fst_matrices<-list()
for (n in 1:100)  {
    #Bootstrap  
#  sample_loci<-loci[sample(1:length(loci),length(loci),replace=TRUE)]
    #Jackknife
    sample_loci<-loci[sample(1:length(loci),(prop*length(loci)),replace=FALSE)]

  csv_GenosNF_genepop_sample<-csv_GenosNF_genepop[,sample_loci]
#  csv_GenosNF_genepop[1:5,1:10]
#  csv_GenosNF_genepop_sample[1:5,1:10]
#  loci<-gsub("\\.","_",loci)
  colnames(csv_GenosNF_genepop_sample)<-sample_loci
  mydata1<-df2genind(csv_GenosNF_genepop_sample,ploidy=4,sep=",", NA.char = "0")
  mydata1@pop<-as.factor(pops_ass)
  Fst_matrices[[n]] <- pairwise.fst(mydata1)
}
Jackknife<-Fst_matrices
#Bootstrap<-Fst_matrices
Jackknife.mat<-lapply(Jackknife,function(x) as.matrix(x))
#Bootstrap.mat<-lapply(Bootstrap,function(x) as.matrix(x))
write.table(Jackknife.mat,file="Fst_jackknife.txt",append=FALSE,quote=FALSE,sep="\t",col.names=T,row.names=T)
#write.table(Bootstrap.mat,file="Fst_bootstrap.txt",append=FALSE,quote=FALSE,sep="\t",col.names=T,row.names=T)




#######----------------------------------------------------------------------------------------------------------------
##minor allele frequency (MAF) by locus by population

table(csv_keep$Reach)
csv_keep$pop<-csv_keep$Reach
head(csv_keep)
#keep only 4N/8N individuals and those related below 0.2
csv_GenosNF<-as.data.frame(csv_keep[(!csv_keep$samples%in%exclude0.2)&csv_keep$ploidy<5,c("samples","pop","ploidy",loci)])
pops<-unique(csv_GenosNF$pop)
maf_by_locus_by_pop<-as.data.frame(matrix(nrow=0,ncol=(length(loci)+2)))
for(p in 1:length(pops)){
  df<-csv_GenosNF[csv_GenosNF$pop==pops[p],]
  df[,loci]<-apply(df[,loci],2,function(x) as.character(x))
  maf_by_locus_by_pop[p,]<-c(as.character(pops[p]),nrow(df),unlist(apply(df[,loci],2,function(x) sum(unlist(strsplit(x,","))=="2")/length(unlist(strsplit(x,","))))))
}
colnames(maf_by_locus_by_pop)<-c("pop","N",loci)
maf_by_locus_by_pop[,loci]<-apply(maf_by_locus_by_pop[,loci],2,function(x) as.numeric(as.character(x)))
maf_by_locus_by_pop$N<-as.numeric(maf_by_locus_by_pop$N)
#only keep localities with samples size 8+
maf_by_locus_by_pop_8<-droplevels(maf_by_locus_by_pop[maf_by_locus_by_pop$N>7,])
#calculate mean MAF across Columbia-only localities
mean_MAF<-colMeans(maf_by_locus_by_pop_8[!maf_by_locus_by_pop_8$pop%in%c("Sacramento","Fraser"),loci])
maf_by_locus_by_pop_8<-maf_by_locus_by_pop_8[,c(1,2,(order(mean_MAF)+2))]
maf_by_locus_by_pop_8$inside_COL<-"Columbia"
#if mean MAF >0.5, switch which is "minor allele" in Columbia Basin
mean_MAF_corr<-mean_MAF
mean_MAF_corr[mean_MAF_corr>0.5]<-1-mean_MAF_corr[mean_MAF_corr>0.5]
boxplot(mean_MAF_corr)
#order loci by corrected mean MAF
mean_MAF_corrord<-mean_MAF_corr[order(mean_MAF_corr)]

#create a stacked dataframe for plotting
maf_stack<-stack(maf_by_locus_by_pop_8[,names(mean_MAF_corrord)])
colnames(maf_stack)<-c("MAF","locus")
#head(maf_stack)
maf_stack$pop<-rep(maf_by_locus_by_pop_8$pop,times=length(loci))
maf_stack$num<-rep(seq(1:length(loci)),each=nrow(maf_by_locus_by_pop_8))
maf_stack$inside_COL<-rep(maf_by_locus_by_pop_8$inside_COL,times=length(loci))
maf_stack$mean_MAF_corr<-as.numeric(as.character(rep(mean_MAF_corrord,each=length(maf_by_locus_by_pop_8$pop))))

library(ggplot2)
ggplot(data=maf_stack,aes(x=num,y=MAF))+geom_point(aes(color=inside_COL))+theme_classic()+
  geom_line(aes(x=num,y=mean_MAF_corr,color="red"),size=1)+scale_color_manual(values=c("darkgoldenrod2", "green4", "black"))+
  xlab("Locus Rank by Mean MAF")+ylab("Minor Allele Frequency")+theme(legend.position = "none")+
  scale_y_continuous(breaks=c(seq(0,0.8,by=0.1)))+
  theme(axis.text=element_text(size=22,color="black"),
        axis.title=element_text(size=28,face="bold"),
        axis.title.y = element_text(margin = margin(t = 0, r = 6, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))



######----------------------------------------------------------------------------------------------------------------
##calculate Mendelian incompatibilities ("percent mismatch") by locus and by trio

require(tidyr)
require(gtools)
require(ggplot2)

head(csv_GenosNF[,1:10])
colnames(csv_GenosNF)<-csv_GenosNF[1,]
csv_GenosNF<-csv_GenosNF[2:nrow(csv_GenosNF),]

pop<-grep("SIRE",csv_GenosNF$Sample,value=T)
mom<-grep("DAM",csv_GenosNF$Sample,value=T)
kids<-grep("x",csv_GenosNF$Sample,value=T)
genos_both_unique<-rbind(csv_GenosNF[csv_GenosNF$Sample%in%pop,],
                         csv_GenosNF[csv_GenosNF$Sample%in%mom,],
                         csv_GenosNF[csv_GenosNF$Sample%in%kids,])

#genos_both_unique<-read.csv("compiled_genos_test.csv",header=TRUE)
#genos_both_unique$Sample<-gsub("\\.genos","",genos_both_unique$Sample)
head(genos_both_unique[,1:10])
nrow(genos_both_unique)
names(genos_both_unique)
loci<-grep("Atr",names(genos_both_unique),value=T)
names(genos_both_unique)[!names(genos_both_unique)%in%loci]
rownames(genos_both_unique)<-as.character(genos_both_unique$Sample)
colnames(genos_both_unique)<-gsub(" ","\\.",colnames(genos_both_unique))
colnames(genos_both_unique)<-gsub("%","X\\.",colnames(genos_both_unique))
colnames(genos_both_unique)[!colnames(genos_both_unique)%in%loci]<-gsub("-","\\.",colnames(genos_both_unique)[!colnames(genos_both_unique)%in%loci])

meta<-as.data.frame(cbind(genos_both_unique$Sample,c(rep("adult",times=(length(pop)+length(mom))),rep("offspring",times=length(kids))),c(rep("male",times=length(pop)),rep("female",times=length(mom)),rep("offspring",times=length(kids)))))
colnames(meta)<-c("Sample","lifestage","sex")
#meta<-read.csv("compiled_genos_test_metadata.csv",header=TRUE)
head(meta)
rownames(meta)<-as.character(meta$Sample)
sum(rownames(meta)%in%rownames(genos_both_unique))
sum(rownames(genos_both_unique)%in%rownames(meta))
table(meta$lifestage)
table(meta$sex)

##EXECUTE THIS LINE TO RUN WITHOUT SEX (all by all parent trio comparison)
#meta$sex<-NULL

genos_both_unique<-genos_both_unique[rownames(genos_both_unique)%in%rownames(meta),]
nrow(genos_both_unique)
genos_both_unique<-merge(meta,genos_both_unique,by="row.names")
names(genos_both_unique)
genos_both_unique$Sample<-genos_both_unique$Sample.x
genos_both_unique$ploidy<-as.numeric(as.character(genos_both_unique$Ploidy.Result))
genos_both_unique<-genos_both_unique[!colnames(genos_both_unique) %in% c("Row.names","Sample.x","Sample.y","Ploidy.Result")]
names(genos_both_unique)
loci<-grep("Atr",colnames(genos_both_unique),value=TRUE)
genos_both_unique[,loci] = apply(genos_both_unique[,loci], 2, function(x) as.character(x));

##was this imported (opened/saved) from Excel? If so the 0's may have been botched so the sum will be >0
sum(genos_both_unique[,loci]=="0")

#missing.genos<-c("0","00","000","0000","00000","000000","0000000","00000000","000000000","0000000000","00000000000","000000000000")
maxP<-max(genos_both_unique$ploidy)
missing.genos<-NULL
for(N in 1:maxP){
  missing.genos<-c(missing.genos,paste(rep("0",times=N),collapse=""))
}
for(r in 1:nrow(genos_both_unique)){
  genos_both_unique[r,loci] = gsub("^0$",missing.genos[as.numeric(as.character(genos_both_unique[r,"ploidy"]))],genos_both_unique[r,loci])
}
sum(genos_both_unique[,loci]=="0")#should now be zero



##we have to translate the genotypes so they are all in alphab. order, regardless of ref/alt allele
##we do this because we predict mismatches based on potential genotypes from each parent pairs' gamete combinations (whether offspring is included)
##pairings of gametes could great all different orders of genotypes; thus we translate them to a standard format (alphab.)

##check classes of heterozygote genotypes, before translation (substitution)
genotypes<-list()
genotypes_all<-NULL
for (l in 1:length(loci)){
  df<-as.data.frame(table(genos_both_unique[,loci[l]]))
  genotypes[[l]]<-as.character(df$Var1)
  genotypes_all<-c(genotypes_all,genos_both_unique[,loci[l]])
}
#check how many genotype classes are currently included
table(genotypes_all)

##create translation table so all heterozygotes of the same nucleotides are the same (alphabetical); now for range of ploidies in input dataset
ploidy<-unique(as.numeric(genos_both_unique$ploidy))
ploidy<-ploidy[!is.na(ploidy)]
ploidy<-ploidy[order(ploidy)]
nucleotides<-rbind(c("A","C"),c("A","G"),c("A","T"),c("C","G"),c("C","T"),c("G","T"))
het.translation<-as.data.frame(matrix(nrow=0,ncol=2))
for (p in ploidy) {
  gamA<-round(p/2-0.001,0)
  gamB<-round(p/2+0.001,0)
  for (n in 1:nrow(nucleotides)){
    AA<-as.data.frame(permutations(n=2,r=gamA,v=nucleotides[n,],repeats.allowed=TRUE))
    two.nuc<-unlist(lapply(as.list(data.frame(t(AA))), function(x) paste(x, collapse="")))
    AB<-as.data.frame(permutations(n=2,r=gamB,v=nucleotides[n,],repeats.allowed=TRUE))
    three.nuc<-unlist(lapply(as.list(data.frame(t(AB))), function(x) paste(x, collapse="")))
    combo.nuc<-unique(c(two.nuc,three.nuc))
    B<-as.data.frame(permutations(n=length(combo.nuc),r=2,v=combo.nuc,repeats.allowed=TRUE))
    paste.nuc<-unique(unlist(lapply(as.list(data.frame(t(B))), function(x) paste(x, collapse=""))))
    keep<-as.character(lapply(as.list(paste.nuc), function(x) length(unlist(strsplit(as.character(x),"")))))==p
    paste.nuc<-paste.nuc[keep]
    paste.nuc<-paste.nuc[!paste.nuc %in% c(paste(rep(nucleotides[n,1],times=p),collapse=""),paste(rep(nucleotides[n,2],times=p),collapse=""))]
    het.translation.sub<-as.data.frame(matrix(nrow=0,ncol=2))
    for (t in 1:length(paste.nuc)){
      df<-as.data.frame(table(strsplit(as.character(paste.nuc[t]),"")))
      het.translation.sub<-rbind(het.translation.sub,cbind(paste.nuc[t],paste( c( as.character(rep(df$Var1[1],times=df$Freq[1])) , as.character(rep(df$Var1[2],times=df$Freq[2])) ),collapse="" ) ))
      
    }
    het.translation.sub<-het.translation.sub[order(het.translation.sub$V2),]
    het.translation<-rbind(het.translation,het.translation.sub)
    print(nrow(het.translation))
  }
}
colnames(het.translation)<-c("orig","named")
head(het.translation)
tail(het.translation)

##replace alternative heterozygotes with alphabetical heterozygotes (e.g. AAAT for TAAA or AATA)
het.translation[,c(1,2)] = apply(het.translation[,c(1,2)], 2, function(x) as.character(x));
genos_both_unique[,loci] = apply(genos_both_unique[,loci], 2, function(x) as.character(x));
#g=1
for(g in 1:nrow(het.translation)){
  print(paste(c("finding ",as.character(het.translation[g,"orig"])," and replacing with ",as.character(het.translation[g,"named"])),collapse=""))
  genos_both_unique[genos_both_unique==as.character(het.translation[g,"orig"])]<-as.character(het.translation[g,"named"])
}

##check classes of heterozygote genotypes are alphabetical
genotypes<-list()
genotypes_all<-NULL
for (l in 1:length(loci)){
  df<-as.data.frame(table(genos_both_unique[,loci[l]]))
  genotypes[[l]]<-as.character(df$Var1)
  genotypes_all<-c(genotypes_all,genos_both_unique[,loci[l]])
}
#should now be many fewer classes than before; with the alleles in alphab. order
table(genotypes_all)





##subset offspring genotypes by missingness
kids_genos<-droplevels(genos_both_unique[genos_both_unique$lifestage=="offspring",])
nrow(kids_genos)
names(kids_genos)
kids_genos$X.GT<-as.numeric(as.character(kids_genos$X.GT))
boxplot(kids_genos$X.GT)
# ggplot(data=kids_genos,aes(y=as.numeric(as.character(X.GT))))+geom_boxplot()
# ggplot(data=kids_genos,aes(x=as.numeric(as.character(OTreads))))+geom_density()
# ggplot(data=kids_genos,aes(x=as.numeric(as.character(OTreads)),y=as.numeric(as.character(X.GT))))+geom_point()
# ggplot(data=kids_genos,aes(x=as.numeric(as.character(OTreads)),y=as.numeric(as.character(X.GT)),color=mom_sample))+geom_point()

kids_genos<-kids_genos[kids_genos$X.GT>=80,]
#kids_genos<-kids_genos[kids_genos$On.Target.Reads>20000,]
nrow(kids_genos)




#genos_both_unique$sex<-NULL

###parentage: count MI among trios, with sex (females x males) or without sex (all x all adult)
if("sex" %in% colnames(genos_both_unique)){print("SEX included; Making Trio Comparisons")
  
  #subset genotypes by sex and missing data
  mom_genos<-droplevels(genos_both_unique[genos_both_unique$sex=="female",])
  #mom_genos<-droplevels(genos_both_unique[genos_both_unique$sex=="female"&genos_both_unique$X.GT>90,])
  pop_genos<-droplevels(genos_both_unique[genos_both_unique$sex=="male",])
  #pop_genos<-droplevels(genos_both_unique[genos_both_unique$sex=="male"&genos_both_unique$X.GT>90,])
  nrow(mom_genos)
  nrow(pop_genos)
  
  
  
  
  
  ##compare parental combo genotypes with each offspring, make ploidy screen first step, make gametes half of ploidy
  #double reduction not currently considered
  missing.genos<-NULL
  for(N in 1:maxP){
    missing.genos<-c(missing.genos,paste(rep("0",times=N),collapse=""))
  }
  # m=1
  # p=1
  # o=1
  locus_mismatches<-as.data.frame(matrix(nrow=0,ncol=(length(loci)+5)))
  summed_mismatches<-as.data.frame(matrix(nrow=0,ncol=9))
  for (m in 1:nrow(mom_genos)){
    malleles<-strsplit(as.character(t(mom_genos[m,loci])),"")#returns list with vector for each locus
    ##if ploidy is even, half this produces integer gametes; otherwise produce gametes of both integers above/below half-ploidy
    ##i.e. one half of 5N makes 2.5, so produce 2N and 3N gametes
    if((mom_genos[m,"ploidy"]/2)%%1==0){
      malleles.combos<-lapply(malleles, function(x) combn(x,(mom_genos[m,"ploidy"]/2),simplify=FALSE))
      mom.gametes<-lapply(lapply(malleles.combos, function(x) lapply(x, function(x) paste(x, collapse=""))), function(x) unique(unlist(x)))
    } else if(!(mom_genos[m,"ploidy"]/2)%%1==0){
      malleles.combos.1<-lapply(malleles, function(x) combn(x,floor((mom_genos[m,"ploidy"]/2)),simplify=FALSE))
      mom.gametes.1<-lapply(lapply(malleles.combos.1, function(x) lapply(x, function(x) paste(x, collapse=""))), function(x) unique(unlist(x)))
      malleles.combos.2<-lapply(malleles, function(x) combn(x,floor((mom_genos[m,"ploidy"]/2)+1),simplify=FALSE))
      mom.gametes.2<-lapply(lapply(malleles.combos.2, function(x) lapply(x, function(x) paste(x, collapse=""))), function(x) unique(unlist(x)))
      mom.gametes<-mapply(c, mom.gametes.1,mom.gametes.2, SIMPLIFY=FALSE)
    }
    mom.gametes
    
    for (p in 1:nrow(pop_genos)){
      palleles<-strsplit(as.character(t(pop_genos[p,loci])),"")#returns list with vector for each locus
      
      if((pop_genos[p,"ploidy"]/2)%%1==0){
        palleles.combos<-lapply(palleles, function(x) combn(x,(pop_genos[p,"ploidy"]/2),simplify=FALSE))
        pop.gametes<-lapply(lapply(palleles.combos, function(x) lapply(x, function(x) paste(x, collapse=""))), function(x) unique(unlist(x)))
      } else if(!(pop_genos[p,"ploidy"]/2)%%1==0){
        palleles.combos.1<-lapply(palleles, function(x) combn(x,floor((pop_genos[p,"ploidy"]/2)),simplify=FALSE))
        pop.gametes.1<-lapply(lapply(palleles.combos.1, function(x) lapply(x, function(x) paste(x, collapse=""))), function(x) unique(unlist(x)))
        palleles.combos.2<-lapply(palleles, function(x) combn(x,floor((pop_genos[p,"ploidy"]/2)+1),simplify=FALSE))
        pop.gametes.2<-lapply(lapply(palleles.combos.2, function(x) lapply(x, function(x) paste(x, collapse=""))), function(x) unique(unlist(x)))
        pop.gametes<-mapply(c, pop.gametes.1,pop.gametes.2, SIMPLIFY=FALSE)
      }
      pop.gametes
      
      pot.offspr.genos<-list()
      # l=1
      for (l in 1:length(pop.gametes)) {
         raw.pot.offspr.genos<-as.character(unlist(as.list(unite(expand.grid(mom.gametes[[l]],pop.gametes[[l]]),"geno",sep="", remove=TRUE))))
        raw.pot.offspr.genos[raw.pot.offspr.genos %in% het.translation[het.translation$orig %in% raw.pot.offspr.genos,"orig"]] <- as.character(het.translation[het.translation$orig %in% raw.pot.offspr.genos,"named"][match(raw.pot.offspr.genos, het.translation[het.translation$orig %in% raw.pot.offspr.genos,"orig"], nomatch = 0)])
        par.alleles<-unique(c(unique(unlist(strsplit(mom.gametes[[l]],""))),unique(unlist(strsplit(pop.gametes[[l]],"")))))
        if(sum(par.alleles %in% "0")>0) {
          raw.pot.offspr.genos<-missing.genos
        }
        pot.offspr.genos[[l]]<-raw.pot.offspr.genos
      }
      names(pot.offspr.genos)<-loci
      pot.offspr.genos
      
      ###collect mismatches from all; mark those non-parents that also differ by ploidy
      ###check only set of kid_genos consistent with ploidies of parent gamete combinations (allow for 4N, 5N and 6N gametes from 5N parents...)?
      ##alternatively: set all (non-missing?) loci to mismatches if ploidy is mismatched
      #       o=1
      #    o=90
      if(mean(c(mom_genos[m,"ploidy"],pop_genos[p,"ploidy"]))%%1==0){
        poss.ploidy<-mean(c(mom_genos[m,"ploidy"],pop_genos[p,"ploidy"]))
      } else if (!mean(c(mom_genos[m,"ploidy"],pop_genos[p,"ploidy"]))%%1==0){
        poss.ploidy<-c( floor(mean(c(mom_genos[m,"ploidy"],pop_genos[p,"ploidy"]))) , floor(mean(c(mom_genos[m,"ploidy"],pop_genos[p,"ploidy"])))+1 )
      }
      
      for (o in 1:nrow(kids_genos)){
        print(paste(c("Now comparing dam #",m,", ",as.character(mom_genos$Sample[m]),", and sire #",p,", ",as.character(pop_genos$Sample[p]),", with offspring #",o,", ",as.character(kids_genos$Sample[o])),collapse=""))
        genos<-as.list(kids_genos[o,loci])
        #parent.cross<-paste(c(as.character(mom_genos[m,"Sample"]),as.character(pop_genos[p,"Sample"])),collapse=" ")
        
        if(!kids_genos[o,"ploidy"] %in% poss.ploidy) {
          match<-as.numeric(rep("1",times=length(loci)))
          #exclude combinations with 0
          match[as.logical(as.character(unlist(lapply(pot.offspr.genos, function(x) unique(ifelse(x %in% missing.genos,TRUE,FALSE))))))|genos %in% missing.genos]<-NA
          class<-"ploidy_mismatch"
        } else if (kids_genos[o,"ploidy"] %in% poss.ploidy) {
          #        t=2
          #        genos[[t]]
          #        pot.offspr.genos[[t]]
          match<-mapply(function(x, y) ifelse(x %in% y,0,1), genos, pot.offspr.genos)
          #exclude combinations with 0
          match[as.logical(as.character(unlist(lapply(pot.offspr.genos, function(x) unique(ifelse(x %in% missing.genos,TRUE,FALSE))))))|genos %in% missing.genos]<-NA
          class<-"ploidy_match"
        }
        mismatches<-sum(na.omit(match))
        no.surveyed<-sum(!is.na(match))
        perc.mismatches<-mismatches/no.surveyed
        
        locus_mismatches<-rbind(as.matrix(locus_mismatches), c(as.character(kids_genos[o,c("Sample")]) ,as.character(kids_genos[o,c("On.Target.Reads")]),as.character(kids_genos[o,c("X.GT")]),as.character(mom_genos[m,"Sample"]),as.character(pop_genos[p,"Sample"]),match))
        summed_mismatches<-rbind(summed_mismatches,cbind( as.character(kids_genos[o,c("Sample")]) ,as.character(kids_genos[o,c("On.Target.Reads")]),as.character(kids_genos[o,c("X.GT")]),as.character(mom_genos[m,"Sample"]),as.character(pop_genos[p,"Sample"]),mismatches,no.surveyed,perc.mismatches,class))
      }
    }
  }
  colnames(summed_mismatches)<-c("offspring","OTreads","X.GT","compared_adult1","compared_adult2","mismatches","loci.compared","percent.mismatch","ploidy_match_class")
  #head(summed_mismatches)
  locus_mismatches<-as.data.frame(locus_mismatches)
  colnames(locus_mismatches)<-c("offspring","OTreads","X.GT","compared_adult1","compared_adult2",loci)
#  head(locus_mismatches)
#  table(locus_mismatches[,loci])
#  sum(as.numeric(as.character(summed_mismatches$mismatches)))
  
  
  
}else if(!"sex" %in% colnames(genos_both_unique)){print("NO SEX; Making Single Parent Comparisons")
  
  #create two dataframes; one with all adults, one with second through last adults
  #the latter dataframe will be diminished by one with each pass through first dataframe to avoid self-on-self comparisons
  mom_genos<-droplevels(genos_both_unique[genos_both_unique$lifestage=="adult",])
  nrow(mom_genos)
  
  mom_genos$X.GT<-as.numeric(as.character(mom_genos$X.GT))
  mom_genos<-mom_genos[mom_genos$X.GT>=80,]
  nrow(mom_genos)
  mom_genos$Sample
  pop_genos<-tail(mom_genos,nrow(mom_genos)-1)
  nrow(pop_genos)
  pop_genos$Sample
  
  ##compare parental combo genotypes with each offspring, make ploidy screen first step, make gametes half of ploidy
  missing.genos<-NULL
  for(N in 1:maxP){
    missing.genos<-c(missing.genos,paste(rep("0",times=N),collapse=""))
  }
  locus_mismatches<-as.data.frame(matrix(nrow=0,ncol=(length(loci)+1)))
  summed_mismatches<-as.data.frame(matrix(nrow=0,ncol=9))
  for (m in 1:(nrow(mom_genos)-1)){
    malleles<-strsplit(as.character(t(mom_genos[m,loci])),"")#returns list with vector for each locus
    ##if ploidy is even, half this produces integer gametes; otherwise produce gametes of both integers above/below half-ploidy
    ##i.e. one half of 5N makes 2.5, so produce 2N and 3N gametes
    if((mom_genos[m,"ploidy"]/2)%%1==0){
      malleles.combos<-lapply(malleles, function(x) combn(x,(mom_genos[m,"ploidy"]/2),simplify=FALSE))
      mom.gametes<-lapply(lapply(malleles.combos, function(x) lapply(x, function(x) paste(x, collapse=""))), function(x) unique(unlist(x)))
    } else if(!(mom_genos[m,"ploidy"]/2)%%1==0){
      malleles.combos.1<-lapply(malleles, function(x) combn(x,floor((mom_genos[m,"ploidy"]/2)),simplify=FALSE))
      mom.gametes.1<-lapply(lapply(malleles.combos.1, function(x) lapply(x, function(x) paste(x, collapse=""))), function(x) unique(unlist(x)))
      malleles.combos.2<-lapply(malleles, function(x) combn(x,floor((mom_genos[m,"ploidy"]/2)+1),simplify=FALSE))
      mom.gametes.2<-lapply(lapply(malleles.combos.2, function(x) lapply(x, function(x) paste(x, collapse=""))), function(x) unique(unlist(x)))
      mom.gametes<-mapply(c, mom.gametes.1,mom.gametes.2, SIMPLIFY=FALSE)
    }
    mom.gametes
    
    for (p in 1:nrow(pop_genos)){
      palleles<-strsplit(as.character(t(pop_genos[p,loci])),"")#returns list with vector for each locus
      
      if((pop_genos[p,"ploidy"]/2)%%1==0){
        palleles.combos<-lapply(palleles, function(x) combn(x,(pop_genos[p,"ploidy"]/2),simplify=FALSE))
        pop.gametes<-lapply(lapply(palleles.combos, function(x) lapply(x, function(x) paste(x, collapse=""))), function(x) unique(unlist(x)))
      } else if(!(pop_genos[p,"ploidy"]/2)%%1==0){
        palleles.combos.1<-lapply(palleles, function(x) combn(x,floor((pop_genos[p,"ploidy"]/2)),simplify=FALSE))
        pop.gametes.1<-lapply(lapply(palleles.combos.1, function(x) lapply(x, function(x) paste(x, collapse=""))), function(x) unique(unlist(x)))
        palleles.combos.2<-lapply(palleles, function(x) combn(x,floor((pop_genos[p,"ploidy"]/2)+1),simplify=FALSE))
        pop.gametes.2<-lapply(lapply(palleles.combos.2, function(x) lapply(x, function(x) paste(x, collapse=""))), function(x) unique(unlist(x)))
        pop.gametes<-mapply(c, pop.gametes.1,pop.gametes.2, SIMPLIFY=FALSE)
      }
      pop.gametes
      
      pot.offspr.genos<-list()
      # l=1
      for (l in 1:length(pop.gametes)) {
        raw.pot.offspr.genos<-as.character(unlist(as.list(unite(expand.grid(mom.gametes[[l]],pop.gametes[[l]]),"geno",sep="", remove=TRUE))))
        raw.pot.offspr.genos[raw.pot.offspr.genos %in% het.translation[het.translation$orig %in% raw.pot.offspr.genos,"orig"]] <- as.character(het.translation[het.translation$orig %in% raw.pot.offspr.genos,"named"][match(raw.pot.offspr.genos, het.translation[het.translation$orig %in% raw.pot.offspr.genos,"orig"], nomatch = 0)])
        par.alleles<-unique(c(unique(unlist(strsplit(mom.gametes[[l]],""))),unique(unlist(strsplit(pop.gametes[[l]],"")))))
        if(sum(par.alleles %in% "0")>0) {
          raw.pot.offspr.genos<-missing.genos
        }
        pot.offspr.genos[[l]]<-raw.pot.offspr.genos
      }
      names(pot.offspr.genos)<-loci
      pot.offspr.genos
      
      ###collect mismatches from all; mark those non-parents that also differ by ploidy
      ###check only set of kid_genos consistent with ploidies of parent gamete combinations (allow for 4N, 5N and 6N gametes from 5N parents...)?
      ##set all (non-missing) loci to mismatches if ploidy is mismatched
      if(mean(c(mom_genos[m,"ploidy"],pop_genos[p,"ploidy"]))%%1==0){
        poss.ploidy<-mean(c(mom_genos[m,"ploidy"],pop_genos[p,"ploidy"]))
      } else if (!mean(c(mom_genos[m,"ploidy"],pop_genos[p,"ploidy"]))%%1==0){
        poss.ploidy<-c( floor(mean(c(mom_genos[m,"ploidy"],pop_genos[p,"ploidy"]))) , floor(mean(c(mom_genos[m,"ploidy"],pop_genos[p,"ploidy"])))+1 )
      }
      poss.ploidy
      
      for (o in 1:nrow(kids_genos)){
        print(paste(c("Now comparing dam #",m,", ",as.character(mom_genos$Sample[m]),", and sire #",p,", ",as.character(pop_genos$Sample[p]),", with offspring #",o,", ",as.character(kids_genos$Sample[o])),collapse=""))
        genos<-as.list(kids_genos[o,loci])
        
        if(!kids_genos[o,"ploidy"] %in% poss.ploidy) {
          match<-as.numeric(rep("1",times=length(loci)))
          #exclude combinations with 0
          match[as.logical(as.character(unlist(lapply(pot.offspr.genos, function(x) unique(ifelse(x %in% missing.genos,TRUE,FALSE))))))|genos %in% missing.genos]<-NA
          class<-"ploidy_mismatch"
        } else if (kids_genos[o,"ploidy"] %in% poss.ploidy) {
           match<-mapply(function(x, y) ifelse(x %in% y,0,1), genos, pot.offspr.genos)
          #exclude combinations with 0
          match[as.logical(as.character(unlist(lapply(pot.offspr.genos, function(x) unique(ifelse(x %in% missing.genos,TRUE,FALSE))))))|genos %in% missing.genos]<-NA
          class<-"ploidy_match"
        }
        mismatches<-sum(na.omit(match))
        no.surveyed<-sum(!is.na(match))
        perc.mismatches<-mismatches/no.surveyed
        mismatches
        no.surveyed
        perc.mismatches
        
        summed_mismatches<-rbind(summed_mismatches,cbind( as.character(kids_genos[o,c("Sample")]) ,as.character(kids_genos[o,c("On.Target.Reads")]),as.character(kids_genos[o,c("X.GT")]),as.character(mom_genos[m,"Sample"]),as.character(pop_genos[p,"Sample"]),mismatches,no.surveyed,perc.mismatches,class))
        locus_mismatches<-rbind(as.matrix(locus_mismatches), c(as.character(kids_genos[o,c("Sample")]) ,as.character(kids_genos[o,c("On.Target.Reads")]),as.character(kids_genos[o,c("X.GT")]),as.character(mom_genos[m,"Sample"]),as.character(pop_genos[p,"Sample"]),match))
        
      }
    }
    pop_genos<-tail(pop_genos,nrow(pop_genos)-1)
    nrow(pop_genos)
    pop_genos$Sample
  }
  colnames(summed_mismatches)<-c("offspring","OTreads","X.GT","compared_adult1","compared_adult2","mismatches","loci.compared","percent.mismatch","ploidy_match_class")
  locus_mismatches<-as.data.frame(locus_mismatches)
  colnames(locus_mismatches)<-c("offspring","OTreads","X.GT","compared_adult1","compared_adult2",loci)
  
}
head(summed_mismatches)
head(locus_mismatches)

#write.table(summed_mismatches,file="Mendelian-mismatches.txt",append=FALSE,col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t");

summed_mismatches$percent.mismatch<-as.numeric(as.character(summed_mismatches$percent.mismatch))
ggplot(data=summed_mismatches,aes(x=percent.mismatch))+geom_density(adjust=0.25)

df<-summed_mismatches[!summed_mismatches$ploidy_match_class=="ploidy_mismatch",]
nrow(df)
ggplot(data=df,aes(x=percent.mismatch))+geom_density(adjust=0.25)


#####Calculate percent MI (mismatches) per locus based on known parentage

summed_mismatches$adult1_short<-as.character(unlist(as.data.frame(strsplit(as.character(summed_mismatches$compared_adult1),"-"))[3,]))
summed_mismatches$adult2_short<-as.character(unlist(as.data.frame(strsplit(as.character(summed_mismatches$compared_adult2),"-"))[3,]))
summed_mismatches$adult1_short[summed_mismatches$adult1_short=="0E5E"]<-"OE5E"
summed_mismatches$adult2_short[summed_mismatches$adult2_short=="0E5E"]<-"OE5E"
summed_mismatches$parent1<-as.character(unlist(as.data.frame(strsplit(as.character(unlist(as.data.frame(strsplit(as.character(summed_mismatches$offspring),"_"))[2,])),"x"))[1,]))
summed_mismatches$parent2<-as.character(unlist(as.data.frame(strsplit(as.character(unlist(as.data.frame(strsplit(as.character(summed_mismatches$offspring),"_"))[2,])),"x"))[2,]))
head(summed_mismatches)
#which((summed_mismatches$adult1_short==summed_mismatches$parent1&summed_mismatches$adult2_short==summed_mismatches$parent2)|(summed_mismatches$adult1_short==summed_mismatches$parent2&summed_mismatches$adult2_short==summed_mismatches$parent1))
locus_mismatches_parents<-locus_mismatches[which((summed_mismatches$adult1_short==summed_mismatches$parent1&summed_mismatches$adult2_short==summed_mismatches$parent2)|(summed_mismatches$adult1_short==summed_mismatches$parent2&summed_mismatches$adult2_short==summed_mismatches$parent1)),]
head(locus_mismatches_parents)
locus.percent.mismatches<-apply(locus_mismatches_parents[,loci],2,function(x) sum(as.numeric(as.character(x)),na.rm=T)/sum(!is.na(x)) )

#write.table(locus.percent.mismatches,file="Mendelian-mismatches-by-locus.txt",append=FALSE,col.names = FALSE,row.names = TRUE,quote=FALSE,sep="\t");

###single parent mismatches

parent_genos<-droplevels(genos_both_unique[genos_both_unique$lifestage=="adult",])
loci<-grep("Atr",colnames(parent_genos),value=TRUE)
parent_genos$X.GT<-as.numeric(as.character(parent_genos$X.GT))
boxplot(parent_genos$X.GT)
#parent_genos<-parent_genos[parent_genos$X.GT>=80,]
nrow(parent_genos)

missing.genos<-NULL
for(N in 1:maxP){
  missing.genos<-c(missing.genos,paste(rep("0",times=N),collapse=""))
}
summed_mismatches_single<-as.data.frame(matrix(nrow=0,ncol=8))
for (p in 1:nrow(parent_genos)){
  palleles<-strsplit(as.character(t(parent_genos[p,loci])),"")#returns list with vector for each locus

  if((parent_genos[p,"ploidy"]/2)%%1==0){
    palleles.combos<-lapply(palleles, function(x) combn(x,(parent_genos[p,"ploidy"]/2),simplify=FALSE))
    parent.gametes<-lapply(lapply(palleles.combos, function(x) lapply(x, function(x) paste(x, collapse=""))), function(x) unique(unlist(x)))
  } else if(!(parent_genos[p,"ploidy"]/2)%%1==0){
    palleles.combos.1<-lapply(palleles, function(x) combn(x,floor((pop_genos[p,"ploidy"]/2)),simplify=FALSE))
    parent.gametes.1<-lapply(lapply(palleles.combos.1, function(x) lapply(x, function(x) paste(x, collapse=""))), function(x) unique(unlist(x)))
    palleles.combos.2<-lapply(palleles, function(x) combn(x,floor((pop_genos[p,"ploidy"]/2)+1),simplify=FALSE))
    parent.gametes.2<-lapply(lapply(palleles.combos.2, function(x) lapply(x, function(x) paste(x, collapse=""))), function(x) unique(unlist(x)))
    parent.gametes<-mapply(c, parent.gametes.1,parent.gametes.2, SIMPLIFY=FALSE)
  }
  parent.gametes
  
  #  o=1
  for (o in 1:nrow(kids_genos)){
    print(paste(c("Now comparing parent #",p,", ",as.character(parent_genos$Sample[p]),", with offspring #",o,", ",as.character(kids_genos$Sample[o])),collapse=""))
    
    oalleles<-strsplit(as.character(t(kids_genos[o,loci])),"")#returns list with vector for each locus
    
    if((kids_genos[o,"ploidy"]/2)%%1==0){
      oalleles.combos<-lapply(oalleles, function(x) combn(x,kids_genos[o,"ploidy"]/2,simplify=FALSE))
      offspring.gametes<-lapply(lapply(oalleles.combos, function(x) lapply(x, function(x) paste(x, collapse=""))), function(x) unique(unlist(x)))
      
     } else if(!(kids_genos[o,"ploidy"]/2)%%1==0){
      oalleles.combos.1<-lapply(oalleles, function(x) combn(x,floor(kids_genos[o,"ploidy"]/2),simplify=FALSE))
      offspring.gametes.1<-lapply(lapply(oalleles.combos.1, function(x) lapply(x, function(x) paste(x, collapse=""))), function(x) unique(unlist(x)))

      oalleles.combos.2<-lapply(oalleles, function(x) combn(x,floor((kids_genos[o,"ploidy"]/2)+1),simplify=FALSE))
      offspring.gametes.2<-lapply(lapply(oalleles.combos.2, function(x) lapply(x, function(x) paste(x, collapse=""))), function(x) unique(unlist(x)))

      offspring.gametes<-mapply(c, offspring.gametes.1,offspring.gametes.2, SIMPLIFY=FALSE)
    }
    offspring.gametes
    
    #exclude combinations with 0
    match<-mapply(function(x, y) ifelse(x %in% y,1,0), offspring.gametes, parent.gametes)
    match<-unlist(lapply(lapply(match, sum), function(x) ifelse(x>0,0,1)))
    match[as.logical(as.character(unlist(lapply(parent.gametes, function(x) unique(ifelse(x %in% missing.genos,TRUE,FALSE))))))|as.logical(as.character(unlist(lapply(offspring.gametes, function(x) unique(ifelse(x %in% missing.genos,TRUE,FALSE))))))]<-NA
    
    if(kids_genos[o,"ploidy"]==parent_genos[p,"ploidy"]){class<-"ploidy_match"}
    if(!kids_genos[o,"ploidy"]==parent_genos[p,"ploidy"]){class<-"ploidy_mismatch_large"}
    if(abs(kids_genos[o,"ploidy"]-parent_genos[p,"ploidy"])==1){class<-"ploidy_mismatch_small"}
    mismatches<-sum(na.omit(match))
    no.surveyed<-sum(!is.na(match))
    perc.mismatches<-mismatches/no.surveyed
    summed_mismatches_single<-rbind(summed_mismatches_single,cbind( as.character(kids_genos[o,c("Sample")]) ,as.character(kids_genos[o,c("On.Target.Reads")]),as.character(kids_genos[o,c("X.GT")]),as.character(parent_genos[p,"Sample"]),mismatches,no.surveyed,perc.mismatches,class))
  }
}
colnames(summed_mismatches_single)<-c("offspring","OTreads","X.GT","compared_adult","mismatches","loci.compared","percent.mismatch","ploidy_match_class")
head(summed_mismatches_single)

summed_mismatches_single$X.GT<-as.numeric(as.character(summed_mismatches_single$X.GT))
summed_mismatches_single$percent.mismatch<-as.numeric(as.character(summed_mismatches_single$percent.mismatch))
ggplot(data=summed_mismatches_single,aes(x=percent.mismatch,fill=ploidy_match_class))+geom_density(adjust=0.25)

df<-summed_mismatches_single[!summed_mismatches_single$ploidy_match_class=="ploidy_mismatch_large",]
nrow(df)
ggplot(data=df,aes(x=percent.mismatch,fill=ploidy_match_class))+geom_density(adjust=0.25)

df<-summed_mismatches_single[summed_mismatches_single$ploidy_match_class=="ploidy_match",]
nrow(df)
ggplot(data=df,aes(x=percent.mismatch,fill=ploidy_match_class))+geom_density(adjust=0.1)

df<-summed_mismatches_single[summed_mismatches_single$ploidy_match_class=="ploidy_mismatch_small",]
nrow(df)
ggplot(data=df,aes(x=percent.mismatch,fill=ploidy_match_class))+geom_density(adjust=0.1)

#write.table(summed_mismatches_single,file="single-parent-mismatches.txt",append=FALSE,col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")
