#creation of circosplots
#Tjitske Los


# you need to specify the Path to file directories and the suffix of the files.
# or you can specify filenames, copynumbersbed, callsbed, name, translocation file, mutationfile, circos plot file name

# probalby you need to make some alterations to make the script work for you, but you can use some snippits from this script
# more options and extensive manual can be found at:
# https://jokergoo.github.io/circlize_book/book/

# install.packages("Biobase")
# install.packages("circlize")

suppressMessages(library(Biobase))
suppressMessages(library(circlize))

# create directory for the plots
circosdir<-"circos/"
dir.create(circosdir)

# read in a document with the samplenames of SR50 and PE150 data - to find data of multiple samples
Samplenamesfile<-"samplenames.txt"
samplenames<-read.delim(Samplenamesfile, stringsAsFactors=F, sep=" ")

#paths to file directories
Path_to_bedfiles_cn<-"../SR50/100kbp/BED/"
Path_to_translocation_file<-"merged/"
Path_to_mutationfile<-"vcf/CNA_added/"

for(i in 1:length(samplenames[,1])){

	#find files or specify them yourself
		copynumbersbed<-paste(Path_to_bedfiles_cn,samplenames[i,"sampleCNA"], "-copynumbers.bed", sep="")
	callsbed<-paste0(Path_to_bedfiles_cn, samplenames[i,"sampleCNA"], "_CNAs.bed")
	name<-samplenames[i,2]
	summary<-paste0(Path_to_translocation_file,samplenames[i,1],"-trl_summary_brkpt_freq.csv")

	#mutationfile<-paste0(Path_to_mutationfile, list.files(Path_to_mutationfile, pattern=paste0(samplenames[i,1], ".HIGHeffidence_CNA_filtered.csv")))
	mutationfile<-paste0(Path_to_mutationfile, samplenames[i,1], ".CCG_134_Functional_Mutations_CNA_filt.csv")
	circlize<-paste0(circosdir,samplenames[i,1],"_circlize.png")

	#read in data
	copynumbers<-read.delim(copynumbersbed, skip=1, header=FALSE)
	bedcalls<-read.delim(callsbed)
	translocations<-read.delim(summary,stringsAsFactors = FALSE, header = TRUE)
	mutations<-read.delim(mutationfile, stringsAsFactors=FALSE, header=TRUE)

	#add chr prefix to copynumber data to be equeal to circlize levels
	copynumbers[,1]<-paste("chr",copynumbers[,1],sep="")
	bedcalls<-bedcalls[,c(1:3,6)]  #chromosome, start, end and segmentvalue

	if(nrow(bedcalls) != 0){
		bedcalls[,1]<-paste("chr",bedcalls[,1],sep="")
		bedcalls[bedcalls$chromo=="chr23","chromo"]<-"chrX"
	}

	#create bedfiles of translocations
	side_one<-translocations[translocations$TPorFP=="TP",c("CHROM","POS","POS","GENE")]
	side_one[,3]<-side_one[,3]+1

	side_two<-translocations[translocations$TPorFP=="TP",c("CHROM2","POS2","POS2","GENE2")]
	side_two[,3]<-side_two[,3]+1
	names(side_two)<-names(side_one)

	genes<-rbind(side_one,side_two)

	if(nrow(genes) != 0){genes[,5]<-"translocation"
		#remove additional info of genes behind dash "-"
		genes$GENE<-sub("-.*","",genes$GENE)
		#genes[,1]<-paste0("chr",genes[,1])
	}else{genes[,5]<-character()}

	#remove duplicate genes from dataframe
	genes<-genes[!duplicated(genes[,4]),]

	# reformat mutation data
	mut_genes<-mutations[,c("CHROM", "POS", "POS", "ANN.0..GENE", "ANN.0..IMPACT")]
	mut_genes[,3]<-mut_genes[,3]+1
	names(mut_genes)<-names(genes)

	#remove duplicate genes from dataframe
	mut_genes<-mut_genes[!duplicated(mut_genes[,4]),]

	#combine with translocation genes
	genes<-rbind(genes,mut_genes)

	#order the genes
	chrOrder<-c(paste("chr",1:22,sep=""),"chrX","chrY")
	genes$CHROM<-factor(genes$CHROM, levels=chrOrder)
	genes<-genes[order(genes$CHROM, genes$POS),]

	png(circlize,width=10, height=10, unit="in" ,res=500)
	#create circosplot
	circos.par("start.degree" = 90, track.margin=c(0,0),points.overflow.warning=FALSE)
	par(cex=0.9, mar=c(5.1,4.1,4.1,2.1), oma=c(1,0,1,0)) #,oma outer margin bottom, left, top, right   pin=width height plot dimensions origional 5.749583 5.159377
	#initialize empty circos plot
	circos.initializeWithIdeogram(plotType = "NULL")

	#place genes around ideogram
	if(nrow(genes)!=0){
	    circos.genomicLabels(genes, labels.column = 4, side="outside", padding=0, cex=1.3, col=ifelse(genes[,5]=="HIGH","#231F20",
	        ifelse(genes[,5]=="translocation","#009F4D",ifelse(genes[,5]=="MODERATE","#231F20", "#EDA366")))) #"MODERATE", "A8AEB5"
	}
	#create ideogram with labels
	circos.trackPlotRegion(ylim = c(0, 1.8), bg.border = NA, track.height = 0.05,
	    panel.fun = function(x, y) {
	        xlim = get.cell.meta.data("xlim")
	        chr = get.cell.meta.data("sector.index")
		circos.rect(xlim[1],-0.55,xlim[2],2)
	        circos.text(mean(xlim),1, labels = gsub(pattern = "chr", x=chr,replacement = ""), facing = "bending.inside", niceFacing = TRUE, cex=0.8)
	    })
	circos.genomicIdeogram(track.height=0.03)

	# The copynumber track
	circos.genomicTrack(copynumbers,numeric.column=5,bg.border=NA,ylim=c(-2,3),track.height=0.3, panel.fun=function(region, value,...){
		circos.genomicPoints(region, value, pch = 20, cex = 0.2, col="#00000010",  ... )
	})

	#copynumber calls - plot over copynumber profile
	if(nrow(bedcalls) != 0){
		circos.genomicTrack(bedcalls, numeric.column=4,bg.border=NA, ylim=c(-2,3), panel.fun=function(region, value,...){
		circos.genomicRect(region, ytop=3.4, ybottom=-2, col=ifelse(value[[1]]>1,"#FF000033","#0000FF4D"), border=NA,track.index=max(get.all.track.index())-1 )
		})
	}

	#plot translocations
	if(nrow(side_one)!=0){circos.genomicLink(side_one, side_two , rou=0.5)}

	#create legends

	legend("bottomleft", legend = c("Mutation", "Translocation"), fill=c("#231F20",  "#009F4D"),box.lty = 0, border=NULL, title="Genes", title.adj=0, xpd=TRUE) #"bottomleft"x=par("usr")[1], y=par("usr")[3]#"#A8AEB5","#EDA366",
	legend("bottomright", legend = c("Copy number gain", "Copy number loss"),fill=c( "#FF000033","#0000FF4D"),box.lty =  0, border=NULL, title="Copy number profile" ,title.adj=0, xpd=TRUE) #"bottomright"x=par("usr")[2], y=par("usr")[3] "#D6A3A27F","#4392CC7F"
	#legend("bottomleft", legend = c("High impact mutations", "Moderate impact mutation"),
	#        fill=c("red", "blue"),box.lty = 0, border=NULL, title="Genes")
	title(name) #, line=6)outer=TRUE
	dev.off()
}
