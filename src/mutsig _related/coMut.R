## Nicolas Stransky
## The Broad Institute of MIT and Harvard / Cancer program.
## stransky@broadinstitute.org

## the file is edited by Qing Zhang to be adapted for new version of mutsig


## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.

## You should have received a copy of the GNU Lesser General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################

library(optparse)
library(R.matlab)
library(RColorBrewer)

print("packages done")


option_list <- list(                                                                                                                         
		make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
				help="Print more output [default %default]"),
		make_option(c("--full.verbose"), action="store_true", default=FALSE,
				help="Print even more output [default %default]"),
		make_option(c("-o", "--outputdir"), action="store", type="character", default=".",
				help="Output directory", metavar="dir.path"),
		make_option(c("--firehose.mode"), action="store_true", default=FALSE,
				help="Run coMut as a FH task [default FALSE]"),
		make_option(c("-f", "--firehose.output"), action="store", type="character",
				help="Firehose output directory [default %default]", default="/local/cga-fh/cga/", metavar="dir.path"),
		make_option(c("-w", "--firehose.workspace"), action="store", type="character", default=NULL,
				help="Firehose workspace"),
		make_option(c("-r", "--firehose.mutsigrun"), action="store", type="character",
				help="MutSig version. e.g. mutsig1.0 mutsig1.5 etc. [default %default]", default="mutsig1.5"),
		make_option(c("--firehose.mutsig.mutcategs"), action="store", type="character", default=NULL,
				help="Full path to the mutcategs.txt file from a MutSig run. Only required in firehose.mode. [default FALSE]"),
		make_option(c("-p", "--mutsig.full.path"), action="store", type="character",
				help="Full path of a MutSig output directory (overrides automated guessing of Firehose directory)", default=NULL, metavar="dir.path"),
		make_option(c("-a", "--analysis.set"), action="store", type="character",
				help="Name of the analysis sample set"),
		make_option(c("-q", "--qthresh"), action="store", type="numeric", default=.1,
				help="q-value threshold of genes to display [default %default]", metavar="threshold"),
		make_option(c("-b", "--blacklist.file"), action="store", type="character",
				help="List of genes to ignore", default=NULL, metavar="file.path"),
		make_option(c("-l", "--whitelist.file"), action="store", type="character",
				help="Force the list of genes to appear on the plot", default=NULL, metavar="file.path"),
		make_option(c("-s", "--significance.file"), action="store", type="character",
				help="Provide the list of significant genes and q-values. Tab-delimited genes and q-values, no header. (overrides MutSig sig_genes)", 
				default=NULL, metavar="file.path"),
		make_option(c("-c", "--coverage.file"), action="store", type="character",
				help="Provide coverage for all samples. Tab-delimited sample and number of bases covered, no header. (overrides MutSig patients.counts_and_rates)", 
				default=NULL, metavar="file.path"),
		make_option(c("-m", "--maf"), action="store", type="character",
				help="MAF file containing all the required fields. (overrides MutSig maf)", 
				default=NULL, metavar="file.path", dest="maf.file"),
		make_option(c("--reference.genome"), action="store", type="character",
				help="Path to an alternate reference genome [default %default]", 
				default="/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta", metavar="file.path"),
		make_option(c("--no.allelic.fraction.boxplot"), action="store_false", default=TRUE,
				help="Suppress allelic fraction boxplots [default FALSE]", dest="allelic.fraction.boxplot"),
		make_option(c("--no.mutation.spectrum.plot"), action="store_false", default=TRUE,
				help="Suppress mutation spectrum plots [default FALSE]", dest="mutation.spectrum.plot"),
		make_option(c("--sort.by.mutation.status"), action="store_true", default=FALSE,
				help="Sort by mutation status of each gene instead of mutation frequency of each sample [default %default]"),
		make_option(c("--sort.genes.by.prevalence"), action="store_true", default=FALSE,
				help="Sort genes by gene mutation prevalence instead of q-value [default %default]"),
		make_option(c("--png"), action="store_true", default=FALSE,
				help="Output a PNG image (rasterized format) instead of a PDF (vectorial, editable) [default %default]")
)

opt <- parse_args(OptionParser(option_list=option_list, usage = "Rscript %prog [options]"), print_help_and_exit=FALSE)

if (opt[["help"]]) {
	print_help(OptionParser(option_list=option_list, usage = "Rscript %prog [options]"))
	cat("There are three modes to run coMut.R:
					- MutSig mode: coMut needs to know where MutSig results are located in the Firehose directory structure. This is if you want to run coMut on the 
					command line, to visualize MutSig results. To do so, the FH workspace (-w) and the FH analysis set (-a) are needed.
					- Firehose mode: Like the above, but don't guess anything, use all paths provided.
					- Standalone mode: coMut does not assume that some results have been precomputed in MutSig and recalculates them. 
					For that mode, you need to tell coMut where the maf file (-m), the coverage file (-c) and the list of q-values (-s) files are located, 
					as well as the name of your analysis set (-a, can be anything). If you don't provide a workspace or the complete path to the MutSig output, 
					it is assumed that you are running the standalone mode.
					
					Examples
					- Firehose mode
					Rscript coMut.R -v -o <output_dir> --firehose.mode -a <analysis_set> -s sig_genes.txt -c patients.counts_and_rates.txt -m final_analysis_set.maf --firehose.mutsig.mutcategs mutcategs.txt
					
					- MutSig mode
					Rscript coMut.R -v -o <output_dir> -a <analysis_set> -w <workspace>
					
					- Standalone mode
					Rscript coMut.R -v -o <output_dir> -a <analysis_set> -s sig_genes.txt -c patients.counts_and_rates.txt -m final_analysis_set.maf
					
					")
	quit(status = 1)
}


## Options parsing
if (FALSE) {
    mutsig.output = ""
	mutsigrun <- "mutsig1.5"
	verbose <- TRUE
	full.verbose <- TRUE
	outputdir <- "."
	firehose.output <- ""
	firehose.mode <- TRUE
	mutsig.full.path <- "."
	workspace <- NULL
	qthresh <- .1
	allelic.fraction.boxplot <- TRUE
	mutation.spectrum.plot <- TRUE
	sort.by.mutation.status <- TRUE
	sort.genes.by.prevalence <- FALSE
	reference.genome <- ""
	png.output <- FALSE
	analysis.set <- "test" # to name the output graph
	significance.file <- "sig_genes.txt"
	maf.file <- "final_analysis_set.maf"
	coverage.file <- "patient_counts_and_rates.txt"
	blacklist.file <- ""
	whitelist.file <- ""
	firehose.mutsig.mutcategs <- "mutcategs.txt"
}

min.reads.to.count.allelic.fraction <- 20
max.mutation.rate <- 100

### read directly from directory





##### parse the options input

verbose <- opt$verbose
full.verbose <- opt$full.verbose
outputdir <- opt$outputdir
firehose.output <- opt$firehose.output
workspace <- opt$firehose.workspace
mutsigrun <- opt$firehose.mutsigrun
firehose.mode <- opt$firehose.mode
firehose.mutsig.mutcategs <- opt$firehose.mutsig.mutcategs
qthresh <- opt$qthresh



analysis.set <- opt$analysis.set
# append the full mutsig path to the each of the input files
mutsig.full.path <- opt$mutsig.full.path


maf.file <- ifelse(is.null(opt$maf.file), "", opt$maf.file)
maf.file <- paste(mutsig.full.path, maf.file, sep = "/")

blacklist.file <- ifelse(is.null(opt$blacklist.file), "", opt$blacklist.file) # file path
if (!blacklist.file == "") blacklist.file <- paste(mutsig.full.path, blacklist.file, sep = "/")

whitelist.file <- ifelse(is.null(opt$whitelist.file), "", opt$whitelist.file) # file path
whitelist.file <- paste(mutsig.full.path, whitelist.file, sep = "/")
if (!whitelist.file == "") whitelist.file <- paste(mutsig.full.path, whitelist.file, sep = "/")

significance.file <- ifelse(is.null(opt$significance.file), "", opt$significance.file) # flie path
significance.file <- paste(mutsig.full.path, significance.file, sep = "/")

coverage.file <- ifelse(is.null(opt$coverage.file), "", opt$coverage.file) # file path
coverage.file <- paste(mutsig.full.path, coverage.file, sep = "/")

firehose.mutsig.mutcategs <- paste(mutsig.full.path, firehose.mutsig.mutcategs, sep = "/")


print(list(maf.file, significance.file, coverage.file, firehose.mutsig.mutcategs))


reference.genome <- opt$reference.genome
allelic.fraction.boxplot <- opt$allelic.fraction.boxplot
mutation.spectrum.plot <- opt$mutation.spectrum.plot
sort.by.mutation.status <- opt$sort.by.mutation.status
sort.genes.by.prevalence <- opt$sort.genes.by.prevalence
png.output <- opt$png

if(full.verbose) verbose <- TRUE
mutcategs <- read.delim(firehose.mutsig.mutcategs)
if (verbose) cat("coMut v2.2 -- co-ocurrence of Mutations\n")

suppressPackageStartupMessages(require(RColorBrewer))
suppressPackageStartupMessages(require(gplots))

## Check arguments
if (is.null(analysis.set)) stop("analysis.set is missing")
if (firehose.mode & is.null(firehose.mutsig.mutcategs)) stop("firehose.mutsig.mutcategs is required with --firehose.mode")

## usuallty not anything from white and blacklist

if (file.exists(blacklist.file)) {
	if (verbose) cat("Reading blacklist")
	blacklist <- read.table(blacklist.file)[,1]
	if (verbose) cat(".\n")
} else {
	if (blacklist.file!="") warning("blacklist not found!")
	blacklist <- c()
}

if (file.exists(whitelist.file)) {
	if (verbose) cat("Reading whitelist")
	forceSelectedGenes <- read.table(whitelist.file)[,1]
	if (verbose) cat(".\n")
} else {
	if (whitelist.file!="") warning("whitelist not found!")
}



####### sig input ##########

mutsig.output = ""
if (verbose) cat("Reading significance values")
if (file.exists(significance.file)) {
	if(any(c("gene", "name", "q", "q-value", "qvalue")%in%read.delim(significance.file, nrow=1, header=FALSE, colClasses="character")[1,])) {
		significant.gene.list <- read.delim(significance.file, header=TRUE, quote="")
	} else {
		significant.gene.list <- read.delim(significance.file, header=FALSE, quote="")
		colnames(significant.gene.list) <- c("gene", "q")
	}
} else {
	if (significance.file!="") stop("significance list not found!")
	significant.gene.list <- read.delim(paste(mutsig.output, "/", analysis.set,".sig_genes.txt", sep=""), header=TRUE, quote="")
}
if (verbose) cat(".\n")

# the significance gene list has 18144 rows (genes)
# the first is PIK3CA / TP53

######## coverage ########## done
if (verbose) cat("Reading coverage values")
if (file.exists(coverage.file)) {
	if(all(c("name", "N_tot")%in%read.delim(coverage.file, nrow=1, header=FALSE, colClasses="character")[1,])) {
		patients.counts_and_rates <- read.delim(coverage.file, header=TRUE, quote="")
	} else {
		patients.counts_and_rates <- read.delim(coverage.file, header=FALSE, quote="")
		colnames(patients.counts_and_rates) <- c("name", "N_tot")
	}
} else {
	if (coverage.file!="") stop("coverage values not found!")
	patients.counts_and_rates <- read.delim(paste(mutsig.output, "/", analysis.set,".patients.counts_and_rates.txt", sep=""), header=TRUE, quote="")
}
rownames(patients.counts_and_rates) <- patients.counts_and_rates$name
if (verbose) cat(".", nlevels(patients.counts_and_rates$name), "samples.\n")


######### maf ############
if (verbose) cat("Reading maf file")
if (file.exists(maf.file)) {
	final_analysis_set <- read.delim(maf.file, header=TRUE, quote="", comment.char = "#")
} else {
	if (maf.file!="") stop("MAF file not found!")
	final_analysis_set <- read.delim(paste(mutsig.output, "/", analysis.set,".final_analysis_set.maf", sep=""), header=TRUE, quote="")
}
# change patidx
if (verbose) cat(".", nrow(final_analysis_set), "mutations in", nlevels(factor(final_analysis_set$pat_idx)), "samples.")
# tmp <- lapply(list(c("_Mutation", ""), c("Silent", "Synonymous"), c("Splice_Site.*", "Splice_Site"), 
# 				c("Nonstop|De_novo_Start.*|Start_.*|Translation_.*|Read\\-through.*", "Other_non_syn."), c("In_frame.*", "In_frame_Indel"), c("Frame_Shift.*", "Frame_Shift"),
# 				c("3'\\-?UTR|5'\\-?UTR|3'\\-?Flank|5'\\-?Flank|IGR|Intron", "Silent")), 
# 				# categ change
# 		function(x) final_analysis_set$categ <<- sub(x[1], x[2], as.character(final_analysis_set$categ), ignore.case=TRUE))



effchar <- c("flank", "syn", "missense", "nonsense", "splice")
final_analysis_set$effchar <-  factor(final_analysis_set$effect_idx,labels=c("flank", "symm", "missense", "nonsense", "splice"))
#final_analysis_set$chr <- factor(sub("^M$", "MT", final_analysis_set$chr))
if (verbose) cat("\n")

final_analysis_set$pat_idx <- factor(final_analysis_set$pat_idx)
final_analysis_set$patient_name <- factor(final_analysis_set$pat_idx, labels =patients.counts_and_rates$name)



## get the hugo gene symbol
geneid <- final_analysis_set$gene_idx
res <- R.matlab::readMat(paste0(mutsig.full.path,"/input_data.M/gene/name.mat"))
syms_db <- unlist(res$tmp)

final_analysis_set$symbol <- syms_db[geneid]
all_sig_genes <- significant.gene.list




if (!exists("forceSelectedGenes")) {
	selectedGenes <- subset(all_sig_genes, q <= qthresh)$gene
} else {
	selectedGenes <- intersect(all_sig_genes$gene,forceSelectedGenes)
}

# To prevent failure, select top 15 if there are no selected genes
if (length(selectedGenes) < 1) {
	selectedGenes <- all_sig_genes$gene[1:15]
}

selectedGenes <- setdiff(selectedGenes, blacklist)
sig_genes <- subset(all_sig_genes, gene%in%selectedGenes)
if ('n' %in% names(sig_genes)) {
	nonsilentColumnName = "n"
	print(nonsilentColumnName)
} else if ('nnon' %in% names(sig_genes)) {
	nonsilentColumnName = "nnon"
	print(nonsilentColumnName)
} else {
	stop("Nonsilent mutations column missing from significant genes file. Require 'n' or 'nnon' column.")
}


if(sort.genes.by.prevalence) {
	if ('n' %in% names(sig_genes)) {
		sig_genes <- sig_genes[order(sig_genes$n, -log10(sig_genes$q), sig_genes$gene, decreasing=TRUE),]
	} else {
		sig_genes <- sig_genes[order(sig_genes$nnon, -log10(sig_genes$q), sig_genes$gene, decreasing=TRUE),]
	}
}

if (verbose) cat(length(selectedGenes),"genes to display...\n q-values:\n")
if (verbose) print(summary(sig_genes$q))

## Build mutation matrix
final_analysis_subset <- subset(final_analysis_set, symbol%in%sig_genes$gene)
# if ("miRNA"%in%(final_analysis_subset$Variant_Classification)) {
# 	warning("coMut does not handle miRNA mutations in significant genes yet. Please report this!")
# }
# final_analysis_subset <- subset(final_analysis_subset, Variant_Classification_Num<=7)
final_analysis_subset$symbol <- factor(final_analysis_subset$symbol, levels=sig_genes$gene)
print(final_analysis_subset$effect_idx)
print(final_analysis_subset$symbol)
print(final_analysis_subset$patient_name)
final_analysis_subset <- aggregate(effect_idx ~ symbol + patient_name, max, na.rm=TRUE, data=final_analysis_subset)

genematrix <- xtabs(effect_idx ~ symbol + patient_name, data=final_analysis_subset, drop.unused.levels = FALSE)

# calculate the appearance of each mutation by category

bkgd_N <- lapply(split(mutcategs, mutcategs$categ),function(df) data.frame(categ = unique(df$categ), sum_categ = sum(as.numeric(df$N))))
bkgd_categ <- do.call(rbind, bkgd_N)

## build spectrum index legend
mutcategs_mat <- xtabs(rate ~ categ + name, data = mutcategs)
spectrum_legends <- apply(mutcategs_mat, 1, function(v) {
    len <- length(v[v>0])
    # print(len)
    legends <- names(v)[order(v, decreasing = T)[1:min(len, 2)]]
    #lapply(legends, function(l) paste0(l, collapse = "."))
    paste0(legends, collapse = " or ")
})

## Build mutation categories matrix
if (mutation.spectrum.plot && is.null(final_analysis_set$categ) == FALSE) {
	final_analysis_set$count <- 1
	total_mut <- aggregate(count ~ patient_name, data=final_analysis_set, sum)
	mutational_signatures <- aggregate(count ~ patient_name + categ, data=final_analysis_set, sum)
	mutational_signatures <- merge(mutational_signatures, total_mut, by="patient_name")
	mutational_signatures <- merge(mutational_signatures, bkgd_categ, by = "categ")
	mutational_signatures$rate <- mutational_signatures$count.x/(mutational_signatures$count.y*mutational_signatures$sum_categ)
	mutational_signatures <- xtabs(rate ~ categ + patient_name, mutational_signatures)
	# make sure column sum is one
	mutational_signatures <- apply(mutational_signatures, 2, function(x) x/sum(x))
	
	
	if (nrow(mutational_signatures) == 0) {
		stop("No mutational signatures found")
	}
	new_rownames <- spectrum_legends
	if (length(new_rownames) != nrow(mutational_signatures)) {
		stop("Mutational categories do not match the categories provided in the MAF")
	}
	rownames(mutational_signatures) <- new_rownames
	# mutational_signatures <- mutational_signatures[, order(mutational_signatures[rownames(mutational_signatures)[1],], decreasing=T), drop=FALSE]
	present_mut_categs <- rownames(mutational_signatures)
}

if (sort.by.mutation.status) {
	bin.genematrix <- rbind((genematrix[,as.character(patients.counts_and_rates$name), drop=FALSE]>1)+0, rate=rowSums(patients.counts_and_rates[,c("rate_non", "rate_sil"), drop=FALSE]))
	clustlabels <- as.character(patients.counts_and_rates$name)[do.call("order", c(lapply(1:nrow(bin.genematrix), function(x) bin.genematrix[x,]), decreasing=TRUE))]
} else {
	clustlabels <- names(sort(rowSums(patients.counts_and_rates[,c("rate_non", "rate_sil"), drop=FALSE]), decreasing=TRUE, na.last = TRUE))
}

## Calculate allelic fractions
if (allelic.fraction.boxplot) {
	if (length(unlist(lapply(c("t_alt_count$", "t_ref_count$"), grep, colnames(final_analysis_set))))==2) {
		final_analysis_set$t_alt_count <- suppressWarnings(as.numeric(as.character(final_analysis_set[,grep("t_alt_count$", colnames(final_analysis_set))])))
		final_analysis_set$t_ref_count <- suppressWarnings(as.numeric(as.character(final_analysis_set[,grep("t_ref_count$", colnames(final_analysis_set))])))
		allelic.fraction.subset <- subset(final_analysis_set, (t_ref_count+t_alt_count)>=min.reads.to.count.allelic.fraction)
		allelic.fraction.subset$patient_name <- factor(allelic.fraction.subset$patient_name, levels=clustlabels)
	} else {
		allelic.fraction.boxplot <- FALSE
		warning("t_alt_count and t_ref_count not found in the maf file, or ambiguous, disabling the allelic fraction bloxplots!")
	}
}

mutperc <- rev(sig_genes$npat)/nrow(patients.counts_and_rates)*100

# remove silent
genematrix[which(genematrix==1, arr.ind=TRUE)] <- 0


vsize <- nrow(sig_genes)/10+3

if (allelic.fraction.boxplot) {
    extra.vsize <- .8
} else {
    extra.vsize <- 0
}

mutsigrunsuffix <- paste(ifelse(firehose.mode | mutsig.output=="", "", "_"), ifelse(firehose.mode | mutsig.output=="", "", mutsigrun), sep="")

plot.coMut <- function() {
    # X11("", 8.5, vsize+1+extra.vsize)
    
    if (verbose) cat("  initial image plot\n")
    ## Mutation matrix
    layout(matrix(c(6, 4, 7, 2, 1, 3, 0, 9, 0, 8, 5, 0), ncol=3, byrow=T), widths=c(1,3,1), heights=c(1, vsize-2, extra.vsize, 2))
    par(mar=c(0,0,0,0), las=1)
    image(1:ncol(genematrix), 1:nrow(genematrix), t(genematrix[rev(as.character(sig_genes$gene)), clustlabels, drop=FALSE]), col=c("grey94", brewer.pal(5, "Set1")), zlim=c(0,5), xlab="", ylab="", axes=FALSE)
    segments(.5 + 1:(ncol(genematrix)-1), .5, .5 + 1:(ncol(genematrix)-1), .5 + nrow(genematrix), col="grey96", lwd=ifelse(ncol(genematrix)>200, .2, .5))
    segments(.5, .5 + 1:(nrow(genematrix)-1), .5 + ncol(genematrix), .5 + 1:(nrow(genematrix)-1), col="grey96", lwd=ifelse(ncol(genematrix)>200, .2, .5))
    
    ## Mutation counts
    par(mar=c(0,.5, 0, 2))
    
    if (max(rowSums(sig_genes[,c(nonsilentColumnName, "nsil")]))>60) {
        maxCountSigGenes <- as.integer(10^(quantile(log10(rowSums(sig_genes[,c(nonsilentColumnName, "nsil")])), 3/4) + 1.5*IQR(log10(rowSums(sig_genes[,c(nonsilentColumnName, "nsil")])))))
    } else {
        maxCountSigGenes <- max(rowSums(sig_genes[,c(nonsilentColumnName, "nsil")]))
    }
    if (verbose) cat("  counts\n")
    plot(0, xlim=c(maxCountSigGenes, 0), type="n", axes=FALSE, frame.plot=FALSE, xlab="", ylab="", main = analysis.set)
    par(usr=c(par("usr")[1:2], 0, nrow(sig_genes)), lwd=.8)
    mypos <- barplot(t(sig_genes[nrow(sig_genes):1,c(nonsilentColumnName, "nsil")]), col=c("dodgerblue4", "#4DAF4A"), 
                     horiz=TRUE, names.arg=rep("", nrow(sig_genes)), add=TRUE, border="grey94", space=0, axes=FALSE)
    axis(4, at=mypos, labels=paste(round(mutperc), "%", sep=""), line=-1, cex.axis=1.1, tick=FALSE)
    if (maxCountSigGenes < max(rowSums(sig_genes[,c(nonsilentColumnName, "nsil")])))
        text(rep(maxCountSigGenes*.9, 2), which(rev(rowSums(sig_genes[,c(nonsilentColumnName, "nsil")]))>maxCountSigGenes)-.5, 
             cex=.8, col="white", labels=rev(rowSums(sig_genes[,c(nonsilentColumnName, "nsil")]))[which(rev(rowSums(sig_genes[,c(nonsilentColumnName, "nsil")]))>maxCountSigGenes)])
    par(new=TRUE, lwd=1)
    plot(0, xlim=c(maxCountSigGenes, 0), type="n", yaxt="n", frame.plot=FALSE, xlab="", ylab="")
    mtext("# mutations", side=1, cex=.7, line=2)
    
    ## q-values
    par(mar=c(0,6,0,.5))
    if (verbose) cat("  q values\n")
    plot(0, xlim=c(min(.4, min(-log10(sig_genes$q+0.0001))), max(-log10(sig_genes$q+0.0001))), type="n", yaxt="n", frame.plot=FALSE, xlab="", ylab="")
    par(usr=c(par("usr")[1:2], 0, nrow(sig_genes)), lwd=.8)
    mypos <- barplot(rev(-log10(sig_genes$q+0.001)), horiz=TRUE, axes=FALSE, add=TRUE, names.arg=rep("", nrow(sig_genes)), border="grey94", space=0, cex.names=.9,
                     col=c("grey70", "grey55")[factor((rev(sig_genes$q)<=.1)+1, levels=c(1,2))])
    axis(2, at=mypos, labels=rev(sig_genes$gene), lwd=0, , cex.axis=1.1, line=-.2)
    mtext("-log10(q-value)", side=1, cex=.7, line=2.2, adj=1.1)
    abline(v=-log10(.1), col="red", lwd=1)
    abline(v=-log10(.25), col="purple", lty=2, lwd=1)
    
    ## Mutation rates
    par(mar=c(.5,0,2,0), las=1)
    #plot(1, ylim=range(rowSums(patients.counts_and_rates[clustlabels,c("rate_non", "rate_sil")], na.rm=TRUE)*1e6), type="n", axes=F, xlab="", ylab="")
    if (verbose) cat("  mutations rates\n")
    plot(1, ylim=c(0,min(max.mutation.rate, max(rowSums(patients.counts_and_rates[clustlabels,c("rate_non", "rate_sil")], na.rm=TRUE)*1e6))), type="n", axes=F, xlab="", ylab="", main = analysis.set)
    par(usr=c(0, nrow(patients.counts_and_rates), par("usr")[3:4]), lwd=ifelse(ncol(genematrix)>200, .2, .5))
    barplot(t(patients.counts_and_rates[clustlabels,c("rate_non", "rate_sil")])*1e6, col=c("dodgerblue4", "#4DAF4A"), axes=FALSE, add=TRUE, names.arg=rep("", nrow(patients.counts_and_rates)), border="grey94", space=0)
    par(lwd=1)
    axis(2, cex.axis=1, line=.3)
    if (any(rowSums(patients.counts_and_rates[clustlabels,c("rate_non", "rate_sil")]*1e6)>max.mutation.rate))
        text(which(rowSums(patients.counts_and_rates[clustlabels,c("rate_non", "rate_sil")]*1e6)>max.mutation.rate)-.5, max.mutation.rate*.9, 
             cex=.8, col="white", labels=round(rowSums(patients.counts_and_rates[clustlabels,c("rate_non", "rate_sil")]*1e6))[which(rowSums(patients.counts_and_rates[clustlabels,c("rate_non", "rate_sil")]*1e6)>max.mutation.rate)],
             srt=90)
    mtext("# mutations/Mb", side=2, cex=.7, line=2.7, las=0)
    
    ## Mutation categories
    if (verbose) cat("  mutations categories\n")
    plot(0, ylim=c(0,100), type="n", axes=FALSE, xlab="", ylab="")
    if (mutation.spectrum.plot && is.null(final_analysis_set$categ) == FALSE) {
        par(mar=c(7,0,.5,0), las=1)
        par(usr=c(0, ncol(mutational_signatures), 0,100))
        box()
        mutation_spectrum_plot_brewer <- brewer.pal(length(present_mut_categs) +1, "Spectral")[seq_along(present_mut_categs)]
        mypos <- barplot(mutational_signatures[,clustlabels]*100, beside=FALSE, col=mutation_spectrum_plot_brewer, names.arg=rep("", ncol(mutational_signatures)), add=TRUE, border=NA, space=0, axes=FALSE)
        axis(4, cex.axis=1, line=.3)
        mtext("%", side=4, cex=.8, line=2.5, las=1)
        par(las=2)
        axis(1, at=mypos, labels=clustlabels, tick=FALSE, cex.axis=min((ifelse(length(clustlabels)<=130, .2, 0) + .7/log10(length(clustlabels))), 1), line=-.7)
    }
    
    ## Display legends
    if (verbose) cat("  legends\n")
    par(mar=c(.3,0,.5,0), las=1)
    plot(0, type="n", axes=FALSE, xlab="", ylab="")
    legend("center", "center", c("Syn.", "Non syn."), fill=c("#4DAF4A", "dodgerblue4"), cex=.9, bty = "n")
    
    plot(0, type="n", axes=FALSE, xlab="", ylab="")
    par(xpd=NA)
    legend("right", "top", effchar,
           fill=brewer.pal(length(effchar), "Set1"), cex=.9, ncol=2, bty = "n", inset=0)
    par(xpd=FALSE)
    
    par(mar=c(.5, 0, .5, 0.5), las=1)
    plot(0, type="n", axes=FALSE, xlab="", ylab="")
    if (mutation.spectrum.plot && is.null(final_analysis_set$categ) == FALSE) {
        mutation.spectrum.legend <- rownames(mutational_signatures)
        legend("right", ifelse(allelic.fraction.boxplot,"top", "center"), 
               rev(mutation.spectrum.legend),
               fill=rev(mutation_spectrum_plot_brewer), cex=.9, ncol=1, inset=0, bty = "n")
    }
    
    ## Allelic fraction boxplot
    if (allelic.fraction.boxplot) {
        if (verbose) cat("  allelic fraction boxplot\n")
        plot(0, ylim=c(0,100), type="n", axes=FALSE, xlab="", ylab="")
        par(mar=c(0.1, 0, .5, 0), las=1)
        par(usr=c(0, length(clustlabels), 0,100))
        box()
        b <- boxplot(allelic.fraction.subset$t_alt_count/(allelic.fraction.subset$t_alt_count+allelic.fraction.subset$t_ref_count)*100 ~ allelic.fraction.subset$patient_name, 
                     pch=20, cex=.3, range=0, lty=1, col=brewer.pal(3, "Set1")[2], ylim=c(0,1), cex.axis=.6, axes=FALSE, add=TRUE, at=(1:length(clustlabels))-.5, plot=FALSE)
        segments(x0=(1:length(clustlabels))-.5, y0=b$stats[1,], y1=b$stats[5,], col=brewer.pal(4, "Purples")[3], lwd=ifelse(ncol(genematrix)>200, .8, 1.5))
        segments(x0=(1:length(clustlabels))-.5, y0=b$stats[2,], y1=b$stats[4,], col=brewer.pal(4, "Purples")[4], lwd=ifelse(ncol(genematrix)>200, 1.5, 3))
        abline(h=median(allelic.fraction.subset$t_alt_count/(allelic.fraction.subset$t_alt_count+allelic.fraction.subset$t_ref_count)*100), col="red", lwd=1, lty=2)
        axis(4, cex.axis=1, line=.3)
        mtext("Allelic\nfraction", side=4, cex=.7, line=3.6, las=0)
    }
    
}


if (verbose) cat("Plotting...\n")
if (!png.output | firehose.mode) {
    if (as.numeric(sub("(.*)\\..*","\\1",R.Version()$minor))>=14) {
        if (full.verbose) cat("  Using Cairo\n")
        cairo_pdf(paste(outputdir, "/", analysis.set, mutsigrunsuffix, "_coMut.pdf", sep=""), 8.5, vsize+1+extra.vsize) ## Only in R 2.14 and up
        plot.coMut()
        dev.off()
    } else {
        pdf(paste(outputdir, "/", analysis.set, mutsigrunsuffix, "_coMut.pdf", sep=""), 8.5, vsize+1+extra.vsize)
        plot.coMut()
        dev.off()
    }
} 

cat(paste(clustlabels, collapse = '\n'), '\n')

if (png.output | firehose.mode) {
    png(paste(outputdir, "/", analysis.set, mutsigrunsuffix, "_coMut.png", sep=""), 8.5*120, (vsize+1+extra.vsize)*120, type="cairo", pointsize=20)
    plot.coMut()
    dev.off()
}

