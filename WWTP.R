library(devtools)
install_github('https://github.com/Mass23/MetaPoly')
library(MetaPoly)

setwd('~/Documents/PhD/Others/MetaPoly/MetaPoly')
vcf = read.vcfR('data/WWTP/Bio17-1_NCBI_filtered.bcf.gz')
genome = read.dna('data/WWTP/Bio17-1_NCBI.fa', format = "fasta")
gff <- read.delim("data/WWTP/Bio17-1_NCBI_prokka.gff", header=F, comment.char="#", sep='\t', quote = '')
gff = gff[gff$V3 == 'CDS',]

metadata = read.csv('data/WWTP/WWTP_samples.csv')
metadata = metadata[grepl('D',metadata$Sample),]
samples_vec = metadata$Sample
metadata$Date = as.Date(metadata$Date, format = "%d-%m-%Y")
metadata$time_diff = abs(difftime(metadata$Date, as.Date('2011-11-23', tryFormats = "%Y-%m-%d"), units='days'))
names(samples_vec) = log(as.numeric(metadata$time_diff)+1)

# Load data, get polymorphism summary
data_mt = GetGenesData(gff, vcf)
mt_poly = PolySummary(data_mt, samples_vec[grepl('D',samples_vec)])

mt_poly_f = mt_poly[mt_poly$DEPTH > 0,]

library(pscl)
library(zoib)
snp_n_mod <- zeroinfl(SNP_N ~ log(DEPTH) + gene_length + as.factor(sample) | log(DEPTH), data = mt_poly_f)
summary(snp_n_mod)

even_mod <- zoib(model= EVENNESS ~ log(DEPTH)+as.factor(sample)|1|1|log(DEPTH)|log(SNP_N) ,data=mt_poly_f)
summary(mod)

coefs = snp_n_mod$coefficients$count
coefs = coefs[startsWith(names(coefs),'as.factor(sample)')]
names(coefs) = vapply(names(coefs), function(x) strsplit(as.character(x), split=')')[[1]][2], FUN.VALUE = character(1))
metadata = metadata[metadata$Sample %in% names(coefs),]
metadata$SNP_N_COEFS = vapply(metadata$Sample, function(x) coefs[names(coefs) == x], FUN.VALUE = numeric(1))

ggplot(metadata, aes(x=Date,y=SNP_N_COEFS,color=Season)) + geom_point()
