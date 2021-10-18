# MetaPoly
## An R framework to analyse metagenomic polymorphism data
### 1. The framework
- Analysis of GFF and VCF formatted data in R
- In these format, columns are usually considered as individuals, here we consider the columns as the SFS for a given metagenomic sample (=population)

### 2. PiCorr: A correlation based approach to identify selective sweeps by comparing polymorphism patterns across samples

Example:

### 1. Load data
```
vcf = read.vcfR('Polymorphism/metabat_res.836_filtered.bcf.gz')
genome = read.dna('Genomes/metabat_res.836.fa', format = "fasta")
gff <- read.delim("Genomes/metabat_res.836/PROKKA_09242021.gff", header=F, comment.char="#", sep='\t', quote = '')
gff = gff[gff$V3 == 'CDS',]
```

### 2. Load the GFF and VCF file using MetaPoly's "GetGenesData" function
```
# Load the VCF and GFF data per gene
data = GetGenesData(gff, vcf)
```

### 3. Run the PiCorr algorithm on the formatted data, plotting the results
```
# We create a binomial {0,1} variabl to correlate with the SNP density
samples = colnames(data[[1]]$gene_data)[2:62]
names(samples) = ifelse(startsWith(samples, 'GL_R'), 1, 0)
rock_sed_res = PiCorr(data, 4, samples)

PlotPiCorr(rock_sed_res$pi_corr_res, rock_sed_res$coefs, 'Rock_sed', binomial_var = TRUE)
```

### 4. Functional enrichment analysis
```
# load the COG categories data
cog_func = read.table('cog-20.def.tab', sep='\t')
pos_genes = rock_sed_res$pi_corr_res$gene_id[(rock_sed_res$pi_corr_res$padj < 0.05) & (rock_sed_res$pi_corr_res$cor > 0)]
neg_genes = rock_sed_res$pi_corr_res$gene_id[(rock_sed_res$pi_corr_res$padj < 0.05) & (rock_sed_res$pi_corr_res$cor < 0)]

# Enrichment of functional categories for positive, and negative correlations separately
pos_enrich_df = CalcEnrichment(gff, pos_genes, cog_func)
pos_enrich_df[pos_enrich_df$padj < 0.05,]
#              Function            p       OR   low_CI   high_CI         padj
# odds ratio8         E 1.985170e-10 5.566774 3.211261  9.752151 9.925852e-09
# odds ratio12        F 2.182167e-04 6.841413 2.241347 23.007305 1.003797e-02
# odds ratio16        H 2.634289e-05 4.150955 2.063813  8.349346 1.264459e-03
# odds ratio22        I 8.729403e-04 2.473693 1.431056  4.196891 3.928231e-02
# odds ratio24        J 3.974258e-05 3.160563 1.785382  5.532929 1.867901e-03
# odds ratio39        O 2.193319e-06 4.589560 2.370398  8.936059 1.074727e-04

neg_enrich_df = CalcEnrichment(gff, neg_genes, cog_func)
neg_enrich_df[neg_enrich_df$padj < 0.05,]
#              Function            p        OR     low_CI   high_CI         padj
# odds ratio8         E 5.826731e-07 0.2189155 0.10218237 0.4290036 2.855098e-05
# odds ratio22        I 1.629519e-04 0.3599734 0.19632710 0.6329195 7.658739e-03
# odds ratio24        J 2.159340e-10 0.1118691 0.03910910 0.2615571 1.079670e-08
# odds ratio39        O 2.410850e-06 0.1635080 0.05613308 0.3935008 1.157208e-04
```
