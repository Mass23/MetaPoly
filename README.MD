# MetaPoly
## An R framework to analyse metagenomic polymorphism data
- Analysis of GFF and VCF formatted data in R
- In these format, columns are usually considered as individuals, here we consider the columns as the SFS for a given metagenomic sample (=population)

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

### 3. Summarise the polymorphism using MetaPoly's "PolySummary" function
