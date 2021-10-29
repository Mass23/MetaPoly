cog_functions = c('J'='Translation, ribosomal structure and biogenesis',
                  'A'='RNA processing and modification',
                  'K'='Transcription',
                  'L'='Replication, recombination and repair',
                  'B'='Chromatin structure and dynamics',
                  'D'='Cell cycle control, cell division, chromosome partitioning',
                  'Y'='Nuclear structure',
                  'V'='Defense mechanisms',
                  'T'='Signal transduction mechanisms',
                  'M'='Cell wall/membrane/envelope biogenesis',
                  'N'='Cell motility',
                  'Z'='Cytoskeleton',
                  'W'='Extracellular structures',
                  'U'='Intracellular trafficking, secretion, and vesicular transport',
                  'O'='Posttranslational modification, protein turnover, chaperones',
                  'X'='Mobilome: prophages, transposons',
                  'C'='Energy production and conversion',
                  'G'='Carbohydrate transport and metabolism',
                  'E'='Amino acid transport and metabolism',
                  'F'='Nucleotide transport and metabolism',
                  'H'='Coenzyme transport and metabolism',
                  'I'='Lipid transport and metabolism',
                  'P'='Inorganic ion transport and metabolism',                     
                  'Q'='Secondary metabolites biosynthesis, transport and catabolism',
                  'R'='General function prediction only',
                  'S'='Function unknown')

################# MKTEST ################# 
calc_dn_ds_pn_ps <- function(gene_data, gene_seq, gene_strand, gene_start, gene_end, min_depth, samp_vec, alleles_data){
  for (i in 1:nrow(gene_data)){
    snp_data = gene_data[i,]
    
    # get the codon position and snp pos in codon
    codon_n = floor((snp_data[['POS']] - gene_start) / 3)
    codon_start = codon_n*3
    codon_end = (codon_n*3)+2
    pos_snp = snp_data[['POS']] - gene_start
    pos_in_codon = (3*((pos_snp/3) - floor(pos_snp/3)) ) + 1
    ref_codon = gene_seq[codon_start:codon_end]
    if (gene_strand == '-'){ref_codon = rev(ref_codon)}
    
    # get the amino acids for each allele
    alleles = lapply(alleles_data$alleles[(alleles_data$CHROM == snp_data$CONTIG) & (alleles_data$POS == snp_data$POS)], tolower)
    amino_acids = unlist(lapply(alleles[[1]], function(x)translate(replace(ref_codon,pos_in_codon,x))))
    print(amino_acids)
    
  }
}

MKTest <- function(data, vcf, min_depth = 10, samp_vec, fasta){
  print('Launching - MetaPoly MKTest: a McDonald-Kreitman test tool for metagenomic data')
  t0 = Sys.time()
  
  print(' - Loading fasta and GFF data')
  fasta_loaded = read.fasta(fasta)
  alleles_data = as.data.frame(vcf@fix)[,c('CHROM','POS','REF','ALT')]
  alleles_data$POS = as.integer(alleles_data$POS)
  
  print(' - Computing dn, ds, pn ,ps...')
  mk_df = data.frame()
  cols = as.vector(samp_vec)
  count=0
  for (i in 1:10){
    count = count +  1
    cat("\r",count)
    if (nrow(data[[i]]$gene_data) > 0){
      # Subset snps in the gene
      alleles_data_sub = alleles_data[(alleles_data$CHROM == data[[i]]$gene_contig) & (alleles_data$POS > data[[i]]$gene_start) & (alleles_data$POS < data[[i]]$gene_end),]
      alleles_data_sub$snp = vapply(1:nrow(alleles_data_sub), function(i) paste(alleles_data_sub$CHROM[i], alleles_data_sub$POS[i], sep = '_'), FUN.VALUE = character(1))
      alleles_data_sub$alleles = strsplit(paste0(alleles_data_sub$REF, alleles_data_sub$ALT),'')
      alleles_data_sub = alleles_data_sub[alleles_data_sub$snp %in% data[[i]]$gene_data[['rn']],]

      gene_seq = getSequence(fasta_loaded[data[[i]]$gene_contig])[[1]][data[[i]]$gene_start:data[[i]]$gene_end]
      gene_res = calc_dn_ds_pn_ps(data[[i]]$gene_data, gene_seq, data[[i]]$gene_strand, data[[i]]$gene_start, data[[i]]$gene_end, min_depth, samp_vec, alleles_data_sub)
      mk_df = rbind(mk_df , data.frame(gene_id = data[[i]]$gene_id,
                                       dn = gene_res$dn,
                                       ds = gene_res$ds,
                                       pn = gene_res$pn,
                                       ps = gene_res$ps,
                                       mk = (gene_res$dn / gene_res$ds) / (gene_res$pn / gene_res$ps),
                                       gene_length = data[[i]]$gene_length))}}
  cat(' genes done\n')
  print(' - Analysis done!')
  return(mk_df)}


