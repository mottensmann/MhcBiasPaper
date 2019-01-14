###############################################################
### Functions used to conduct analyses shown in R_code.Rmd  ###
### ------------------------------------------------------  ###
### A header on top of each functions definitions describes ###
### the main operations and the used parameters in roxygen2 ###
### notation. Type source("R/functions.R") to make each     ###
### function executable                                     ### 
###############################################################


#####
# B #
#####

#' @description  convert barcode counts to data frame
#' @param file barcode counts file returned by mhc_cluster2.py
#' 
barcode.counts2df <- function(file) {
  header <- strsplit(readLines(file, n = 1), "\t")[[1]]
  lines <- R.utils::countLines(file)
  steps <- seq(2,lines,2)
  out <- lapply(steps, function(x) {
    mat <- matrix(data = NA, nrow = 1, ncol = 3)
    mat[1, 1] <- as.character(read.table(file, skip = x - 1, nrows = 1, header = F)[[1]])  
    mat[1, 2] <- as.character(read.table(file, skip = x, nrows = 1, header = F)[[1]]) 
    mat[1, 3] <- as.character(read.table(file, skip = x, nrows = 1, header = F)[[2]]) 
    return(as.data.frame(mat))
  })
  df <- do.call("rbind", out)
  names(df) <- header
  return(df)
}

#' @description Computes critical abundance skew (Edgar 2016)
#' @param d nucleotide differences
#' @param a alpha parameter
#' @return critical skew(M,C)
beta <- function(d, a = 2) {
  1/(2^(a*d + 1))
}

#####
# C #
#####

#' @description complement DNA sequence
#' @param x DNA sequence as character vector
comp_seq <- function(x) {
  bases = c("A","C","G","T")
  xx <- unlist(strsplit(toupper(x),NULL))
  paste(unlist(lapply(xx,function(bbb) {
    if (bbb == "A") compString <- "T"
    if (bbb == "C") compString <- "G"
    if (bbb == "G") compString <- "C"
    if (bbb == "T") compString <- "A"
    if (!bbb %in% bases) compString <- "N"
    return(compString)
  })),collapse = "")
}


#####
# D #
#####

#' @description split fastq files with pooled data into amplicons
#' @param reads path to fastq file
#' @param out path to outputfolder
#'
demultiplex_fastq <- function(reads = "miseq_reads/DQB-Pool/merged_reads/merged_filtered.fastq", out = "miseq_reads/DQB-Pool/merged_reads/demultiplexed", outname = "seq.fastq") {
  suppressPackageStartupMessages(library(ShortRead))
  reads_path <- reads
  
  ## check file size > 2GB
  seq_num <- as.numeric(countLines(reads_path)/4)
  steps <- floor(seq(0, seq_num, length.out = 10))
  
  for (i in 1:(length(steps) - 1)) {
    cat("\n", "Step", i, "of", length(steps) - 1)
    reads <- readFastq(reads_path)[(steps[i] + 1):steps[i + 1]]
    
    # extract header from fastq files
    header <- as.character(id(reads)) 
    
    # get barcodes
    num_cols <- stringr::str_count(header[1], pattern = ":")
    labels <- lapply(header, function(x) strsplit(x, ":")[[1]][num_cols + 1]) 
    labels <- unlist(stringr::str_replace(labels,"barcodelabel=","")) 
    
    # get uniques
    barcodes <- unique(labels)
    
    if (i == 1) {
      ## erase files if they exist
      silent <- lapply(barcodes, function(x) {
        if (file.exists(file.path(out, x, outname))) file.remove(file.path(out, x, outname))
        
      })
    }
    
    # write seqs to file 
    write_fa <- function(barcode) {
      ind <- which(labels %in% barcode)
      sub <- reads[ind]
      ifelse(!dir.exists(file.path(out, barcode)), 
             dir.create(file.path(out, barcode), showWarnings = F),
             FALSE)
      
      writeFastq(sub, file.path(out, barcode, outname), "a", compress = F)  
      # fastq2disk(data = sub, file = file.path(out, barcode, outname))
    }
    x <- lapply(barcodes, write_fa)
    
  }
}


#' @description DensityAmplicon (Sommer et al. 2013, S4)
#' @param reads observed number of reads within amplicon
#' @param proba efficiency values for alleles with non-null read numbers
#' @param logOutput: if TRUE return the logdensity, if FALSE return the density.
DensityAmplicon <- function(reads, proba, logOutput) {
  allele.indexes <- which(reads != 0L)
  return(dmultinom(x = reads[allele.indexes], prob = proba[allele.indexes],
                   log = logOutput))
}

#' Compute the loglikelihood of the dataset (Sommer et al. 2013, S4)
#' @param proba: the efficiency values for all alleles.
#' @param data: the name of the dataset.
#' @return the loglikelihood of the dataset
LoglikData <- function(proba, data) {
  densities <- apply(data, 1,
                     function(x) DensityAmplicon(x, proba, logOutput = T))
  return(sum(densities))
}


#####
# E #
#####

#' @description Extract sequences 
#' @param file BLAST output (-outfmt 6)
#' @param flanking num of additional bases 
#' @param file.out output file
#' @param rows select rows of the BLAST output
#' 
extract_hits <- function(file, file.out = NULL,  flanking = NULL, rows = NULL) {
  # read data
  data <- read.table(file, header = F)[,c(2,9,10,1,3,4)]
  if (!is.null(rows)) data <- data[rows,]
  data <- data[!duplicated(data),]
  # extract start and end 
  out <- lapply(1:nrow(data), function(x) {
    temp <- cbind(data[x,], "plus")
    if (temp[3] < temp[2]) {
      temp <- temp[,c(1,3,2,4:6)]
      temp[7] <- "minus"
    }
    if (!is.null(flanking)) {
      temp[,2:3] <- temp[,2:3] + c(-flanking, flanking)
    }
    names(temp) <- paste0("V", 1:7)
    return(temp)
    })
  out <- do.call("rbind", out)
  
  if (is.null(file.out)) file.out <- stringr::str_replace(file, ".txt", ".bed")
  
  # export to disk
  write.table(x = out[,1:3], 
              file = file.out,
              row.names = F,
              col.names = F,
              quote = FALSE,
              sep = "\t",
              fileEncoding = "UTF-8")
  return(out)
}


#####
# F #
#####

#' Filters merged reads prior to (i) clustering Zotus and (ii) assigning raw reads to them
#' 
#' @description 
#' This script does the following tasks:
#' 1.) Ensure that reads contain expected primer pair
#' 2.) Discard reads with unexpected length
#' 3.) Adds barcodes to headers of the fastq file
#' 4.) Removes reads that contain unexpected barcodes
#' 5.) Export subsetted to new FASTQ files.
#' 
#' @param reads
#' path to a file containing merged reads in fastq format with file ending '.fastq'
#' 
#' @param barcodes
#' path to a file containing barcodes for \code{reads} with file ending '.fastq'
#' 
#' @param mapping_file
#' mapping file that contains barcode sequences for all samples
#' 
#' @param forward_primer
#' character giving the forward primer sequence in 5'-3' orientation
#' 
#' @param reverse_primer
#' character giving the reverse primer sequence in 3'-5' orientation 
#' 
#' @param max.mismatch
#' Interger giving the maximum number of allowed mismatches between sequences and primer. By default 1
#' 
#' @param with.indels
#' Boolean. By default no indels in primer sequences are allowed.
#' 
#' @param suffix 
#' string added to file names for generating the output. Note, existing files may be overwritten.
#' 
#' @param illumina_adapter
#' illumina adapter sequence written at the end of the identifier line. Same for all samples. This string will be replaced by 'barcodelabel=' to allow demultiplexing by usearch
#' 
#' @param pcr_size
#' expected fragment length. If specified all reads of differing size will be filtered out.
#' 
#' @param min_size
#' minimum fragment length. If specifed all reads that are shorter will be discarded
#' 
filter_merged_reads <- function(reads = NULL, barcodes = NULL, mapping_file = NULL, forward_primer = NULL, reverse_primer = NULL, max.mismatch = 1, with.indels = F, suffix = "_filtered", illumina_adapter = "GCCAAT", pcr_size = NULL, min_size = NULL) {
  
  suppressPackageStartupMessages(library(ShortRead))
  source("R/fasta_fastq_functions.R")
  
  # checks  
  if (!file.exists(reads)) stop("Specify a valid path to the reads file")
  if (!file.exists(barcodes)) stop("Specify a valid path to the barcodes file")
  if (!file.exists(mapping_file)) stop("Specify a valid path to the mapping file")
  if (!is.character(forward_primer)) stop("Give a forward primer sequence")
  if (!is.character(reverse_primer)) stop("Give a reverse primer sequece")
  
  ## keep path information for later use
  reads_path <- reads
  barcodes_path <- barcodes
  
  ## get reads and associated barcodes
  cat("Read files ... ")
  reads <- readFastq(reads)
  initial_size <- length(reads)
  barcodes <- readFastq(barcodes)
  cat("done\n")
  
  # extract sequences from fastq files
  cat("Extract sequences ... ")
  reads_seq = sread(reads)
  barcodes_seq <- sread(barcodes)
  cat("done\n")
  
  # search for forward primer
  cat("Extract forward primer region ... ")
  primer_f = DNAStringSet(substr(reads_seq, 1, nchar(forward_primer)))  
  cat("done\nDetect matches with the primer sequence ... ")
  hits_f = vcountPattern(forward_primer, 
                         primer_f, 
                         max.mismatch = max.mismatch,
                         with.indels = with.indels)
  cat("done\n")
  perc_forward <- round(sum(hits_f/length(hits_f))*100,2)
  size_forward <- sum(hits_f)
  cat(perc_forward, "% of all reads contain the expected forward primer\n")
  
  # search for reverse primer in the reverse complement
  cat("Extract reverse primer region ... ")
  primer_r = DNAStringSet(substr(reverseComplement(reads_seq), 1, nchar(reverse_primer))) 
  cat("done\n")
  cat("Detect matches with the primer sequence ... ")
  hits_r = vcountPattern(reverse_primer,
                         primer_r, 
                         max.mismatch = max.mismatch,
                         with.indels = with.indels)
  cat("done\n")
  perc_reverse <- round(sum(hits_r/length(hits_r))*100,2) 
  size_reverse <- sum(hits_r)
  cat(perc_reverse, "% of all raw reads contain the expected reverser primer\n")
  
  # merge hits and turn into a logical
  hits <- hits_f + hits_r 
  hits <- ifelse(hits == 2, 1, 0) 
  hits <- as.logical(hits)
  
  # keep reads with bother primers present
  reads <- reads[hits]
  primers_truncated <- length(reads)
  size_primers <- length(reads)
  perc_primers <- round(size_primers/initial_size*100,2)
  cat(perc_primers, "% of all raw reads contain both expected primers\n")
  
  # size selection. Retain only reads of the expected length
  if (!is.null(pcr_size)) {
    cat("Size selection with expected n =", pcr_size," ... ")
    width <- reads@sread@ranges@width
    width.matched <- which(width %in% seq(from = floor(pcr_size*0.99), to = floor(pcr_size*1.01)))
    # apply selection
    reads <- reads[width.matched]
    cat("done\n")
    cat("Retained", round(100*(length(width.matched)/initial_size),2),"% of all reads")
  }
  
  if (!is.null(min_size)) {
    cat("Size selection with expected length >= ", min_size," ... ")
    width <- reads@sread@ranges@width
    width.matched <- which(width >= min_size)
    # apply selection
    reads <- reads[width.matched]
    #barcodes <- barcodes[width.matched]
    cat("done\n")
    cat("Retained", round(100*(length(width.matched)/initial_size),2),"% of all reads")
  }
  
  
  # get sequence header
  reads_id <- id(reads)
  barcodes_id <- id(barcodes)
  # subset barcodes based on id of reads
  barcodes <- barcodes[barcodes_id %in% reads_id]
  
  ## define output names
  make_name <- function(x, suffix = suffix) {
    path <- strsplit(x, "/")[[1]]
    name <- strsplit(path[length(path)], ".fastq")[[1]]
    prefix <- paste0(path[-length(path)], collapse = "/")
    return(paste0(prefix, "/", name, suffix, ".fastq"))
  }
  
  out_names <- unlist(lapply(X = c(reads_path, barcodes_path), make_name, suffix))
  
  cat("\nAdd barcodes to sequence headers ... ")
  
  # preapre labels
  x <- as.character(reads@id)
  # remove whitespaces
  x <- sub(pattern = " ", replacement = "", x = x)
  # subsitute illumina label by 'barcodelabel='
  x <- sub(pattern = illumina_adapter, replacement = "barcodelabel=", x = x)
  x <- paste0(x, as.character(barcodes@sread))
  
  ## polish header
  reads <- ShortReadQ(sread(reads), quality(reads), BStringSet(x))
  barcodes <- ShortReadQ(sread(barcodes), quality(barcodes), BStringSet(x))
  cat("done\n")
  
  #  extract sequences
  cat("Filter for expected barcodes ... ")
  barcodes_seq <- sread(barcodes)
  
  # read mapping file
  mapping <- read.table(file = mapping_file, header = T)[,2]
  
  ## match barcodes to samples
  indices <- which(barcodes_seq %in% mapping)
  
  reads <- reads[indices]
  barcodes_selected <- length(reads)
  
  ## take a look at unmatched barcodes
  unexp_barcode <- barcodes[-indices]
  
  barcodes <- barcodes[indices]
  cat("done\n")
  
  
  cat("Write to disk ... ")
  fastq2disk(data = reads, file = out_names[1])
  fastq2disk(data = barcodes, file = out_names[2])
  cat("done\n")
  cat(round(100 * barcodes_selected/initial_size,2), " % of all reads are retained\n")
  
  ## write log file
  log_out <- strsplit(reads_path, "/")[[1]]
  log_out <- paste0(log_out[-length(log_out)], collapse = "/")
  log_out <- paste0(log_out, "/log", suffix, ".txt")
  if (file.exists(log_out)) {
    sink(log_out, append = F)
    sink()
  }
  
  sink(log_out)
  cat("Input file:\t", reads_path, "\n\t", initial_size,"raw reads\n")
  cat("Filter for primers:","\n\t",
      size_forward, paste0("(", perc_forward,"%) matches with forward primer sequence\n\t"),
      size_reverse, paste0("(", perc_reverse,"%) matches with reverse primer sequence\n\t"),
      size_primers, paste0("(", perc_primers,"%) matches with both primer sequences\n"))
  if (!is.null(pcr_size)) {
    cat("Filter for expected size of", pcr_size,"bp:\n\t",
        length(width.matched), paste0("(", round(length(width.matched)/initial_size*100,2),"%) sequences are retained\n"))
  }
  if (!is.null(min_size)) {
    cat("Filter for expected size >", min_size,"bp:\n\t",
        length(width.matched), paste0("(", round(length(width.matched)/initial_size*100,2),"%) sequences are retained\n"))
    
  }
  cat("Filter for expected barcodes:\n\t",
      barcodes_selected, paste0("(", round(100 * barcodes_selected/initial_size,2),"%) of all raw reads are retained"))
  
  sink()
  
}

#' Filters unpaired FASTQ files prior to clustering sequencing into OTUs
#' 
#' @param reads
#' path to a file containing raw reads in fastq format with file ending .fastq
#' 
#' @param barcodes
#' path to a file containing barcodes for \code{reads} with file ending .fastq
#' 
#' @param mapping_file
#' mapping file that contains barcode sequences for all samples
#' 
#' @param forward_primer
#' character giving the forward primer sequence in 5'-3' orientation
#' 
#' @param reverse_primer
#' character giving the reverse primer sequence in 3'-5' orientation 
#' 
#' @param max.mismatch
#' Interger giving the maximum number of allowed mismatches between sequences and primer. By default 1.
#' 
#' @param with.indels
#' Boolean. By default no indels in primer sequences are allowed.
#' 
#' @param suffix 
#' string added to file names for generating the output. Note, existing files may be overwritten.
#' 
#' @param illumina_adapter
#' illumina adapter sequence written at the end of the identifier line. Same for all samples.
#' 
#' @param pcr_size
#' expected fragment length. If specified all reads of differing size will be filtered out.
#' 
#' @param min_size
#' minimum fragment length
#' 
#' @export
#' 
filter_single_reads <- function(reads = NULL, barcodes = NULL, mapping_file = NULL, forward_primer = NULL, reverse_primer = NULL, max.mismatch = 1, with.indels = F, suffix = "_filtered", illumina_adapter = "GCCAAT", pcr_size = NULL, min_size = NULL, splits = 100) {
  
  suppressPackageStartupMessages(library(ShortRead))
  #source("R/fasta_fastq_functions.R")
  
  reads_path <- reads
  barcodes_path <- barcodes
  
  ## check file size > 2GB
  cat("Load data ...")
  seq_num <- as.numeric(countLines(reads_path)/4)
  seq_num2 <- as.numeric(countLines(barcodes_path)/4)
  cat(seq_num, "reads\n")
  if (seq_num2 > seq_num) {
    cat("Subset barcodes based on reads ...")
    reads <- readFastq(reads_path)
    reads <- id(reads)
    barcodes <- readFastq(barcodes_path)
    barcodes <- barcodes[id(barcodes) %in% reads]
    rm("reads")
    barcodes_path <- strsplit(barcodes_path, "/")[[1]]
    barcodes_path <- 
      file.path(paste0(barcodes_path[1:length(barcodes_path) - 1], collapse = "/"), "barcodes_filt.fastq")
    if (file.exists(barcodes_path)) file.remove(barcodes_path)
    
    writeFastq(object = barcodes,
               file = barcodes_path,
               compress = F)
    rm("barcodes")
    rm("seq_num2")
  }
  
  steps <- floor(seq(0, seq_num, length.out = splits + 1))
  
  # read mapping file
  mapping <- read.table(file = mapping_file, header = T)[,2]
  
  ## define output names
  make_name <- function(x, suffix = suffix) {
    path <- strsplit(x, "/")[[1]]
    name <- strsplit(path[length(path)], ".fastq")[[1]]
    prefix <- paste0(path[-length(path)], collapse = "/")
    return(paste0(prefix, "/", name, suffix, ".fastq"))
  }
  
  out_names <- unlist(lapply(X = c(reads_path, barcodes_path), make_name, suffix))
  out <- lapply(out_names, function(x) {
    if (file.exists(x)) file.remove(x)
  })
  
  rm(list = c("seq_num", "mapping_file", "suffix"))
  
  for (i in 1:(length(steps) - 1)) {
    cat("\n", "Step", i, "of", length(steps) - 1)
    
    cat("\nRead files ... ")
    reads <- readFastq(reads_path)[(steps[i] + 1):steps[i + 1]]
    barcodes <- readFastq(barcodes_path)[(steps[i] + 1):steps[i + 1]]
    cat("done\n")
    
    # extract sequences from fastq files
    cat("Extract sequences ... ")
    reads_seq = sread(reads)
    barcodes_seq <- sread(barcodes)
    cat("done\n")
    
    # search for primer
    if (!is.null(forward_primer)) {
      cat("Extract forward primer region ... ")
      primer_f = DNAStringSet(substr(reads_seq, 1, nchar(forward_primer)))  
      cat("done\nDetect matches with the primer sequence ... ")
      hits_f = vcountPattern(forward_primer, 
                             primer_f, 
                             max.mismatch = max.mismatch,
                             with.indels = with.indels)
      hits <- as.logical(hits_f)
      rm(list = c("primer_f", "hits_f"))
      cat("\n kept", round(length(hits[hits == TRUE])/length(hits)*100, 2), "%\n")
    } else {
      cat("Extract reverse primer region ... ")
      # search for reverse primer in the reverse complement
      cat("Extract reverse primer region ... ")
      primer_r = DNAStringSet(substr(reverseComplement(reads_seq), 1, nchar(reverse_primer))) 
      cat("done\n")
      cat("Detect matches with the primer sequence ... ")
      hits_r = vcountPattern(reverse_primer,
                             primer_r, 
                             max.mismatch = max.mismatch,
                             with.indels = with.indels)
      cat("done\n")
      hits <- as.logical(hits_r)
    }
    
    # filter the reads
    reads <- reads[hits]
    
    # subset barcodes file based on ids matching the reads
    barcodes <- barcodes[id(barcodes) %in% id(reads)]
    
    cat("\nAdd barcodes to sequence headers ... ")
    
    # clean workspace
    
    # for (i in 1:(length(steps) - 1)) {
    #   cat("\n", "Step", i, "of", length(steps) - 1)
    # prepare lablels
    x <- as.character(reads@id)
    
    # remove whitespaces
    x <- sub(pattern = " ", replacement = "", x = x)
    
    # subsitute illumina label by 'barcodelabel='
    x <- sub(pattern = illumina_adapter, replacement = "barcodelabel=", x = x)
    x <- paste0(x, as.character(barcodes@sread))
    
    ## polish header
    reads <- 
      ShortReadQ(sread(reads), quality(reads), BStringSet(x))
    
    barcodes <- 
      ShortReadQ(sread(barcodes), quality(barcodes), BStringSet(x))
    
    #  extract sequences
    barcodes_seq <- sread(barcodes)
    
    ## match barcodes to samples
    indices <- which(barcodes_seq %in% mapping)
    
    reads <- reads[indices]
    barcodes <- barcodes[indices]
    
    writeFastq(reads, out_names[1] , "a", compress = F)  
    writeFastq(barcodes, out_names[2] , "a", compress = F)    
    
    rm(list = c("barcodes", "barcodes_seq",
                "reads", "reads_seq",
                "x", "indices", "hits"))
    
  }
  
}






#' @description convert fasta file to matrix
#' @param fasta file contaning otus in fasta format
#' @return matrix with aligned sequences
fasta2mat <- function(fasta = "miseq_reads/DQB-Pool/clustered_reads/dqb_pct_1.0_a_1.0.fixed.otus.fa") { 
  
  library(magrittr)
  library(ShortRead)
  library(stringr)
  
  # read fasta sequences
  data <- readFasta(fasta)
  
  # shorten names
  short_name <- 
    as.character(id(data)) %>%
    lapply(., function(x) str_replace(x, pattern = "otu", replacement = "")) %>%
    unlist()
  
  
  # Extract sequence and split
  seqs <- 
    sread(data) %>%
    clustal(.) %>%
    set_rownames(.,short_name) %>%
    as.matrix()
  
  return(seqs)
}

#' @description write to fastq file
#' @param data input fastq data
#' @param file name of local file
fastq2disk = function(data = NULL, file = NULL) {
  if (file.exists(file)) {
    sink(file, append = F)
    sink()
  }
  ShortRead::writeFastq(object = data,
                        file = file,
                        full = T,
                        compress = F,
                        mode = "a")
}

#' @description write to fasta file
#' @param data input fasta data
#' @param file name of local file
fasta2disk = function(data = NULL, file = NULL) {
  if (file.exists(file)) {
    sink(file, append = F)
    sink()
  }
  ShortRead::writeFasta(object = data,
                        file = file)
}


#' @description subset fastq file 
#' @param read path to read
#' @param n number of reads to keep
#' @param out suffix of output
#' 
fastq_subset <- function(read = "miseq_reads/DQB-Pool/parsed_barcodes/reads1.fastq", n = 1000000, out = "sub") {
  
  loadR <- function(x, n = n, out = out) {
    source("R/fastq2disk.R")
    fastq <- ShortRead::readFastq(x)
    cat("\nWrite to file...\n")
    out_name <- paste0(strsplit(x, ".fastq")[[1]][1],"_",out,".fastq")
    fastq2disk(data = fastq[1:n], file = out_name)
  }
  
  cat("Loading...\n\t", read)
  loadR(x = read, n = n, out = out)
  cat("done")  
  
} 

#####
# G #
#####

#' Genotype estimation using DOC approach (modified after Lighten et al. 2014)
#' 
#' @description 
#' (I) Calculate cumulative sequencing depth among n variants
#' (II) Compute rate of change (ROC, first derivative)
#' (III) Compute degree of change (DOC) to estimate inflection point
#' 
#' @details 
#' Mathematically, inflections points are identical for 1 or 2 alleles, therefore
#' the two are discriminated based on the gain in sequencing depth from 1 to 2 
#' (i.e. stay with 1 allele, if step to 2 is below 'gain', default 0.05)
#' 
#' @param x vector of reads per variant
#' @param n numbers of most abundant variants to consider. 
#' @param names names of variants
#' @param gain relative increase in sequencing depth standardised over total depth
#' @param doc_min minimal required degree of change giving the significane level of inflection points
#' @param depth minimum relative sequencing depth of alleles. Value between 0 and 1
#' @param .drop drop variants that are zero

get_genotypes <- function(x = NULL,
                          n = 8, 
                          names = NULL,
                          plot = F, 
                          gain = 0.05,
                          doc_min = 40,
                          depth_min = 0.7,
                          .drop = T) {
  # 1. Calculate degree of changes per for n variants
  # 1.1 Cumulative sequence depth
  
  # ensure that n <= nmax
  if (n > length(x)) n <- length(x)
  cum_dep <- cumsum(sort(x, decreasing = T)[1:n])
  
  # ensure that n <= nmax per sample
  if (any(diff(cum_dep) == 0) & isTRUE(.drop)) {
  cum_dep <- cum_dep[-which(diff(cum_dep) == 0)]
  n <- length(cum_dep)
  }
  
  # 1.2. Rate of change i.e. first derivative
  # roc <- diff(cum_dep)
  roc <- c(cum_dep[1], diff(cum_dep))
  
  # 1.3 Degree of change i.e. second derivative
  # if there is only a single variant:
  if (length(roc) == 1) {
    doc <- 1
    
  } else {
  doc <- rep(0, (length(roc) - 1))
  for (i in 1:length(doc)) {
    doc[i] <-  roc[i]/roc[i + 1]
  } 
  }
  
  # 1.4 Standardise doc
  doc <- doc/sum(doc)*100
  
  # 2. Get inflection point of curve
  # infl <- which(doc == max(doc))[1] + 1
  infl <- which(doc == max(doc))[1]
  
  # 3. Distinguish between 1 and 2 alleles
  # Assume one allele, when raise in cumsum is less than 5 %
  
  # if (infl == 2) {
  #   if ((cum_dep[2]/max(cum_dep)) - (cum_dep[1]/max(cum_dep)) < gain)  {
  #     infl <- 1 
  #   }
  # }
  loop <- "Start"
  if (infl > 2) {
    while (loop != "Stop") {
      if ((cum_dep[infl]/max(cum_dep)) - (cum_dep[infl - 1]/max(cum_dep)) < gain) {
        infl <- infl - 1 
      } else {
        loop <- "Stop"
      }
      if (infl <= 1) {
        loop <- "Stop"
      }
    }
  }
  # 3.1 Set quality score to high initially
  quality_score <- "High"
  
  # data frame for plotting
  df2 <- data.frame(x = 0:length(cum_dep),
                    y = c(0, cum_dep/max(cum_dep)*100),
                    quality = "High")
  
  # Linear models
  lm_fit <- with(df2, lm(y~x))
  slope <- lm_fit$coefficients[[2]]
  intercept <- lm_fit$coefficients[[1]]
  
  # correlation between depth and order of variants
  cor_val <- with(df2, cor(x,y))
  if (cor_val > 0.90) {
    df2$quality <- "Low" 
    quality_score <- "Low"
  }
  
  # Check Depth of variants
  if (cum_dep[[infl]]/max(cum_dep) < depth_min) {
    df2$quality <- "Low"
    quality_score <- "Low"
  }
  
  # Require doc to be significant i.e. > doc_min
  if (max(doc) < doc_min) {
    df2$quality <- "Low"
    quality_score <- "Low"
  }
  
  
  # If amplicon quality is low, then:
  # a) Cumsum is straight line i.e. r > 0.85
  # b) Gain to n-1 is small AND
  # c) increase in depth marginal i.e < 10% of n-1
  
  # 4. Score alleles
  alleles <- names(x)[order(x, decreasing = T)][1:infl]
  df <- 
    data.frame(alleles = paste0(alleles, collapse = ";"),
               total_depth = max(cum_dep),
               allele_depth = cum_dep[[infl]]/max(cum_dep),
              # doc = ifelse(infl == 1, NA, max(doc)),
               doc = max(doc),
               n_alleles = infl,
               quality = quality_score)
  
  if (plot == TRUE) {
    
    
    plot <- ggplot(df2, aes(x = x,y = y, shape = quality)) +
      geom_point() +
      geom_line() +
      theme_classic() +
      xlab("Number of variants") +
      ylab("Cumulative sequencing depth [%]") +
      scale_y_continuous(expand = c(0,0),
                         breaks = seq(0,100,10)) +
      scale_x_continuous(expand = c(0,0),
                         breaks = 0:16) +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 14)) +
      geom_vline(xintercept = infl,
                 linetype = "dashed", 
                 color = "darkorange",
                 size = 1.2) +
      geom_abline(slope = slope,
                  intercept = intercept,
                  linetype = "dotted",
                  color = "blue") +
      annotate("text", x = n - 1, y = 20,
               label = paste0("R: ", round(cor_val, 2)),
               color = "blue")
    
    return(plot)
  } else {
    ## create output in form of a list object
    return(list(df = df,
                alleles = alleles,
                coord = 
                  data.frame(x = 1:n,
                             y = cum_dep/max(cum_dep)*100,
                             group = infl,
                             quality = df$quality)))
  }
  
  
}

#####
# H #
#####

#' @description  hamming distance among pairs of sequences
#' @param  char character vector of DNA or Protein sequences
hamming.distance <- function(char) {
  # helper function
  pair.diff <- function(v1,v2) {
    diff <- 0
    for (i in 1:length(v1))  diff <- diff + ifelse(v1[i] == v2[i], 0, 1)
    diff
  }
  # split strings and convert to matrix  
  x <- lapply(char, function(x) strsplit(x, "")[[1]]) 
  x <- t(as.data.frame(x))
  mat <- matrix(0, nrow = nrow(x), ncol = nrow(x))
  rownames(mat) <- names(char)
  colnames(mat) <- names(char)
  
  for (k in 1:(nrow(mat) - 1)) {
    for (l in (k + 1):nrow(mat)) {
      mat[k,l] <- pair.diff(x[k,], x[l,])
      mat[l, k] <- mat[k, l]
    }
  }
  # set above diagonal to NA
  
  for (i in 1:nrow(mat)) mat[i, i:ncol(mat)] <- NA
  
  return(mat)
}

## calculate pairwise difference to primer, account for variable alignment length
Hamming.dist <- function(seq, ref, method = c("rel", "abs")) {
  method <- match.arg(method)
  # discard gaps and binding N
  gaps <- which(seq %in% c("-", "N"))
  seqx <- seq[-gaps]
  refx <- ref[-gaps]
  
  # estimate diff
  diff <- 0
  for (i in 1:length(seqx))  diff <- diff + ifelse(seqx[i] == refx[i], 0, 1)
  # correct for sequence length
  if (method == "rel") {
    diff <- 
      ifelse(length(diff) > 0,diff/length(seqx), NA)
  } 
  return(diff)
}


#####
# P #
#####

#' @description pool unique otus generated for individual amplicons into a single file
#' @param parentfolder path
#' @param filen filename 
#' @param out output filename
pool_zotus <- function(parentfolder = NULL, filen = NULL, out = NULL) {
  suppressPackageStartupMessages(library(ShortRead))
  library(magrittr)
  
  ## list directories  
  dirs <- list.dirs(parentfolder, full.names = F, recursive = F)
  
  ## read zotus
  zotus <- lapply(dirs, function(x) {
    readFasta(file.path(parentfolder, x, filen)) %>%
      sread() %>%
      as.character()})
  variants <- unlist(zotus) %>%
    unique()
  
  ## output 
  if (is.null(out)) {
    out <- stringr::str_replace(string = filen,
                                                pattern = ".fixed.otus.fa", 
                                                replacement = ".pooled.otus.fa")
  }
  sink(file = file.path(parentfolder, out))
  for (i in 1:length(variants)) {
    cat(paste0(">", "Zotu", i, "\n"))
    cat(variants[i], "\n")
  }
  sink()
}

#' @description analyse BLAST hits to ensure primer are target specific
#' @param Blats output of blastn in format -outfmt 6
#' @param threshold defines maximum 'in silico' fragment size to consider
#' 
primer_check <- function(Blast = NULL , threshold = NULL) {
  ## read results
  data <- read.table(file = Blast)[,c(1,2,4,9,10)]
  names(data) <- c("primer","contig","length","start","end")
  data[["primer"]] <- as.character(data[["primer"]])
  # define small helper functions
  #' remove '_F' or '_R' of a string
  remove_F_or_R <- function(primer_list){
    for (i in 1:length(primer_list)) {
      temp <- strsplit(primer_list[i],split = "_")[[1]]
      primer_list[i] <- stringr::str_c(temp[1:(length(temp) - 1)],collapse = "_")
    }
    primer_list
  }
  
  #' get type of a primer by extracting last letter of name
  primer_type <- function(primer_list) {
    for (i in 1:length(primer_list)) { 
      temp <- strsplit(primer_list[i],split = "_")[[1]]
      primer_list[i] <- temp[length(temp)]
    }
    primer_list
  } 
  
  data[["id"]] <- as.character(remove_F_or_R(data[["primer"]]))
  data[["type"]] <- as.character(primer_type(data[["primer"]]))
  
  pairs <- unique(data[["id"]])
  out <- matrix(data = NA, nrow = 1, ncol = ncol(data))
  colnames(out) <- names(data)
  for (n in 1:length(pairs)) {
    temp <- subset(data, data[["id"]] == pairs[n]) 
    candidates <- unique(temp[["contig"]][which(duplicated(temp[["contig"]]))])
    for (i in 1:length(candidates)) {
      temp2 <- temp[temp[["contig"]] == candidates[i],] 
      if (length(unique(temp2[["type"]])) == 2) {
        if (is.null(threshold)) {
          out <- rbind(out, temp2)
        }
        if (length(dist(temp2[["start"]])[abs(dist(temp2[["start"]]) < threshold)] != 0)) { 
          out <- rbind(out, temp2)
        }
      }
    }
  }
  out <- out[-1,]
  rownames(out) <- NULL
  return(out)
}

#' @description This function reads and formats clustering results and exports them as RData
#' @param locus targeted gene in order to load mapping file
#' @param fname folder of the data
process_otus <- function(locus = c('dqb', 'drb'), fname = NULL) {
  library(magrittr)
  
  locus <- match.arg(locus)
  
  ## load mapping file
  if (locus == 'dqb') {
    load("MiSeq/RData/dqb_mapping_file.RData")
    mapping_file <- dqb_mapping_file
  } 
  # else {
  #   load("MiSeq/DRB/RData/drb_mapping_file.RData")  
  #   mapping_file <- drb_mapping_file
  # }
  
  mapping_file$pos <- as.character(mapping_file$pos)
  mapping_file$BarcodeSequence <- as.character(mapping_file$BarcodeSequence)
  
  ## load otus
  otu_tab <-
    read.table(paste0(fname,".otu_table.txt"), header = T)
  
  ## extract barcodes
  otu_barcode <-
    as.character(names(otu_tab)[-1])
  
  ## sort barcodes with respect to the OTU table
  mapping_file <-
    mapping_file[match(otu_barcode, mapping_file$BarcodeSequence),]
  
  ## count reads per barcode
  x <- 
    apply(otu_tab[,-1], 2, sum) %>%
    as.data.frame() %>%
    cbind(rownames(.),.) %>%
    set_colnames(., c("BarcodeSequence", "mapped")) 
  
  ## format otu table replace OTU barcodes by rack location
  names(otu_tab)[-1] <- 
    as.character(mapping_file$pos)
  
  # sort rows by allele number
  otu_number <- lapply(otu_tab$OTUId, function(x) strsplit(as.character(x), split = "Zotu")[[1]][2])
  otu_number <- as.numeric(unlist(otu_number))
  
  otu_tab <- 
    otu_tab[order(otu_number),] # match(paste0("Zotu", as.character(1:nrow(otu_tab))), otu_tab[,1]),
  rownames(otu_tab) <- as.character(otu_tab$OTUId)
  otu_tab <- otu_tab[,-1]
  
  ## summarise barcode frequency
  bc_counts <-
    barcode.counts2df(file = paste0(fname,".barcode.counts.txt"))
  bc_counts[["BarcodeSequence"]] <- as.character(bc_counts[["BarcodeSequence"]])
  x[["BarcodeSequence"]] <- as.character(x[["BarcodeSequence"]])
  
  bc_counts <-  dplyr::left_join(x = bc_counts, y = mapping_file, by = "BarcodeSequence") %>%
    dplyr::left_join(x = ., y = x, by = "BarcodeSequence")
     
  bc_counts[,c(2:3,9)] <-  
    apply(bc_counts[,c(2:3,9)], MARGIN =  2, FUN =  function(x) as.numeric(as.character(x)))
  
  otu_freq <- 
    data.frame(allele = rownames(otu_tab),
               freq = apply(otu_tab, 1, sum))
  
  otu_freq <- otu_freq[order(otu_freq$freq, decreasing = T),]
  otu_freq$allele <- factor(otu_freq$allele, levels = as.character(otu_freq$allele))
  
  
  # normalise otu
  otu_tab2 <- apply(otu_tab, 2, function(x) x/sum(x))
  
  ## create list
  out <- list(otu_tab = otu_tab,
              otu_freq = otu_freq,
              bc_counts = bc_counts)
}

#' @description purify otu table
#' @param x otu table
#' @param y 
purify_otus <- function(x, y) {
  sample <- rownames(x)
  for (i in 1:length(sample)) {
    purified <- rep(0, ncol(x))
    purified[which(colnames(x) %in% y[[sample[i]]])] <- 
      x[which(rownames(x) == sample[i]),y[[sample[i]]]]
    x[i,] <- purified
  }
  x
}

#####
# R #
#####

#' @description subsample amplicon call alleles
#' @param m observed otus
#' @param n sample size
#' @param gain minimum increase in depth from 1 to 2 alleles
#' @param doc_min significance level of inflection point
#' @param depth_min minimum relative sequencing depth over called alleles
rarefaction <- function(m = NULL, n = NULL, gain = 0.05, doc_min = 40, depth_min = 0.7) {
  out <- 
    rrarefy(m, n) %>%
    apply(., 1, function(x) {
      out <- get_genotypes(x, gain = gain, doc_min = doc_min, depth_min = depth_min)
      length(out[["alleles"]])
    })
  out
}

#' @description reverse complement DNA
#' @param x DNA sequence as character vector
revcomp_seq <- function(x) {
  xx <- comp_seq(x)
  xx <- unlist(strsplit(toupper(xx),NULL))
  xx <- xx[seq(length(xx),1,-1)]
  xx <- paste(xx, collapse = "")
  return(xx)
}

#' remove files created by mhc_cluster2.py for demultiplexed samples and only keep
#' raw sequences and list of fixed otus
remove_mhc_cluster_files <- function(
  parentfolder = NULL,
  keep = c("\\.fixed.otus.fa", "\\.otu_table.txt"),
  subfolder = T) {
  deleteR <- function(dir, keep_pattern) {
    files_keep <- lapply(keep, function(x) {
      list.files(dir, pattern = x)
    }) %>%
      unlist()
    files_all <- list.files(dir) 
    files_all <- files_all[!files_all %in% list.dirs(dir, full.names = F)[-1]]
    files_remove <- subset(files_all, !(files_all %in% files_keep))
    file.remove(file.path(dir, files_remove))
  }
  
  if (isTRUE(subfolder)) {
    dirs <- list.dirs(parentfolder, full.names = T, recursive = F)  
    out <- lapply(dirs, deleteR, keep_pattern = keep)
  } else {
    out <- deleteR(parentfolder, keep_pattern = keep)
  }
}  

#' @description execute genotyping by DOC pipeline
#' @param fname folder
#' @param locus locus
#' @param gain minimum increase in depth among scored alleles
#' @param doc_min minimum relative change at inflection point
#' @param depth_min minimum total sequencing over scored alleles
run_genotyping <- function(fname = NULL, 
                           locus = "dqb",
                           gain = 0.05,
                           doc_min = 40,
                           depth_min = 0.7) {
  
  ## get data
  mhc_data <- process_otus(locus = locus, fname = fname)
  
  # View(mhc_data$otu_tab)
  
  ## call genotypes
  genotypes_list <- apply(mhc_data$otu_tab, 2,
                          get_genotypes,
                          names = rownames(mhc_data$otu_tab),
                          gain = gain,
                          doc_min = doc_min, 
                          depth_min = depth_min) 
  
  ## count alleles and discard bad amplicons
  allele_num_df <- lapply(genotypes_list, function(x) x[["df"]]) %>%
    do.call("rbind",.) %>%
    subset(., quality == "High")
  allele_num_df$row <- rownames(allele_num_df)
  
  ## subset otu table using amplicon quality classification
  otu_table <- 
    mhc_data$otu_tab[rownames(allele_num_df[allele_num_df$quality == "High",])] %>%
    t()
  
  ## check that good quality amplicon exists
  if (nrow(otu_table) >= 1) {
    
  
  ## call alleles
  genotypes <- apply(otu_table, 1, function(x) {
    out <- get_genotypes(x)
    out[["alleles"]]
  })
  
  ## unlist and summarise by observed abundance
  df <- unlist(genotypes) %>%
    as.factor() %>%
    summary()
  
  ## order alleles 
  allele_order <- lapply(names(df), function(x) strsplit(x, split = "Zotu")[[1]][2]) %>%
    unlist()  %>%
    as.numeric() %>%
    order(., decreasing = T)
  
  df <- data.frame(x = names(df), y = df)
  df$x <- factor(df$x, levels = df$x[allele_order])
  df$z <- ifelse(df$y == 1,  "Putative Artefact", "Putative Allele")
  
  zotus <- df
  temp <- summary(as.factor(zotus$z))
  
  ## check if there is at least one putative allele
  if (!any(names(temp) == "Putative Allele")) {
    message("Bad data: No putative alleles detected.\nSkip dataset ...\n")
  
    output <- list(zotu_summary = data.frame(
      putative_artefact = 0,
      putative_allele = 0,
      Artefact = 0))  
    return(output)
  } else {
  
  ## check if there is no artefact
  if (length(temp) == 1 && names(temp) == "Putative Allele") {
    temp <- c(temp[1], 0)
  } 
  
  zotu_summary <- data.frame(putative_artefact = temp[2],
                             putative_allele = temp[1],
                             Artefact  = nrow(mhc_data$otu_tab) - sum(temp))
  
  output <- list(called_alleles = genotypes,
                 zotus = zotus,
                 zotu_summary = zotu_summary)
  
  return(output)
  }
  } else {
    message("bad data: No high quality amplicon!", "\nSkip dataset ...")
    output <- list(zotu_summary = data.frame(
      putative_artefact = 0,
      putative_allele = 0,
      Artefact = 0))
    return(output) 
  }
}

#####
# S #
#####

#' @description calculate polymorphism per site (Reche & Reinherz 2003)
#' @param vec_x vector of amino acids across sequences per site
shannon_entropy <- function(vec_x) {
  # count amino acid occurrences
  aminos <- summary(as.factor(vec_x[!is.na(vec_x)]))
  
  # estimate proportions
  prop <- aminos/sum(aminos)
  
  # estimate variation
  var <- sum(prop*log2(prop))*-1
  return(var)
} 

#' @description Summarizes data 
#' @param data a data frame
#' @param measurevar character giving column name of data to summarise
#' @param groupvars character giving column names of grouping variables
#' @param na.rm boolean
#' @param conf.interval confidence interval (default 0.95)
#' @param .drop boolean
#'
#' @source
#' Taken from the R cookbook (cookbook-r.com/Manipulating_data/Summarizing_data/)
#'
summary_stats <- function(data = NULL, measurevar = NULL, groupvars = NULL, na.rm = TRUE, conf.interval = 0.95, .drop = TRUE) {
  
  length2 <- function(x, na.rm = FALSE) {
    if (na.rm) {
      sum(!is.na(x))
    } else {
      length(x)
    }
  }
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- plyr::ddply(data, groupvars, .drop = .drop,
                       .fun = function(xx, col) {
                         c(N = length2(xx[[col]], na.rm = na.rm),
                           mean = mean(xx[[col]], na.rm = na.rm),
                           sd = sd(xx[[col]], na.rm = na.rm)
                         )
                       },
                       measurevar
  )
  
  # Rename the "mean" column
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N - 1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
