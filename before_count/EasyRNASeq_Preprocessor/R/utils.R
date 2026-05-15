#' Utility Functions for EasyRNA-Seq Preprocessor
#' 
#' Includes GTF parsing and count aggregation logic.

library(dplyr)
library(readr)
library(tximport)

#' Parse GTF to create tx2gene table
#' 
#' Extracts transcript_id and gene_id pairs.
#' @param gtf_path Path to GTF file.
#' @return data.frame with columns TXNAME and GENEID
make_tx2gene_from_gtf <- function(gtf_path) {
  # Robust parser using base R and regex
  # We read the file line by line or using fread
  # Filter for lines containing "transcript_id" and "gene_id"
  
  message("Reading GTF file: ", gtf_path)
  
  # Read as safe vector of lines (avoiding tabular issues with comments)
  # GTF lines start with comment # usually at top
  # Read as safe vector of lines (avoiding tabular issues with comments)
  # GTF lines start with comment # usually at top
  
  if (grepl("\\.gz$", gtf_path)) {
      con <- gzfile(gtf_path, "rt")
      lines <- readLines(con)
      close(con)
  } else {
      lines <- readLines(gtf_path)
  }
  lines <- lines[!grepl("^#", lines)]
  
  # Filter for 'transcript' features if possible, but actually we just need lines that HAVE both IDs.
  # Often 'exon' lines or 'transcript' lines have both.
  # Let's simple grep for lines with both keys.
  target_lines <- lines[grepl("transcript_id", lines) & grepl("gene_id", lines)]
  
  if (length(target_lines) == 0) {
    stop("No lines with both transcript_id and gene_id found in GTF.")
  }
  
  # Helper to extract value by key
  extract_attribute <- function(lines, key) {
    # Regex: key "value";  or key "value"
    # match key, space?, quote, capture, quote
    regex <- paste0(key, '\\s+"([^"]+)"')
    matches <- regexpr(regex, lines)
    
    # Extract
    res <- regmatches(lines, matches)
    # Remove key and quotes
    res <- gsub(paste0(key, '\\s+"'), "", res)
    res <- gsub('"$', "", res)
    return(res)
  }
  
  tx_ids <- extract_attribute(target_lines, "transcript_id")
  gene_ids <- extract_attribute(target_lines, "gene_id")
  
  # Create map
  tx2gene <- data.frame(
    TXNAME = tx_ids,
    GENEID = gene_ids,
    stringsAsFactors = FALSE
  )
  
  # Remove duplicates
  tx2gene <- distinct(tx2gene)
  
  return(tx2gene)
}

#' Aggregate Salmon Quant.sf files to Gene Counts
#' 
#' @param quant_files Named vector of paths to quant.sf files. Names should be sample names.
#' @param gtf_path Path to GTF file
#' @param output_csv Path to save counts_matrix.csv
#' @param output_len_csv Path to save gene_lengths.csv
aggregate_salmon_counts <- function(quant_files, gtf_path, output_csv, output_len_csv) {
  
  # 1. Create tx2gene
  tx2gene <- make_tx2gene_from_gtf(gtf_path)
  
  # Check first quant file IDs to align with GTF & Force Intersection
  if (length(quant_files) > 0) {
      q1 <- read.delim(quant_files[1], stringsAsFactors=FALSE)
      
      # 1. Try exact match
      common_strict <- intersect(tx2gene$TXNAME, q1$Name)
      
      # 2. Try adding "rna-" prefix (RefSeq fix)
      tx2gene_rna <- tx2gene
      tx2gene_rna$TXNAME <- paste0("rna-", tx2gene$TXNAME)
      common_rna <- intersect(tx2gene_rna$TXNAME, q1$Name)
      
      # Decide which is better
      if (length(common_rna) > length(common_strict)) {
          message(paste("INFO: Better match with 'rna-' prefix. (Strict:", length(common_strict), " vs RNA-prefix:", length(common_rna), ")"))
          tx2gene <- tx2gene_rna
          common_final <- common_rna
      } else {
          message(paste("INFO: Using strict match. (Common IDs:", length(common_strict), ")"))
          common_final <- common_strict
      }
      
      if (length(common_final) == 0) {
          # Extreme Debugging
          message("ERROR: NO COMMON TRANSCRIPT IDS FOUND!")
          message("Head GTF IDs: ", paste(head(tx2gene$TXNAME), collapse=", "))
          message("Head Quant IDs: ", paste(head(q1$Name), collapse=", "))
          stop("Transcript IDs in GTF and Salmon Output do not match at all. Please check if Reference and GTF are from the same source.")
      }
      
      # 3. Filter tx2gene to ONLY keep transcripts in the quant file
      # This prevents tximport error "different number of columns" if sets don't match
      tx2gene <- tx2gene[tx2gene$TXNAME %in% q1$Name, ]
  }

  # 2. Run tximport
  # ensure file existence
  if (!all(file.exists(quant_files))) {
    stop("Some quant.sf files are missing.")
  }
  
  # ignoreTxVersion tries to strip .vX from IDs to match Gencode/Ensembl logic
  # Since Salmon index was built with the FASTA which might have versions, 
  # and GTF might have them too, mismatch is common if one has it and other doesn't or formatted differently.
  txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
  
  # 3. Counts
  counts_df <- as.data.frame(txi$counts)
  # Add Geneid column
  counts_df <- cbind(Geneid = rownames(counts_df), counts_df)
  
  # 4. Save Counts
  write_csv(counts_df, output_csv)
  
  # 5. Lengths
  # tximport calculates average transcript length per gene (weighted by coverage usually)
  # txi$length contains matrix of lengths. We can take average or just output the matrix.
  # Spec says: "gene_lengths.csv (Geneid, Length)" - implying single length per gene?
  # Usually length varies by sample due to isoform usage.
  # Existing app `module_data_upload_metadata.R` likely expects a single column or matrix.
  # Let's verify spec: "Gene Length: gene_lengths.csv (Geneid, Length) also output separately."
  # It implies a static file. But with Salmon/tximport, length is dynamic.
  # We will output the median or mean length across samples for simple single-column output,
  # OR output the whole matrix if compatible. 
  # "Geneid, Length" singular suggests simple 2-column CSV.
  
  # We'll take the row means of length matrix
  avg_lengths <- rowMeans(txi$length)
  len_df <- data.frame(
    Geneid = names(avg_lengths),
    Length = avg_lengths
  )
  write_csv(len_df, output_len_csv)
  
  return(list(counts = counts_df, lengths = len_df))
}

#' Aggregate featureCounts outputs to Gene Counts
#' 
#' @param quant_files Named vector of paths to _featurecounts.txt files. Names should be sample names.
#' @param output_csv Path to save counts_matrix.csv
#' @param output_len_csv Path to save gene_lengths.csv
aggregate_featurecounts <- function(quant_files, output_csv, output_len_csv) {
  
  if (length(quant_files) == 0 || !all(file.exists(quant_files))) {
    stop("Some featureCounts files are missing.")
  }
  
  # Read the first file to get Geneid and Length
  first_file <- read.delim(quant_files[1], comment.char = "#", stringsAsFactors = FALSE)
  counts_df <- data.frame(Geneid = first_file$Geneid, stringsAsFactors = FALSE)
  
  # Extract Length
  len_df <- data.frame(Geneid = first_file$Geneid, Length = first_file$Length, stringsAsFactors = FALSE)
  write_csv(len_df, output_len_csv)
  
  # Loop through all files and bind counts
  for (s_name in names(quant_files)) {
    cf <- read.delim(quant_files[s_name], comment.char = "#", stringsAsFactors = FALSE)
    # Ensure order matches
    if (!all(cf$Geneid == counts_df$Geneid)) {
        stop("Geneid order mismatch between featureCounts files.")
    }
    # Count is usually the 7th column
    counts_df[[s_name]] <- cf[[7]]
  }
  
  write_csv(counts_df, output_csv)
  
  return(list(counts = counts_df, lengths = len_df))
}
