#' Docker Wrapper Functions for EasyRNA-Seq Preprocessor
#'
#' This file contains functions to interact with Docker Desktop.
#' It handles path mounting and platform compatibility (Apple Silicon).

library(processx)
library(tools)

# --- Configuration ---

# Docker Images
DOCKER_IMG_FASTP <- "quay.io/biocontainers/fastp:0.23.4--h5f740d0_0"
DOCKER_IMG_SALMON <- "combinelab/salmon:1.10.0"
DOCKER_IMG_TRINITY <- "trinityrnaseq/trinityrnaseq:2.15.1"
DOCKER_IMG_STAR <- "quay.io/biocontainers/star:2.7.11b--h5ca1c30_8"
DOCKER_IMG_SUBREAD <- "quay.io/biocontainers/subread:2.1.1--h577a1d6_0"
DOCKER_IMG_CORSET <- "quay.io/biocontainers/corset:1.09--h9f5acd7_3"
DOCKER_IMG_BUSCO <- "quay.io/biocontainers/busco:5.5.0--pyhdfd78af_0"
DOCKER_IMG_MULTIQC <- "quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0"
DOCKER_IMG_FALCO   <- "quay.io/biocontainers/falco:1.2.5--h077b44d_0"

# --- Helper Functions ---

# Check if Docker is running
# Check if Docker is running
# Check if Docker is running
check_docker_status <- function() {
  tryCatch(
    {
      # Check for docker executable in path
      docker_path <- Sys.which("docker")
      message("DEBUG: Docker path found: '", docker_path, "'")

      # If not found, check common locations on macOS
      if (docker_path == "") {
        common_paths <- c(
          "/usr/local/bin/docker",
          "/opt/homebrew/bin/docker",
          "/usr/bin/docker",
          "/Applications/Docker.app/Contents/Resources/bin/docker"
        )

        for (p in common_paths) {
          if (file.exists(p)) {
            docker_path <- p
            message("DEBUG: Found docker at: ", p)
            # Add to PATH for the session so processx can find it simply by "docker"
            old_path <- Sys.getenv("PATH")
            Sys.setenv(PATH = paste(dirname(p), old_path, sep = ":"))
            break
          }
        }
      }

      if (docker_path == "") {
        message("DEBUG: Docker executable not found in PATH even after searching common locations.")
        return(FALSE)
      }

      res <- processx::run("docker", c("info"), error_on_status = FALSE)
      message("DEBUG: Exit Code: ", res$status)
      if (res$status != 0) {
        message("DEBUG: Stderr: ", res$stderr)
      }
      return(res$status == 0)
    },
    error = function(e) {
      message("DEBUG: Error checking docker status: ", e$message)
      return(FALSE)
    }
  )
}

# Normalize path for mounting (Mac specific mainly)
# Ensures path is absolute and clean.
get_abs_path <- function(path) {
  if (is.null(path) || path == "") {
    return(NULL)
  }
  normalizePath(path, mustWork = FALSE, winslash = "/")
}

# Construct -v mount arguments and translate paths
# Returns a list: list(mounts = c("-v", "host:container", ...), mapped_paths = list(name = "/container/path"))
map_files_to_container <- function(file_list) {
  # file_list: list(key = absolute_host_path)
  # Strategy: Mount the dirname of each file to a unique mount point in container

  mounts <- c()
  mapped_paths <- list()

  # Group files by directory to minimize mounts
  dir_map <- list() # host_dir -> container_mnt
  mnt_counter <- 1

  for (key in names(file_list)) {
    host_path <- file_list[[key]]
    if (is.null(host_path)) next

    host_dir <- dirname(host_path)
    base_name <- basename(host_path)

    if (is.null(dir_map[[host_dir]])) {
      # New directory to mount
      container_mnt <- paste0("/mnt/vol_", mnt_counter)
      dir_map[[host_dir]] <- container_mnt
      mounts <- c(mounts, "-v", paste0(host_dir, ":", container_mnt))
      mnt_counter <- mnt_counter + 1
    }

    mapped_paths[[key]] <- file.path(dir_map[[host_dir]], base_name)
  }

  return(list(mounts = mounts, mapped_paths = mapped_paths))
}

# Generic Docker Run Command Builder
build_docker_cmd <- function(image, cmd_args, mounts, environment = NULL, platform = NULL, workdir = NULL) {
  args <- c("run", "--rm")

  # Platform (essential for Trinity on M1/M2)
  if (!is.null(platform)) {
    args <- c(args, "--platform", platform)
  }

  # Workdir
  if (!is.null(workdir)) {
    args <- c(args, "-w", workdir)
  }

  # Mounts
  args <- c(args, mounts)

  # Environment vars
  if (!is.null(environment)) {
    for (key in names(environment)) {
      args <- c(args, "-e", paste0(key, "=", environment[[key]]))
    }
  }

  # Image and Command
  args <- c(args, image)
  args <- c(args, cmd_args)

  return(args)
}

# --- Module Specific Functions ---

# 1. Salmon Index
# input_fasta: Host path to FASTA
# output_idx_dir: Host path to output directory
run_salmon_index <- function(input_fasta, output_idx_dir, threads = 2) {
  # Create output dir if not exists
  if (!dir.exists(output_idx_dir)) dir.create(output_idx_dir, recursive = TRUE)

  # Parse paths
  files <- list(
    fasta = get_abs_path(input_fasta),
    out_dir = get_abs_path(output_idx_dir)
  )

  mapping <- map_files_to_container(files)

  # Salmon command: salmon index -t transcripts.fa -i index_name
  # Note: Salmon expects the index DIRECTORY to be created or specified.
  # We mapped the parent of index dir as out_dir? No, we mapped the actual dir safely.
  # Wait, if we map `output_idx_dir` to `/mnt/vol_2`, then inside we write to `/mnt/vol_2`.

  # CMD: salmon index -t <fasta> -i <index_path_in_container>
  # Since mapped_paths$out_dir IS the folder, we point index output to that folder

  cmd_args <- c(
    "salmon",
    "index",
    "-t", mapping$mapped_paths$fasta,
    "-i", mapping$mapped_paths$out_dir,
    "-p", as.character(threads)
  )

  # Salmon supports ARM64 usually, but if not, use linux/amd64
  # Using default platform (should be auto-detected or native if available)

  run_args <- build_docker_cmd(
    image = DOCKER_IMG_SALMON,
    cmd_args = cmd_args,
    mounts = mapping$mounts,
    platform = "linux/amd64" # Enforce amd64 for stability/compatibility if image is multi-arch or x86 only
  )

  return(run_args)
}

# 2. Fastp (QC)
# r1, r2: Host paths (r2 is NULL for single end)
# output_dir: Host path for reports and clean fq
# sample_name: Prefix for output
run_fastp <- function(r1, r2 = NULL, output_dir, sample_name) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  files <- list(r1 = get_abs_path(r1), out = get_abs_path(output_dir))
  if (!is.null(r2)) files$r2 <- get_abs_path(r2)

  mapping <- map_files_to_container(files)

  # Container paths
  c_r1 <- mapping$mapped_paths$r1
  c_r2 <- if (!is.null(r2)) mapping$mapped_paths$r2 else NULL

  c_out_html <- file.path(mapping$mapped_paths$out, paste0(sample_name, "_fastp.html"))
  c_out_json <- file.path(mapping$mapped_paths$out, paste0(sample_name, "_fastp.json"))
  c_out_r1 <- file.path(mapping$mapped_paths$out, paste0(sample_name, "_clean_1.fastq.gz"))
  c_out_r2 <- if (!is.null(r2)) file.path(mapping$mapped_paths$out, paste0(sample_name, "_clean_2.fastq.gz")) else NULL

  cmd_args <- c(
    "fastp",
    "-i", c_r1,
    "-h", c_out_html,
    "-j", c_out_json,
    "-o", c_out_r1
  )

  if (!is.null(r2)) {
    cmd_args <- c(cmd_args, "-I", c_r2, "-O", c_out_r2)
  }

  # Fastp is usually x86_64, use emulation
  run_args <- build_docker_cmd(
    image = DOCKER_IMG_FASTP,
    cmd_args = cmd_args,
    mounts = mapping$mounts,
    platform = "linux/amd64"
  )

  return(list(cmd = run_args, clean_r1 = paste0(sample_name, "_clean_1.fastq.gz"), clean_r2 = if (!is.null(r2)) paste0(sample_name, "_clean_2.fastq.gz") else NULL))
}

# 3. Salmon Quant
# index_dir: Host path (the directory generated by salmon index)
# r1, r2: Host paths (clean fastq)
# output_dir: Host path specific for this sample
run_salmon_quant <- function(index_dir, r1, r2 = NULL, output_dir, threads = 4) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  files <- list(
    idx = get_abs_path(index_dir),
    r1 = get_abs_path(r1),
    out = get_abs_path(output_dir)
  )
  if (!is.null(r2)) files$r2 <- get_abs_path(r2)

  mapping <- map_files_to_container(files)

  c_idx <- mapping$mapped_paths$idx
  c_r1 <- mapping$mapped_paths$r1
  c_out <- mapping$mapped_paths$out

  cmd_args <- c(
    "salmon",
    "quant",
    "-i", c_idx,
    "-l", "A", # Automatic library type
    "-p", as.character(threads),
    "--validateMappings",
    "-o", c_out,
    "--gcBias" # Good practice
  )

  if (!is.null(r2)) {
    c_r2 <- mapping$mapped_paths$r2
    cmd_args <- c(cmd_args, "-1", c_r1, "-2", c_r2)
  } else {
    cmd_args <- c(cmd_args, "-r", c_r1)
  }

  run_args <- build_docker_cmd(
    image = DOCKER_IMG_SALMON,
    cmd_args = cmd_args,
    mounts = mapping$mounts,
    platform = "linux/amd64"
  )

  return(run_args)
}

# 4. Trinity
# r1_list, r2_list: comma separated list of files or vector
# output_dir: Host path
# max_memory: e.g. "32G"
# cpu: e.g. 8
run_trinity <- function(r1_list, r2_list = NULL, output_dir, max_memory = "30G", cpu = 4) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # r1_list could be multiple files
  files <- list(out = get_abs_path(output_dir))

  # Flatten inputs to list
  for (i in seq_along(r1_list)) {
    files[[paste0("r1_", i)]] <- get_abs_path(r1_list[i])
  }
  if (!is.null(r2_list)) {
    for (i in seq_along(r2_list)) {
      files[[paste0("r2_", i)]] <- get_abs_path(r2_list[i])
    }
  }

  mapping <- map_files_to_container(files)

  # Reconstruct comma separated lists inside container
  c_r1s <- c()
  for (i in seq_along(r1_list)) c_r1s <- c(c_r1s, mapping$mapped_paths[[paste0("r1_", i)]])

  c_r2s <- c()
  if (!is.null(r2_list)) {
    for (i in seq_along(r2_list)) c_r2s <- c(c_r2s, mapping$mapped_paths[[paste0("r2_", i)]])
  }

  cmd_args <- c(
    "Trinity",
    "--seqType", "fq",
    "--max_memory", max_memory,
    "--CPU", as.character(cpu),
    "--output", mapping$mapped_paths$out
  )

  if (!is.null(r2_list)) {
    cmd_args <- c(cmd_args, "--left", paste(c_r1s, collapse = ","), "--right", paste(c_r2s, collapse = ","))
  } else {
    cmd_args <- c(cmd_args, "--single", paste(c_r1s, collapse = ","))
  }

  # Trinity MUST be linux/amd64 on Apple Silicon
  run_args <- build_docker_cmd(
    image = DOCKER_IMG_TRINITY,
    cmd_args = cmd_args,
    mounts = mapping$mounts,
    platform = "linux/amd64"
  )

  return(run_args)
}

# 5. STAR Index
run_star_index <- function(input_fasta, input_gtf = NULL, output_idx_dir, threads = 4) {
  if (!dir.exists(output_idx_dir)) dir.create(output_idx_dir, recursive = TRUE)

  files <- list(
    fasta = get_abs_path(input_fasta),
    out_dir = get_abs_path(output_idx_dir)
  )
  if (!is.null(input_gtf) && input_gtf != "") {
      files$gtf <- get_abs_path(input_gtf)
  }

  mapping <- map_files_to_container(files)

  cmd_args <- c(
    "STAR",
    "--runMode", "genomeGenerate",
    "--genomeDir", mapping$mapped_paths$out_dir,
    "--genomeFastaFiles", mapping$mapped_paths$fasta,
    "--runThreadN", as.character(threads)
  )

  if (!is.null(input_gtf) && input_gtf != "") {
      cmd_args <- c(cmd_args, "--sjdbGTFfile", mapping$mapped_paths$gtf)
  }

  run_args <- build_docker_cmd(
    image = DOCKER_IMG_STAR,
    cmd_args = cmd_args,
    mounts = mapping$mounts,
    platform = "linux/amd64"
  )

  return(run_args)
}

# 6. STAR Align
run_star_align <- function(index_dir, r1, r2 = NULL, gtf, output_dir, sample_name, threads = 4) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  files <- list(
    idx = get_abs_path(index_dir),
    r1 = get_abs_path(r1),
    out = get_abs_path(output_dir),
    gtf = get_abs_path(gtf)
  )
  if (!is.null(r2)) files$r2 <- get_abs_path(r2)

  mapping <- map_files_to_container(files)

  c_idx <- mapping$mapped_paths$idx
  c_r1 <- mapping$mapped_paths$r1
  c_gtf <- mapping$mapped_paths$gtf
  c_out_prefix <- file.path(mapping$mapped_paths$out, paste0(sample_name, "_"))

  cmd_args <- c(
    "STAR",
    "--genomeDir", c_idx,
    "--runThreadN", as.character(threads),
    "--outFileNamePrefix", c_out_prefix,
    "--sjdbGTFfile", c_gtf,
    "--outSAMtype", "BAM", "SortedByCoordinate"
  )

  if (grepl("\\.gz$", r1)) {
    cmd_args <- c(cmd_args, "--readFilesCommand", "zcat")
  }

  if (!is.null(r2)) {
    c_r2 <- mapping$mapped_paths$r2
    cmd_args <- c(cmd_args, "--readFilesIn", c_r1, c_r2)
  } else {
    cmd_args <- c(cmd_args, "--readFilesIn", c_r1)
  }

  run_args <- build_docker_cmd(
    image = DOCKER_IMG_STAR,
    cmd_args = cmd_args,
    mounts = mapping$mounts,
    platform = "linux/amd64"
  )
  
  return(run_args)
}

# 7. featureCounts
run_featurecounts <- function(bam, gtf, output_dir, sample_name, is_paired = FALSE, threads = 4) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  files <- list(
    bam = get_abs_path(bam),
    gtf = get_abs_path(gtf),
    out = get_abs_path(output_dir)
  )

  mapping <- map_files_to_container(files)

  c_out_file <- file.path(mapping$mapped_paths$out, paste0(sample_name, "_featurecounts.txt"))

  cmd_args <- c(
    "featureCounts",
    "-t", "exon",
    "-g", "gene_id",
    "-a", mapping$mapped_paths$gtf,
    "-o", c_out_file,
    "-T", as.character(threads)
  )

  if (is_paired) {
    cmd_args <- c(cmd_args, "-p", "--countReadPairs")
  }
  
  cmd_args <- c(cmd_args, mapping$mapped_paths$bam)

  run_args <- build_docker_cmd(
    image = DOCKER_IMG_SUBREAD,
    cmd_args = cmd_args,
    mounts = mapping$mounts,
    platform = "linux/amd64"
  )
  return(run_args)
}

# 8. Corset
run_corset <- function(eq_classes_files, output_dir, threads = 4) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  files <- list()
  for (i in seq_along(eq_classes_files)) {
    files[[paste0("eq_", i)]] <- get_abs_path(eq_classes_files[i])
  }
  files$out <- get_abs_path(output_dir)

  mapping <- map_files_to_container(files)

  c_eq_files <- c()
  for (i in seq_along(eq_classes_files)) {
    c_eq_files <- c(c_eq_files, mapping$mapped_paths[[paste0("eq_", i)]])
  }

  cmd_args <- c(
    "corset",
    "-p", as.character(threads),
    "-i", "salmon"
  )
  cmd_args <- c(cmd_args, c_eq_files)

  run_args <- build_docker_cmd(
    image = DOCKER_IMG_CORSET,
    cmd_args = cmd_args,
    mounts = mapping$mounts,
    workdir = mapping$mapped_paths$out,
    platform = "linux/amd64"
  )
  return(run_args)
}

# 9. BUSCO
run_busco <- function(fasta, lineage, output_dir, threads = 4) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  files <- list(
    fasta = get_abs_path(fasta),
    out = get_abs_path(output_dir)
  )
  
  # For busco downloads dataset automatically if it's missing, but it needs a download path.
  # We will map a busco_downloads folder inside output_dir so it persists if needed, or just let it download to output_dir
  
  mapping <- map_files_to_container(files)

  cmd_args <- c(
    "busco",
    "-i", mapping$mapped_paths$fasta,
    "-l", lineage,
    "-o", "busco_output",
    "-m", "transcriptome",
    "-c", as.character(threads),
    "--download_path", file.path(mapping$mapped_paths$out, "busco_downloads")
  )

  run_args <- build_docker_cmd(
    image = DOCKER_IMG_BUSCO,
    cmd_args = cmd_args,
    mounts = mapping$mounts,
    workdir = mapping$mapped_paths$out,
    platform = "linux/amd64"
  )
  return(run_args)
}

# 10. Falco (FastQC-compatible QC)
# reads: single file path (call separately for R1 and R2)
# output_dir: host dir to write results into
# sample_name, read_num ("R1"/"R2"), suffix ("pre"/"post") compose the subdir name
run_falco <- function(reads, output_dir, sample_name, read_num = "R1", suffix = "pre", threads = 2) {
  subdir <- file.path(output_dir, paste0(sample_name, "_", read_num, "_", suffix))
  dir.create(subdir, recursive = TRUE, showWarnings = FALSE)

  files <- list(
    read    = get_abs_path(reads),
    out_sub = get_abs_path(subdir)
  )

  mapping <- map_files_to_container(files)

  cmd_args <- c(
    "falco",
    "--outdir", mapping$mapped_paths$out_sub,
    "-t", as.character(threads),
    mapping$mapped_paths$read
  )

  run_args <- build_docker_cmd(
    image    = DOCKER_IMG_FALCO,
    cmd_args = cmd_args,
    mounts   = mapping$mounts,
    platform = "linux/amd64"
  )
  return(run_args)
}

# 11. MultiQC
run_multiqc <- function(target_dir, output_dir) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  files <- list(
    tgt = get_abs_path(target_dir),
    out = get_abs_path(output_dir)
  )

  mapping <- map_files_to_container(files)

  cmd_args <- c(
    "multiqc", 
    mapping$mapped_paths$tgt,
    "-o", "multiqc_report",
    "-f"
  )

  run_args <- build_docker_cmd(
    image = DOCKER_IMG_MULTIQC,
    cmd_args = cmd_args,
    mounts = mapping$mounts,
    workdir = mapping$mapped_paths$out,
    platform = "linux/amd64"
  )
  return(run_args)
}
