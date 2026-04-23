#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# gAIRR Read Extraction Tool  -  R wrapper (port of gAIRR_suite/gAIRR_extract.sh)
#
# Takes a sample list TSV and, per sample, extracts reads that either
#   (a) fall inside a BED describing AIRR loci on a reference, OR
#   (b) are unmapped,
# then writes paired on_AIRR_unmapped_R1/R2.fastq to an output directory.
#
# Supports BAM / CRAM inputs directly, or FASTQ inputs (which are first aligned with bwa). 
# The real work is still done by samtools / bwa - this script only wires them together, exactly like the original bash version.
#
# Usage (identical to the bash version):
#     Rscript gAIRR_extract.R \
#         -l  <samples.tsv>  \
#         -o  <output_dir>   \
#         [-b <bed_file>]    \
#         [-r <reference.fa>] \
#         [-t <threads>]
#
# Only base R is required.
# ------------------------------------------------------------------------------

options(stringsAsFactors = FALSE)

# ---- defaults (same as the bash original) ------------------------------------
DEFAULT_BED     <- "./material/extract_materials/gAIRR_allele.bed"
DEFAULT_REF     <- "./material/references/hg38_mod.fa"
DEFAULT_THREADS <- 20L

# ---- logging infrastructure --------------------------------------------------
# The bash version uses  `exec > >(tee -a "$LOGFILE") 2>&1`  which ties stdout and stderr of the whole script (including its subprocesses) to tee.
# R cannot reassign its own process file-descriptors that cleanly, so we do two things:
#   1. All R-side messages go through log_echo(), which writes to both console and the log file.
#   2. Every external command is wrapped so its stdout+stderr is piped through `tee -a <logfile>`.  
#      This preserves the user-visible output while still capturing everything in the log.
log_file <- NULL  # set in main()

log_echo <- function(...) {
    msg <- paste0(..., collapse = "")
    cat(msg, "\n", sep = "")
    if (!is.null(log_file)) {
        cat(msg, "\n", sep = "", file = log_file, append = TRUE)
    }
    invisible(NULL)
}

# Run an external pipeline, mirroring the bash behaviour of tee-ing into the logfile.
# Returns the shell exit status.
run_with_log <- function(cmd) {
    if (is.null(log_file)) {
        rc <- suppressWarnings(system(cmd))
    } else {
        wrapped <- paste0("( ", cmd, " ) 2>&1 | tee -a ", shQuote(log_file))
        # NB: tee always exits 0; the producer exit status is captured via
        # bash's PIPESTATUS only if we invoke bash explicitly.
        wrapped <- paste0("bash -c ", shQuote(paste0(wrapped, "; exit ${PIPESTATUS[0]}")))
        rc <- suppressWarnings(system(wrapped))
    }
    invisible(rc)
}

# ----- CLI -----
print_usage <- function() {
    cat(
        "Usage: Rscript gAIRR_extract.R [OPTIONS]\n\n",
        "Required arguments:\n",
        "  -l, --list FILE          Sample list TSV file with columns: sample_ID, path1, [path2]\n",
        "  -o, --outdir DIR         Output directory\n\n",
        "Optional arguments:\n",
        "  -b, --bed FILE           BED file defining AIRR loci regions\n",
        "                           (default: ", DEFAULT_BED, ")\n",
        "  -r, --ref FILE           Reference genome FASTA file\n",
        "                           (default: ", DEFAULT_REF, ")\n",
        "  -t, --threads INT        Number of threads (default: ", DEFAULT_THREADS, ")\n",
        "  -h, --help              Show this help message\n\n",
        "Sample list format:\n",
        "  - For BAM/CRAM: sample_ID<TAB>path_to_bam_or_cram\n",
        "  - For FASTQ:    sample_ID<TAB>R1.fastq<TAB>R2.fastq\n\n",
        "  The tool auto-detects input format based on file extension.\n\n",
        "Example:\n",
        "  # Using default bed and reference\n",
        "  Rscript gAIRR_extract.R -l samples.tsv -o output_dir\n\n",
        "  # Using custom bed and reference\n",
        "  Rscript gAIRR_extract.R -l samples.tsv -b custom.bed -r custom.fa -o output_dir -t 30\n",
        sep = ""
    )
}

parse_args <- function(argv) {
    opts <- list(
        sample_list = NULL,
        bed         = DEFAULT_BED,
        ref         = DEFAULT_REF,
        outdir      = NULL,
        threads     = DEFAULT_THREADS
    )
    i <- 1L; n <- length(argv)
    while (i <= n) {
        a <- argv[i]
        if (a %in% c("-h","--help")) {
            print_usage(); quit(save = "no", status = 0L)
        } else if (a %in% c("-l","--list")) {
            opts$sample_list <- argv[i+1L]; i <- i + 2L
        } else if (a %in% c("-b","--bed")) {
            opts$bed <- argv[i+1L]; i <- i + 2L
        } else if (a %in% c("-r","--ref")) {
            opts$ref <- argv[i+1L]; i <- i + 2L
        } else if (a %in% c("-o","--outdir")) {
            opts$outdir <- argv[i+1L]; i <- i + 2L
        } else if (a %in% c("-t","--threads")) {
            v <- suppressWarnings(as.integer(argv[i+1L]))
            if (is.na(v) || v < 1L) {
                cat("Error: --threads must be a positive integer\n"); quit(save="no", status=1L)
            }
            opts$threads <- v; i <- i + 2L
        } else {
            cat("Unknown option: ", a, "\n", sep = "")
            print_usage()
            quit(save = "no", status = 1L)
        }
    }
    opts
}

# ----- utilities -----

# Detect file format the same way the bash version does, by extension.
detect_format <- function(filepath) {
    fp <- tolower(filepath)
    if (grepl("\\.bam$",          fp)) return("bam")
    if (grepl("\\.cram$",         fp)) return("cram")
    if (grepl("\\.(fastq|fq)\\.gz$", fp)) return("fastq.gz")
    if (grepl("\\.(fastq|fq)$",   fp)) return("fastq")
    "unknown"
}

fmt_hms <- function(seconds) {
    seconds <- as.integer(seconds)
    hours   <- seconds %/% 3600L
    minutes <- (seconds %% 3600L) %/% 60L
    secs    <- seconds %% 60L
    sprintf("%02d:%02d:%02d", hours, minutes, secs)
}

# ----- per-file processing -----

process_alignment <- function(sample_id, input_path, output_dir, bed_file, ref_fasta, threads) {
    log_echo("[", sample_id, "] Processing alignment file: ", input_path)

    # --- Extract gAIRR reads ---
    log_echo("[", sample_id, "] Extracting reads on AIRR loci...")
    extract_start <- Sys.time()
    cmd_extract <- sprintf(
        paste(
            "samtools view -@ %d -bh -L %s -T %s %s |",
            "samtools collate -@ %d -O -u - |",
            "samtools fastq -1 %s -2 %s -0 /dev/null -s /dev/null -n -"
        ),
        threads,
        shQuote(bed_file), shQuote(ref_fasta), shQuote(input_path),
        threads,
        shQuote(file.path(output_dir, "on_AIRR_R1.fastq")),
        shQuote(file.path(output_dir, "on_AIRR_R2.fastq"))
    )
    rc1 <- run_with_log(cmd_extract)
    extract_end <- Sys.time()

    # --- Extract unmapped reads ---
    log_echo("[", sample_id, "] Extracting unmapped reads...")
    unmapped_start <- Sys.time()
    cmd_unmapped <- sprintf(
        paste(
            "samtools view -@ %d -bh -f 4 -T %s %s |",
            "samtools collate -@ %d -O -u - |",
            "samtools fastq -1 %s -2 %s -0 /dev/null -s /dev/null -n -"
        ),
        threads,
        shQuote(ref_fasta), shQuote(input_path),
        threads,
        shQuote(file.path(output_dir, "unmapped_R1.fastq")),
        shQuote(file.path(output_dir, "unmapped_R2.fastq"))
    )
    rc2 <- run_with_log(cmd_unmapped)
    unmapped_end <- Sys.time()

    # --- Merge ---
    log_echo("[", sample_id, "] Merging reads...")
    merge_start <- Sys.time()
    cmd_merge_r1 <- sprintf("cat %s %s > %s",
        shQuote(file.path(output_dir, "on_AIRR_R1.fastq")),
        shQuote(file.path(output_dir, "unmapped_R1.fastq")),
        shQuote(file.path(output_dir, "on_AIRR_unmapped_R1.fastq"))
    )
    cmd_merge_r2 <- sprintf("cat %s %s > %s",
        shQuote(file.path(output_dir, "on_AIRR_R2.fastq")),
        shQuote(file.path(output_dir, "unmapped_R2.fastq")),
        shQuote(file.path(output_dir, "on_AIRR_unmapped_R2.fastq"))
    )
    run_with_log(cmd_merge_r1)
    run_with_log(cmd_merge_r2)

    # Cleanup intermediates.
    inter <- file.path(output_dir,
        c("on_AIRR_R1.fastq","on_AIRR_R2.fastq","unmapped_R1.fastq","unmapped_R2.fastq")
    )
    suppressWarnings(file.remove(inter))
    merge_end <- Sys.time()

    log_echo("  - Read extraction time: ",
             as.integer(difftime(extract_end,  extract_start,  units = "secs")), " seconds")
    log_echo("  - Unmapped read extraction time: ",
             as.integer(difftime(unmapped_end, unmapped_start, units = "secs")), " seconds")
    log_echo("  - Merge time: ",
             as.integer(difftime(merge_end,    merge_start,    units = "secs")), " seconds")

    invisible(list(rc1 = rc1, rc2 = rc2))
}

process_fastq <- function(sample_id, r1_path, r2_path, output_dir, bed_file, ref_fasta, threads) {
    log_echo("[", sample_id, "] Processing FASTQ files")
    log_echo("[", sample_id, "] R1: ", r1_path)
    log_echo("[", sample_id, "] R2: ", r2_path)

    log_echo("[", sample_id, "] Mapping reads to reference...")
    aligned_bam <- file.path(output_dir, "aligned.bam")
    cmd_map <- sprintf(
        paste(
            "bwa mem -t %d -T 10 %s %s %s |",
            "samtools view -@ %d -bh - |",
            "samtools sort -@ %d -o %s -"
        ),
        threads,
        shQuote(ref_fasta), shQuote(r1_path), shQuote(r2_path),
        threads, threads,
        shQuote(aligned_bam)
    )
    run_with_log(cmd_map)

    cmd_idx <- sprintf("samtools index -@ %d %s", threads, shQuote(aligned_bam))
    run_with_log(cmd_idx)

    process_alignment(sample_id, aligned_bam, output_dir, bed_file, ref_fasta, threads)

    suppressWarnings(file.remove(c(aligned_bam, paste0(aligned_bam, ".bai"))))
}

# ----- main -----
main <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    opts <- parse_args(args)

    if (is.null(opts$sample_list) || is.null(opts$outdir)) {
        cat("Error: Missing required arguments (-l and -o are required)\n")
        print_usage()
        quit(save = "no", status = 1L)
    }

    # Validate input files (same three the bash version checks).
    for (f in c(opts$sample_list, opts$bed, opts$ref)) {
        if (!file.exists(f)) {
            cat("Error: File not found: ", f, "\n", sep = "")
            quit(save = "no", status = 1L)
        }
    }

    dir.create(opts$outdir, recursive = TRUE, showWarnings = FALSE)

    date_tag <- format(Sys.time(), "%Y%m%d%H%M%S")
    log_file <<- file.path(opts$outdir, paste0("gAIRR_extraction_", date_tag, ".log"))
    # Create empty log file so tee -a always has a target.
    file.create(log_file, showWarnings = FALSE)

    # Read sample list as a TSV, preserving the bash semantics:
    #   IFS=$'\t' read -r sample_id path1 path2
    # Skip header whose first cell is literally "sample_ID".  Skip empty lines.
    lines <- readLines(opts$sample_list, warn = FALSE)
    line_num <- 0L
    for (ln in lines) {
        line_num <- line_num + 1L
        if (!nzchar(ln)) next

        fields <- strsplit(ln, "\t", fixed = TRUE)[[1]]
        sample_id <- if (length(fields) >= 1L) fields[1] else ""
        path1     <- if (length(fields) >= 2L) fields[2] else ""
        path2     <- if (length(fields) >= 3L) fields[3] else ""

        if (identical(sample_id, "sample_ID") || !nzchar(sample_id)) next

        log_echo("")
        log_echo("==================================================")
        log_echo("Processing sample ", line_num, ": ", sample_id)
        log_echo("==================================================")

        sample_outdir <- file.path(opts$outdir, sample_id)
        dir.create(sample_outdir, recursive = TRUE, showWarnings = FALSE)

        start_time <- Sys.time()

        format_ <- detect_format(path1)
        if (format_ %in% c("bam","cram")) {
            if (!file.exists(path1)) {
                log_echo("Error: File not found: ", path1); next
            }
            process_alignment(sample_id, path1, sample_outdir,
                              opts$bed, opts$ref, opts$threads)
        } else if (format_ %in% c("fastq","fastq.gz")) {
            if (!nzchar(path2)) {
                log_echo("Error: R2 path required for FASTQ input"); next
            }
            if (!file.exists(path1) || !file.exists(path2)) {
                log_echo("Error: FASTQ file(s) not found"); next
            }
            process_fastq(sample_id, path1, path2, sample_outdir,
                          opts$bed, opts$ref, opts$threads)
        } else {
            log_echo("Error: Unknown file format for ", path1); next
        }

        end_time  <- Sys.time()
        total_sec <- as.integer(difftime(end_time, start_time, units = "secs"))
        log_echo("------------------------------------------------")
        log_echo(sprintf("[%s] Total processing time: %s (hh:mm:ss)",
                         sample_id, fmt_hms(total_sec)))
        log_echo("------------------------------------------------")
    }

    log_echo("==================================================")
    log_echo("All samples processed successfully!")
    log_echo("Completed at: ", format(Sys.time(), "%a %b %d %H:%M:%S %Z %Y"))
    log_echo("==================================================")

    invisible(NULL)
}

if (sys.nframe() == 0L) {
    main()
}
