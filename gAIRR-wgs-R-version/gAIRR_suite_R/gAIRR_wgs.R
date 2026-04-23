#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# gAIRR-wgs pipeline of gAIRR-suite  -  R wrapper (port of gAIRR_wgs.py)
#
# This is a direct 1:1 port of gAIRR_suite/gAIRR_wgs.py.
# It preserves the original command-line interface, the locus mapping, 
# the order of the three underlying bash stages (novel_allele.sh, allele_calling.sh, flanking_sequence.sh), 
# and the exact argument order each stage expects.
#
# Under the hood we still call bwa / samtools / SPAdes and the 50+ Python helper scripts that ship inside gAIRR_suite/scripts/.  
# Those scripts do the real sequence analysis and are intentionally NOT re-implemented in R:
# re-writing them would almost certainly introduce silent numerical / parsing differences relative to the published pipeline.
#
# Usage (identical to the Python version):
#     Rscript gAIRR_wgs.R \
#         -wd  <work_dir> \
#         -id  <sample_ID> \
#         -rd1 <read.R1.fastq.gz> \
#         -rd2 <read.R2.fastq.gz> \
#         -lc  TRV TRJ TRD
#
# Only base R is required - no external R packages.
# ------------------------------------------------------------------------------

# Never save/restore an .RData file when invoked via Rscript
options(stringsAsFactors = FALSE)

# ----- helper: find the directory this script lives in -----
#
# Need this to locate scripts/ and material/, exactly like the Python version's os.path.dirname(__file__).
# When run with Rscript, the path is passed as --file=...; inside an interactive R session we fall back to cwd.
script_dir_of_this_file <- function() {
    ca <- commandArgs(trailingOnly = FALSE)
    tag <- "--file="
    idx <- grep(tag, ca, fixed = TRUE)
    if (length(idx) > 0) {
        f <- substring(ca[idx[1]], nchar(tag) + 1L)
        return(normalizePath(dirname(f), winslash = "/", mustWork = FALSE))
    }
    normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

# ----- helper: port of get_max_thread() -----
get_max_thread <- function() {
    sysname <- Sys.info()[["sysname"]]
    cmd <- if (identical(sysname, "Darwin")) "sysctl -n hw.ncpu" else "nproc"
    out <- tryCatch(
        suppressWarnings(system(cmd, intern = TRUE, ignore.stderr = TRUE)),
        error = function(e) character(0)
    )
    n <- suppressWarnings(as.integer(trimws(out[1])))
    if (length(n) == 0L || is.na(n) || n < 1L) 1L else n
}

# ----- helper: port of is_tool() / check_program_install() -----
is_tool <- function(name) {
    w <- Sys.which(name)[[1]]
    nzchar(w)
}

check_program_install <- function(names) {
    bad <- FALSE
    for (n in names) {
        if (!is_tool(n)) {
            cat(n, "is a prerequisite program, please install it before running biastools\n", sep = " ")
            bad <- TRUE
        }
    }
    if (bad) quit(save = "no", status = 2L)
}

# ----- helper: run a shell command identically to subprocess.call(shell=TRUE) -----
#
# Go through a shell (system(), not system2()) because several cleanup commands in the original Python rely on glob expansion:
#       rm .../asm_check/*.sam
#       rm -rf .../asm_contigs/*/
#       rm .../asm_contigs/*.fasta.*
# system2() does NOT invoke a shell and would try to literally delete a file
# named "*.sam".
run_shell <- function(cmd) {
    rc <- suppressWarnings(system(cmd))
    invisible(rc)
}

# ----- argparse replacement (base R) -----
#
# Python's argparse accepts `nargs='+'` for -lc, which optparse/argparser in R cannot express directly.
# We parse by hand so we match the Python CLI character-for-character.
print_usage <- function() {
    cat(
        "usage: Rscript gAIRR_wgs.R [-h] [-wd WORK_DIR] [-lc LOCUS [LOCUS ...]]\n",
        "                           [-id SAMPLE_ID] -rd1 READ1 -rd2 READ2\n",
        "                           [--flanking] [--keep] [-t THREAD] [-REF REFERENCE]\n",
        "\n",
        "gAIRR-wgs pipeline of gAIRR-suite (R version).\n",
        sep = ""
    )
}

print_help <- function() {
    print_usage()
    cat(
        "\n",
        "optional arguments:\n",
        "  -h, --help                        show this help message and exit\n",
        "  -wd,  --work_dir   WORK_DIR       Path to output directory [target_call/].\n",
        "  -lc,  --locus      LOCUS ...      Target Locus [TRV TRJ TRD IGV IGJ IGD]\n",
        "  -id,  --sample_id  SAMPLE_ID      Sample ID ['sample']\n",
        "  -rd1, --read1      READ1          Path to gAIRR-seq. Pair-end 1   (required)\n",
        "  -rd2, --read2      READ2          Path to gAIRR-seq. Pair-end 2   (required)\n",
        "        --flanking                  Option to do flanking sequence analysis.\n",
        "        --keep                      Option specify to keep to sam files.\n",
        "  -t,   --thread     THREAD         Number of threads to use [max].\n",
        "  -REF, --reference  REFERENCE      Path to material directory [default: material/ next to this script].\n",
        sep = ""
    )
}

parse_args <- function(argv, default_material) {
    opts <- list(
        work_dir   = "target_call/",
        locus      = c("TRV","TRJ","TRD","IGV","IGJ","IGD"),
        sample_id  = NULL,
        read1      = NULL,
        read2      = NULL,
        flanking   = FALSE,
        keep       = FALSE,
        thread     = NA_integer_,
        reference  = default_material
    )

    i <- 1L
    n <- length(argv)
    while (i <= n) {
        a <- argv[i]
        if (a %in% c("-h", "--help")) {
            print_help(); quit(save = "no", status = 0L)

        } else if (a %in% c("-wd", "--work_dir")) {
            if (i + 1L > n) { cat("Error: ", a, " requires a value\n", sep=""); quit(save="no", status=2L) }
            opts$work_dir <- argv[i + 1L]; i <- i + 2L

        } else if (a %in% c("-id", "--sample_id")) {
            if (i + 1L > n) { cat("Error: ", a, " requires a value\n", sep=""); quit(save="no", status=2L) }
            opts$sample_id <- argv[i + 1L]; i <- i + 2L

        } else if (a %in% c("-rd1", "--read1")) {
            if (i + 1L > n) { cat("Error: ", a, " requires a value\n", sep=""); quit(save="no", status=2L) }
            opts$read1 <- argv[i + 1L]; i <- i + 2L

        } else if (a %in% c("-rd2", "--read2")) {
            if (i + 1L > n) { cat("Error: ", a, " requires a value\n", sep=""); quit(save="no", status=2L) }
            opts$read2 <- argv[i + 1L]; i <- i + 2L

        } else if (identical(a, "--flanking")) {
            opts$flanking <- TRUE; i <- i + 1L

        } else if (identical(a, "--keep")) {
            opts$keep <- TRUE; i <- i + 1L

        } else if (a %in% c("-t", "--thread")) {
            if (i + 1L > n) { cat("Error: ", a, " requires a value\n", sep=""); quit(save="no", status=2L) }
            v <- suppressWarnings(as.integer(argv[i + 1L]))
            if (is.na(v)) { cat("Error: ", a, " must be an integer\n", sep=""); quit(save="no", status=2L) }
            opts$thread <- v; i <- i + 2L

        } else if (a %in% c("-REF", "--reference")) {
            if (i + 1L > n) { cat("Error: ", a, " requires a value\n", sep=""); quit(save="no", status=2L) }
            opts$reference <- argv[i + 1L]; i <- i + 2L

        } else if (a %in% c("-lc", "--locus")) {
            # Consume every subsequent token until we hit the next flag or EOA.
            j <- i + 1L
            vals <- character(0)
            while (j <= n && !startsWith(argv[j], "-")) {
                vals <- c(vals, argv[j]); j <- j + 1L
            }
            if (length(vals) == 0L) {
                cat("Error: -lc requires at least one locus name\n"); quit(save="no", status=2L)
            }
            opts$locus <- vals
            i <- j

        } else {
            cat("Unknown argument: ", a, "\n", sep="")
            print_usage()
            quit(save = "no", status = 2L)
        }
    }

    # Mirror argparse's `required=True` for rd1/rd2.
    if (is.null(opts$read1) || is.null(opts$read2)) {
        cat("Error: -rd1/--read1 and -rd2/--read2 are required\n")
        print_usage()
        quit(save = "no", status = 2L)
    }

    opts
}

# ----- main -----
main <- function() {
    script_dir <- script_dir_of_this_file()
    # Python: os.path.dirname(__file__) + '/material/'
    default_material <- paste0(script_dir, "/material/")

    args <- commandArgs(trailingOnly = TRUE)
    opts <- parse_args(args, default_material)

    workspace     <- opts$work_dir
    target_locus  <- opts$locus
    person_name   <- opts$sample_id
    path_read1    <- opts$read1
    path_read2    <- opts$read2
    flag_flanking <- isTRUE(opts$flanking)
    flag_keeping  <- isTRUE(opts$keep)
    thread        <- opts$thread
    if (is.na(thread)) thread <- get_max_thread()
    path_material <- opts$reference
    if (!endsWith(path_material, "/")) path_material <- paste0(path_material, "/")

    # Python: os.path.dirname(__file__) + '/scripts/'
    path_module <- paste0(script_dir, "/scripts/")

    # Same prerequisite checks as the Python version.
    if (flag_flanking) {
        check_program_install(c("bwa", "spades.py"))
    } else {
        check_program_install(c("bwa"))
    }

    # Same locus map.
    dict_locus_map <- c(
        TRV = "TCRV",
        TRJ = "TCRJ",
        TRD = "TCRD_plusHep",
        IGV = "BCRV",
        IGJ = "BCRJ",
        IGD = "BCRD_plusHep"
    )
    for (loc in target_locus) {
        if (!(loc %in% names(dict_locus_map))) {
            cat("Only TRV, TRJ, TRD, IGV, IGJ, IGD are allowed for target locus call.\n")
            print_usage()
            quit(save = "no", status = 2L)
        }
    }
    list_allele_names <- unname(dict_locus_map[target_locus])

    # Create directories (preserves the Python "mkdir -p" semantics, including the fact that Python concatenates with '/', leading to harmless '//').
    run_shell(paste0("mkdir -p ", shQuote(workspace)))
    run_shell(paste0("mkdir -p ", shQuote(paste0(workspace, "/", person_name))))

    # Helper to build the bash invocation used by the three stages.
    # The two main stages (novel_allele.sh, allele_calling.sh) take 7 args;
    # flanking_sequence.sh takes 8 args (inserts a SPAdes path between the read2 path and the thread count).
    # We pass `extra` only for the latter.
    build_stage_cmd <- function(stage_sh, allele_name, allele_path, extra = NULL) {
        parts <- c(
            "bash",
            shQuote(paste0(path_module, stage_sh)),
            shQuote(paste0(workspace, "/", person_name)),
            shQuote(allele_name),
            shQuote(allele_path),
            shQuote(person_name),
            shQuote(path_read1),
            shQuote(path_read2)
        )
        if (!is.null(extra)) parts <- c(parts, extra)
        parts <- c(parts, shQuote(as.character(thread)))
        paste(parts, collapse = " ")
    }

    for (allele_name in list_allele_names) {
        # --- novel allele calling ---
        allele_path <- paste0(path_material, allele_name, "_alleles_parsed.fasta")
        cat("[gAIRRCall: gAIRR-wgs]", person_name, allele_name, "novel allele calling...\n")
        run_shell(build_stage_cmd("novel_allele.sh", allele_name, allele_path))
        if (!flag_keeping) {
            run_shell(paste0(
                "rm ", workspace, "/", person_name, "/",
                person_name, "_", allele_name, "_novel/bwa_read_to_allele.sam"
            ))
        }

        # --- allele calling (using the novel-augmented reference) ---
        allele_path <- paste0(
            workspace, "/", person_name, "/",
            person_name, "_", allele_name, "_novel/",
            allele_name, "_with_novel.fasta"
        )
        cat("[gAIRRCall: gAIRR-wgs]", person_name, allele_name, "allele calling...\n")
        run_shell(build_stage_cmd("allele_calling.sh", allele_name, allele_path))
        if (!flag_keeping) {
            run_shell(paste0(
                "rm ", workspace, "/", person_name, "/",
                person_name, "_", allele_name, "/bwa_read_to_allele_all.sam"
            ))
        }

        # --- optional flanking-sequence calling ---
        if (flag_flanking) {
            cat("[gAIRRCall: gAIRR-wgs]", person_name, allele_name, "flanking sequence calling...\n")
            # flanking_sequence.sh takes 8 args: $1..$6 same as other stages,
            # $7 = path to SPAdes executable, $8 = thread count.
            run_shell(build_stage_cmd("flanking_sequence.sh", allele_name, allele_path,
                                      extra = "spades.py"))
            if (!flag_keeping) {
                # These require shell glob expansion; go through run_shell.
                run_shell(paste0("rm ", workspace, "/", person_name, "/",
                                 person_name, "_", allele_name,
                                 "_flanking/asm_check/*.sam"))
                run_shell(paste0("rm ", workspace, "/", person_name, "/",
                                 person_name, "_", allele_name,
                                 "_flanking/haplotype_sam/bwa_reads_to_flanking.sam"))
                run_shell(paste0("rm -rf ", workspace, "/", person_name, "/",
                                 person_name, "_", allele_name,
                                 "_flanking/asm_contigs/*/"))
                run_shell(paste0("rm ", workspace, "/", person_name, "/",
                                 person_name, "_", allele_name,
                                 "_flanking/asm_contigs/*.fasta.*"))
                run_shell(paste0("rm ", workspace, "/", person_name, "/",
                                 person_name, "_", allele_name,
                                 "_flanking/group_allele_reads/*.fasta"))
            }
        }
    }

    invisible(NULL)
}

# Only run main() when invoked as a script (not when source()'d).
if (sys.nframe() == 0L) {
    main()
}
