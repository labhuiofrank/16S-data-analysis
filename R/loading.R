#' Check file extension and throw a warning/error
#' if it's not correct
#'
#' @param file Path to file to check
#' @param expected Expected file extension
#' @param abort whether to abort if the extension
#' doesn't match
#'
check_ext <- function(file, expected, abort=FALSE) {
    pattern <- sprintf("%s$", expected)
    if(!grepl(pattern, file)) {
        msg <- sprintf(
            "Expected extension %s for file %s",
            expected, file
        )
        if(abort) { stop(msg) } else { warning(msg) }
    }
}

#' Find a file with a given pattern in folder
#' - Throws an error (or warning if abort is FALSE)
#' if it's not there
#' - Throws a warning if multiple matches are found
#'
#' @param folder Folder to check
#' @param pattern Regex pattern
#' @param abort Whether to abort if no matches are found
#' @return file path
#'
find_file <- function(folder, pattern, abort=FALSE, warning=TRUE) {
    matches <- list.files(
        path=folder, pattern=pattern,
        recursive=TRUE, full.names=TRUE
    )
    if(length(matches) == 0) {
        msg <- sprintf(
            "Could not find any file with pattern %s",
            pattern
        )
        if(abort) { stop(msg) } else if(warning) { warning(msg) }
    }
    if(length(matches) > 1) {
        warning(sprintf(
            "Found multiple files matching pattern : %s. Picking the first one (%s)",
            pattern, matches[1]
        ))
    }
    return(matches[1])
}

#' Load OTU abundance file from mothur (shared format)
#'
#' @param file Path to shared file
#' @return Abundance dataframe (sample x OTU)
#' @examples
#' load_shared("abundance.100.shared")
load_shared <- function(file) {
    check_ext(file, "shared", abort=TRUE)
    abundance <- data.table::fread(file, drop=c("label", "numOtus"), header=T, blank.lines.skip=T) %>%
        tibble::tibble() %>%
        tibble::column_to_rownames("Group")
    return(abundance)
}

#' Load OTU abundance file from mothur (count table format)
#'
#' @param file Path to count_table file
#' @return Abundance dataframe (OTU x sample)
#' @examples
#' load_count_table("OTUs.100.rep.count_table")
load_count_table <- function(file) {
    check_ext(file, "count_table", abort=TRUE)
    abundance <- data.table::fread(file, drop=c("total"), header=T, blank.lines.skip=T) %>%
        tibble::tibble() %>%
        tibble::column_to_rownames("Representative_Sequence")
    return(abundance)
}

#' Load OTU abundance file from mothur
#' (for either count_table or shared format)
#'
#' @param file Path to shared file
#' @return Abundance dataframe (sample x OTU)
#' @export
#' @examples
#' load_abundance("abundance.100.shared")
#' load_abundance("abundance.100.rep.count_table")
load_abund <- function(file) {
    if (grepl("\\.shared$", file)) {
        abundance <- load_shared(file)
    } else {
        abundance <- load_count_table(file)
    }
    return(abundance)
}

#' Load consensus taxonomy from mothur
#' 
#' @param file Path to cons.taxonomy file
#' @return Taxonomy dataframe (OTU x taxonomic ranks)
#' @export
#' @examples
#' load_constaxonomy("annotations.100.cons.taxonomy")
load_constaxonomy <- function(file) {
    check_ext(file, "taxonomy", abort=FALSE)
    tax_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    taxonomy <- read.table(file, header=T, row.names=1, sep="\t") %>%
        tidyr::separate(Taxonomy, tax_ranks, sep=";") %>%
        dplyr::select(-c(Size, Species))
    return(taxonomy)
}

#' Load metadata file
#' 
#' @param file Path to metadata file
#' @return Data Frame
#' @export
#' @examples
#' load_metadata("metadata.csv")
load_metadata <- function(file) {
    check_ext(file, "[ct]sv", abort=FALSE)

    metadata <- data.table::fread(file, header=T, blank.lines.skip=T)  %>% tibble::tibble()
    metadata <- metadata %>% 
        tibble::column_to_rownames(colnames(metadata)[1])

    return(metadata)
}

#' Check if the output folder from the cmaiki pipeline
#' contains all the necessary files (abund, tax) for a
#' given OTU identity theshold and returns
#' the paths in a list
#' 
#' @param folder Path to pipeline output folder
#' @param otu_id OTU Identity threshold
#' @return list of file paths
#' @export
#' @examples
#' find_pipeline_files("2022-07-01/")
find_pipeline_files <- function(folder, otu_id=100) {
    ## Special case for abundance file
    ## since there are 2 valid patterns (rep.count_table or .shared)
    abund_file <- find_file(folder, sprintf(".%s.rep.count_table$", otu_id), abort=FALSE)
    if(length(abund_file) == 0) {
        abund_file <- find_file(folder, sprintf(".%s.shared$", otu_id), abort=TRUE)
    }

    return(list(
        abund=abund_file,
        tax=find_file(folder, sprintf(".%s.cons.taxonomy$", otu_id), abort=TRUE),
        metadata=find_file(folder, ".[ct]sv$", warning=FALSE),
        tree=find_file(folder, sprintf(".%s.nwk$", otu_id), warning=FALSE)
    ))
}


#' Load outputs from C-MAIKI pipeline (MetaFlow|mics)
#' into phyloseq. You can either provide the output
#' folder from the pipeline run or the path to all
#' the individual files.
#'
#' @param folder Output folder from MetaFlow|mics
#' @param otu_id OTU clustering threshold in case
#' multiple values where used
#' @param abund_file Path to abundance file
#' (either .count_table or .shared)
#' @param tax_file Path to .cons.taxonomy file
#' @param metadata_file Path to sample metadata (csv formatted).
#' @param tree_file Path to phylogenetic tree (otpional)
#' Assumes there is a header and the first column are the sample names
#' @return phyloseq object
#' @export
#' @examples
#' load_cmaiki("pipeline-2022-07-9/", )
#' load_cmaiki(abund_file="OTUs.100.rep.count_table", tax_file="OTUs.100.cons.taxonomy", metadata="proj/metadata.csv")
load_cmaiki <- function(folder=NULL, abund_file=NULL, tax_file=NULL, metadata_file=NULL, tree_file=NULL, otu_id=100){

    if(is.null(folder) && (is.null(abund_file) | is.null(tax_file) | is.null(metadata_file))) {
        stop("You must either provide the folder or all the individual files")
    } else if(!is.null(folder)) {
        if(!is.null(abund_file) | !is.null(tax_file)) {
            stop("You can't provide both the folder and all the individual files")
        }
        # Find all files in provided folder
        files <- find_pipeline_files(folder, otu_id)
        abund_file <- files$abund
        tax_file <- files$tax
        metadata_file <- ifelse(is.null(metadata_file), files$metadata, metadata_file)
        tree_file <- files$tree
    }
    if(any(is.na(metadata_file) || is.null(metadata_file))) {
        stop("Could not find metadata file.")
    }
    ps_data <- list(
        phyloseq::otu_table(load_abund(abund_file),
                            taxa_are_rows=endsWith(abund_file, "count_table")),
        phyloseq::tax_table(load_constaxonomy(tax_file) %>% as.matrix()),
        phyloseq::sample_data(load_metadata(metadata_file))
    )
    # Add tree if provided
    if(any(!is.null(tree_file) & !is.na(tree_file))) {
        check_ext(tree_file, "nwk", abort=FALSE)
        ps_data <- append(ps_data, phyloseq::read_tree(tree_file))
    }
    # Compile into phyloseq object
    ps <- do.call(phyloseq::phyloseq, ps_data)

    return(ps)
}
