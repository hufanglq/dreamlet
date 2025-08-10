#' Processing expression data from assay
#'
#' For raw counts, filter genes and samples, then estimate precision weights using linear mixed model weighting by number of cells observed for each sample.  For normalized data, only weight by number of cells.
#'
#' @param y matrix of counts or log2 CPM
#' @param formula regression formula for differential expression analysis
#' @param data metadata used in regression formula
#' @param n.cells array of cell count for each sample
#' @param min.cells minimum number of observed cells for a sample to be included in the analysis
#' @param min.count used to compute a CPM threshold of \code{CPM.cutoff = min.count/median(lib.size)*1e6}.  Passed to \code{edgeR::filterByExpr()}
#' @param min.samples minimum number of samples passing cutoffs for cell cluster to be retained
#' @param min.prop minimum proportion of retained samples with \code{CPM > CPM.cutoff}
#' @param min.total.count minimum total count required per gene for inclusion
#' @param isCounts logical, indicating if data is raw counts
#' @param normalize.method normalization method to be used by \code{calcNormFactors}
#' @param span Lowess smoothing parameter using by \code{variancePartition::voomWithDreamWeights()}
#' @param quiet show messages
#' @param weights matrix of precision weights
#' @param rescaleWeightsAfter default = FALSE, should the output weights be scaled by the input weights
#' @param BPPARAM parameters for parallel evaluation
#' @param cache logical, whether to use formula caching for repeated operations
#' @param ... other arguments passed to \code{dream}
#'
#' @return \code{EList} object storing log2 CPM and precision weights
#'
#' @importFrom BiocParallel SerialParam bplapply
#' @importClassesFrom limma EList
#' @importFrom variancePartition voomWithDreamWeights
#' @importFrom edgeR calcNormFactors filterByExpr DGEList
#' @importFrom methods is new
#' @importFrom stats model.matrix var
#' @importFrom SummarizedExperiment colData assays
#' @importFrom S4Vectors as.data.frame
#' @importFrom lme4 subbars
#' @importFrom MatrixGenerics colMeans2
#' @importFrom Matrix Matrix sparse.model.matrix
#' @importFrom data.table as.data.table setkey setDT
#' @importFrom digest digest
#'
processOneAssay <- function(y, formula, data, n.cells, min.cells = 5, min.count = 2, 
                          min.samples = 4, min.prop = .4, min.total.count = 15, 
                          isCounts = TRUE, normalize.method = "TMM", span = "auto", 
                          quiet = TRUE, weights = NULL, rescaleWeightsAfter = FALSE, 
                          BPPARAM = SerialParam(), cache = TRUE, ...) {
    
    # Create environment for formula caching if enabled
    if(cache) {
        if(!exists("formula_cache")) {
            formula_cache <- new.env(parent = emptyenv())
        }
    }
    
    # Early validation
    if (is.null(n.cells)) {
        stop("n_cells must not be NULL")
    }
    
    # Convert to sparse matrix if possible
    if (!is(y, "dgCMatrix") && requireNamespace("Matrix", quietly = TRUE)) {
        y <- Matrix::Matrix(as.matrix(y), sparse = TRUE)
    }
    
    # Efficient sample filtering
    include <- n.cells >= min.cells
    if (sum(include) == 0) return(NULL)
    
    # Efficient subsetting
    y <- y[, colnames(y)[include], drop = FALSE]
    data <- droplevels(data[include, , drop = FALSE])
    
    if (nrow(data) < min.samples || nrow(y) == 0) return(NULL)
    if (!isCounts) stop("isCounts = FALSE is not currently supported")
    
    # Efficient DGEList creation and normalization
    y <- suppressMessages({
        dge <- DGEList(y, remove.zeros = TRUE)
        calcNormFactors(dge, method = normalize.method)
    })
    
    # Cache formula processing
    if(cache) {
        cache_key <- digest::digest(list(formula, names(data)))
        if (exists(cache_key, envir = formula_cache)) {
            formula <- get(cache_key, envir = formula_cache)
        } else {
            formula <- removeConstantTerms(formula, data)
            formula <- dropRedundantTerms(formula, data)
            assign(cache_key, formula, envir = formula_cache)
        }
    } else {
        formula <- removeConstantTerms(formula, data)
        formula <- dropRedundantTerms(formula, data)
    }
    
    # Efficient gene filtering
    keep <- suppressWarnings({
        filterByExpr(y, 
                    min.count = min.count, 
                    min.prop = min.prop, 
                    min.total.count = min.total.count)
    })
    
    if (sum(keep) == 0) return(NULL)
    
    # Efficient weight processing
    if (!is.null(weights) && !is(weights, "function")) {
        if (!all(rownames(y)[keep] %in% rownames(weights))) {
            stop("All genes retained in count matrix must be present in weights matrix.\n",
                 "Make sure getExprGeneNames() and processAssays() use same parameter values.")
        }
        precWeights <- weights[rownames(y)[keep], colnames(y)]
    } else {
        precWeights <- rep(1, ncol(y))
    }
    
    # Efficient voom processing
    geneExpr <- voomWithDreamWeights(y[keep, ], formula, data,
                                    weights = precWeights,
                                    rescaleWeightsAfter = rescaleWeightsAfter,
                                    BPPARAM = BPPARAM, ...,
                                    save.plot = TRUE,
                                    quiet = quiet,
                                    span = span,
                                    hideErrorsInBackend = TRUE)
    
    if (is.null(geneExpr) || nrow(geneExpr) == 0) return(NULL)
    
    geneExpr$formula <- formula
    geneExpr$isCounts <- isCounts
    
    return(geneExpr)
}

#' Processing SingleCellExperiment to dreamletProcessedData
#'
#' @inheritParams processOneAssay
#' @param sceObj SingleCellExperiment object
#' @param assays array of assay names to include in analysis
#' @param weightsList list storing matrix of precision weights for each cell type
#'
#' @return Object of class \code{dreamletProcessedData}
#'
#' @export
processAssays <- function(sceObj, formula, assays = assayNames(sceObj), 
                         min.cells = 5, min.count = 5, min.samples = 4, 
                         min.prop = .4, isCounts = TRUE, normalize.method = "TMM", 
                         span = "auto", quiet = FALSE, weightsList = NULL, 
                         BPPARAM = SerialParam(), num_workers = 4L, cache = TRUE, ...) {
    
    # Input validation
    stopifnot(is(sceObj, "SingleCellExperiment"))
    stopifnot(is(formula, "formula"))
    
    if (is.null(colnames(sceObj))) {
        stop("colnames(sceObj) is NULL. Column names are needed for internal filtering")
    }
    
    # Extract and prepare metadata
    data_constant <- droplevels(as.data.table(as.data.frame(colData(sceObj)), keep.rownames = TRUE))
    row.names(data_constant) <- data_constant$rn
    
    # Validate assays
    invalid_assays <- setdiff(assays, assayNames(sceObj))
    if (length(invalid_assays) > 0) {
        stop("Assays not found in dataset: ", paste(head(invalid_assays), collapse = ", "))
    }
    
    # Validate weightsList
    if (!is.null(weightsList)) {
        missing_assays <- setdiff(assays, names(weightsList))
        if (length(missing_assays) > 0) {
            stop("Assays not found in weightsList: ", paste(missing_assays, collapse = ", "))
        }
    }
    
    # Extract cell counts efficiently
    n.cells_full <- cellCounts(sceObj)
    colNamesAll <- unique(unlist(lapply(assayNames(sceObj), colnames)))
    
    if (any(!colNamesAll %in% rownames(n.cells_full))) {
        stop("Cell counts extraction failed. Check that colnames(sceObj) or rownames(colData(sceObj)) are intact")
    }
    
    # Process assays in parallel
    if (!quiet) {
        message("Processing ", length(assays), " assays...")
        pb <- txtProgressBar(min = 0, max = length(assays), style = 3)
    }

    suppressMessages(require(future))
    suppressMessages(require(furrr))
    
    if (Sys.info()['sysname'] == "Windows") {
      plan(multisession, workers = num_workers)
    } else {
      plan(multiprocess, workers = num_workers)
    }
    
    resList <- future_map(seq_along(assays), function(i) {
        suppressMessages(require(stats))
        k <- assays[i]
        if (!quiet) setTxtProgressBar(pb, i)
        
        startTime <- Sys.time()
        y <- assay(sceObj, k)
        n.cells <- n.cells_full[colnames(y), k, drop = FALSE]
        
        # Merge metadata efficiently using data.table
        data <- merge_metadata_dt(
            data_constant,
            get_metadata_aggr_means(sceObj),
            k,
            metadata(sceObj)$agg_pars$by
        )
        
        # Process weights
        weights <- if (!is.null(weightsList)) {
            weightsList[[k]][, rownames(data), drop = FALSE]
        } else {
            matrix(1, nrow = nrow(y), ncol = nrow(data))
        }
        
        if (!is.null(weights)) {
            colnames(weights) <- rownames(data)
            rownames(weights) <- rownames(y)
        }
        
        # Process assay
        result <- processOneAssay(
            y[, rownames(data), drop = FALSE],
            formula = formula,
            data = data,
            n.cells = n.cells[rownames(data), , drop = FALSE],
            min.cells = min.cells,
            min.count = min.count,
            min.samples = min.samples,
            min.prop = min.prop,
            isCounts = isCounts,
            normalize.method = normalize.method,
            span = span,
            weights = weights,
            BPPARAM = BPPARAM,
            cache = cache,
            ...
        )
        
        if (!quiet) {
            message(k, " processed in ", format(Sys.time() - startTime, digits = 2))
        }
        
        return(result)
    })
    
    if (!quiet) close(pb)
    
    names(resList) <- assays
    
    # Handle empty results
    exclude <- vapply(resList, is.null, logical(1))
    if (any(exclude)) {
        warning("Not enough samples retained or model fit fails: ",
                paste(names(resList)[exclude], collapse = ", "))
    }
    resList <- resList[!exclude]
    
    # Process errors
    error.initial <- lapply(resList, `[[`, "error.initial")
    errors <- lapply(resList, attr, "errors")
    names(error.initial) <- names(errors) <- names(resList)
    
    # Create details dataframe efficiently
    df_details <- data.table(
        assay = names(resList),
        n_retain = vapply(resList, ncol, integer(1)),
        formula = vapply(resList, function(x) Reduce(paste, deparse(x$formula)), character(1)),
        formDropsTerms = vapply(resList, function(x) !equalFormulas(x$formula, formula), logical(1)),
        n_genes = vapply(resList, nrow, integer(1)),
        n_errors = vapply(resList, function(x) length(attr(x, "errors")), integer(1)),
        error_initial = vapply(resList, function(x) !is.null(x$error.initial), logical(1))
    )
    
    # Report warnings and errors
    ndrop <- sum(df_details$formDropsTerms)
    if (ndrop > 0) {
        warning("Terms dropped from formulas for ", ndrop, " assays.\n",
                "Run details() on result for more information")
    }
    
    failure_frac <- sum(df_details$n_errors) / sum(df_details$n_genes)
    if (is.nan(failure_frac)) {
        stop("All models failed. Consider changing formula")
    }
    
    if (failure_frac > 0) {
        message("\nOf ", format(sum(df_details$n_genes), big.mark = ","),
                " models fit across all assays, ",
                format(failure_frac * 100, digits = 3), "% failed\n")
    }
    
    # Return results
    new("dreamletProcessedData",
        resList,
        data = as.data.frame(data_constant),
        metadata = get_metadata_aggr_means(sceObj),
        by = metadata(sceObj)$agg_pars$by,
        df_details = as.data.frame(df_details),
        errors = errors,
        error.initial = error.initial
    )
}

# merge data_constant (data constant for all cell types)
# with metadata(sceObj)$aggr_means (data that varies)
#' @importFrom dplyr filter_at
merge_metadata <- function(dataIn, md, cellType, by) {
  # PASS R CMD check
  cell <- NULL

  data <- merge(dataIn,
    dplyr::filter_at(md, by[1], ~ . == cellType),
    by.x = "row.names",
    by.y = by[2]
  )
  rownames(data) <- data$Row.names
  id <- rownames(dataIn)[rownames(dataIn) %in% rownames(data)]
  data <- data[id, , drop = FALSE]
  droplevels(data)
}

# Optimized metadata merging using data.table
merge_metadata_dt <- function(dataIn, md, cellType, by) {
    suppressMessages(require(data.table))
    if(!inherits(dataIn, 'data.table'))
      dt1 <- as.data.table(dataIn, keep.rownames = TRUE)
    else
      dt1 <- dataIn
    dt2 <- as.data.table(md)
    
    setkeyv(dt2, by[1])
    dt2 <- dt2[get(by[1]) == cellType]
    
    result <- dt1[dt2, on = c("rn" = by[2])]
    setDF(result)
    rownames(result) <- result$rn
    result$rn <- NULL
    
    return(droplevels(result))
}
 
