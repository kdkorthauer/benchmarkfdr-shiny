#' Scatterplot of nominal alpha vs rejections
#' 
#' This function inputs a SummarizedBenchmark object build by the
#' common benchDesign in this repository. It calculates the number 
#' of rejections for several nominal alpha values and plots them. 
#'   
#' @param sb A SummarizedBenchmark object
#' @param asFraction Logical. Whether to plot the fraction of hypotheses
#' that were rejected (TRUE) or the absolute number of rejections (FALSE).
#' @param supplementary Logical. Is the figure a supplementary figure?
#'        
#' @return a ggplot2 object
#' 
#' 
#' @author Alejandro Reyes
rejections_scatter <- function( sb, as_fraction=FALSE, supplementary=TRUE,
                                palette = candycols, alpha = 0.05){
  stopifnot( is(sb, "SummarizedBenchmark" ) )
  alphas <- unique( as.numeric( as.character( colData( sb )$param.alpha ) ) )
  alphas <- alphas[!is.na(alphas)]
  if( as_fraction ){
    deno <- nrow(sb)
    yl <- "Fraction of hypotheses rejected"
  }else{
    deno <- 1
    yl <- "Number of rejections"
  }
  if(sum(colSums(is.na(assays(sb)$qvalue)) == nrow(sb)) > 0){
    miss <- which(colSums(is.na(assays(sb)$qvalue)) == nrow(sb)) 
    sb <- sb[,-miss]
  }
  
  plotDF <- estimatePerformanceMetrics( sb, alpha=alphas, tidy=TRUE ) %>%
    dplyr:::filter( !(grepl("ihw", blabel) & param.alpha != alpha ) ) %>%
    dplyr:::mutate( blabel=gsub("(ihw)-a\\d+", "\\1", blabel ) ) %>%
    #dplyr:::select( blabel, key, value, assay, performanceMetric, alpha ) %>%
    dplyr:::filter( blabel!="unadjusted", performanceMetric == "rejections" ) 
  if( !supplementary ){
    plotDF <- plotDF %>%
      dplyr:::mutate( param.smooth.df=gsub("L", "", param.smooth.df ) ) %>%
      dplyr:::filter( !( grepl("bl-df", blabel) & 
                           as.numeric(as.character(param.smooth.df) != 3 ) ) ) %>%
      dplyr:::mutate( Method=gsub("-df03", "", blabel))
  }
  
  # add color palette
  plotDF <- dplyr::left_join(plotDF, palette, by="Method") 
  
  col <- as.character(plotDF$col)
  names(col) <- as.character(plotDF$Method)
  
  lty <- as.character(plotDF$lty)
  names(lty) <- as.character(plotDF$Method)
  
  plotDF %>%
    ggplot( aes(alpha, value/deno, col=Method) ) +
    geom_line(alpha = 3/4, aes(linetype=Method)) + 
    geom_point(alpha = 3/4, show.legend = FALSE) +
    xlab(expression(paste("Nominal"~alpha))) +
    scale_color_manual(values=col) +
    scale_linetype_manual(values = lty) +
    ylab(yl)
}

#' Plot overlaps between FDR methods
#'   
#' @param object A SummarizedBenchmark object
#' @param alpha An alpha value.
#' @param supplementary Logical. Is the figure a supplementary figure?
#'        
#' @return an upsetr plot
#' 
#' @author Alejandro Reyes
plotFDRMethodsOverlap <- function( object, supplementary=TRUE, alpha=0.1, ... ){
  stopifnot( is( object, "SummarizedBenchmark" ) )
  stopifnot( any( colData(object)$param.alpha == alpha, na.rm=TRUE) )
  object <- object[,!( grepl("^ihw", as.character( colData( object )$blabel ) ) & colData( object )$param.alpha != alpha )]
  colData(object)$blabel <- gsub("(ihw)-.*", "\\1", colData( object )$blabel)
  qvals <- assays( object )[["qvalue"]]
  object <- object[,!apply( is.na( assays( object )[["qvalue"]] ), 2, all )]
  object <- object[,colData(object)$blabel != "unadjusted"]
  if( !supplementary ){
    object <- object[,!( grepl("bl", colData( object )$blabel) & as.numeric( gsub( "L", "", colData( object )$param.smooth.df ) ) != 3 )]
    colData(object)$blabel <- gsub( "(bl)-.*", "\\1", colData(object)$blabel )
  }
  colnames(object) <- colData( object )$blabel
  plotMethodsOverlap( object, alpha=alpha, ... )
}

#' Standardize Simulation Results
#'
#' Given a list of SummarizedBenchmark objects with a "qvalue" assay,
#' this function calculates performance metrics and returns a single
#' data.table that can be used for plotting with the average_plots function.
#' 
#' @param res list of SummarizedBenchmark objects with "qvalue" assay that
#'        correspond to repilications of the same simulation setting. 
#' @param alpha sequence of alpha cutoffs at which metrics should be
#'        calculated. (default = seq(0.01, 0.10, 0.01))
#' 
#' @return
#' a tibble of the following performance metrics:
#' 'rejections', 'TPR', 'TNR', 'FPR', 'FNR', 'FWER'.
#' 
#' @author Patrick Kimes & Keegan Korthauer
plotsim_standardize <- function(res, alpha = seq(0.01, 0.10, 0.01)) {
  
  sbl <- lapply(res, addDefaultMetrics)
  sbl <- lapply(sbl, addPerformanceMetric,
                evalMetric = "FWER", assay = "qvalue",
                evalFunction =
                  function( query, truth, alpha = 0.1) {
                    any((query < alpha) & !as.logical(truth), na.rm = TRUE)
                  })
  sbl <- lapply(sbl, addPerformanceMetric,
                evalMetric = "rejectprop", assay = "qvalue",
                evalFunction =
                  function( query, truth, alpha = 0.1) {
                    mean(query < alpha, na.rm = TRUE)
                  })
  
  # add FPR metric for ROC curve
  sbl <- lapply(sbl, addPerformanceMetric,
                assay="qvalue", evalMetric="FPR",
                evalFunction = 
                  function( query, truth, alpha=0.1 ){
                    sum( query <= alpha & truth == 0, na.rm = TRUE ) / 
                      sum( truth == 0, na.rm = TRUE )
                  })
  
  # override default FDR perf. metric st fdr is zero (instead of NA) when
  # there are zero rejections
  sbl <- lapply(sbl, addPerformanceMetric,
                evalMetric = "FDR", assay = "qvalue",
                evalFunction =
                  function( query, truth, alpha = 0.1) {
                    rejections <- sum( query <= alpha, na.rm = TRUE )
                    fdr <- sum( query <= alpha & truth == 0, na.rm = TRUE ) / 
                      sum( query <= alpha, na.rm = TRUE )
                    fdr[rejections == 0] <- 0
                    return(fdr)
                  })
  
  tsb <- lapply(sbl, estimatePerformanceMetrics,
                alpha = alpha, tidy = TRUE)
  
  tsb <- bind_rows(tsb, .id = "rep")
  as_tibble(tsb)
}


#' Plot Mean Performance Metrics Across Replications
#'
#' Given a standardized metric table created using "standardize_results",
#' this function calculates the mean metric across replications and
#' generates a standard plot.
#' 
#' @param tsb standardized metric data.frame generated using
#'        standardize_results.
#' @param met name(s) of metric to plot - must be one of the performance
#'        metrics added by default with "addDefaultMetrics" or FWER or FPR.
#'        May be a vector of names of length 2 - in this case the first will
#'        be used for the x-axis and the second for the y-axis. 
#' @param filter_set character vector of "blabel" IDs of methods which
#'        should be included in the plot. Alternatively, a subset of
#'        methods can be chosen by filtering the data.table before passing
#'        to this function. (default = NULL)
#' @param merge_ihw logical whether IHW results should be merged across
#'        alpha values by matching "param.alpha" and "alpha" columns.
#'        (default = TRUE)
#' @param clean_names logical whether to clean-up method names in 'blabel'
#'        column. The tsb table must only contain one column with each of the
#'        following labels, else method labels will not be changed: 'ashq',
#'        'bh', 'bl', 'ihw', 'lfdr', 'qvalue', 'fdrreg-e',
#'        'fdrreg-t'. (default = FALSE)
#' @param errorBars logical indicating whether to include error bars 
#' @param palette data.frame containing the color palette - should contain 
#'        three columns: 'Method', 'col', and 'lty'
#' @param facetMethodType logical indicating whether to facet the plot into 
#'        two panels - methods that use only one piece of information, and 
#'        those that use more than just the p-value (default = FALSE)
#' @param diffplot logical indicating whether 'value' in tsb is a difference
#'        between informative and uninformative covariates. (default = FALSE)
#' @param grpVars vector of character names of additional columns of tsb to keep 
#' @param type character indicating type of comparison; if not "de", then throw error
#' @param rocstyle logical whether plotting roc style curve
#'
#' @return
#' a ggplot object.
#' 
#' @author Patrick Kimes & Keegan Korthauer
plotsim_average <- function(tsb, met, 
                            filter_set = c("unadjusted", "bl-df02", "bl-df04", "bl-df05"),
                            merge_ihw = TRUE,
                            clean_names = FALSE, errorBars=FALSE,
                            palette = candycols, facetMethodType = FALSE,
                            diffplot = FALSE, grpVars = NULL, type = "de",
                            rocstyle = FALSE){
  if (type != "de"){
    #Can't calculate TPR for null comparison
    return(NULL)
  }
  if (length(met)>2)
    stop("Can only plot 2 metrics at a time")
  
  ## cacluate mean per replication
  if(!is.null(grpVars)){
    tsba <- tsb %>%
      group_by(blabel, performanceMetric, alpha, param.alpha, key, 
               .dots = grpVars) %>%
      summarize(n = sum(!is.na(value)),
                se = sd(value, na.rm = TRUE) / sqrt(n),
                value = mean(value, na.rm = TRUE))
  }else{
    tsba <- tsb %>%
      group_by(blabel, performanceMetric, alpha, param.alpha, key) %>%
      summarize(n = sum(!is.na(value)),
                se = sd(value, na.rm = TRUE) / sqrt(n),
                value = mean(value, na.rm = TRUE))
  }
  
  tsba$value[tsba$n == 0] <- NA
  tsba$se[tsba$n == 0] <- NA
  
  ## group IHW methods if specified
  if (merge_ihw) {
    tsba <- filter(tsba, is.na(param.alpha) | (param.alpha == alpha))
    tsba$blabel[grepl("^ihw-", tsba$blabel)] <- "ihw"
  }
  
  # filter by performance metric
  tsba_m <- NULL 
  for (m in seq_along(met)){
    tmp <- tsba %>% 
      filter(performanceMetric == met[m]) %>%
      dplyr::mutate(Method = gsub("-df03", "", blabel)) 
    if (m > 1){
      tsba_m <- left_join(tsba_m, tmp, 
                          by = c("Method", "alpha", "n", 
                                 "param.alpha", "blabel"))
    }else{
      tsba_m <- tmp
    }
  }
  tsba <- tsba_m
  
  if (clean_names) {
    ulabs <- unique(tsba$blabel)
    vlabs <- c('ashq', 'bh', 'bl', 'ihw', 'lfdr', 'qvalue')
    clabs <- c("ASH q-value", "Benjamini-Hochberg", "Boca-Leek", "IHW", "local FDR", "Storey's q-value")
    if (any(grepl("fdrreg", ulabs))) {
      vlabs <- c(vlabs, 'fdrreg-e', 'fdrreg-t')
      clabs <- c(clabs, "FDRreg (emp)", "FDRreg (theor)")
    }
    if (any(grepl("adapt", ulabs))) {
      vlabs <- c(vlabs, 'adapt-glm', 'adapt-gam')
      clabs <- c(clabs, 'AdaPT (GLM)', 'AdaPT (GAM)')
    }
    lcnts <- sapply(vlabs, function(x) { grep(paste0("^", x), ulabs, value=TRUE) })
    if (!is(lcnts, "list")) {
      names(lcnts) <- clabs
      tsba$blabel <- do.call(forcats::fct_recode, c(list("f" = tsba$blabel), as.list(lcnts)))
    }
  }
  
  # add color palette
  tsba <- dplyr::left_join(tsba, palette, by="Method") 
  
  ## remove methods if any specified
  if (!is.null(filter_set)) {
    tsba <- filter(tsba, (blabel %in% filter_set))
  }
  
  col <- as.character(tsba$col)
  names(col) <- as.character(tsba$Method)
  
  lty <- as.character(tsba$lty)
  names(lty) <- as.character(tsba$Method)
  
  # add method type
  tsba <- tsba %>%
    mutate(Type = ifelse(Method %in% c("unadjusted", "bonf", "bh", "qvalue"),
                         "Univariate (p-value only)", "Multivariate"))
  if(length(met)==1){
    gp <- tsba %>%
      ggplot(aes(x = alpha, y = value, color = Method)) +
      geom_line(alpha = 0.85, aes(linetype=Method)) +
      theme_classic() +
      theme(axis.title = element_text(face="bold"),
            plot.title = element_text(face="bold")) +
      expand_limits(x = 0) +
      scale_x_continuous("alpha cutoff", breaks=seq(0, 1, by=0.01)) +
      ylab(met) +
      ggtitle(ifelse(diffplot,
                     paste0("Mean Difference ", met, " Over ", max(tsba$n), " Replications"),
                     paste0("Mean ", met, " Over ", max(tsba$n), " Replications"))) +
      scale_color_manual(values = col) +
      scale_linetype_manual(values = lty)
  }else{
    gp <- tsba %>% 
      ggplot(aes(x = value.x, y = value.y, color = Method)) +
      geom_line(alpha = 0.85, aes(linetype=Method)) +
      theme_classic() +
      theme(axis.title = element_text(face="bold"),
            plot.title = element_text(face="bold")) +
      expand_limits(x = 0) +
      scale_x_continuous(breaks=seq(0, 1, by=0.01)) +
      ylab(met[2]) +
      xlab(met[1]) +
      ggtitle("Mean ROC Curve over 100 Replications") +
      scale_color_manual(values = col) +
      scale_linetype_manual(values = lty)
  }
  
  if(facetMethodType){
    gp <- gp +
      facet_wrap( ~ Type)
  }
  
  if(errorBars & length(met)==1){
    gp <- gp + geom_errorbar(width=0.0025, alpha=0.5,
                             aes(ymin=value-se, ymax=value+se))
  }
  
  ## use percentage on y-axis labels when appropriate
  if (sum(met %in% c("FDR", "FNR", "TPR", "TNR", "FWER", "rejectprop")) > 0){
    if (length(met)==1){
      gp <- gp +
        scale_y_continuous(ifelse(diffplot, paste(met, "(informative - uninformative)"), met),
                           labels=scales::percent)
    }else{
      gp <- gp +
        scale_x_continuous(met[1], labels=scales::percent) +
        scale_y_continuous(met[2], labels=scales::percent) 
    }
  } else if (diffplot) {
    gp <- gp + ylab(paste(met, "(informative - uninformative)")) 
  }
  
  
  if (diffplot) {
    gp <- gp + expand_limits(y = 0)
    gp <- gp + geom_hline(yintercept = 0, lty = 2, color = "blue", alpha = 1/2)
  } else {
    ## include 0% or 100% in plotting range depending on metric
    if (sum(met %in% c("FDR", "FNR", "FWER", "rejectprop")) > 0){
      if (length(met)==1)
        gp <- gp + expand_limits(y = 0)
    }
    if (sum(met %in% c("TNR")) > 0){
      if (length(met)==1)
        gp <- gp + expand_limits(y = 1)
    }
    
    ## add identity line for FPR/FDR plotting
    if (sum(met == "FDR")>0) {
      if (length(met)==1)
        gp <- gp +
          geom_abline(intercept = 0, slope = 1, lty = 2, color = "blue", alpha = 1/2)
    }
  }
  
  if(rocstyle){
    # assumes plotsim_average(tsb, met=c("FDR", "TPR"))
    sim_res <- tsb
    sim_res_sum <- sim_res %>% 
      group_by(blabel, performanceMetric, alpha, param.alpha, key) %>%
      summarize(n = sum(!is.na(value)),
                se = sd(value, na.rm = TRUE) / sqrt(n),
                value = mean(value, na.rm = TRUE)) %>%
      filter(round(alpha,2) %in% c(0.01, 0.05, 0.1),
             performanceMetric %in% c("FDR", "TPR"),
             blabel %in% filter_set,
             is.na(param.alpha) | (param.alpha == alpha)) 
    tsba_m <- NULL 
    met = c("FDR", "TPR")
    for (m in seq_along(met)){
      tmp <- sim_res_sum %>% 
        filter(performanceMetric == met[m]) %>%
        dplyr::mutate(Method = gsub("-df03", "", blabel)) 
      if (m > 1){
        tsba_m <- left_join(tsba_m, tmp, 
                            by = c("Method", "alpha", "n", 
                                   "param.alpha", "blabel"))
      }else{
        tsba_m <- tmp
      }
    }
    sim_res_sum <- tsba_m
    sim_res_sum <- sim_res_sum %>%
      mutate(control = value.x < alpha) 
    # add color palette
    sim_res_sum$Method[grepl("^ihw-", sim_res_sum$Method)] <- "ihw"
    sim_res_sum <- dplyr::left_join(sim_res_sum, candycols, by="Method") 
    
    gp <- gp +
      geom_point(data=sim_res_sum, aes(x = value.x, y = value.y, 
                                       shape = control), size = 3.5) +
      geom_point(data=sim_res_sum %>% filter(!control), aes(x = value.x, y = value.y), 
                 shape = 19, size = 2.5, color = "white") +
      geom_point(data=sim_res_sum %>% filter(!control), aes(x = value.x, y = value.y), 
                 shape = 1, size = 3.5) +
      scale_shape_manual(values = c(1, 19)) +
      scale_x_continuous(breaks=seq(0, 1, by=0.02)) +
      labs(shape = "FDR control") +
      ggtitle("")
    
  }
  
  return(gp)
}


#' Aggregated UpSet Plot
#'
#' Generate an UpSet plot showing the (rounded) median overlap between methods
#' across a collection of SummarizedBenchmark objects, e.g. corresponding
#' to simulation replicates.
#' 
#' @param res list of SummarizedBenchmark objects to be combined in the plot.
#' @param alpha significance threshold to use for distinguishing significant
#'        and non-significant tests.
#' @param supplementary logical whether plot is for supplementary materials.
#'        (default = FALSE)
#' @param return_list logical whether frequency list should be returned instead
#'        of upset plot. The returned list can be used to generate the upset plot
#'        using `upset(fromExpression(freq_list)`. This can be useful if the user
#'        wants to experiment with upset plot styles. (default = FALSE)
#' @param nintersects scalar value representing number of sets to look at.  
#'        Default is 40 (same as default in UpSetR package).
#' @param filter_set optional character vector of methods to restrict plotting to.
#'        Should be in "clean" format (IHW/BL simplified).
#' 
#' @return
#' an upset plot if `return_list = FALSE` (default), else a list of frequencies that
#' can be used to generate the same upset plot.
#'
#' @details
#' Note: this can get incredibly slow if the number of methods being compared is
#' large since the possible number of overlaps grows exponentially. Anecdotally,
#' in simulations comparing 9 methods (+ truth) with 20,000 tests takes approximately
#' 15 to 20 seconds for 20 replications, and 40 to 45 seconds with 100 replications.  
#' 
#' @import dplyr magrittr
#' @author Patrick Kimes and Keegan Korthauer
aggupset <- function(res, alpha, supplementary = FALSE, return_list = FALSE,
                     nintersects = 40, filter_set = NULL) { 
  
  ## find significant hits at alpha cutoff for all replicates
  hits_tabs <- lapply(res, sb2hits, a = alpha, s = supplementary)
  
  # check if enough methods with rejections to compute overlaps
  if(any(sapply(lapply(hits_tabs, colSums), function(x) sum(x > 0) <= 1))){
    #message("Not enough methods reject anything")
    return(NULL) 
  }
  
  if(!is.null(filter_set))
    hits_tabs <- suppressWarnings(lapply(hits_tabs, 
                                         function(x) select(x, c(truth, one_of(filter_set)))))
  
  ## replace NAs with 0s (not called significant)
  fails <- lapply(hits_tabs, sapply, function(x) { all(is.na(x)) })
  hits_tabs <- lapply(hits_tabs, function(x) { x[is.na(x)] <- 0; x })
  
  ## count up frequencies in each intersection
  n_cols <- unique(sapply(hits_tabs, ncol))
  if (length(n_cols) > 1) {
    stop("not all SummarizedBenchmarks have the same set of methods")
  }
  freq_tabs <- lapply(hits_tabs, hits2freq, nm = n_cols)
  
  ## convert anything that failed completely to NAs
  freq_tabs <- mapply(function(x, y) {
    if (!any(x)) { return(y) }
    failid <- make.names(names(x))[x]
    failid <- match(failid, names(y))
    y$freq[rowSums(y[, failid]) > 0] <- NA
    y
  }, x = fails, y = freq_tabs, SIMPLIFY = FALSE) 
  
  ## merge all freqs into single table - first rename 'freq' columns to 'freq.i' (i = 1..100)
  method_names <- setdiff(names(freq_tabs[[1]]), "freq")
  freq_tab <- mapply(function(itab, idx) { dplyr::rename(itab, !!(paste0("freq.", idx)) := freq) },
                     itab = freq_tabs, idx = 1:length(freq_tabs),
                     SIMPLIFY = FALSE) %>%
    purrr::reduce(dplyr::left_join, by = method_names)
  
  ## summarize across 100 replications of each setting
  freq_tab <- freq_tab %>%
    gather(repl, cnt, starts_with("freq")) %>%
    group_by_at(method_names) %>%
    summarize(freq_mean = round(mean(cnt, na.rm = TRUE))) %>%
    ungroup() %>%
    mutate(freq_mean = ifelse(is.nan(freq_mean), 0, freq_mean)) 
  
  ## convert binary design matrix to UpSetR format (method names separated by "&")
  freq_tab <- freq_tab %>%
    unite("design", method_names, sep = "&", remove = FALSE) %>%
    gather(method, val, -design, -freq_mean) %>%
    mutate(val = ifelse(val, method, "")) %>%
    spread(method, val) %>%
    select(-design) %>%
    unite("setname", method_names, sep = "&") %>%
    mutate(setname = setname %>% gsub("&{2,}", "&", .) %>%
             gsub("^&", "", .) %>%
             gsub("&$", "", .)) 
  
  ## convert to vector to pass to UpSetR package
  freq_list <- freq_tab$freq_mean
  names(freq_list) <- freq_tab$setname
  
  ## return frequency list if requested
  if (return_list) {
    return(freq_list)
  }
  
  ## draw upset plot if frequency list not returned
  upset(fromExpression(freq_list),
        nsets = n_cols,
        nintersects = nintersects,
        mb.ratio = c(0.55, 0.45),
        order.by = "freq",
        decreasing = TRUE,
        set.metadata = list(data = data.frame(sets = method_names,
                                              isTruth = grepl("truth", method_names)),
                            plots = list(
                              list(type = "matrix_rows", 
                                   column = "isTruth",
                                   colors = c("TRUE" = "blue", "FALSE" = "gray"), 
                                   alpha = 0.2))))
}


#' Helper to Parse Significant Hits for Specified Alpha
#'
#' Determines which tests are significant for each method based
#' on a specified alpha threshold and returns as a binary data.frame
#' for easier downstream parsing.
#' 
#' @param x SummarizedBenchmark w/ qvalue assay.
#' @param a alpha cutoff.
#' @param s logical whether for supplementary materials or not.
#'
#' @return
#' data.frame of 0/1s; rows are test, columns are methods.
#'
#' @import dplyr magrittr
#' @author Patrick Kimes (modified by Keegan Korthauer)
sb2hits <- function(x, a, s) {
  ## make quick table of significant tests w/ groundTruth
  ht <- as_tibble(cbind((assay(x, "qvalue") < a) + 0,
                        truth = rowData(x)$qvalue))
  ## keep only IHW matching "alpha" parameter
  ihw_keep <- paste0("ihw-a", sprintf("%02i", 100 * a ))
  if (ihw_keep %in% names(ht)) {
    ## note - using mutate instead of rename so next 'select' call to drop
    ## extra "ihw-*" columns doesn't throw an error if correct alpha was only alpha
    ht <- dplyr::mutate(ht, ihw = get(ihw_keep))   }
  ht <- dplyr::select(ht, -dplyr::contains("ihw-"))
  ## if not plotting for supplementary materials, remove BL w/ multiple DoF 
  if (!s) {
    suppressWarnings({
      ht <- ht %>%
        dplyr::select(-one_of("bl-df02", "bl-df04", "bl-df05")) %>%
        dplyr::rename(bl = `bl-df03`)
    })
  }
  suppressWarnings({
    ht <- dplyr::select(ht, -one_of("unadjusted"))
  })
  as.data.frame(ht)
}


#' Helper to Count Intersection Frequencies Across Methods
#'
#' Counts overlaps/intersections between methods based on the
#' binary data.frame generated using the `sb2hits()` function.
#' This function is just a wrapper to the `Counter()` function
#' in the `UpSetR` package.
#' 
#' @param x data.frame returned by sb2hits
#' @param nm integer number of methods in comparison
#'
#' @return
#' tibble with one (binary) column per method, and a `freq` column, with
#' each row corresponding to a single overlap of methods - methods
#' contained in the overlap are set to 1, those not in the overlap
#' are set to 0 - with the `freq` column containing the number of test
#' statistics in the overlap.
#' 
#' @import dplyr magrittr
#' @importFrom UpSetR Counter
#' @author Patrick Kimes
hits2freq <- function(x, nm) {
  UpSetR:::Counter(x, nm, 1, names(x), nintersections = 2^nm,
                   mbar_color = "gray23",
                   order_mat = "degree", aggregate = "degree",
                   cut = NULL, empty_intersects = TRUE,
                   decrease = TRUE) %>%
    as_tibble() %>%
    select(-x, -color)
}


# Custom color palette to be consistent across plots and case studies
# Specific color choices adapted from Alyssa Frazee's RSkittleBrewer 
# package (https://github.com/alyssafrazee/RSkittleBrewer)
# Also including ideas from discussion in meeting on 4/6/2018

# Bonferroni - dashed black
# BH - solid black
# qvalue - gray45
# ashq - red3
# ashs - hotpink
# ihw - green3
# bl - darkorange1
# locfdr - dodgerblue3
# scott empirical - purple4
# scott theoretical - mediumpurple3
# AdaPT GLM - darkgoldenrod
# AdaPT GAM - darkgoldenrod3

candycols <- data.frame(Method = c("bonf", "bh", "qvalue", "ashs", "ashq", 
                                   "ihw", "bl", "lfdr", 
                                   "fdrreg-t", "fdrreg-e",
                                   "adapt-glm", "adapt-gam"),
                        col = c("black", "black", "gray45", "hotpink", "red3",
                                "green3", "darkorange1", "dodgerblue3", 
                                "purple4", "mediumpurple3",
                                "darkgoldenrod", "darkgoldenrod3"),
                        lty = c("dashed", rep("solid", 11)),
                        stringsAsFactors = FALSE)


addDefaultMetrics <- function(sb) {
  sb <- addPerformanceMetric(sb, evalMetric = c("TPR", "FDR", "TNR", "FNR", "rejections"),
                             assay = "qvalue")
  return(sb)
}

