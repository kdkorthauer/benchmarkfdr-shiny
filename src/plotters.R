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


