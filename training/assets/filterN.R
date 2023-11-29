#' Date: 2023-11-29
#' This function will be added to enchantr. Look for the most recent
#' version there
#' 
#' filterN: filter sequences with Ns in the whole sequences or in a region
#'
#' Find sequences which have Ns in a particular region and are scattered as opposed
#' to being concentrated in a long stretch. Being in a long stretch could mean these N where
#' introduced during alignment to a reference while being scattered is more
#' likely to be a sign of poor quality.
#'
#' @param   db    data.frame containing sequence data
#' @param   sequence    column containing the sequence to be analyzed
#' @param   start   1   last position to subset the region in the sequence
#'                      to be analyzed or name of the column with this information
#' @param   end     last position to subset the region in the sequence
#'                  to be analyzed or the name of the column with this
#'                  information. Defaults to the length of \code{sequence}
#' @param   max_n   Threshold for the maximum total number of N characters allowed.
#' @param   max_n_stretch   Threshold for the maximum number of stretches of N
#'                          characters allowed.
#' @param  color   name of the field(s) in \code{db} use to color the plot
#' @param  label   prefix added to the name of the output fields.
#' @param  plot
#' @return  If \code{plot=TRUE}, it will return a list with
#'          \itemize{
#'            \item \code{db_pass}: a T/F vector showing if the rows in \code{db} pass
#'                                  meet the filtering criteria
#'            \item \code{p}:  a \code{ggplot} object
#'          }
#'          If \code{plot=FALSE}, it will return the original \code{db} with additional
#'          columns with names ending in:
#'          \itemize{
#'               \item n  Number of Ns in the region
#'               \item n_stretches Number of stretches of N in the region
#'               \item n_pass Whether the sequence meets the filtering criteria
#'          }
#' @examples
#' library(ggplot2)
#' data(ExampleDb, package="alakazam")
#' filterN(ExampleDb, color=NULL)
#' filterN(ExampleDb, color="rev_comp")
#' @export
filterN <- function(db, sequence="sequence_alignment", start=1, end=NULL,
                    max_n=20, max_n_stretch=15,
                    color=NULL, label=NULL, plot=TRUE) {
    
    check <- alakazam::checkColumns(db, sequence)
    if (any(check != TRUE)) { stop(check)}
    
    if (length(color)>1) {
        stop("The use of more than one variable to color the plot is not implemented. You can create a feature request in https://bitbucket.org/kleinstein/enchantr/issues. ")
    }
    
    if (is.character(start)) {
        check <- alakazam::checkColumns(db, start)
        if (any(check != TRUE)) { stop(check)}
        start <- db[[start]]
    }
    
    if (is.character(end)) {
        check <- alakazam::checkColumns(db, end)
        if (any(check != TRUE)) { stop(check)}
        end <- db[[end]]
    } else if (is.null(end)) {
        end <- nchar(db[[sequence]])
    }
    
    filter_n_region <- substr(
        db[[sequence]],
        start,
        end)
    
    if (is.null(label)) {
        sep <- ""
    } else {
        sep="_"
    }
    n_label <- paste(label,"n", sep=sep)
    n_stretches_label <- paste(label,"n_stretches",sep=sep)
    n_pass_label <- paste(label, "n_pass", sep=sep)
    
    db[[ n_label ]] <- nchar(gsub("[^Nn]","",filter_n_region))
    db[[ n_stretches_label ]] <- sapply(filter_n_region, function(seq) {
        matches <- gregexpr("[Nn]+", gsub("\\.","",seq))[[1]]
        n_stretches <- 0
        if (matches[1]>0) {
            n_stretches <- length(matches)
        }
        n_stretches
    }, USE.NAMES = F)
    
    db[[n_pass_label]] <- db[[ n_label ]] <= max_n & db[[ n_stretches_label ]] <= max_n_stretch
    
    if (plot) {
        p <- ggplot(db, aes(x=!!rlang::sym(n_label),
                            y=!!rlang::sym(n_stretches_label),
                            shape=!!rlang::sym(n_pass_label)))
        if (!is.null(color)) {
            p <- p + aes(color=!!rlang::sym(color))
        }
        p <- p +
            geom_abline(slope = 1, intercept = 0, color="grey80") +
            geom_vline(xintercept = max_n, color="grey80") +
            geom_hline(yintercept = max_n_stretch, color="grey80") +
            geom_point(alpha=0.5) +
            xlab("Number of N") +
            ylab("Number of N stretches")
        list("db_pass"=db[[n_pass_label]],
             "p"=p)
    } else {
        db
    }
}

