#' q
#'
#' Shortcut to quot
#'
#' @param save Whether to save the workspace
#' @param ... Passed on to \code{quit}
#'
#' @export
q <- function (save="no", ...) {
  quit(save=save, ...)
}

#' lsos
#'
#' Lists loaded R objects, their dimensions and memory usage
#'
#' @param order.by Which object parameter to sort by
#' @param decreasing Sort order
#' @param head Subset list of objects?
#' @param n Number of objects to print if \code{head} is \code{TRUE}
#'
#' @export
#'
#' @author stolen from somewhere on stackoverflow
lsos <- function (order.by="Size", decreasing=TRUE, head=TRUE, n=10) {
    napply <- function(names, fn) sapply(names, function(x)
                                         fn(get(x, pos=1)))
    names <- base::ls(pos=1)
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.size <- napply(names, object.size)
    obj.prettysize <- sapply(obj.size, function(r) prettyNum(r, big.mark = ",") )
    obj.dim <- t(napply(names, function(x)
                        as.numeric(dim(x))[1:2]))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, length)[vec]
    out <- data.frame(obj.type, obj.size,obj.prettysize, obj.dim)
    names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
    if (!missing(order.by))
        out <- out[order(out[[order.by]], decreasing=decreasing), ]
        out <- out[c("Type", "PrettySize", "Rows", "Columns")]
        names(out) <- c("Type", "Size", "Rows", "Columns")
    if (head)
        out <- head(out, n)
    out
}

