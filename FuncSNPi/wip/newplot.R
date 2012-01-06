LDtableUPdate <- function (x, colorcut = c(0, 0.01, 0.025, 0.5, 0.1, 1), colors = heat.colors(length(colorcut)), 
    textcol = "black", digits = 3, show.all = FALSE, which = c("D", 
        "D'", "r", "X^2", "P-value", "n"), colorize = "P-value", 
    cex, ...) 
{
    if (!colorize %in% names(x)) 
        stop(colorize, " not an element of ", deparse(substitute(x)))
    datatab <- summary(x)
    missmatch <- which[!(which %in% names(x))]
    if (length(missmatch) > 0) 
        stop(missmatch, " not an element of ", deparse(substitute(x)))
    matform <- function(value, template) {
        dim(value) <- dim(template)
        dimnames(value) <- dimnames(template)
        value
    }
    tmp <- cut(x[[colorize]], colorcut, include.lowest = TRUE)
    colormat <- matform(as.numeric(tmp), x[[colorize]])
    n <- matform(paste("(", x$n, ")", sep = ""), x$n)
    if (!show.all) {
        colormat <- colormat[-nrow(colormat), -1, drop = FALSE]
        n <- n[-nrow(n), -1, drop = FALSE]
    }
    image(x = 1:ncol(colormat), y = 1:ncol(colormat), z = t(colormat[nrow(colormat):1, 
        ]), col = colors, xlab = "", ylab = "", 
        xaxt = "n", yaxt = "n", ...)
    abline(v = -0.5 + 1:(ncol(colormat) + 1))
    abline(h = -0.5 + 1:(nrow(colormat) + 1))
    axis(3, 1:ncol(colormat), colnames(colormat), las = 1, cex.axis=0.4)
    axis(2, 1:nrow(colormat), rev(rownames(colormat)), las = 1, cex.axis=0.4)
    cex.old <- par("cex")
    if (missing(cex)) 
        cex <- min(c(1/10, 1/(length(which) + 1))/c(strwidth("W"), 
            strheight("W") * 1.5))
    par(cex = cex)
    lineheight <- strheight("W") * 1.5
    center <- lineheight * length(which)/2
    for (i in 1:length(which)) {
        displaymat <- x[[which[i]]]
        if (!show.all) 
            displaymat <- displaymat[-nrow(displaymat), -1, drop = FALSE]
        if (which[i] == "P-value") 
            displaymat <- format.pval(displaymat, digits = digits)
        else if (which[i] != "n") 
            displaymat <- format(displaymat, digits = digits)
        displaymat[] <- gsub("NA.*", "", as.character(displaymat))
        text(x = col(colormat), y = nrow(colormat) - row(colormat) + 
            1 + center - lineheight * (i - 1), displaymat, col = textcol, 
            adj = c(0.5, 1))
    }
    text(x = 1, y = 1, paste(which, collapse = "\n"), adj = c(0.5, 
        0.5))
    par(cex = cex.old)
    title(main = "Linkage Disequilibrium\n")
    invisible(colormat)
}
