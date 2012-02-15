\name{TSList-class}
\Rdversion{1.1}
\docType{class}
\alias{TSList-class}
\alias{show,TSList-method}
\alias{summary,TSList-method}

\title{Class \code{"TSList"}}
\description{
ffff%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("TSList", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{snp.data}:}{Object of class \code{"list"} ~~ }
    \item{\code{summary.data}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{elementType}:}{Object of class \code{"character"} ~~ }
    \item{\code{elementMetadata}:}{Object of class \code{"DataTableORNULL"} ~~ }
    \item{\code{metadata}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Extends}{
%Class \code{"\linkS4class{List}"}, directly.
%Class \code{"\linkS4class{Vector}"}, by class "List", distance 2.
%Class \code{"\linkS4class{Annotated}"}, by class "List", distance 3.
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "TSList")}: ... }
    \item{summary}{\code{signature(object = "TSList")}: ... }
	 }
}
\references{
    Coetzee SG et al. submitted for review. 2012
}
\author{
Simon G. Coetzee and Houtan Noushmehr
}
\note{
    NA
}


\seealso{
    \code{\link{getFSNPs}}, \code{\link{FunciSNPplot}}, \code{\link{FunciSNPAnnotateSummary}}, \code{\link{FunciSNPtable}}, \code{\link{FunciSNPbed}}
}
\examples{
showClass("TSList")
}
\keyword{classes}
