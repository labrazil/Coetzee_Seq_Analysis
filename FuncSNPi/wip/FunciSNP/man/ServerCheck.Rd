\name{ServerCheck}
\alias{ServerCheck}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ServerCheck(primary.server)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{primary.server}{
%%     ~~Describe \code{primary.server} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (primary.server) 
{
    ncbi <- "ftp://ftp-trace.ncbi.nih.gov/1000genomes/"
    ebi <- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/"
    test.file <- "ftp/release/20110521/phase1_integrated_calls.20101123.ALL.panel"
    if (primary.server == "ebi") {
        primary.server <- ebi
        secondary.server <- ncbi
        primary.server.name <- "ebi"
        secondary.server.name <- "ncbi"
    }
    else {
        primary.server <- ncbi
        secondary.server <- ebi
        primary.server.name <- "ncbi"
        secondary.server.name <- "ebi"
    }
    cat("trying", primary.server.name, "as 1000 genomes server\n")
    server.up <- try(url(paste(primary.server, test.file, sep = ""), 
        open = "rt"), silent = TRUE)
    server.error <- inherits(server.up, "try-error")
    if (server.error) {
        cat(primary.server.name, "failed \ntrying", secondary.server.name, 
            "as 1000 genomes server\n")
        rm(server.error)
        server.up <- try(url(paste(secondary.server, test.file, 
            sep = ""), open = "rt"), silent = TRUE)
        server.error <- inherits(server.up, "try-error")
        if (server.error) {
            cat(secondary.server.name, "failed \n")
            close(server.up)
            stop("Neither EBI nor NCBI mirrors for the 1000 genomes project were found")
        }
        else {
            cat("OK using", secondary.server.name, ":", secondary.server, 
                "\n")
            close(server.up)
            secondary.server
        }
    }
    else {
        cat("OK using", primary.server.name, ":", primary.server, 
            "\n")
        close(server.up)
        primary.server
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
