% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orfquant.R
\name{annotate_splicing}
\alias{annotate_splicing}
\title{Annotate splice features of detected ORFs}
\usage{
annotate_splicing(orf_gen, ref_cds)
}
\arguments{
\item{orf_gen}{Exon structure of a detected ORF}

\item{ref_cds}{Exon structure of a reference ORF}
}
\value{
Exon structure of detected ORF including possible missing exons from reference, together with a \code{spl_type} column
including the annotation for each exon (e.g. alternative acceptors or donor).
}
\description{
This function detects usage of different exons and exonic boundaries of one ORF with respect to a reference ORF.
}
\details{
each exon is aligned to the closest one to match acceptor and donor sites, or to annotate missing exons.
\code{5ss} and \code{3ss} indicate exon 5' and 3', respectively. \code{CDS_spanning} indicates retained intron;
\code{missing_CDS} indicates no overlapping exon (missed or included); \code{monoCDS} indicates a single-exon ORF; 
\code{firstCDS} and \code{lastCDS} indicate first CDS exon or last CDS exon.
}
\seealso{
\code{\link{detect_translated_orfs}}, \code{\link{annotate_ORFs}}
}
\author{
Lorenzo Calviello, \email{calviello.l.bio@gmail.com}
}
\keyword{ORFquant}
