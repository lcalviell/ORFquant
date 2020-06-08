#' @import Rsamtools
NULL


DEFAULT_CIRC_SEQS <- unique(c("chrM","MT","MtDNA","mit","Mito","mitochondrion",
                              "dmel_mitochondrion_genome","Pltd","ChrC","Pt","chloroplast",
                              "Chloro","2micron","2-micron","2uM",
                              "Mt", "NC_001879.2", "NC_006581.1","ChrM"))

   
################################################################################
#' A simple extension to the FaFile class that
#' allows one to include a list of circular
#' ranges, e.g. chrM
#'
#' @field circularRanges A character vector describing which seqnames have circular ranges
#' @importFrom Rsamtools FaFile
#' @examples
#' mytempfile=tempfile()
#' writeXStringSet(setNames(DNAStringSet(c('AAAAAAAAGG','AAAAAAAAGG')),
#'   c('chrM','chr2')),filepath=mytempfile)
#' Rsamtools::indexFa(mytempfile)
#' cREF<-FaFile_Circ(Rsamtools::FaFile(mytempfile),circularRanges='chrM')
#' cREF
#' @export FaFile_Circ
#' @exportClass FaFile_Circ

FaFile_Circ<-setRefClass("FaFile_Circ",
      contains="FaFile",
      fields=representation(circularRanges="character"),
      inheritPackage=TRUE
)


#' Yields a seqinfo object for the FaFile_Circ with it's circular ranges slot
#' set appropriately
#'
#'
#' @param x FaFile_Circ; the object to get the seqinfo for
#'
#' @return A Seqinfo object
#'
#' @seealso \code{\link{create_html_report}}
#' @examples
#'
#' mytempfile=tempfile()
#' writeXStringSet(setNames(DNAStringSet(c('AAAAAAAAGG','AAAAAAAAGG')),
#'   c('chrM','chr2')),
#'   filepath=mytempfile)
#' Rsamtools::indexFa(mytempfile)
#' cREF<-FaFile_Circ(Rsamtools::FaFile(mytempfile),circularRanges='chrM')
#' seqinfo(cREF)
#' @export

setMethod(seqinfo,'FaFile_Circ',function (x){
   gr <- scanFaIndex(x$path)
   Seqinfo(as.character(seqnames(gr)), width(gr),
               isCircular = as.character(seqnames(gr)) %in% x$circularRanges )
})


#' Yields the sequence for a particular range on a circular Fasta File
#' note that
#'
#'
#' @param x FaFile_Circ; the object to get the seqinfo for
#'
#' @return A Seqinfo object
#'
#' @seealso \code{\link{create_html_report}}
#' @examples
#'
#' mytempfile=tempfile()
#' writeXStringSet(setNames(DNAStringSet(c('AAAAAAAAGG','AAAAAAAAGG')),
#'   c('chrM','chr2')),filepath=mytempfile)
#' Rsamtools::indexFa(mytempfile)
#' cREF<-FaFile_Circ(Rsamtools::FaFile(mytempfile),circularRanges='chrM')
#' seqinfo(cREF)
#' @export


setMethod('getSeq','FaFile_Circ', function (x, ...){
   dots<-list(...)
   if(length(dots)==1 && is(dots[[1]],'GRanges')){
      names(dots)[1] <- 'param'
   }
   .local <- function (x, param, ...){
      if (missing(param)) {
         scanFa(x$path, ...)
      }
      else {
         ################################################################################

         if (is(param, "GRanges")) {
            stopifnot(all(unique(seqnames(param)) %in% seqinfo(x)@seqnames))
            seqinfo(param) <- seqinfo(x)[seqlevels(param)]
            idx <- as.logical(strand(param) == "-")
            chrends <- seqlengths(param)[as.character(param@seqnames)]
            is_wrapping <- end(param) > chrends
            is_wrapping <- vapply(is_wrapping,isTRUE,TRUE)
            wrapchrs<-as.character(seqnames(param[is_wrapping]))
            circs_wrap <- isCircular(param)[wrapchrs]

            if(!all( vapply(circs_wrap,isTRUE,TRUE) )) {
               stop('Out of bounds ranges detected which are ',
                      'not on circular chromosomes')
            }
            #
            if(any(vapply(is_wrapping,isTRUE,TRUE))) {
               wrapped_grs <- restrict(param[is_wrapping],chrends+1)
               #get ends for the wrapped ranges
               wrappedchrs<-as.character(wrapped_grs@seqnames)
               chrends_wrap <- seqlengths(wrapped_grs)[wrappedchrs]
               #shift these wrapped ranges back to 1
               wrapped_grs <- shift(wrapped_grs,-chrends_wrap)
               #make sure they don't wrap twice
               if(any(end(wrapped_grs) > chrends_wrap)) stop("Ranges wrapping twice',
                                                                         ' isn't implemented yet...")
               #also get the within bounds ranges
               nonwrapped_grs <- restrict(param,end=chrends)
               #scan seperately
               dna <- scanFa(x$path, nonwrapped_grs,...)
               wrap_dna <- scanFa(x$path, wrapped_grs,...)

               #append the positive wraps to the end of the '+' rnages
               if(any(is_wrapping[!idx])){
                  isposwrap<- is_wrapping & (!idx)
                  posthatwrap<-!idx[is_wrapping]
                  dna[isposwrap] <- xscat( dna[isposwrap],
                                                       wrap_dna[posthatwrap]
                  )
               }
               #reverse our negative ranges nad append their wraps to the start
               if(any(idx)){
                  dna[idx] <- reverseComplement(dna[idx])
                  isnegwrap<-is_wrapping & (idx)
                  if(any(is_wrapping[idx])){
                     negsthatwrap<-idx[is_wrapping]
                     wrap_dna[negsthatwrap] <- reverseComplement(
                        wrap_dna[negsthatwrap]
                     )
                     dna[isnegwrap] <- xscat(
                        wrap_dna[negsthatwrap],
                        dna[isnegwrap]
                     )
                  }
               }
            }else{
               dna <- scanFa(x$path, param, ...)
               idx <- as.logical(strand(param) == "-")
               if (any(idx)) dna[idx] <- reverseComplement(dna[idx])
            }
         }else{
            dna <- scanFa(x$path, param, ...)
         }
         dna
      }
   }
   do.call(.local,c(list(x=x), dots))
})




