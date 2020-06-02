# ORFquant, a software to perform splice-aware 
# quantification of ORF translation using Ribo-seq data
#
# Authors: 
# Lorenzo Calviello (calviello.l.bio@gmail.com)
# Uwe Ohler (Uwe.Ohler@mdc-berlin.de)
#
# This software is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This software is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this software. If not, see
# <http://www.gnu.org/licenses/>.

setMethods(f = "unlist","GRanges",function(x){return(x)})


DEFAULT_CIRC_SEQS <- unique(c("chrM","MT","MtDNA","mit","Mito","mitochondrion",
                              "dmel_mitochondrion_genome","Pltd","ChrC","Pt","chloroplast",
                              "Chloro","2micron","2-micron","2uM",
                              "Mt", "NC_001879.2", "NC_006581.1","ChrM"))


#' Find ATG-starting ORFs in a sequence
#'
#' This function loads the annotation created by the \code{prepare_annotation_files function}
#' @keywords ORFquant
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param tx_name transcript_id
#' @param sequence DNAString object containing the sequence of the transcript
#' @param get_all_starts Output all possible start codons? Defaults to \code{TRUE}
#' @param Stop_Stop Find Stop-Stop pairs (no defined start codon)? Defaults to \code{FALSE}
#' @param scores Deprecated
#' @param genetic_code_table GENETIC_CODE table to use
#' @return \code{GRanges} object containing coordinates for the detected ORFs
#' @seealso \code{\link{detect_translated_orfs}}
#' @export


get_orfs<-function(tx_name,sequence,get_all_starts=T,Stop_Stop=F,scores=c(1,.5),genetic_code_table){
    list_frames<-list()
    length<-nchar(sequence)
    for(u in 0:2){
        pept<-NA
        pept<-unlist(strsplit(as.character(suppressWarnings(translate(subseq(sequence,start=u+1),genetic.code = genetic_code_table,if.fuzzy.codon = "solve"))),split=""))
        
        
        starts<-pept=="M"
        
        stops<-pept=="*"
        
        start_pos<-((1:length(pept))[starts])*3
        if(length(start_pos)>0){
            start_pos<-start_pos+u-2
        } else {start_pos<-NA}
        
        stop_pos<-((1:length(pept))[stops])*3-3
        if(length(stop_pos)>0){
            stop_pos<-stop_pos+u
        } else {stop_pos<-NA}
        
        st2vect<-c()
        for(h in 1:length(start_pos)){
            st1<-start_pos[h]
            diff<-stop_pos-st1
            diff<-diff[diff>0]
            if(length(diff)>0){st2<-st1+min(diff)}
            if(length(diff)==0){st2<-NA}
            st2vect[h]<-st2
            
        }
        st_st<-data.frame(cbind(start_pos,st2vect),stringsAsFactors=F)
        st_st<-st_st[!is.na(st_st[,1]),]
        st_st<-st_st[!is.na(st_st[,2]),]
        if(dim(st_st)[1]==0){list_frames[[paste("frame",u,sep = "_")]]<-GRanges(); next}
        if(get_all_starts==T){
            gra_orf<-GRanges(seqnames = paste(tx_name,"frame",u,sep = "_"),strand = "+",ranges = IRanges(start = st_st[,1],end = st_st[,2]))
        }
        if(get_all_starts==F){
            gra_orf<-reduce(GRanges(seqnames = paste(tx_name,"frame",u,sep = "_"),strand = "*",ranges = IRanges(start = st_st[,1],end = st_st[,2])))
        }
        if(length(gra_orf)>0){
            gra_orf$type="ORF"
            gra_orf$score<-scores[1]
        }
        
        gra_par<-GRanges()
        if(Stop_Stop==T){
            stops<-pept=="*"
            stop_pos<-((1:length(pept))[stops])*3
            if(length(stop_pos)>0){
                stop_pos<-stop_pos+u
            } else {stop_pos<-NA}
            
            
            vector_starts<-c()
            vector_stops<-c()
            
            if(sum(stops)==0){vector_starts<-u+1; vector_stops<-length}
            
            if(sum(stops)==1){
                vector_starts[1]<-u+1
                vector_starts[2]<-stop_pos+1
                vector_stops[1]<-stop_pos-3
                vector_stops[2]<-length
                
            }
            
            if(sum(stops)>0){
                for(stoppo in 0:(length(stop_pos))){
                    start<-stop_pos[stoppo]+1
                    if(stoppo==0){
                        start<-u+1
                    }
                    end<-stop_pos[stoppo+1]-3
                    if(stoppo==length(stop_pos)){
                        end<-length
                    }
                    
                    if((end-start)<5){
                        next
                    }
                    vector_starts[stoppo+1]<-start
                    vector_stops[stoppo+1]<-end
                    
                }
                
                
            }
            
            st_st<-data.frame(cbind(vector_starts,vector_stops))
            st_st<-st_st[!is.na(st_st[,"vector_starts"]),]
            st_st<-st_st[!is.na(st_st[,"vector_stops"]),]
            #st_st<-st_st[(st_st[,2]-st_st[,1]>11),]
            gra_par<-GRanges(seqnames = paste(tx_name,"frame",u,sep = "_"),strand = "+",ranges = IRanges(start = st_st[,1],end = st_st[,2]))
            if(length(gra_par)>0){
                gra_par<-setdiff(gra_par,gra_orf)
                
            }
            if(length(gra_par)>0){
                gra_par$type="frame_ok"
                gra_par$score<-scores[2]                        
            }
        }
        gra<-c(gra_orf,gra_par)
        
        list_frames[[paste("frame",u,sep = "_")]]<-gra
    }
    GRangesList(list_frames)
}

#' Extract output from multitaper analysis of a signal
#'
#' This function uses the multitaper tool to extract F-values and multitaper spectral coefficients
#' @details Values reported correspond to the closest frequency to 1/3 (same parameters as in RiboTaper). \cr
#' Padding to a minimum length of 1024 is performed to increase spectral resolution.
#' @keywords ORFquant
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param x numeric signal to analyze
#' @param n_tapers n of tapers to use
#' @param time_bw time_bw parameter
#' @param slepians_values set of calculated slepian functions to use in the multitaper analysis
#' @return two numeric values representing the F-value for the multitaper test and its corresponding spectral coefficient at the closest frequency to 1/3
#' @seealso \code{\link{detect_translated_orfs}}, \code{\link{spec.mtm}}, \code{\link{dpss}}
#' @export

take_Fvals_spect<-function(x,n_tapers,time_bw,slepians_values){
    if(length(x)<25){
        remain<-50-length(x)
        x<-c(rep(0,as.integer(remain/2)),x,rep(0,remain%%2+as.integer(remain/2)))
    }
    if(length(x)<1024/2){padding<-1024}
    if(length(x)>=1024/2){padding<-"default"}
    resSpec1 <- spec.mtm(as.ts(x), k=n_tapers, nw=time_bw, nFFT = padding, centreWithSlepians = TRUE, Ftest = TRUE, maxAdaptiveIterations = 100,returnZeroFreq=F,plot=F,dpssIN=slepians_values)
    
    resSpec2<-dropFreqs(resSpec1,0.29,0.39)
    
    freq_max_3nt<-resSpec1$freq[which(abs((resSpec1$freq-(1/3)))==min(abs((resSpec1$freq-(1/3)))))]
    
    Fmax_3nt<-resSpec1$mtm$Ftest[which(abs((resSpec1$freq-(1/3)))==min(abs((resSpec1$freq-(1/3)))))]
    spect_3nt<-resSpec1$spec[which(abs((resSpec1$freq-(1/3)))==min(abs((resSpec1$freq-(1/3)))))]
    return(c(Fmax_3nt,spect_3nt))
    
}

#' Select start codon 
#'
#' This function selects the start codon for ORFs in the same transcript
#' @details ORFs are divided based on stop codon and Ribo-seq signal between start codons is used to select one.\cr
#' When more than \code{cutoff_ave} fraction of codons is in-frame between two candidate start codons, the most upstream
#' is selected.
#' @keywords ORFquant
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param ORFs Set of detected ORFs
#' @param P_sites_rle Rle signal of P_sites along the transcript
#' @param cutoff cutoff of total in-frame signal between start codons (sensitive to outliers). Defaults to NA
#' @param cutoff_ave cutoff for frequency of in-frame codons between two start codons (less sensitive to outliers). Defaults to .5
#' @return Set of detected ORFs, including info about the possible longest ORF for that frame.
#' @seealso \code{\link{detect_translated_orfs}}, \code{\link{get_orfs}}
#' @export

select_start<-function(ORFs,P_sites_rle,cutoff=NA,cutoff_ave=.5){
    ORFs<-ORFs[order(end(ORFs),start(ORFs),decreasing = F)]
    names(ORFs)<-NULL
    longest_ORF<-split(ORFs,end(ORFs))
    maxo<-which.max(width(longest_ORF))
    longest_ORF<-unlist(longest_ORF[splitAsList(unname(maxo), names(maxo))])
    P_sites_rle<-RleList(P_sites_rle)
    names(P_sites_rle)<-seqnames(ORFs[1])
    covss<-P_sites_rle[ORFs]
    ok<-sapply(covss,function(x){sum(as.vector(x)>0)>2})
    if(length(ok)==0){return(GRanges())}
    ORFs<-ORFs[ok]
    covss<-covss[ok]
    rrr<-lapply(covss,function(x){
        fr<-suppressWarnings(matrix(as.vector(x),nrow = 3))
        fr<-fr[,colSums(fr)>0,drop=F]
        if(dim(fr)[2]==0){return(GRanges())}
        fra<-apply(fr,2,function(y){y/sum(y)})
        infr_freq<-rowMeans(fra)[1]
        infr<-((rowSums(fr))[1])/sum(fr)
        y<-DataFrame(ave_pct_fr=round(infr_freq*100,digits = 4))
        y$pct_fr<-round(infr*100,digits = 4)
        y$ave_pct_fr_st<-NA
        y$pct_fr_st<-NA
        y
    })
    mcols(ORFs)<-do.call(rrr,what = rbind)
    
    if(!is.na(cutoff)){
        covss<-covss[ORFs$pct_fr>=cutoff]
        ORFs<-ORFs[ORFs$pct_fr>=cutoff]
    }
    if(!is.na(cutoff_ave)){
        covss<-covss[ORFs$ave_pct_fr>=cutoff_ave]
        ORFs<-ORFs[ORFs$ave_pct_fr>=cutoff_ave]
    }
    if(length(ORFs)==0){return(GRanges())}
    ORFs$endorf<-end(ORFs)
    orfs_spl<-split(ORFs,end(ORFs))
    names(covss)<-end(ORFs)
    ORFs<-endoapply(orfs_spl,function(x){
        if(length(x)==1){x$endorf<-NULL;return(x)}
        covo<-covss[names(covss)%in%as.character(x$endorf[1])]
        okorfa<-c()
        for(cnt in 1:length(x)){
            psit<-as.vector(covo[[cnt]])
            
            if(cnt<length(x)){
                psit<-psit[1:(start(x)[cnt+1]-start(x)[cnt])]
            }
            if(cnt==length(x)){okorfa<-cnt}
            
            frames<-suppressWarnings(matrix(as.vector(psit),nrow = 3))
            frames<-frames[,colSums(frames)>0,drop=F]
            if(dim(frames)[2]==0){next}
            frames<-apply(frames,2,function(y){y/sum(y)})
            infr_freq<-rowMeans(frames)[1]
            infr<-sum(psit[seq(1,length(psit),by=3)])/sum(psit)
            if(!is.na(cutoff)){
                if(infr>=cutoff & !is.na(infr)){
                    x$ave_pct_fr_st[cnt]<-round(infr_freq,digits = 4)
                    x$pct_fr_st[cnt]<-round(infr,digits = 4)
                    okorfa<-cnt
                    break
                }
            }
            
            if(!is.na(cutoff_ave)){
                if(infr_freq>=cutoff_ave & !is.na(infr_freq)){
                    x$ave_pct_fr_st[cnt]<-round(infr_freq*100,digits = 4)
                    x$pct_fr_st[cnt]<-round(infr*100,digits = 4)
                    okorfa<-cnt
                    break
                }
                
            }
            
            
        }
        x<-x[okorfa]
        x$endorf<-NULL
        return(x)
    })
    ORFs<-unlist(ORFs)
    names(ORFs)<-NULL
    ORFs$longest_ORF<-longest_ORF[match(end(ORFs),names(longest_ORF))]
    return(ORFs)
    
}

#' Collect ORF Ribo-seq statistics
#'
#' This function calculates statistics for the analysis of P_sites profiles for each ORF
#' @details Number of P_sites (uniquely mapping or all), frame percentage and multitaper test 
#' statistics are collected for each ORF. The parameter space for the multitaper analysis was explored in the RiboTaper paper.
#' @keywords ORFquant
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param ORFs Set of detected ORFs
#' @param P_sites_rle Rle signal of P_sites along the transcript
#' @param P_sites_uniq_rle Rle signal of uniquely mapping P_sites along the transcript
#' @param P_sites_uniq_mm_rle Rle signal of uniquely mapping P_sites with mismatches along the transcript
#' @param cutoff cutoff of average in-frame signal for each codon in the ORF. Defaults to .5
#' @param tapers Number of tapers to use in the multitaper analysis. Defaults to 24
#' @param bw time_bw parameter to use in the multitaper analysis. Defaults to 12
#' @return Set of detected ORFs, including info about the possible longest ORF for that frame.
#' @seealso \code{\link{detect_translated_orfs}}, \code{\link{get_orfs}}, \code{\link{take_Fvals_spect}}
#' @export

calc_orf_pval<-function(ORFs,P_sites_rle,P_sites_uniq_rle,P_sites_uniq_mm_rle,cutoff=.5,tapers=24,bw=12){
    ORFs$pval<-NA
    ORFs$pval_uniq<-NA
    ORFs$P_sites_raw<-NA
    ORFs$P_sites_raw_uniq<-NA
    ORFs$P_sites_raw_uniq_mm<-NA
    ORFs$pct_fr<-NA
    
    
    for(i in 1:length(ORFs)){
        psit<-as.vector(P_sites_rle[ORFs[i]@ranges])
        psit_uniq<-as.vector(P_sites_uniq_rle[ORFs[i]@ranges])
        psit_uniq_mm<-as.vector(P_sites_uniq_mm_rle[ORFs[i]@ranges])
        ORFs$P_sites_raw[i]<-sum(psit)
        ps_unq<-round(sum(psit_uniq)/sum(psit)*100,digits = 2)
        if(is.na(ps_unq)){ps_unq<-0}
        ORFs$P_sites_raw_uniq[i]<-sum(psit_uniq)
        
        ps_unq<-round((sum(psit_uniq)-sum(psit_uniq_mm))/sum(psit)*100,digits = 2)
        if(is.na(ps_unq)){ps_unq<-0}
        
        ORFs$P_sites_raw_uniq_mm[i]<-sum(psit_uniq_mm)
        ORFs$ORF_id_tr[i]<-paste(as.character(seqnames(ORFs[i])[1]),start(ORFs[i]),end(ORFs[i]),sep = "_")
        if(sum(psit)>0){
            infr<-round(sum(psit[seq(1,length(psit),by=3)])/sum(psit),digits = 4)
            ORFs$pct_fr[i]<-infr
        }
        if(sum(psit>0)>2){
            if(infr>cutoff){
                if(length(psit)<25){slepians<-dpss(n=length(psit)+(50-length(psit)),k=tapers,nw=bw)}
                if(length(psit)>=25){slepians<-dpss(n=length(psit),k=tapers,nw=bw)}
                vals<-take_Fvals_spect(x = psit,n_tapers = tapers,time_bw = bw,slepians_values = slepians)
                ORFs$pval[i]<-pf(q=vals[1],df1=2,df2=(2*24)-2,lower.tail=F)
                vals<-take_Fvals_spect(x = psit_uniq,n_tapers = tapers,time_bw = bw,slepians_values = slepians)
                ORFs$pval_uniq[i]<-pf(q=vals[1],df1=2,df2=(2*24)-2,lower.tail=F)
                
                
            }
        }
    }
    return(ORFs)
}


#' Detect actively translated ORFs 
#'
#' This function detects translated ORFs
#' @details A set of transcripts, together with genome sequence and Ribo-signal are analyzed to extract translated ORFs
#' @keywords ORFquant
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param selected_txs set of selected transcripts, output from \code{select_txs}
#' @param genome_sequence BSgenome object
#' @param annotation Rannot object containing annotation of CDS and transcript structures (see \code{prepare_annotation_files})
#' @param P_sites GRanges object with P_sites positions
#' @param P_sites_uniq GRanges object with uniquely mapping P_sites positions
#' @param P_sites_uniq_mm GRanges object with uniquely mapping (with mismatches) P_sites positions
#' @param genomic_region GRanges object with genomic coordinates of the genomic region analyzed
#' @param genetic_code GENETIC_CODE table to use
#' @param all_starts \code{get_all_starts} parameter for the \code{get_orfs} function
#' @param nostarts \code{Stop_Stop} parameter for the \code{get_orfs} function
#' @param start_sel_cutoff \code{cutoff} parameter for the \code{select_start} function
#' @param start_sel_cutoff_ave \code{cutoff_ave} parameter for the \code{select_start} function
#' @param cutoff_fr_ave \code{cutoff} parameter for the  \code{calc_orf_pval} functions
#' @param uniq_signal Use only signal from uniquely mapping reads? Defaults to \code{FALSE}.
#' @return A list with transcript coordinates, exonic coordinates and statistics for each ORF exonic bin and junction(from \code{select_txs}).\cr\cr
#' The value for each column is as follows:\cr\cr
#' \code{ave_pct_fr}: average percentage of in-frame reads for each codon in the ORF
#' \code{pct_fr}: percentage of in-frame reads in the ORF
#' \code{ave_pct_fr}: average percentage of in-frame reads for each codon in the ORF
#' \code{ave_pct_fr_st}: average percentage of in-frame reads per each codon between the selected start codon and the next candidate one
#' \code{pct_fr_st}: percentage of in-frame reads between the selected start codon and the next candidate one
#' \code{longest_ORF}: GRanges coordinates for the longest ORF with the same stop codon
#' \code{pval}: P-value for the multitaper F-test at 1/3 using the ORF P_sites profile
#' \code{pval_uniq}: P-value for the multitaper F-test at 1/3 using the ORF P_sites profile (only uniquely mapping reads)
#' \code{P_sites_raw}: Raw number of P_sites mapping to the ORF
#' \code{P_sites_raw_unique}: Uniquely mapping P_sites mapping to the ORF
#' \code{ORF_id_tr}: ORF id containing <tx_id>_<start>_<end>
#' \code{Protein}: AAString sequence of the translated protein
#' \code{region}: Genomic coordinates of the analyzed region
#' \code{gene_id}: gene_id for the corresponding analyzed transcript
#' \code{gene_biotype}: gene biotype for the corresponding analyzed transcript
#' \code{gene_name}: gene name for the corresponding analyzed transcript
#' \code{transcript_id}: transcript_id for the corresponding analyzed ORF
#' \code{transcript_biotype}: transcript biotype for the corresponding analyzed ORF
#' @seealso \code{\link{select_txs}}, \code{\link{get_orfs}}, \code{\link{take_Fvals_spect}}, \code{\link{select_start}}, \code{\link{prepare_annotation_files}}
#' @export

detect_translated_orfs<-function(selected_txs,genome_sequence,annotation,P_sites,P_sites_uniq,P_sites_uniq_mm,genomic_region,genetic_code,
                                 all_starts=T,nostarts=F,start_sel_cutoff=NA,start_sel_cutoff_ave=.5,
                                 cutoff_fr_ave=.5,uniq_signal=F){
    
    orfs_gr<-GRangesList()
    orfs_gen_gr<-GRangesList()
    annot_tx_cds_gr<-GRangesList()
    cdss<-annotation$cds_txs
    exss<-annotation$exons_txs
    intr_txs<-annotation$introns_txs
    tr_gen<-annotation$trann
    txs_sels<-unique(unlist(selected_txs$txs_selected))
    annot_sels<-annotation$exons_txs[txs_sels]
    mapp<-mapToTranscripts(P_sites,annot_sels)
    mapp$reads<-P_sites$score[mapp$xHits]
    lens_sels<-sum(width(annot_sels))
    seqm<-seqlengths(mapp)
    seqlengths(mapp)<-lens_sels[match(names(seqm),names(lens_sels))]
    cov_txs<-coverage(mapp,weight = mapp$reads)
    
    mapp<-mapToTranscripts(P_sites_uniq,annot_sels)
    if(length(P_sites_uniq)>0){
        mapp$reads<-P_sites_uniq$score[mapp$xHits]
        seqm<-seqlengths(mapp)
        seqlengths(mapp)<-lens_sels[match(names(seqm),names(lens_sels))]
        cov_uniq_txs<-coverage(mapp,weight = mapp$reads)
    }
    
    if(length(P_sites_uniq)==0){cov_uniq_txs<-coverage(mapp)}
    
    mapp<-mapToTranscripts(P_sites_uniq_mm,annot_sels)
    
    if(length(P_sites_uniq_mm)==0){cov_uniq_mm_txs<-coverage(mapp)}
    
    if(length(P_sites_uniq_mm)>0){
        mapp$reads<-P_sites_uniq_mm$score[mapp$xHits]
        seqm<-seqlengths(mapp)
        seqlengths(mapp)<-lens_sels[match(names(seqm),names(lens_sels))]
        cov_uniq_mm_txs<-coverage(mapp,weight = mapp$reads)
    }
    txs_seqs<-extractTranscriptSeqs(genome_sequence,annot_sels)
    
    for(tx in txs_sels){
        
        ex_txs<-exss[tx]
        ex_tx<-ex_txs[[1]]
        intr_tx<-intr_txs[[tx]]
        nm_cds<-which(names(cdss)==tx)
        #map cds in tx space
        
        
        covtx<-cov_txs[[tx]]
        covtx_uniq<-cov_uniq_txs[[tx]]
        covtx_uniq_mm<-cov_uniq_mm_txs[[tx]]
        
        seq_tx<-txs_seqs[[tx]]
        if(length(covtx)==0){next}
        orfs<-unlist(get_orfs(tx_name = tx,sequence = seq_tx,get_all_starts=all_starts,Stop_Stop = F,genetic_code_table=genetic_code))
        #no need to know frame for now
        #names(orfs)<-NULL
        if(length(orfs)==0){next}
        
        orfs<-GRanges(seqnames = tx,ranges = orfs@ranges,strand = strand(orfs))
        orfs<-select_start(ORFs = orfs,P_sites_rle = covtx,cutoff = start_sel_cutoff,cutoff_ave = start_sel_cutoff_ave)
        if(length(orfs)==0){next}
        orfs<-calc_orf_pval(ORFs = orfs,P_sites_rle = covtx,P_sites_uniq_rle = covtx_uniq,P_sites_uniq_mm_rle = covtx_uniq_mm,cutoff = cutoff_fr_ave)
        if(length(orfs)==0){next}
        #uniq_flag
        if(uniq_signal){
            orfs<-subset(orfs,pval_uniq<.05)
        }
        if(!uniq_signal){
            orfs<-subset(orfs,pval<.05)
        }
        if(length(orfs)==0){next}
        #orfs$gene_id<-mapIds(keys = tx,x = annot,column = "GENEID",keytype = "TXNAME")
        orfs$Protein<-AAStringSet(rep("NA",length(orfs)))
        for(h in 1:length(orfs)){
            orfs$Protein[h]<-AAStringSet(as.character(translate(seq_tx[orfs@ranges[h]],genetic.code = genetic_code,if.fuzzy.codon = "solve")))
        }
        #orfs$Protein<-AAStringSet(orfs$Protein)
        orfs_gen<-from_tx_togen(ORFs = orfs,exons = ex_txs,introns = intr_tx)
        #here I should annotate
        
        tr_gen_tx<-tr_gen[as.vector(match(seqnames(orfs)[1],tr_gen[,"transcript_id"])),]
        
        orfs$region<-genomic_region
        orfs$gene_id<-unique(as.character(tr_gen_tx[,"gene_id"]))
        orfs$gene_biotype<-unique(as.character(tr_gen_tx[,"gene_biotype"]))
        orfs$gene_name<-unique(as.character(tr_gen_tx[,"gene_name"]))
        orfs$transcript_id<-unique(as.character(tr_gen_tx[,"transcript_id"]))
        orfs$transcript_biotype<-unique(as.character(tr_gen_tx[,"transcript_biotype"]))
        
        #must add the other compatible txs, to avoid calculating same stuff
        for(w in 1:length(orfs)){
            orf<-orfs[w]
            nam<-orf$ORF_id_tr
            orf$compatible_with<-NA
            orfs_gr[[nam]]<-orf
            orfs_gen_gr[[nam]]<-orfs_gen[[nam]]
            
        }
        
    }
    if(length(orfs_gr)==0){
        return(list())
    }
    tx_orfs<-unique(sapply(strsplit(names(orfs_gr),split="_"),function(x){
        len<-length(x)
        x[-c(len-1,len)]
    }))
    a<-lapply(selected_txs$txs,function(x){x[x%in%tx_orfs]})
    selected_txs$txs_orfs<-CharacterList(a)
    
    check<-sapply(selected_txs$txs_orfs,FUN = length)
    use<-rep("shared",length(check))
    use[check==1]<-"unique"
    use[check==0]<-"absent"
    use[check>1]<-"shared"
    selected_txs$use_ORFs<-use
    
    orfs_unq_gr<-GRangesList()
    
    for(j in names(orfs_gr)){
        
        orf_gen<-orfs_gen_gr[[j]]
        orf_tx<-orfs_gr[[j]]
        featexs<-selected_txs[selected_txs$type=="E"]
        featjuns<-selected_txs[selected_txs$type=="J"]
        
        over_tx<-featexs[featexs%over%orf_gen]
        if(length(featjuns)>0){
            over_tx<-sort(c(over_tx,featjuns[which(featjuns%in%gaps(orf_gen))]))
        }
        
        a<-sapply(over_tx$txs,function(x){length(x[x%in%as.character(seqnames(orf_tx))])})
        over_tx<-over_tx[a>0]
        
        orfs_unq_gr[[j]]<-over_tx
    }
    
    if(length(orfs_gr)>1){
        
        ident_mat<-matrix(FALSE,nrow=length(orfs_gr),ncol=length(orfs_gr))
        for(i in names(orfs_gr)){
            orf<-orfs_gr[[i]]
            orf$compatible_with<-NA
            gen<-orfs_gen_gr[[i]]
            
            ident<-c()
            for(j in 1:length(orfs_gen_gr)){
                ident[j]<-identical(gen,orfs_gen_gr[[j]])
            }
            ident_mat[which(names(orfs_gr)==i),]<-ident
        }
        #diag(ident_mat)<-NA
        ident_mat[lower.tri(ident_mat,diag=T)]<-NA
        ide<-which(ident_mat,arr.ind=T)
        if(length(ide)>0){
            rems<-c()
            ide<-split(ide,ide[,1])
            for(i in 1:length(ide)){
                iden<-ide[[i]]
                iden<-iden[!iden==as.numeric(names(ide)[i])]
                if(sum(iden%in%rems)>0){next}
                ind<-names(orfs_gr)[as.numeric(names(ide)[i])]
                ch<-names(orfs_gr)[iden]
                rems<-unique(c(rems,iden))
                #new_cat<-paste(unlist(orfs_gr[ch])$category_tx,collapse=";")
                new_id<-paste(unlist(orfs_gr[ch])$ORF_id_tr,collapse=";")
                orfs_gr[[ind]]$compatible_with<-new_id
                #if(length(unique(unlist(orfs_gr[ch])$category_tx))==1){next}
                #orfs_gr[[ind]]$compatible_categories<-new_cat
                #orfs_gr[[ind]]$category_tx<-"multiple"
            }
            orfs_gr<-orfs_gr[-rems]
            orfs_gen_gr<-orfs_gen_gr[-rems]
            orfs_unq_gr<-orfs_unq_gr[-rems]
        }
    }
    list_res<-list(orfs_gr,orfs_gen_gr,orfs_unq_gr)
    names(list_res)<-c("ORFs_tx_position","ORFs_genomic_position","ORFs_features")
    return(list_res)
}


#' Map transcript coordinates to genomic coordinates
#'
#' This function uses the \code{mapFromTranscripts} function to switch between transcript
#' and genomic coordinates
#' @keywords ORFquant
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param ORFs Set of detected ORFs from the \code{calc_orf_pval} function
#' @param exons exonic regions of the analyzed transcripts, as a GRangesList object
#' @param introns intronic regions of the analyzed transcripts, as a GRangesList object
#' @return exonic coordinates for each ORF.
#' @seealso \code{\link{mapFromTranscripts}}
#' @export

from_tx_togen<-function(ORFs,exons,introns){
    strand(ORFs)<-rep("*",length(ORFs))
    orfs_gen<-mapFromTranscripts(x = ORFs,transcripts = exons,ignore.strand=F)
    strand(orfs_gen)<-strand(exons[[1]][1])
    list_ma<-GRangesList()
    for(i in 1:length(ORFs)){
        or<-orfs_gen[i]
        or<-setdiff(or,introns)
        list_ma[[ORFs$ORF_id_tr[i]]]<-or
    }
    return(list_ma)
}


#' Select a subset of transcripts with Ribo-seq data
#'
#' This function flattens all annotated transcript structures and uses Ribo-seq to select a subset of transcripts.
#' @details Features (bins and junctions) are divided into shared and unique features, and into with support and without
#' support (with  or without reads mapping). A set of logical rules filters out transcripts with internal features with no support
#' and no unique features with reads. More specific details can be found in the ORFquant manuscript.
#' @keywords ORFquant
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param region genomic region being analyzed
#' @param annotation Rannot object containing annotation of CDS and transcript structures (see \code{prepare_annotation_files})
#' @param P_sites GRanges object with P_sites positions
#' @param P_sites_uniq GRanges object with uniquely mapping P_sites positions
#' @param junction_counts GRanges object containing Ribo-seq counts on the set of annotated junctions
#' @param uniq_signal Use only signal from uniquely mapping reads? Defaults to \code{FALSE}.
#' @return GRanges object with the set of counts on each exonic bin and junctions, together with the list
#' of selected transcripts
#' @seealso \code{\link{prepare_annotation_files}}
#' @export

select_txs<-function(region,annotation,P_sites,P_sites_uniq,junction_counts,uniq_signal=F){
    
    stra<-as.character(region@strand)
    gene_feat <- junction_counts[junction_counts %over% region]
    
    nsns<-annotation$exons_bins
    sel_nsns<-nsns%over%region
    gen_nsns<-nsns[sel_nsns]
    genbin<-gen_nsns
    genbin$exonic_part<-NULL
    genbin$type<-"E"
    genbin$reads<-0
    genbin$unique_reads<-0
    
    d<-genbin$reads
    
    hts<-findOverlaps(genbin,P_sites,ignore.strand=F)
    hts<-cbind(queryHits(hts),P_sites[subjectHits(hts)]$score*width(P_sites[subjectHits(hts)]))
    if(length(hts)>0){
        hts<-aggregate(x = hts[,2],list(hts[,1]),FUN=sum)
        for(i in 1:dim(hts)[1]){
            d[hts[i,1]]<-hts[i,2]
        }
        genbin$reads<-d
        
    }
    d<-genbin$unique_reads
    hts<-findOverlaps(genbin,P_sites_uniq,ignore.strand=F)
    hts<-cbind(queryHits(hts),P_sites_uniq[subjectHits(hts)]$score*width(P_sites_uniq[subjectHits(hts)]))
    #IMPORTANT, HERE THERE WAS A IF NO UNIQ RETURN GRANGESLIST()
    if(length(hts)>0){
        hts<-aggregate(x = hts[,2],list(hts[,1]),FUN=sum)
        for(i in 1:dim(hts)[1]){
            d[hts[i,1]]<-hts[i,2]
        }
        genbin$unique_reads<-d
        
    }
    mcols(gene_feat)<-mcols(gene_feat)[,names(mcols(genbin))]
    
    gene_feat<-sort(c(gene_feat,genbin))
    gene_feat$txs<-gene_feat$tx_name
    gene_feat$tx_name<-NULL
    
    
    rang<-gene_feat
    a<-gene_feat$txs
    b<-gene_feat$gene_id
    c<-gene_feat$type
    
    #HERE uniq_flag
    if(uniq_signal){
        d<-gene_feat$unique_reads
    }
    if(!uniq_signal){
        d<-gene_feat$reads
    }
    
    
    if(sum(d)<4){
        return(GRangesList())
    }
    
    
    txs_gene<-unique(unlist(a))     
    
    
    gen_bins_junct<-gene_feat
    if(length(txs_gene)<2){
        gen_bins_junct$genes<-b
        gen_bins_junct$txs<-a
        gen_bins_junct$genes_selected<-b
        gen_bins_junct$txs_selected<-a
        gen_bins_junct$use<-"unique"
        final_ranges<-sort(gen_bins_junct)
        return(final_ranges)
    }
    
    
    
    #first round
    txs_sofar<-txs_gene
    
    mat<-matrix(data=0,nrow=length(d),ncol=length(txs_sofar))
    colnames(mat)<-txs_sofar
    for(i in 1:length(txs_sofar)){
        mat[,i]<-sapply(a,function(x){sum(x==txs_sofar[i])})
    }
    
    nest<-c()
    ident<-c()
    for(i in 1:dim(mat)[2]){
        yes<-which(mat[,i]==1)
        nesti<-c()
        for(j in (1:dim(mat)[2])[-i]){
            yesj<-which(mat[,j]==1)
            if(identical(yes,yesj)){ident<-c(ident,paste(sort(colnames(mat)[c(i,j)]),collapse=";"))}
            #added if length> otherwise txs with same structure are both deleted
            nesti<-c(nesti,sum(yesj%in%yes)==length(yes) & length(mat[,j])>length(yes))
        }
        
        if(sum(nesti)>0){nest<-c(nest,colnames(mat)[i])}
        
    }
    if(length(ident)>0){nest<-nest[!nest%in%unique(sapply(strsplit(ident,";"),"[[",1))]}
    txs_sofar<-txs_sofar[!txs_sofar%in%nest]
    
    change<-1
    
    while(change>0){
        
        mat<-matrix(data=0,nrow=length(d),ncol=length(txs_sofar))
        colnames(mat)<-txs_sofar
        for(i in 1:length(txs_sofar)){
            mat[,i]<-sapply(a,function(x){sum(x==txs_sofar[i])})
        }
        
        mat_orig<-mat
        d_count<-paste(1:length(d),d,sep="_")
        good<-d_count[which(d>0)]
        bad<-d_count[which(d==0)]
        
        txs_good<-c()
        expl_good<-c()
        for(i in 1:dim(mat)[2]){
            expl_good_old<-expl_good
            #if new good feature
            tx<-d_count[which(mat[,i]>0)]
            tx_good<-tx[which(tx%in%good)]
            tx_bad<-tx[which(tx%in%bad)]
            
            if(length(tx_good)==0){next}
            
            if(sum(!tx_good%in%expl_good_old)>0){
                
                tx_torem<-c()
                tx_toscreen<-which(colnames(mat)%in%txs_good)
                if(length(tx_toscreen)>0){
                    for(j in tx_toscreen){
                        tx_contr<-d_count[which(mat[,j]>0)]
                        tx_contr_good<-tx_contr[which(tx_contr%in%good)]
                        tx_contr_bad<-tx_contr[which(tx_contr%in%bad)]
                        if(sum(tx_contr_good%in%tx_good)==length(tx_contr_good)){
                            tx_torem<-c(tx_torem,colnames(mat)[j])
                            
                        }
                        
                        if(length(tx_torem)>0){
                            txs_good<-txs_good[!txs_good%in%tx_torem]
                            
                        }
                    }
                    
                    
                    
                }
                txs_good<-unique(c(txs_good,colnames(mat)[i]))
                expl_good<-unique(c(expl_good,tx_good))
            }
            #if same good feature, but fewer bad INTERNAL features than others.
            if(sum(!tx_good%in%expl_good_old)==0){
                tx_torem<-c()
                tx_toscreen<-which(colnames(mat)%in%txs_good)
                #here a counter when at least one good feature more than competing
                moref<-c()
                for(j in tx_toscreen){
                    tx_contr<-d_count[which(mat[,j]>0)]
                    tx_contr_good<-tx_contr[which(tx_contr%in%good)]
                    tx_contr_bad<-tx_contr[which(tx_contr%in%bad)]
                    moref<-c(moref,sum(!tx_good%in%tx_contr_good)>0)
                    if(sum(tx_contr_good%in%tx_good)==length(tx_contr_good)){
                        
                        if(length(tx_good)>length(tx_contr_good)){
                            tx_torem<-c(tx_torem,colnames(mat)[j])
                        }
                        #here
                        if(length(tx_good)==length(tx_contr_good)){      
                            fi<-which(tx==tx_good[1])
                            la<-which(tx==tx_good[length(tx_good)])
                            int_tx<-tx[fi:la]
                            int_tx_bad<-tx_bad[tx_bad%in%int_tx]
                            
                            contr_fi<-which(tx_contr==tx_contr_good[1])
                            contr_la<-which(tx_contr==tx_contr_good[length(tx_contr_good)])
                            contr_int_tx<-tx_contr[contr_fi:contr_la]
                            contr_int_tx_bad<-tx_contr_bad[tx_contr_bad%in%contr_int_tx]
                            #SAME INTERNAL, TAKE
                            if(length(int_tx_bad)<=length(contr_int_tx_bad)){
                                
                                txs_good<-unique(c(txs_good,colnames(mat)[i]))
                                expl_good<-unique(c(expl_good,tx_good))
                                #LESS INTERNAL, TAKE AND REMOVE OTHER
                                if(length(int_tx_bad)<length(contr_int_tx_bad)){
                                    
                                    tx_torem<-c(tx_torem,colnames(mat)[j])}
                            }
                            
                            
                        }
                    }
                    
                    
                }
                if(length(tx_torem)>0){
                    txs_good<-txs_good[!txs_good%in%tx_torem]
                    txs_good<-unique(c(txs_good,colnames(mat)[i]))
                    expl_good<-unique(c(expl_good,tx_good))
                }
                
                if(length(tx_torem)==0 & sum(moref)==length(tx_toscreen)){
                    txs_good<-unique(c(txs_good,colnames(mat)[i]))
                    expl_good<-unique(c(expl_good,tx_good))
                }
                
                
            }
            
            
        }
        txs_sofar<-txs_good
        change<-abs(length(txs_good)-length(txs_sofar))
    }
    
    gen_bins_junct$genes<-b
    gen_bins_junct$txs<-a
    a<-lapply(a,function(x){x[x%in%txs_sofar]})
    
    genes_sofar<-unique(subset(annotation$trann,annotation$trann$transcript_id%in%txs_sofar)$gene_id)
    
    b<-lapply(b,function(x){x[x%in%genes_sofar]})
    gen_bins_junct$genes_selected<-CharacterList(b)
    gen_bins_junct$txs_selected<-CharacterList(a)
    
    check<-sapply(gen_bins_junct$txs,FUN = length)
    use<-rep("shared",length(check))
    use[check==1]<-"unique"
    use[check==0]<-"absent"
    use[check>1]<-"shared"
    
    gen_bins_junct$use<-use
    final_ranges<-sort(gen_bins_junct)
    return(final_ranges)
}


#' Extract possible readthrough sequences (beta) 
#'
#' This function extracts readthrough regions for subsequent analysis
#' @details The function looks for stop-stop pairs after the stop codon of the detected ORF
#' @keywords ORFquant
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param tx_name transcript_id
#' @param sequence DNAString object containing the sequence of the transcript
#' @param orf transcript-level ORF coordinates 
#' @param genetic_code GENETIC_CODE table to use
#' @return GRanges object with the set of possible readthrough sequences
#' @seealso \code{\link{detect_translated_orfs}}, \code{\link{select_quantify_ORFs}}
#' @export

get_reathr_seq<-function(tx_name,orf,sequence,genetic_code){
    length<-nchar(sequence)
    u=start(orf)%%3-1
    if(u==-1){u=2}
    pept<-NA
    pept<-unlist(strsplit(as.character(suppressWarnings(translate(subseq(sequence,start=u+1),genetic.code = genetic_code,if.fuzzy.codon = "solve"))),split=""))
    
    
    starts<-pept=="M"
    
    stops<-pept=="*"
    
    start_pos<-((1:length(pept))[starts])*3
    if(length(start_pos)>0){
        start_pos<-start_pos+u-2
    } else {start_pos<-NA}
    
    stop_pos<-((1:length(pept))[stops])*3-3
    if(length(stop_pos)>0){
        #here was -2 at the end, I'll -1, as the GRanges puts 1 more ???
        stop_pos<-stop_pos+u
    } else {stop_pos<-NA}
    stop_pos<-stop_pos[stop_pos>=end(orf)]
    
    vector_starts<-c()
    vector_stops<-c()
    
    if(sum(stops)==0){vector_starts<-u+1; vector_stops<-length}
    
    if(sum(stops)==1){
        vector_starts[1]<-u+1
        vector_starts[2]<-stop_pos+1
        vector_stops[1]<-stop_pos-3
        vector_stops[2]<-length
        
    }
    
    if(sum(stops)>0){
        for(stoppo in 0:(length(stop_pos))){
            start<-stop_pos[stoppo]+1
            if(stoppo==0){
                start<-u+1
            }
            end<-stop_pos[stoppo+1]-3
            if(stoppo==length(stop_pos)){
                end<-length
            }
            
            vector_starts[stoppo+1]<-start
            vector_stops[stoppo+1]<-end
            
        }
        
        
    }
    
    st_st<-data.frame(cbind(vector_starts,vector_stops))
    st_st<-st_st[!is.na(st_st[,"vector_starts"]),]
    st_st<-st_st[!is.na(st_st[,"vector_stops"]),]
    #st_st<-st_st[(st_st[,2]-st_st[,1]>11),]
    gra_par<-GRanges(seqnames = paste(tx_name,sep = "_"),strand = "+",ranges = IRanges(start = st_st[,1],end = st_st[,2]))
    gra_par<-gra_par[end(gra_par)>=end(orf)]
    return(gra_par)
}

#' Analyzed translation on possible readthrough regions (beta)
#'
#' This function uses the multitaper method to look for readthrough translation
#' @details The function looks for stop-stop pairs after the stop codon of the detected ORF
#' @keywords ORFquant
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param results_orf Full list of detected ORFs, from \code{select_quantify_ORFs} and \code{annotate_ORFs}
#' @param genome_sequence BSgenome object
#' @param annotation Rannot object containing annotation of CDS and transcript structures (see \code{prepare_annotation_files})
#' @param P_sites GRanges object with P_sites positions
#' @param P_sites_uniq GRanges object with uniquely mapping P_sites positions
#' @param P_sites_uniq_mm Rle signal of uniquely mapping P_sites with mismatches along the transcript
#' @param cutoff_fr_ave \code{cutoff} parameter for the  \code{calc_orf_pval} functions
#' @param genetic_code_table GENETIC_CODE table to use
#' @param uniq_signal Use only signal from uniquely mapping reads? Defaults to \code{FALSE}.
#' @return GRanges object with the set of translated readthrough regions
#' @seealso \code{\link{detect_translated_orfs}}, \code{\link{select_quantify_ORFs}}, \code{\link{annotate_ORFs}}, \code{\link{get_reathr_seq}}
#' @export

detect_readthrough<-function(results_orf,P_sites,P_sites_uniq,P_sites_uniq_mm,genome_sequence,annotation,genetic_code_table,cutoff_fr_ave=.5,uniq_signal=F){
    P_sites<-P_sites[!P_sites%over%unlist(results_orf$ORFs_genomic_position)]
    P_sites_uniq<-P_sites_uniq[!P_sites_uniq%over%unlist(results_orf$ORFs_genomic_position)]
    P_sites_uniq_mm<-P_sites_uniq_mm[!P_sites_uniq_mm%over%unlist(results_orf$ORFs_genomic_position)]
    readthroughs<-GRanges()
    if(length(P_sites)>4){
        for(i in 1:length(results_orf$ORFs_tx_position)){
            orf_tx<-results_orf$ORFs_tx_position[[i]]
            tx<-orf_tx$transcript_id
            ex_tx<-annotation$exons_txs[[tx]]
            
            #tx_gr<-transcripts_ranges[tx_ok]
            if(as.character(strand(ex_tx)[1])=="+"){
                covtx<-suppressWarnings(unlist(coverage(P_sites,weight = P_sites$score)[ex_tx]))
                if(length(P_sites_uniq)>0){
                    covtx_uniq<-suppressWarnings(unlist(coverage(P_sites_uniq,weight = P_sites_uniq$score)[ex_tx]))
                }
                if(length(P_sites_uniq)==0){
                    covtx_uniq<-suppressWarnings(unlist(coverage(P_sites_uniq)[ex_tx]))
                }
                if(length(P_sites_uniq_mm)>0){
                    covtx_uniq_mm<-suppressWarnings(unlist(coverage(P_sites_uniq_mm,weight = P_sites_uniq_mm$score)[ex_tx]))
                }
                if(length(P_sites_uniq_mm)==0){
                    covtx_uniq_mm<-suppressWarnings(unlist(coverage(P_sites_uniq_mm)[ex_tx]))
                }
                
            }
            if(as.character(strand(ex_tx)[1])=="-"){
                covtx<-suppressWarnings(unlist(unlist(RleList(lapply(coverage(P_sites,weight = P_sites$score)[ex_tx],FUN=rev)))))
                if(length(P_sites_uniq)>0){
                    covtx_uniq<-suppressWarnings(unlist(unlist(RleList(lapply(coverage(P_sites_uniq,weight = P_sites_uniq$score)[ex_tx],FUN=rev)))))
                }
                if(length(P_sites_uniq)==0){
                    covtx_uniq<-suppressWarnings(unlist(unlist(RleList(lapply(coverage(P_sites_uniq)[ex_tx],FUN=rev)))))
                }
                if(length(P_sites_uniq_mm)>0){
                    covtx_uniq_mm<-suppressWarnings(unlist(unlist(RleList(lapply(coverage(P_sites_uniq_mm,weight = P_sites_uniq_mm$score)[ex_tx],FUN=rev)))))
                }
                if(length(P_sites_uniq_mm)==0){
                    covtx_uniq_mm<-suppressWarnings(unlist(unlist(RleList(lapply(coverage(P_sites_uniq_mm)[ex_tx],FUN=rev)))))
                }
                
            }
            seq_tx<-unlist(getSeq(x=genome_sequence,ex_tx))
            
            orfs<-get_reathr_seq(tx_name = tx,orf = orf_tx,sequence = seq_tx,genetic_code = genetic_code_table)
            
            if(length(orfs)==0){next}
            vals1<-calc_orf_pval(ORFs = orfs[1],P_sites_rle = covtx,P_sites_uniq_rle = covtx_uniq,P_sites_uniq_mm_rle = covtx_uniq_mm,cutoff = cutoff_fr_ave)
            #uniq_flag
            if(!uniq_signal){
                if(is.na(vals1$pval) | vals1$pct_fr<.5 | vals1$pval>.05 ){next}
            }
            if(uniq_signal){
                if(is.na(vals1$pval_uniq) | vals1$pct_fr<.5 | vals1$pval_uniq>.05 ){next}
            }
            vals1$Protein<-AAStringSet(as.character(translate(seq_tx[vals1@ranges],genetic.code = genetic_code_table,if.fuzzy.codon = "solve")))
            vals1$ORF_orig_tr<-orf_tx$ORF_id_tr
            vals1$n_stops_readth<-1
            vals1$compatible_id<-CharacterList("")
            vals1$compatible_original_id<-CharacterList("")
            mcsva<-mcols(vals1)
            mcsva[,c("gene_id","gene_biotype","gene_name","transcript_id","transcript_biotype")]<-""
            mcols(vals1)<-mcsva
            vals1$gen_coords<-GRangesList(GRanges())
            
            if(length(orfs)>1){
                keep<-1
                while(keep>0){
                    if(length(orfs)<(keep+1)){break}
                    vals<-calc_orf_pval(ORFs = orfs[keep+1],P_sites_rle = covtx,P_sites_uniq_rle = covtx_uniq,P_sites_uniq_mm_rle = covtx_uniq_mm,cutoff = cutoff_fr_ave)
                    keep<-keep+1
                    #uniq_flag
                    
                    if(!uniq_signal){
                        if(is.na(vals1$pval) | vals1$pct_fr<.5 | vals1$pval>.05 ){keep<-0}
                        if(!is.na(vals$pval) & vals$pct_fr>.5 & vals$pval<.05 ){
                            end(vals1)<-end(vals)
                            vals1<-calc_orf_pval(ORFs = vals1,P_sites_rle = covtx,P_sites_uniq_rle = covtx_uniq,P_sites_uniq_mm_rle = covtx_uniq_mm,cutoff = cutoff_fr_ave)
                            vals1$Protein<-AAStringSet(as.character(translate(seq_tx[vals1@ranges],genetic.code = genetic_code_table,if.fuzzy.codon = "solve")))
                            vals1$ORF_orig_tr<-orf_tx$ORF_id_tr
                            vals1$n_stops_readth<-keep
                            vals1$compatible_id<-CharacterList("")
                            vals1$compatible_original_id<-CharacterList("")
                            mcsva<-mcols(vals1)
                            mcsva[,c("gene_id","gene_biotype","gene_name","transcript_id","transcript_biotype")]<-""
                            mcols(vals1)<-mcsva
                            
                        }
                    }
                    if(uniq_signal){
                        if(is.na(vals1$pval_uniq) | vals1$pct_fr<.5 | vals1$pval_uniq>.05 ){keep<-0}
                        if(!is.na(vals$pval_uniq) & vals$pct_fr>.5 & vals$pval_uniq<.05 ){
                            end(vals1)<-end(vals)
                            vals1<-calc_orf_pval(ORFs = vals1,P_sites_rle = covtx,P_sites_uniq_rle = covtx_uniq,P_sites_uniq_mm_rle = covtx_uniq_mm,cutoff = cutoff_fr_ave)
                            vals1$Protein<-AAStringSet(as.character(translate(seq_tx[vals1@ranges],genetic.code = genetic_code_table,if.fuzzy.codon = "solve")))
                            vals1$ORF_orig_tr<-orf_tx$ORF_id_tr
                            vals1$n_stops_readth<-keep
                            vals1$compatible_id<-CharacterList("")
                            vals1$compatible_original_id<-CharacterList("")
                            mcsva<-mcols(vals1)
                            mcsva[,c("gene_id","gene_biotype","gene_name","transcript_id","transcript_biotype")]<-""
                            mcols(vals1)<-mcsva
                            
                        }
                    }
                    
                    
                    
                }
            }
            readthroughs<-unique(suppressWarnings(c(readthroughs,vals1)))
            mcs<-mcols(readthroughs)
            mcs_orfs_orig<-mcs$ORF_orig_tr
            mcs_orfs<-mcs$ORF_id_tr
            mcs$ORF_orig_tr<-NULL
            mcs$ORF_id_tr<-NULL
            dups<-duplicated(mcs)
            readthroughs<-readthroughs[!dups]
            for(j in 1:length(readthroughs)){
                orig_orf<-readthroughs$ORF_orig_tr[j]
                mcs_orig<-mcols(results_orf$ORFs_tx_position[[which(names(results_orf$ORFs_tx_position)==orig_orf)]])[,c("gene_id","gene_biotype","gene_name","transcript_id","transcript_biotype")]
                mcols(readthroughs)[j,c("gene_id","gene_biotype","gene_name","transcript_id","transcript_biotype")]<-mcs_orig
                readthroughs$gen_coords[j]<-from_tx_togen(readthroughs[j],annotation$exons_txs[as.character(seqnames(readthroughs[j]))],annotation$introns_txs[[as.character(seqnames(readthroughs[j]))]])
                readthroughs$compatible_id[j]<-CharacterList(mcs_orfs[which(mcs$Protein==mcs[j,"Protein"])[-j]])
                readthroughs$compatible_original_id[j]<-CharacterList(mcs_orfs_orig[which(mcs$Protein==mcs[j,"Protein"])[-j]])
                
            }
            dups<-duplicated(readthroughs$Protein) & duplicated(readthroughs$P_sites_raw)
            
            if(sum(dups)>0){readthroughs<-readthroughs[!dups]}
        }
    }
    dups<-duplicated(readthroughs$Protein) & duplicated(readthroughs$P_sites_raw)
    
    if(sum(dups)>0){readthroughs<-readthroughs[!dups]}
    return(readthroughs)
}

#' Select and quantify ORF translation
#'
#' This function selects a subset of detected ORFs and quantifies their translation
#' @details ORFs are first selected using the same method as in the \code{select_txs} function, but 
#' using ORF features (ORF structures are treated as transcript structures).\cr
#' Ribo-seq coverage (reads/length) on bins and junctions (set to a length of 60) is used to derive a scaling factor (0-1) for each ORF,
#' which indicates how much of the ORF coverage can be assigned to such ORF (1 when no other ORF is present). 
#' When no unique features are present on an ORF, an adjusted scaling value is calculated subtracting coverage expected from a ORF with a unique feature. 
#' When no unique features are present on any ORF, scaling values are calculated assuming uniform coverage on each ORF.\cr
#' Scaling values are then further scaled to adjust for average coverage (recommended) or total number of reads in the region.\cr
#' ORFs are then further filtered to exclude lowly translated ORFs and quantification/selection is re-iterated until no ORF is further filtered out. 
#' Percentage of total gene translation and length-adjusted quantification estimates are produced.
#' More details about the quantificatin procedure can be found in the ORFquant manuscript.\cr\cr
#' Additional columns are added to the ORFs_tx object:\cr
#' \code{P_sites}: P_sites_raw value from \code{detect_translated_ORFs} divided by the ORF scaling value.\cr
#' \code{ORF_pct_P_sites}: Percentage of gene translation output for the ORF, derived using P_sites values.\cr
#' \code{ORF_pct_P_sites_pN}: Percentage of gene translation ouptut (adjusted by length) for the ORF, derived using P_sites values.\cr
#' \code{unique_features_reads}:  initial number of reads on each unique ORF feature. \code{NA} when no unique feature is present.\cr
#' \code{adj_unique_features_reads}:  final number of reads on each unique ORF feature after the ORF filtering/quantification procedure. \code{NA} when no unique feature is present.\cr
#' \code{scaling_factors}: Set of 3 scaling factors assigned to the ORF using intial unique ORF features, after adjusting for the presence of ORFs with no unique features, and final scaling factor after correcting for average Ribo-seq coverage (or total number of reads) on the ORFs.
#' @keywords ORFquant
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param results_ORFs Full list of detected ORFs, from \code{detect_translated_ORFs}
#' @param P_sites GRanges object with P_sites positions
#' @param P_sites_uniq GRanges object with uniquely mapping P_sites positions
#' @param cutoff_cums cutoff to select ORFs until <x> percentage of total gene translation. Defaults to 99
#' @param cutoff_pct minimum percentage of total gene translation for an ORF to be selected. Defaults to 1
#' @param cutoff_P_sites minimum number of P_sites assigned to the ORF to be selected. Defaults to 10
#' @param optimiz (Beta) should numerical optimization (minimizing distance between observed coverage and expected coverage) 
#' be used to quantify ORF translation? Defaults to FALSE
#' @param scaling Additional scaling value taking into account average or total signal on the detected ORFs to adjust quantification estimates. 
#' Can be average_coverage or total_Psites. Defaults to total_Psites for consistency.
#' @param uniq_signal Use only signal from uniquely mapping reads? Defaults to \code{FALSE}.
#' @return modified \code{results_ORFs} object with the selected ORFs including quantification estimates.
#' @seealso \code{\link{detect_translated_orfs}}, \code{\link{select_txs}}
#' @export


select_quantify_ORFs<-function(results_ORFs,P_sites,P_sites_uniq,cutoff_cums=NA,cutoff_pct=2,cutoff_P_sites=NA,optimiz=FALSE,scaling="total_Psites",uniq_signal=F){
    
    if(!scaling%in%c("total_Psites","average_coverage")){stop(paste("scaling parameter must be either total_Psites (recommended) or average_coverage"),date())}
    select_feat <- results_ORFs[["ORFs_features"]]
    
    select_feat<-endoapply(select_feat,function(x){
        unqid<-paste(x@ranges,x$type,names(x),sep="_")
        x<-x[!duplicated(unqid)]
        x
    })
    
    select_feats<-unlist(select_feat)
    
    nmss<-c()
    for(nm in names(results_ORFs[["ORFs_features"]])){
        nmss<-c(nmss,rep(nm,length(select_feat[[nm]])))
    }
    names(select_feats)<-nmss
    select_feats_jun<-select_feats[select_feats$type=="J"]
    
    # ran_j<-IRanges()
    # allofthem<-unique(unlist(select_feats_jun$txs_orfs))
    # for(i in)
    if(length(select_feats_jun)>0){
        
        
        ran_j<-select_feats_jun@ranges
        #startend!
        df<-data.frame(stend=paste(ran_j@start,end(ran_j),sep="_"),orfs=names(ran_j),stringsAsFactors=F)
        tab_j<-as.matrix(table(df$stend,df$orfs))
        
        orfs_print<-colnames(tab_j)
        listorf_print<-list()
        for(junso in rownames(tab_j)){
            listorf_print[[junso]]<-orfs_print[which(tab_j[junso,]>0)]
        }
        
        orfs_print2<-listorf_print
        orfs_printWRONG<-apply(tab_j,MARGIN=1,FUN=function(x){colnames(tab_j)[which(x>0)]})
        
        #names(orfs_print2)<-NULL
    }
    orfs <- results_ORFs[["ORFs_genomic_position"]]
    stra<-as.character(select_feats[1]@strand)
    
    df_orfs_ex<-data.frame(orfs)
    names(df_orfs_ex)<-c("tx_id","tx_name","tx_chrom","exon_start","exon_end","width","tx_strand")
    df_orfs_ex$exon_rank<-c(unlist(sapply(as.numeric(table(df_orfs_ex$tx_id)),FUN=function(x){1:x})))
    df_orfs<-data.frame(tx_id=sort(unique(df_orfs_ex$tx_id)),tx_name=as.character(names(orfs)),tx_chrom=as.character(unique(seqnames(orfs[[1]]))),tx_start=min(start(orfs)),tx_end=max(end(orfs)),tx_strand=as.character(unique(strand(orfs[[1]]))),stringsAsFactors=F)
    df_genes<-data.frame(tx_name=as.character(names(orfs)),gene_id="OFF")
    orfann<-suppressWarnings(makeTxDb(transcripts=df_orfs,splicings=df_orfs_ex,genes=df_genes))
    #disjointExons(orfann)
    
    
    exbin<-disjointExons(orfann)
    
    
    d<-rep(0,length(exbin))
    
    hts<-findOverlaps(exbin,P_sites,ignore.strand=F)
    hts<-cbind(queryHits(hts),P_sites[subjectHits(hts)]$score*width(P_sites[subjectHits(hts)]))
    if(length(hts)==0){
        return(GRangesList())
        
    }
    hts<-aggregate(x = hts[,2],list(hts[,1]),FUN=sum)
    for(i in 1:dim(hts)[1]){
        d[hts[i,1]]<-hts[i,2]
    }
    
    d2<-rep(0,length(exbin))
    
    hts<-findOverlaps(exbin,P_sites_uniq,ignore.strand=F)
    hts<-cbind(queryHits(hts),P_sites_uniq[subjectHits(hts)]$score*width(P_sites_uniq[subjectHits(hts)]))
    if(length(hts)>0){
        hts<-aggregate(x = hts[,2],list(hts[,1]),FUN=sum)
        for(i in 1:dim(hts)[1]){
            d2[hts[i,1]]<-hts[i,2]
        }
    }
    
    
    exbin$reads<-d
    exbin$unique_reads<-d2
    
    exbin<-GRanges(exbin)
    
    exbin$gene_id<-NULL
    exbin$exonic_part<-NULL
    exbin$X<-NULL
    mcols(exbin)<-exbin[,c("reads","unique_reads","tx_name")]
    gene_feat<-exbin
    gene_feat$type<-rep("E",length(exbin))
    
    if(length(select_feats_jun)>0){
        junc<-select_feats_jun
        names(junc)<-NULL
        reads_j<-junc$reads
        reads_j_uniq<-junc$unique_reads
        mcols(junc)<-NULL
        junc$X<-junc
        junc$reads<-reads_j
        junc$unique_reads<-reads_j_uniq
        
        orfs_jj<-orfs_print2[match(df$stend,names(orfs_print2))]
        #junc$tx_name<-as(orfs_jj,"CharacterList")
        mcols(junc)<-junc[,c("reads","unique_reads")]
        tx_ex<-exbin$tx_name
        mcols(exbin)<-exbin[,c("reads","unique_reads")]
        
        gene_feat<-c(exbin,junc)
        if(is(orfs_jj,"CompressedList")){
            gene_feat$tx_name<-CharacterList(c(tx_ex,orfs_jj))
        }
        
        if(!is(orfs_jj,"CompressedList")){
            if(is.list(orfs_jj)){
                gene_feat$tx_name<-CharacterList(c(tx_ex,CharacterList(orfs_jj)))
            }
            if(!is.list(orfs_jj)){
                gene_feat$tx_name<-CharacterList(c(tx_ex,CharacterList(as.list(orfs_jj))))
            }
        }
        
        gene_feat$type<-c(rep("E",length(exbin)),rep("J",length(junc)))
    }
    
    #first round
    gene_feat_cp<-gene_feat
    #change!
    gene_feat<-gene_feat_cp[!duplicated(paste(GRanges(gene_feat_cp),gene_feat_cp$type,sep="_"))]
    txs_gene<-unique(unlist(gene_feat$tx_name))
    txs_sofar<-txs_gene
    
    #uniq_flag
    
    if(uniq_signal){
        d<-gene_feat$unique_reads
    }
    if(!uniq_signal){
        d<-gene_feat$reads
    }
    
    a<-gene_feat$tx_name
    
    mat<-matrix(data=0,nrow=length(d),ncol=length(txs_sofar))
    colnames(mat)<-txs_sofar
    for(i in 1:length(txs_sofar)){
        mat[,i]<-sapply(a,function(x){sum(x==txs_sofar[i])})
    }
    
    #TAKE AWAY NESTED TXS
    nest<-c()
    ident<-c()
    for(i in 1:dim(mat)[2]){
        yes<-which(mat[,i]==1)
        nesti<-c()
        for(j in (1:dim(mat)[2])[-i]){
            yesj<-which(mat[,j]==1)
            if(identical(yes,yesj)){ident<-c(ident,paste(sort(colnames(mat)[c(i,j)]),collapse=";"))}
            #added if length> otherwise txs with same structure are both deleted
            nesti<-c(nesti,sum(yesj%in%yes)==length(yes) & length(mat[,j])>length(yes))
        }
        
        if(sum(nesti)>0){nest<-c(nest,colnames(mat)[i])}
        
    }
    if(length(ident)>0){nest<-nest[!nest%in%unique(sapply(strsplit(ident,";"),"[[",1))]}
    txs_sofar<-txs_sofar[!txs_sofar%in%nest]
    change<-1
    
    while(change>0){
        
        mat<-matrix(data=0,nrow=length(d),ncol=length(txs_sofar))
        colnames(mat)<-txs_sofar
        for(i in 1:length(txs_sofar)){
            mat[,i]<-sapply(a,function(x){sum(x==txs_sofar[i])})
        }
        
        mat_orig<-mat
        d_count<-paste(1:length(d),d,sep="_")
        good<-d_count[which(d>0)]
        bad<-d_count[which(d==0)]
        
        txs_good<-c()
        expl_good<-c()
        for(i in 1:dim(mat)[2]){
            expl_good_old<-expl_good
            #if new good feature
            tx<-d_count[which(mat[,i]>0)]
            tx_good<-tx[which(tx%in%good)]
            tx_bad<-tx[which(tx%in%bad)]
            
            if(length(tx_good)==0){next}
            
            if(sum(!tx_good%in%expl_good_old)>0){
                
                tx_torem<-c()
                tx_toscreen<-which(colnames(mat)%in%txs_good)
                if(length(tx_toscreen)>0){
                    for(j in tx_toscreen){
                        tx_contr<-d_count[which(mat[,j]>0)]
                        tx_contr_good<-tx_contr[which(tx_contr%in%good)]
                        tx_contr_bad<-tx_contr[which(tx_contr%in%bad)]
                        if(sum(tx_contr_good%in%tx_good)==length(tx_contr_good)){
                            tx_torem<-c(tx_torem,colnames(mat)[j])
                            
                        }
                        
                        if(length(tx_torem)>0){
                            txs_good<-txs_good[!txs_good%in%tx_torem]
                            
                        }
                    }
                    
                    
                    
                }
                txs_good<-unique(c(txs_good,colnames(mat)[i]))
                expl_good<-unique(c(expl_good,tx_good))
            }
            #if same good feature, but fewer bad INTERNAL features than others.
            if(sum(!tx_good%in%expl_good_old)==0){
                tx_torem<-c()
                tx_toscreen<-which(colnames(mat)%in%txs_good)
                #here a counter when at least one good feature more than competing
                moref<-c()
                for(j in tx_toscreen){
                    tx_contr<-d_count[which(mat[,j]>0)]
                    tx_contr_good<-tx_contr[which(tx_contr%in%good)]
                    tx_contr_bad<-tx_contr[which(tx_contr%in%bad)]
                    moref<-c(moref,sum(!tx_good%in%tx_contr_good)>0)
                    #if same good features in competing tx, check bad internal ones
                    if(sum(tx_contr_good%in%tx_good)==length(tx_contr_good)){
                        
                        if(length(tx_good)>length(tx_contr_good)){
                            tx_torem<-c(tx_torem,colnames(mat)[j])
                        }
                        #here
                        if(length(tx_good)==length(tx_contr_good)){      
                            fi<-which(tx==tx_good[1])
                            la<-which(tx==tx_good[length(tx_good)])
                            int_tx<-tx[fi:la]
                            int_tx_bad<-tx_bad[tx_bad%in%int_tx]
                            
                            contr_fi<-which(tx_contr==tx_contr_good[1])
                            contr_la<-which(tx_contr==tx_contr_good[length(tx_contr_good)])
                            contr_int_tx<-tx_contr[contr_fi:contr_la]
                            contr_int_tx_bad<-tx_contr_bad[tx_contr_bad%in%contr_int_tx]
                            #SAME INTERNAL, TAKE
                            if(length(int_tx_bad)<=length(contr_int_tx_bad)){
                                
                                txs_good<-unique(c(txs_good,colnames(mat)[i]))
                                expl_good<-unique(c(expl_good,tx_good))
                                #LESS INTERNAL, TAKE AND REMOVE OTHER
                                if(length(int_tx_bad)<length(contr_int_tx_bad)){
                                    
                                    tx_torem<-c(tx_torem,colnames(mat)[j])}
                            }
                            
                            
                        }
                    }
                    
                    
                }
                if(length(tx_torem)>0){
                    txs_good<-txs_good[!txs_good%in%tx_torem]
                    txs_good<-unique(c(txs_good,colnames(mat)[i]))
                    expl_good<-unique(c(expl_good,tx_good))
                }
                
                if(length(tx_torem)==0 & sum(moref)==length(tx_toscreen)){
                    txs_good<-unique(c(txs_good,colnames(mat)[i]))
                    expl_good<-unique(c(expl_good,tx_good))
                }
                
                
            }
            
            
        }
        txs_sofar<-txs_good
        change<-abs(length(txs_good)-length(txs_sofar))
    }
    
    gene_feat$ORF_id_tr<-gene_feat$tx_name
    gene_feat$tx_name<-NULL
    a<-gene_feat$ORF_id_tr
    a<-lapply(a,function(x){x[x%in%txs_good]})
    gene_feat$ORF_id_tr_selected<-CharacterList(a)
    
    check<-sapply(gene_feat$ORF_id_tr,FUN = length)
    use<-rep("shared",length(check))
    use[check==1]<-"unique"
    use[check==0]<-"absent"
    use[check>1]<-"shared"
    gene_feat$use_ORF<-use
    
    check<-sapply(gene_feat$ORF_id_tr_selected,FUN = length)
    use<-rep("shared",length(check))
    use[check==1]<-"unique"
    use[check==0]<-"absent"
    use[check>1]<-"shared"
    gene_feat$use_ORF_selected<-use
    
    final_ranges<-sort(gene_feat)
    sel_feats<-list()
    for(i in txs_good){
        featexs<-gene_feat[gene_feat$type=="E"]
        featjuns<-gene_feat[gene_feat$type=="J"]
        
        
        ok<-featexs[featexs%over%orfs[[i]]]
        if(length(featjuns)>0){
            ok<-sort(c(ok,featjuns[which(featjuns%in%gaps(orfs[[i]]))]))
        }
        a<-sapply(ok$ORF_id_tr,function(x){length(x[x%in%i])})
        ok<-ok[a>0]
        sel_feats[[i]]<-unique(ok)
    }
    all_feats<-results_ORFs$selected_ORFs_features
    selected_ORFs<-lapply(results_ORFs,FUN=function(x){x[names(x)%in%txs_good]})
    selected_ORFs$selected_ORFs_features<-sel_feats
    
    #quantif 
    
    
    results_ORFs<-selected_ORFs
    feats<-results_ORFs$selected_ORFs_features
    orfs_tx<-results_ORFs$ORFs_tx_position
    orfs_tx<-lapply(orfs_tx,function(x){
        cols<-mcols(x)
        cols[,c("P_sites","ORF_pct_P_sites","ORF_pct_P_sites_pN")]<-NA
        cols[,"unique_features_reads"]<-NumericList("")
        cols[,"adj_unique_features_reads"]<-NumericList("")
        cols[,"scaling_factors"]<-NumericList("")
        mcols(x)<-cols
        x
    })
    
    #first round of unq
    
    orf_del<-c("")
    counter<-1
    gene_feat$ORF_id_tr_selected_quant<- gene_feat$ORF_id_tr_selected
    gene_feat$use_ORF_selected_quant<-gene_feat$use_ORF_selected
    
    while(length(orf_del)>0){
        
        feats<-feats[!names(feats)%in%orf_del]
        nms<-sapply(feats,length)
        nms<-rep(names(nms),nms)
        feats<-unlist(GRangesList(unlist(feats)))
        names(feats)<-nms
        
        feats$ORF_id_tr_selected<-CharacterList(lapply(feats$ORF_id_tr_selected,function(x){
            unique(x[!x%in%orf_del])
        }))
        
        gene_feat$ORF_id_tr_selected_quant<-CharacterList(lapply(gene_feat$ORF_id_tr_selected_quant,function(x){
            unique(x[!x%in%orf_del])
        }))
        lens<-sapply(feats$ORF_id_tr_selected,length)
        lens2<-sapply(gene_feat$ORF_id_tr_selected_quant,length)
        gene_feat$use_ORF_selected_quant<-"shared"
        gene_feat$use_ORF_selected_quant[lens2==1]<-"unique"
        gene_feat$use_ORF_selected_quant[lens2==0]<-"absent"
        
        feats$use_ORF_selected<-"shared"
        feats$use_ORF_selected[lens==1]<-"unique"
        feats<-split(feats,names(feats))
        
        orfs_tx<-orfs_tx[!names(orfs_tx)%in%orf_del]
        
        
        unqs<-rep(NA,length(feats))
        names(unqs)<-names(feats)
        for(i in names(feats)){
            feat<-feats[[i]]
            orf_tx<-orfs_tx[[i]]
            
            #uniq_flag
            if(!uniq_signal){
                riz<-feat$reads
            }
            
            if(uniq_signal){
                riz<-feat$unique_reads
            }
            
            
            cov_feat<-riz/width(feat)
            js<-feat$type=="J"
            if(sum(js)>0){
                cov_feat[js]<-riz[js]/60
            }
            unq<-feat$use_ORF_selected=="unique"
            if(sum(unq)>0){
                unq_rat<-mean(cov_feat[unq])/mean(cov_feat)
                if(unq_rat>1){unq_rat<-1}
                orfs_tx[[i]]$unique_features_reads<-NumericList(riz[unq])
                orfs_tx[[i]]$adj_unique_features_reads<-NumericList(riz[unq])
                
                
            }
            if(sum(unq)==0){
                orfs_tx[[i]]$unique_features_reads<-NumericList(NA)
                unq_rat<-NA
            }
            unqs[i]<-unq_rat
        }
        unqs_adj<-unqs
        unqnas<-unqs_adj[is.na(unqs_adj)]
        unqs_okk<-unqs_adj[!is.na(unqs_adj)]
        #NA_adjustment
        
        #if all NAs
        if(length(unqs_okk)==0){
            unqs_na<-names(unqs_adj[is.na(unqs_adj)])
            for(i in unqs_na){
                feat<-feats[[i]]
                
                riz<-feat$reads
                
                if(uniq_signal){
                    riz<-feat$unique_reads
                }
                
                cov_feat<-riz/width(feat)
                
                js<-feat$type=="J"
                if(sum(js)>0){
                    cov_feat[js]<-riz[js]/60
                }
                cov_adj<-cov_feat
                for(j in 1:length(feat)){
                    fea<-feat[j]
                    txs_fea<-unlist(fea$ORF_id_tr_selected)
                    nass<-txs_fea%in%unqs_na
                    txs_fea<-txs_fea[!nass]
                    adj<-cov_feat[j]-(cov_feat[j]*sum(unqs_adj[txs_fea]))
                    if(sum(nass)>1){
                        adj<-(cov_feat[j]-(cov_feat[j]*sum(unqs_adj[txs_fea])))/sum(nass)
                    }
                    if(length(txs_fea)==0){
                        adj<-cov_feat[j]/(sum(nass))
                    }
                    if(adj<0){adj=0}
                    cov_adj[j]<-adj
                }
                unq<-mean(cov_adj)/mean(cov_feat)
                if(unq>1){unq<-1}
                unqs_adj[i]<-unq
            }
        }
        unqnas<-unqs_adj[is.na(unqs_adj)]
        unqs_okk<-unqs_adj[!is.na(unqs_adj)]
        
        #if all 0
        
        if(sum(unqs_okk==0)==length(unqs_adj)){
            unqs_zero<-names(unqs_okk)
            for(i in unqs_zero){
                feat<-feats[[i]]
                
                riz<-feat$reads
                
                if(uniq_signal){
                    riz<-feat$unique_reads
                }
                
                cov_feat<-riz/width(feat)
                
                js<-feat$type=="J"
                if(sum(js)>0){
                    cov_feat[js]<-riz[js]/60
                }
                cov_adj<-cov_feat
                for(j in 1:length(feat)){
                    fea<-feat[j]
                    txs_fea<-unlist(fea$ORF_id_tr_selected)
                    nass<-txs_fea%in%unqs_zero
                    txs_fea<-txs_fea[!nass]
                    adj<-cov_feat[j]-(cov_feat[j]*sum(unqs_okk[txs_fea]))
                    if(sum(nass)>1){
                        adj<-(cov_feat[j]-(cov_feat[j]*sum(unqs_okk[txs_fea])))/sum(nass)
                    }
                    if(length(txs_fea)==0){
                        adj<-cov_feat[j]/(sum(nass))
                    }
                    if(adj<0){adj=0}
                    cov_adj[j]<-adj
                }
                unq<-mean(cov_adj)/mean(cov_feat)
                if(unq>1){unq<-1}
                unqs_okk[i]<-unq
                unqs_adj[i]<-unq
                
            }
        }
        unqnas<-unqs_adj[is.na(unqs_adj)]
        unqs_okk<-unqs_adj[!is.na(unqs_adj)]
        
        if(length(unqnas)>0){
            nonas<--1
            
            while(nonas<0){
                
                prev_nas<-length(unqnas)
                if(length(unqs_okk)>0){
                    for(i in names(unqnas)){
                        unqs_okk_noi<-unqs_okk[names(unqs_okk)!=i]
                        feat<-feats[[i]]
                        orf_tx<-orfs_tx[[i]]
                        
                        riz<-feat$reads
                        
                        if(uniq_signal){
                            riz<-feat$unique_reads
                        }
                        
                        cov_feat<-riz/width(feat)
                        
                        js<-feat$type=="J"
                        if(sum(js)>0){
                            cov_feat[js]<-riz[js]/60
                        }
                        adj_use<-feat$use_ORF_selected
                        adj_cov_feat<-cov_feat
                        for(j in 1:length(feat)){
                            txs_fea<-feat$ORF_id_tr_selected[[j]]
                            unqs_okk_noi_feat<-unqs_okk_noi[names(unqs_okk_noi)%in%txs_fea]
                            txs_fea_adj<-txs_fea[!txs_fea%in%names(unqs_okk_noi)]
                            adj<-cov_feat[j]
                            usef<-"shared"
                            if(length(txs_fea_adj)==1){usef<-"unique"}
                            
                            if(length(unqs_okk_noi_feat[!is.na(unqs_okk_noi_feat)])>0){
                                adj<-adj-(adj*sum(unqs_okk_noi_feat[!is.na(unqs_okk_noi_feat)]))     
                                if(adj<0){adj=0}
                            }
                            adj_use[j]<-usef
                            adj_cov_feat[j]<-adj
                        }
                        
                        unq<-adj_use=="unique"
                        
                        if(sum(unq)>0){
                            unq_rat<-mean(adj_cov_feat[adj_use=="unique"])/mean(adj_cov_feat)
                            if(mean(adj_cov_feat)==0){unq_rat<-0}
                            if(unq_rat>1){unq_rat<-1}
                            orfs_tx[[i]]$adj_unique_features_reads<-NumericList(riz[unq])
                            
                        }
                        if(sum(unq)==0){
                            orfs_tx[[i]]$adj_unique_features_reads<-NumericList(NA)
                            unq_rat<-NA
                        }
                        unqs_adj[i]<-unq_rat
                        
                        
                    }
                    unqnas<-unqs_adj[is.na(unqs_adj)]
                    unqs_okk<-unqs_adj[!is.na(unqs_adj)]
                    nonas<-length(unqnas)-prev_nas
                    
                }
                
                
            }
            #if still NAs
            if(sum(is.na(unqs_adj))>0){
                unqs_na<-names(unqs_adj[is.na(unqs_adj)])
                for(i in unqs_na){
                    feat<-feats[[i]]
                    
                    riz<-feat$reads
                    
                    if(uniq_signal){
                        riz<-feat$unique_reads
                    }
                    
                    cov_feat<-riz/width(feat)
                    
                    
                    js<-feat$type=="J"
                    if(sum(js)>0){
                        cov_feat[js]<-riz[js]/60
                    }
                    cov_adj<-cov_feat
                    for(j in 1:length(feat)){
                        fea<-feat[j]
                        txs_fea<-unlist(fea$ORF_id_tr_selected)
                        nass<-txs_fea%in%unqs_na
                        txs_fea<-txs_fea[!nass]
                        adj<-cov_feat[j]-(cov_feat[j]*sum(unqs_adj[txs_fea]))
                        if(sum(nass)>1){
                            adj<-(cov_feat[j]-(cov_feat[j]*sum(unqs_adj[txs_fea])))/sum(nass)
                        }
                        if(length(txs_fea)==0){
                            adj<-cov_feat[j]/(sum(nass))
                        }
                        if(adj<0){adj=0}
                        cov_adj[j]<-adj
                    }
                    unq<-mean(cov_adj)/mean(cov_feat)
                    if(unq>1){unq<-1}
                    unqs_adj[i]<-unq
                }
                
            }
            
            
        }
        
        
        if(sum(unqs_adj>0)==0){
            nmmm<-names(unqs_adj)
            unqs_adj<-rep(1/length(nmmm),length(nmmm))
            names(unqs_adj)<-nmmm
        }
        
        orftxs<-unlist(results_ORFs[["ORFs_tx_position"]])
        unqs_optim<-unqs_adj
        
        if(!is(feats,"GRangesList")){feats<-GRangesList(feats)}
        featsall<-sort(unlist(feats))
        featsall<-featsall[!duplicated(paste(GRanges(featsall),featsall$type,sep="_"))]
        
        #optimization option
        
        if(optimiz==TRUE){
            unqs_optim_noopt<-unqs_optim
            #uniq_flag
            
            if(!uniq_signal){
                cov_ps_unqs<-orftxs$P_sites_raw/width(orftxs)
            }
            
            if(uniq_signal){
                cov_ps_unqs<-orftxs$P_sites_raw_uniq/width(orftxs)
            }
            
            
            names(cov_ps_unqs)<-orftxs$ORF_id_tr
            cov_ps_unqs<-cov_ps_unqs[names(unqs_adj)]
            #uniq_flag
            
            if(!uniq_signal){
                featsall$coverage<-featsall$reads/width(featsall)
                jxs<-featsall$type=="J"
                featsall$coverage[jxs]<-featsall$reads[jxs]/60
            }
            
            if(uniq_signal){
                featsall$coverage<-featsall$unique_reads/width(featsall)
                jxs<-featsall$type=="J"
                featsall$coverage[jxs]<-featsall$unique_reads[jxs]/60
            }
            
            vals<-cov_ps_unqs*unqs_optim
            expect<-unlist(lapply(featsall$ORF_id_tr_selected,function(x){sum(vals[x])}))
            maxdist<-sum(abs(expect-featsall$coverage))
            
            calc_dist_exp<-function(unqs_optim,maxval=maxdist){
                vals<-cov_ps_unqs*unqs_optim
                #sizzs<-width(featsall)
                #sizzs[featsall$type=="J"]<-60
                
                expect<-unlist(lapply(featsall$ORF_id_tr_selected,function(x){sum(vals[x])}))
                if(sum(expect)>0){
                    expect<-expect/sum(expect)
                }
                uniq_feats<-featsall$use_ORF_selected=="unique"
                trucov<-featsall$coverage
                if(sum(expect)>0){
                    trucov<-trucov/sum(trucov)
                }
                #or do something different, like minimize % error per each feature after adding pseudocount
                if(sum(trucov>0 & expect==0)>0){return(maxval+1)}
                
                #sum(abs(expect-featsall$coverage)*log(sizzs+1))
                #-cor.test(trucov,expect,method = "p")$estimate
                mean(abs(expect-trucov))
            }
            #inspired by Alpine: https://github.com/mikelove/alpine/blob/master/R/estimate_abundance.R
            optimm<-optim(par = unqs_optim,fn = calc_dist_exp,lower = rep(0,length(unqs_optim)),upper = rep(1,length(unqs_optim)),method = "L-BFGS-B")
            unqs_optim<-optimm$par
            if(sum(unqs_optim>0)==0){unqs_optim<-unqs_optim_noopt}
        }
        
        if(scaling=="average_coverage"){
            
            #scale to adjust coverage?
            
            #uniq_flag
            
            if(!uniq_signal){
                cov_ps_unqs<-orftxs$P_sites_raw/width(orftxs)
            }
            
            if(uniq_signal){
                cov_ps_unqs<-orftxs$P_sites_raw_uniq/width(orftxs)
            }
            
            names(cov_ps_unqs)<-orftxs$ORF_id_tr
            cov_ps_unqs<-cov_ps_unqs[names(unqs_optim)]
            
            #uniq_flag
            if(uniq_signal){
                covo<-featsall$unique_reads/width(featsall)
            }
            
            if(!uniq_signal){
                covo<-featsall$reads/width(featsall)
            }
            featsall$coverage<-covo
            jxs<-featsall$type=="J"
            featsall$coverage[jxs]<-covo[jxs]/60
            
            cov_ps_unqs_adj<-cov_ps_unqs*unqs_optim
            cov_unpcts<-NumericList(lapply(featsall$ORF_id_tr_selected,function(x){cov_ps_unqs_adj[x]}))
            cov_unpcts<-sum(cov_unpcts)
            
            scale_f<-sum(featsall$coverage)/sum(cov_unpcts)
            if(sum(scale_f>1)>0){scale_f<-scale_f/max(scale_f)}
            unqs_optim<-unqs_optim*scale_f
        }
        
        if(scaling=="total_Psites"){
            
            unqs_tott<-unqs_optim
            orfgrps<-list()
            for(grpo in 1:length(featsall)){
                fezzo=featsall$ORF_id_tr_selected[[grpo]]
                presalready<-fezzo%in%unlist(orfgrps)
                if(sum(presalready)==0){orfgrps<-append(orfgrps,fezzo);next()}
                if(sum(presalready)>0){
                    presalready<-which(orfgrps%in%fezzo)
                    for(pres in presalready){
                        orfgrps[[pres]]<-sort(unique(c(orfgrps[[pres]],fezzo)))
                    }
                }
                
            }
            
            orfgrps<-unique(orfgrps)
            if(length(orfgrps)>1){
                orfgrpsdello<-rep(FALSE,length(orfgrps))
                for(grpo in 1:length(orfgrps)){
                    orfes<-orfgrps[[grpo]]
                    otherse<-orfgrps[-grpo]
                    remo<-FALSE
                    for(garpo in 1:length(otherse)){
                        contest<-otherse[[garpo]]
                        if(sum(orfes%in%contest)==length(orfes)){remo=T;break}
                    }
                    orfgrpsdello[grpo]=remo
                }
                
                orfgrps<-orfgrps[!orfgrpsdello]
            }
            
            ps_feats_tott<-c()
            fetse<-featsall[featsall$type=="E"]
            rizz<-fetse$reads
            ps_tott<-orftxs$P_sites_raw
            names(ps_tott)<-orftxs$ORF_id_tr
            ps_tott<-ps_tott[names(unqs_tott)]
            
            if(uniq_signal){
                ps_tott<-orftxs$P_sites_raw_uniq
                names(ps_tott)<-orftxs$ORF_id_tr
                ps_tott<-ps_tott[names(unqs_tott)]
                rizz<-fetse$unique_reads
            }
            
            for(grpo in 1:length(orfgrps)){
                ps_feats_tott<-c(ps_feats_tott,sum(rizz[sum(fetse$ORF_id_tr_selected%in%orfgrps[[grpo]])>0]))
            }
            
            adj_psts<-unqs_tott*ps_tott
            final_scls<-unqs_tott
            
            for(grpo in 1:length(orfgrps)){
                adj_grp<-adj_psts[orfgrps[[grpo]]]
                final_scls[orfgrps[[grpo]]]<-ps_feats_tott[grpo]/sum(adj_grp)
            }
            final_scls[is.infinite(final_scls)]<-0
            
            # to be added?
            # final_scls[is.infinite(final_scls)]<-0
            # final_scls[is.na(final_scls)]<-0
            # final_scls[is.nan(final_scls)]<-0
            # unqs_optim<-final_scls*unqs_tott
            # 
            # unqs_optim[is.infinite(unqs_optim)]<-0
            # unqs_optim[is.na(unqs_optim)]<-0
            # unqs_optim[is.nan(unqs_optim)]<-0
            
            
            unqs_optim<-final_scls*unqs_tott
            
            unqs_optim[unqs_optim>1]<-1
            
            
        }
        #apply quantification factor
        
        for(i in names(feats)){
            feat<-feats[[i]]
            orf_tx<-orfs_tx[[i]]
            
            #uniq_flag
            if(!uniq_signal){
                ps<-orf_tx$P_sites_raw
            }
            
            if(uniq_signal){
                ps<-orf_tx$P_sites_raw_uniq
            }
            unq_rat<-unqs_optim[i]
            
            ps_norm<-ps*unq_rat
            
            orfs_tx[[i]]$P_sites<-ps_norm
            scalss<-c(unqs[i],unqs_adj[i],unqs_optim[i])
            names(scalss)<-c("unq_feats","adj_feats","optim_feats")
            orfs_tx[[i]]$scaling_factors<-NumericList(round(scalss,digits = 4))
            
        }
        
        genes<-unlist(GRangesList(orfs_tx))$gene_id
        genes_unq<-unique(genes)
        orfs_genes<-split(unlist(GRangesList(orfs_tx)),f=unlist(GRangesList(orfs_tx))$gene_id)
        orfs_genes<-GRangesList(lapply(orfs_genes,FUN=function(x){
            xnot<-x[is.na(x$P_sites)]
            x<-x[!is.na(x$P_sites)]
            x$ORF_pct_P_sites<-x$P_sites*100/sum(x$P_sites)
            x$ORF_pct_P_sites_pN<-(x$P_sites/width(x))*100/sum(x$P_sites/width(x))
            #added this when genes, mostly overlapping ones, get no reads
            
            if(sum(x$P_sites)==0){x$ORF_pct_P_sites<-0;x$ORF_pct_P_sites_pN<-0}
            
            c(x[order(x$ORF_pct_P_sites,decreasing=T)],xnot)
            
            
        }))
        
        orf_del_cums<-c()
        orf_del_iso<-c()
        orf_del_ps<-c()
        
        list_genes<-list()
        list_genes_feats<-list()
        
        cums_list<-list()
        for(h in names(orfs_genes)){
            x<-orfs_genes[[h]]
            cums<-cumsum(x$ORF_pct_P_sites)
            names(cums)<-x$ORF_id_tr
            cums_list[[h]]<-cums
        }
        if(is.numeric(cutoff_cums)){
            orf_del_cums<-lapply(cums_list,function(x){
                del<-which(x>cutoff_cums)
                delna<-which(is.na(x))
                if(length(del)==1){del<-c()}
                if(length(del)>1){del<-del[-1]}
                if(length(delna)==1){del<-c(del,delna)}
                del
                
            })
            names(orf_del_cums)<-NULL
            orf_del_cums<-names(unlist(orf_del_cums))
        }
        
        if(is.numeric(cutoff_pct)){
            orfgrlun<-unlist(orfs_genes)
            orf_del_iso<-unique(c(orf_del_iso,orfgrlun$ORF_id_tr[which(orfgrlun$ORF_pct_P_sites<cutoff_pct)]))
        }
        
        if(is.numeric(cutoff_P_sites)){
            orfgrlun<-unlist(orfs_genes)
            orf_del_ps<-orfgrlun$ORF_id_tr[which(orfgrlun$P_sites<cutoff_P_sites)]
            
        }
        
        orf_del<-unique(c(orf_del_cums,orf_del_iso,orf_del_ps))
        
        gene_feat$ORF_id_tr_selected_quant<-CharacterList(lapply(gene_feat$ORF_id_tr_selected_quant,function(x){
            unique(x[!x%in%orf_del])
        }))
        lens2<-sapply(gene_feat$ORF_id_tr_selected_quant,length)
        gene_feat$use_ORF_selected_quant<-"shared"
        gene_feat$use_ORF_selected_quant[lens2==1]<-"unique"
        gene_feat$use_ORF_selected_quant[lens2==0]<-"absent"
        
        fs<-unlist(orfs_genes)
        #put isovalues
        for(g in names(orfs_tx)){
            orfs_tx[[g]]$ORF_pct_P_sites<-round(fs$ORF_pct_P_sites[fs$ORF_id_tr==g],digits = 4)
            orfs_tx[[g]]$ORF_pct_P_sites_pN<-round(fs$ORF_pct_P_sites_pN[fs$ORF_id_tr==g],digits = 4)
        }
        counter<-counter+1
    }
    
    results_ORFs$ORFs_tx_position<-orfs_tx
    results_ORFs$ORFs_genomic_position<-results_ORFs$ORFs_genomic_position[names(results_ORFs$ORFs_tx_position)]
    results_ORFs$ORFs_features<-results_ORFs$ORFs_features[names(results_ORFs$ORFs_tx_position)]
    results_ORFs$tx_annotated_ORFs<-results_ORFs$tx_annotated_ORFs[names(results_ORFs$ORFs_tx_position)]
    #adjust feats!
    
    results_ORFs$selected_ORFs_features<-gene_feat
    return(results_ORFs)
    
}

#' Annotate splice features of detected ORFs
#'
#' This function detects usage of different exons and exonic boundaries of one ORF with respect to a reference ORF.
#' @details each exon is aligned to the closest one to match acceptor and donor sites, or to annotate missing exons.
#' \code{5ss} and \code{3ss} indicate exon 5' and 3', respectively. \code{CDS_spanning} indicates retained intron;
#' \code{missing_CDS} indicates no overlapping exon (missed or included); \code{monoCDS} indicates a single-exon ORF; 
#' \code{firstCDS} and \code{lastCDS} indicate first CDS exon or last CDS exon.
#' @keywords ORFquant
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param orf_gen Exon structure of a detected ORF
#' @param ref_cds Exon structure of a reference ORF
#' @return Exon structure of detected ORF including possible missing exons from reference, together with a \code{spl_type} column
#' including the annotation for each exon (e.g. alternative acceptors or donor).
#' @seealso \code{\link{detect_translated_orfs}}, \code{\link{annotate_ORFs}}
#' @export

annotate_splicing<-function(orf_gen,ref_cds){
    
    overref<-orf_gen%over%ref_cds
    refover<-ref_cds%over%orf_gen
    spl_ran<-GRanges()
    
    if(sum(!refover)>0){
        spl_ran<-c(spl_ran,ref_cds[!refover])
        grliss<-GRangesList()
        refgrl<-ref_cds[!refover]
        for(gri in 1:length(refgrl)){
            grliss[[gri]]<-refgrl[gri]
        }
        spl_ran$ref<-grliss
        spl_ran$spl_type<-"missing_CDS"
        spl_ran$cds_id<-NULL
        spl_ran$cds_name<-NULL
        spl_ran$exon_rank<-NULL
        
    }
    
    orf_gen<-sort(orf_gen)
    if(length(orf_gen)>0){
        for(f in 1:length(orf_gen)){
            ran<-orf_gen[f]
            last_ex<-length(orf_gen)
            if(overref[f]==T){
                ref_over<-ref_cds[ref_cds%over%ran]
                #annotate for 5' and 3'; porcoddio
                
                if(length(ref_over)>1){
                    ran$ref<-GRangesList(ref_over)
                    ran$spl_type<-"CDS_spanning"
                    if(as.vector(strand(orf_gen[1]))=="+"){
                        if(start(ran)==min(start(ref_over))){
                            if(end(ran)==max(end(ref_over))){
                                ran$spl_type<-"CDS_spanning;same_5ss;same_3ss"
                                if(f==last_ex){ran$spl_type<-"CDS_spanning;same_5ss;same_lastCDS"}
                                if(f==1){ran$spl_type<-"CDS_spanning;same_firstCDS;same_3ss"}
                                if(1==last_ex){ran$spl_type<-"CDS_spanning;same_5monoCDS;same_3monoCDS"}
                            }
                            if(end(ran)>max(end(ref_over))){
                                ran$spl_type<-"CDS_spanning;same_5ss;down_3ss"
                                if(f==last_ex){ran$spl_type<-"CDS_spanning;same_5ss;down_lastCDS"}
                                if(f==1){ran$spl_type<-"CDS_spanning;same_firstCDS;down_3ss"}
                                if(1==last_ex){ran$spl_type<-"CDS_spanning;same_5monoCDS;down_3monoCDS"}
                            }
                            if(end(ran)<max(end(ref_over))){
                                ran$spl_type<-"CDS_spanning;same_5ss;up_3ss"
                                if(f==last_ex){ran$spl_type<-"CDS_spanning;same_5ss;up_lastCDS"}
                                if(f==1){ran$spl_type<-"CDS_spanning;same_firstCDS;up_3ss"}
                                if(1==last_ex){ran$spl_type<-"CDS_spanning;same_5monoCDS;up_3monoCDS"}
                            }
                            
                        }
                        if(end(ran)==max(end(ref_over))){
                            if(start(ran)>min(start(ref_over))){
                                ran$spl_type<-"CDS_spanning;down_5ss;same_3ss"
                                if(f==last_ex){ran$spl_type<-"CDS_spanning;down_5ss;same_lastCDS"}
                                if(f==1){ran$spl_type<-"CDS_spanning;down_firstCDS;same_3ss"}
                                if(1==last_ex){ran$spl_type<-"CDS_spanning;down_5monoCDS;same_3monoCDS"}
                            }
                            if(start(ran)<min(start(ref_over))){
                                ran$spl_type<-"CDS_spanning;up_5ss;same_3ss"
                                if(f==last_ex){ran$spl_type<-"CDS_spanning;up_5ss;same_lastCDS"}
                                if(f==1){ran$spl_type<-"CDS_spanning;up_firstCDS;same_3ss"}
                                if(1==last_ex){ran$spl_type<-"CDS_spanning;up_5monoCDS;same_3monoCDS"}
                            }
                            
                        }
                        
                        if(end(ran)>max(end(ref_over))){
                            if(start(ran)>min(start(ref_over))){
                                ran$spl_type<-"CDS_spanning;down_5ss;down_3ss"
                                if(f==last_ex){ran$spl_type<-"CDS_spanning;down_5ss;down_lastCDS"}
                                if(f==1){ran$spl_type<-"CDS_spanning;down_firstCDS;down_3ss"}
                                if(1==last_ex){ran$spl_type<-"CDS_spanning;down_5monoCDS;down_3monoCDS"}
                            }
                            if(start(ran)<min(start(ref_over))){
                                ran$spl_type<-"CDS_spanning;up_5ss;down_3ss"
                                if(f==last_ex){ran$spl_type<-"CDS_spanning;up_5ss;down_lastCDS"}
                                if(f==1){ran$spl_type<-"CDS_spanning;up_firstCDS;down_3ss"}
                                if(1==last_ex){ran$spl_type<-"CDS_spanning;up_5monoCDS;down_3monoCDS"}
                            }
                            
                        }
                        
                        if(end(ran)<max(end(ref_over))){
                            if(start(ran)>min(start(ref_over))){
                                ran$spl_type<-"CDS_spanning;down_5ss;up_3ss"
                                if(f==last_ex){ran$spl_type<-"CDS_spanning;down_5ss;up_lastCDS"}
                                if(f==1){ran$spl_type<-"CDS_spanning;down_firstCDS;up_3ss"}
                                if(1==last_ex){ran$spl_type<-"CDS_spanning;down_5monoCDS;up_3monoCDS"}
                            }
                            if(start(ran)<min(start(ref_over))){
                                ran$spl_type<-"CDS_spanning;up_5ss;up_3ss"
                                if(f==last_ex){ran$spl_type<-"CDS_spanning;up_5ss;up_lastCDS"}
                                if(f==1){ran$spl_type<-"CDS_spanning;up_firstCDS;up_3ss"}
                                if(1==last_ex){ran$spl_type<-"CDS_spanning;up_5monoCDS;up_3monoCDS"}
                            }
                            
                        }
                        
                        
                        
                        
                    }
                    
                    #if - and spanning
                    
                    if(as.vector(strand(orf_gen[1]))=="-"){
                        if(start(ran)==min(start(ref_over))){
                            if(end(ran)==max(end(ref_over))){
                                ran$spl_type<-"CDS_spanning;same_5ss;same_3ss"
                                if(f==1){ran$spl_type<-"CDS_spanning;same_5ss;same_lastCDS"}
                                if(f==last_ex){ran$spl_type<-"CDS_spanning;same_firstCDS;same_3ss"}
                                if(1==last_ex){ran$spl_type<-"CDS_spanning;same_5monoCDS;same_3monoCDS"}
                            }
                            if(end(ran)>max(end(ref_over))){
                                ran$spl_type<-"CDS_spanning;up_5ss;same_3ss"
                                if(f==1){ran$spl_type<-"CDS_spanning;up_5ss;same_lastCDS"}
                                if(f==last_ex){ran$spl_type<-"CDS_spanning;up_firstCDS;same_3ss"}
                                if(1==last_ex){ran$spl_type<-"CDS_spanning;up_5monoCDS;same_3monoCDS"}
                            }
                            if(end(ran)<max(end(ref_over))){
                                ran$spl_type<-"CDS_spanning;down_5ss;same_3ss"
                                if(f==1){ran$spl_type<-"CDS_spanning;down_5ss;same_lastCDS"}
                                if(f==last_ex){ran$spl_type<-"CDS_spanning;down_firstCDS;same_3ss"}
                                if(1==last_ex){ran$spl_type<-"CDS_spanning;down_5monoCDS;same_3monoCDS"}
                            }
                            
                        }
                        if(end(ran)==max(end(ref_over))){
                            if(start(ran)>min(start(ref_over))){
                                ran$spl_type<-"CDS_spanning;same_5ss;up_3ss"
                                if(f==1){ran$spl_type<-"CDS_spanning;same_5ss;up_lastCDS"}
                                if(f==last_ex){ran$spl_type<-"CDS_spanning;same_firstCDS;up_3ss"}
                                if(1==last_ex){ran$spl_type<-"CDS_spanning;same_5monoCDS;up_3monoCDS"}
                            }
                            if(start(ran)<min(start(ref_over))){
                                ran$spl_type<-"CDS_spanning;same_5ss;down_3ss"
                                if(f==1){ran$spl_type<-"CDS_spanning;same_5ss;down_lastCDS"}
                                if(f==last_ex){ran$spl_type<-"CDS_spanning;same_firstCDS;down_3ss"}
                                if(1==last_ex){ran$spl_type<-"CDS_spanning;same_5monoCDS;down_3monoCDS"}
                            }
                            
                        }
                        
                        if(end(ran)>max(end(ref_over))){
                            if(start(ran)>min(start(ref_over))){
                                ran$spl_type<-"CDS_spanning;up_5ss;up_3ss"
                                if(f==1){ran$spl_type<-"CDS_spanning;up_5ss;up_lastCDS"}
                                if(f==last_ex){ran$spl_type<-"CDS_spanning;up_firstCDS;up_3ss"}
                                if(1==last_ex){ran$spl_type<-"CDS_spanning;up_5monoCDS;up_3monoCDS"}
                            }
                            if(start(ran)<min(start(ref_over))){
                                ran$spl_type<-"CDS_spanning;up_5ss;down_3ss"
                                if(f==1){ran$spl_type<-"CDS_spanning;up_5ss;down_lastCDS"}
                                if(f==last_ex){ran$spl_type<-"CDS_spanning;up_firstCDS;down_3ss"}
                                if(1==last_ex){ran$spl_type<-"CDS_spanning;up_5monoCDS;down_3monoCDS"}
                            }
                            
                        }
                        
                        if(end(ran)<max(end(ref_over))){
                            if(start(ran)>min(start(ref_over))){
                                ran$spl_type<-"CDS_spanning;down_5ss;up_3ss"
                                if(f==1){ran$spl_type<-"CDS_spanning;down_5ss;up_lastCDS"}
                                if(f==last_ex){ran$spl_type<-"CDS_spanning;down_firstCDS;up_3ss"}
                                if(1==last_ex){ran$spl_type<-"CDS_spanning;down_5monoCDS;up_3monoCDS"}
                            }
                            if(start(ran)<min(start(ref_over))){
                                ran$spl_type<-"CDS_spanning;down_5ss;down_3ss"
                                if(f==1){ran$spl_type<-"CDS_spanning;down_5ss;down_lastCDS"}
                                if(f==last_ex){ran$spl_type<-"CDS_spanning;down_firstCDS;down_3ss"}
                                if(1==last_ex){ran$spl_type<-"CDS_spanning;down_5monoCDS;down_3monoCDS"}
                            }
                            
                        }
                    }
                    
                    
                }
                if(length(ref_over)==1){
                    ran$ref<-GRangesList(ref_over)
                    ran$spl_type<-NA
                    if(as.vector(strand(orf_gen[1]))=="+"){
                        if(start(ran)==(start(ref_over))){
                            if(end(ran)==(end(ref_over))){
                                ran$spl_type<-"same_5ss;same_3ss"
                                if(f==last_ex){ran$spl_type<-"same_5ss;same_lastCDS"}
                                if(f==1){ran$spl_type<-"same_firstCDS;same_3ss"}
                                if(1==last_ex){ran$spl_type<-"same_5monoCDS;same_3monoCDS"}
                            }
                            if(end(ran)>(end(ref_over))){
                                ran$spl_type<-"same_5ss;down_3ss"
                                if(f==last_ex){ran$spl_type<-"same_5ss;down_lastCDS"}
                                if(f==1){ran$spl_type<-"same_firstCDS;down_3ss"}
                                if(1==last_ex){ran$spl_type<-"same_5monoCDS;down_3monoCDS"}
                            }
                            if(end(ran)<(end(ref_over))){
                                ran$spl_type<-"same_5ss;up_3ss"
                                if(f==last_ex){ran$spl_type<-"same_5ss;up_lastCDS"}
                                if(f==1){ran$spl_type<-"same_firstCDS;up_3ss"}
                                if(1==last_ex){ran$spl_type<-"same_5monoCDS;up_3monoCDS"}
                            }
                            
                        }
                        if(end(ran)==(end(ref_over))){
                            if(start(ran)>(start(ref_over))){
                                ran$spl_type<-"down_5ss;same_3ss"
                                if(f==last_ex){ran$spl_type<-"down_5ss;same_lastCDS"}
                                if(f==1){ran$spl_type<-"down_firstCDS;same_3ss"}
                                if(1==last_ex){ran$spl_type<-"down_5monoCDS;same_3monoCDS"}
                            }
                            if(start(ran)<(start(ref_over))){
                                ran$spl_type<-"up_5ss;same_3ss"
                                if(f==last_ex){ran$spl_type<-"up_5ss;same_lastCDS"}
                                if(f==1){ran$spl_type<-"up_firstCDS;same_3ss"}
                                if(1==last_ex){ran$spl_type<-"up_5monoCDS;same_3monoCDS"}
                            }
                            
                        }
                        
                        if(end(ran)>(end(ref_over))){
                            if(start(ran)>(start(ref_over))){
                                ran$spl_type<-"down_5ss;down_3ss"
                                if(f==last_ex){ran$spl_type<-"down_5ss;down_lastCDS"}
                                if(f==1){ran$spl_type<-"down_firstCDS;down_3ss"}
                                if(1==last_ex){ran$spl_type<-"down_5monoCDS;down_3monoCDS"}
                            }
                            if(start(ran)<(start(ref_over))){
                                ran$spl_type<-"up_5ss;down_3ss"
                                if(f==last_ex){ran$spl_type<-"up_5ss;down_lastCDS"}
                                if(f==1){ran$spl_type<-"up_firstCDS;down_3ss"}
                                if(1==last_ex){ran$spl_type<-"up_5monoCDS;down_3monoCDS"}
                            }
                            
                        }
                        
                        if(end(ran)<(end(ref_over))){
                            if(start(ran)>(start(ref_over))){
                                ran$spl_type<-"down_5ss;up_3ss"
                                if(f==last_ex){ran$spl_type<-"down_5ss;up_lastCDS"}
                                if(f==1){ran$spl_type<-"down_firstCDS;up_3ss"}
                                if(1==last_ex){ran$spl_type<-"down_5monoCDS;up_3monoCDS"}
                            }
                            if(start(ran)<(start(ref_over))){
                                ran$spl_type<-"up_5ss;up_3ss"
                                if(f==last_ex){ran$spl_type<-"up_5ss;up_lastCDS"}
                                if(f==1){ran$spl_type<-"up_firstCDS;up_3ss"}
                                if(1==last_ex){ran$spl_type<-"up_5monoCDS;up_3monoCDS"}
                            }
                            
                        }
                        
                        
                        
                        
                    }
                    
                    if(as.vector(strand(orf_gen[1]))=="-"){
                        if(start(ran)==(start(ref_over))){
                            if(end(ran)==(end(ref_over))){
                                ran$spl_type<-"same_5ss;same_3ss"
                                if(f==1){ran$spl_type<-"same_5ss;same_lastCDS"}
                                if(f==last_ex){ran$spl_type<-"same_firstCDS;same_3ss"}
                                if(1==last_ex){ran$spl_type<-"same_5monoCDS;same_3monoCDS"}
                            }
                            if(end(ran)>(end(ref_over))){
                                ran$spl_type<-"up_5ss;same_3ss"
                                if(f==1){ran$spl_type<-"up_5ss;same_lastCDS"}
                                if(f==last_ex){ran$spl_type<-"up_firstCDS;same_3ss"}
                                if(1==last_ex){ran$spl_type<-"up_5monoCDS;same_3monoCDS"}
                            }
                            if(end(ran)<(end(ref_over))){
                                ran$spl_type<-"down_5ss;same_3ss"
                                if(f==1){ran$spl_type<-"down_5ss;same_lastCDS"}
                                if(f==last_ex){ran$spl_type<-"down_firstCDS;same_3ss"}
                                if(1==last_ex){ran$spl_type<-"down_5monoCDS;same_3monoCDS"}
                            }
                            
                        }
                        if(end(ran)==(end(ref_over))){
                            if(start(ran)>(start(ref_over))){
                                ran$spl_type<-"same_5ss;up_3ss"
                                if(f==1){ran$spl_type<-"same_5ss;up_lastCDS"}
                                if(f==last_ex){ran$spl_type<-"same_firstCDS;up_3ss"}
                                if(1==last_ex){ran$spl_type<-"same_5monoCDS;up_3monoCDS"}
                            }
                            if(start(ran)<(start(ref_over))){
                                ran$spl_type<-"same_5ss;down_3ss"
                                if(f==1){ran$spl_type<-"same_5ss;down_lastCDS"}
                                if(f==last_ex){ran$spl_type<-"same_firstCDS;down_3ss"}
                                if(1==last_ex){ran$spl_type<-"same_5monoCDS;down_3monoCDS"}
                            }
                            
                        }
                        
                        if(end(ran)>(end(ref_over))){
                            if(start(ran)>(start(ref_over))){
                                ran$spl_type<-"up_5ss;up_3ss"
                                if(f==1){ran$spl_type<-"up_5ss;up_lastCDS"}
                                if(f==last_ex){ran$spl_type<-"up_firstCDS;up_3ss"}
                                if(1==last_ex){ran$spl_type<-"up_5monoCDS;up_3monoCDS"}
                            }
                            if(start(ran)<(start(ref_over))){
                                ran$spl_type<-"up_5ss;down_3ss"
                                if(f==1){ran$spl_type<-"up_5ss;down_lastCDS"}
                                if(f==last_ex){ran$spl_type<-"up_firstCDS;down_3ss"}
                                if(1==last_ex){ran$spl_type<-"up_5monoCDS;down_3monoCDS"}
                            }
                            
                        }
                        
                        if(end(ran)<(end(ref_over))){
                            if(start(ran)>(start(ref_over))){
                                ran$spl_type<-"down_5ss;up_3ss"
                                if(f==1){ran$spl_type<-"down_5ss;up_lastCDS"}
                                if(f==last_ex){ran$spl_type<-"down_firstCDS;up_3ss"}
                                if(1==last_ex){ran$spl_type<-"down_5monoCDS;up_3monoCDS"}
                            }
                            if(start(ran)<(start(ref_over))){
                                ran$spl_type<-"down_5ss;down_3ss"
                                if(f==1){ran$spl_type<-"down_5ss;down_lastCDS"}
                                if(f==last_ex){ran$spl_type<-"down_firstCDS;down_3ss"}
                                if(1==last_ex){ran$spl_type<-"down_5monoCDS;down_3monoCDS"}
                            }
                            
                        }
                    }
                    
                    
                    
                    
                }
                
            }
            if(overref[f]==F){
                if(length(ref_cds)>0){ran$ref<-GRangesList(ref_cds[nearest(x=ran,subject=ref_cds)])}
                if(length(ref_cds)==0){ran$ref<-GRangesList(GRanges())}
                ran$spl_type<-"new_CDS"
                if(f==last_ex){
                    if(as.vector(strand(orf_gen[1]))=="+"){
                        ran$spl_type<-"new_lastCDS"
                    }
                    if(as.vector(strand(orf_gen[1]))=="-"){
                        ran$spl_type<-"new_firstCDS"
                    }
                }
                if(f==1){
                    
                    if(as.vector(strand(orf_gen[1]))=="-"){
                        ran$spl_type<-"new_lastCDS"
                    }
                    if(as.vector(strand(orf_gen[1]))=="+"){
                        ran$spl_type<-"new_firstCDS"
                    }                                       
                }
                if(1==last_ex){
                    ran$spl_type<-"new_monoCDS"
                    
                }
            }
            spl_ran<-sort(c(spl_ran,ran))
            
            
        }
        
    }
    
    newspl<-spl_ran$spl_type
    newspl[grep(spl_ran$spl_type,pattern = "new")]<-"new_miss"
    newspl[grep(spl_ran$spl_type,pattern = "missing")]<-"new_miss"
    rlesp<-Rle(newspl)
    
    #change new and missing
    if(runValue(rlesp)[1]=="new_miss"){
        if(as.character(strand(spl_ran)[1])=="+"){spl_ran$spl_type[1:runLength(rlesp)[1]]<-paste(spl_ran$spl_type[1:runLength(rlesp)[1]],"_5prime",sep = "")}
        if(as.character(strand(spl_ran)[1])=="-"){spl_ran$spl_type[1:runLength(rlesp)[1]]<-paste(spl_ran$spl_type[1:runLength(rlesp)[1]],"_3prime",sep = "")}
    }
    lenna<-length(runValue(rlesp))
    if(runValue(rlesp)[lenna]=="new_miss"){
        if(as.character(strand(spl_ran)[1])=="+"){spl_ran$spl_type[(length(spl_ran)-(runLength(rlesp)[lenna]-1)):length(spl_ran)]<-paste(spl_ran$spl_type[(length(spl_ran)-(runLength(rlesp)[lenna]-1)):length(spl_ran)],"_3prime",sep = "")}
        if(as.character(strand(spl_ran)[1])=="-"){paste(spl_ran$spl_type[(length(spl_ran)-(runLength(rlesp)[lenna]-1)):length(spl_ran)],"_5prime",sep = "")}
    }
    spl_ran$spl_type<-gsub(spl_ran$spl_type,pattern = "_5prime_3prime",replacement = "_notoverl")
    spl_ran
}

#' Annotate detected ORFs in transcript and genome space
#'
#' This function annotates quantified ORFs with respect to other detected ORFs and annotated ones, in both genome and transcript space.
#' @details As multiple transcripts can contain the same ORF,
#' all the transcript and transcript biotypes are indicated, with a preference for protein_coding transcripts in the "compatible" 
#' columns (to be conservative when assessing translation of non-protein coding transcripts). Such compatibility is also output considering the most upstream
#' start codon for that ORF. \cr
#' Splice features of each orf is annotated with respect to the longest coding transcripts and to the highest translated ORF in that gene.\cr
#' Variants in N or C terminus of the translated proteins are also indicated (Beta).\cr
#' ORF annotation with respect to the annotated transcript is also indicated, as follows:\cr\cr
#' \code{novel}: no ORF annotated in the transcript.\cr
#' \code{ORF_annotated}: same exact ORF as annotated.\cr
#' \code{N_extension}: N terminal extension.\cr
#' \code{N_truncation}: N terminal extension.\cr
#' \code{uORF}: upstream ORF.\cr
#' \code{overl_uORF}: upstream overlappin uORF.\cr
#' \code{NC_extension}: N and C termini extension.\cr
#' \code{dORF}: downstream ORF.\cr
#' \code{overl_dORF}: downstream overlapping ORF.\cr
#' \code{nested_ORF}: nested ORF.\cr
#' \code{C_truncation}: C terminal truncation.\cr
#' \code{C_extension}: C terminal extension.\cr\cr
#' As transcipt-specific annotation can be misleading due to a plethora of different transcripts, it is important to distinguish ORFs
#' also on the basis of their overlap with know CDS regions.
#' ORF annotation with respect to the entire set of CDS exon for the analyzed genomic regions is  indicated as follows:\cr\cr
#' \code{novel}: No CDS region is annotated in the entire region.\cr
#' \code{novel_Upstream}: ORF is upstream of annotated CDS regions (does not overlap).\cr
#' \code{novel_Downstream}: ORF is downstream of annotated CDS regions (does not overlap).\cr
#' \code{novel_Internal}: genomic location of the ORF is present between the start of the first,
#'  and the end of the last CDS region (does not overlap).\cr
#' \code{exact_start_stop}: Same start and end locations.\cr
#' \code{Alt5_start}: Different start region, upstream.\cr
#' \code{Alt3_start}: Different start region, downstream.\cr
#' \code{Alt5_stop}: Different end region, upstream.\cr
#' \code{Alt3_stop}: Different end region, downstream.\cr\cr
#' Another layer of annotation is performed by checking the position of the ORF stop codon 
#' with respect to the last exon-exon junction.
#' @keywords ORFquant
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param results_ORFs Full list of detected ORFs, from \code{select_quantify_ORFs}
#' @param Annotation Rannot object containing annotation of CDS and transcript structures (see \code{prepare_annotation_files}
#' @param genome_sequence BSgenome object
#' @param region genomic region being analyzed
#' @param genetic_code GENETIC_CODE table to use
#' @return Exon structure of detected ORF including possible missing exons from reference, together with a \code{spl_type} column
#' including the annotation for each exon (e.g. alternative acceptors or donor).\cr\cr
#' Additional columns are added to the ORFs_tx object:\cr
#' \code{compatible_with}: Set of transcript ids possibly containing the entire ORF structure.\cr
#' \code{compatible_biotype}: Compatible transcript biotype; if a protein coding transcript can contain 
#' the ORF, this is set to protein_coding.\cr
#' \code{compatible_tx}: One selected compatible transcript (preference if protein_coding).\cr
#' \code{compatible_ORF_id_tr}: ORF_id_tr id if selecting the compatible transcript.\cr
#' \code{compatible_with_longest}: Same as \code{compatible_with} but using the most upstream start codon.\cr
#' \code{compatible_ORF_id_tr_longest}: Same as \code{compatible_ORF_id_tr} but using the most upstream start codon .\cr
#' \code{ref_id}: transcript_id of the transcript used to annotate splicing (longest) .\cr
#' \code{ref_id_maxORF}: ORF_id_tr of the ORF used to annotated splicing (most translated of the gene).\cr
#' \code{NC_protein_isoform}: Annotation of possible N or C termini variant (when transcript is protein_coding) .\cr
#' \code{ORF_category_Tx}: ORF annotation with respect to ORF position in the transcript .\cr
#' \code{ORF_category_Tx_compatible}: ORF annotation with respect to ORF position in the transcript, using the \code{compatible_ORF_id_tr} .\cr
#' \code{ORF_category_Gen}: ORF annotation with respect to its genomic position .\cr
#' \code{NMD_candidate}: TRUE or FALSE, depending on the presence of an additional exon-exon junction downstream the stop codon.\cr
#' \code{NMD_candidate_compatible_txs}: same as NMD_candidate, but for all transcripts compatible with the ORF structure.\cr
#' \code{Distance_to_lastExEx}: Distance (in nt) between the last exon-exon junction and the stop codon.\cr
#' \code{Distance_to_lastExEx_compatible_txs}: same as  Distance_to_lastExEx, but for all transcripts compatible with the ORF structure.
#' @seealso \code{\link{select_quantify_ORFs}}, \code{\link{annotate_splicing}}
#' @export

annotate_ORFs<-function(results_ORFs,Annotation,genome_sequence,region,genetic_code){
    
    annotated_cds<-Annotation$cds_genes[Annotation$cds_genes%over%region]
    annotated_cds_tx<-Annotation$cds_txs[Annotation$cds_txs%over%region]
    annotated_cds_tx_genes<-Annotation$trann$gene_id[match(names(annotated_cds_tx),Annotation$trann$transcript_id)]
    annotated_exons_tx<-Annotation$exons_txs[Annotation$exons_txs%over%region]
    
    ORFs_tx<-results_ORFs$ORFs_tx_position
    ORFs_gen<-results_ORFs$ORFs_genomic_position
    ORFs_splice_feats<-GRangesList()
    ORFs_splice_feats_tomaxORF<-GRangesList()
    orfssss_tx<-unlist(GRangesList(ORFs_tx))
    orfssss_tx<-orfssss_tx[!is.na(orfssss_tx$ORF_pct_P_sites)]
    maxORF_orf<-c()
    max_cds<-c()
    max_cdsok<-GRanges()
    
    for(i in unique(orfssss_tx$gene_id)){
        okorfsss<-orfssss_tx[orfssss_tx$gene_id==i]
        maxORF_orf[i]<-okorfsss$ORF_id_tr[which.max(okorfsss$ORF_pct_P_sites)]
        anncdss<-annotated_cds_tx[annotated_cds_tx_genes==i]
        if(length(anncdss)>0){
            max_cds[i]<-names(which.max(sum(width(anncdss))))
        }
        
        if(length(anncdss)==0){
            orf_gen<-reduce(unlist(ORFs_gen[okorfsss$ORF_id_tr]))
            
            moret<-sapply(annotated_cds_tx,FUN=function(x){
                sum(width(setdiff(orf_gen,x)))
            })
            moret<-t(data.frame(moret,stringsAsFactors=F))
            max_cds[i]<-colnames(moret)[which.min(colSums(moret))]
        }
        
        
    }
    
    #compatibilities
    
    for(i in names(ORFs_gen)){
        x<-ORFs_gen[[i]]
        comp<-c()
        for(jjj in names(annotated_exons_tx)){
            mapp<-reduce(mapToTranscripts(x,transcripts = annotated_exons_tx[jjj]))
            if(sum(width(mapp))==sum(width(x)) & length(mapp)==1){
                comp<-c(comp,paste(seqnames(mapp),start(mapp),end(mapp),sep="_"))
            }
        }
        names(comp)<-NULL
        ORFs_tx[[i]]$compatible_with<-NULL
        ORFs_tx[[i]]$compatible_with<-CharacterList(unlist(comp))
        
        
        ORFs_tx[[i]]$compatible_biotype<-ORFs_tx[[i]]$transcript_biotype
        ORFs_tx[[i]]$compatible_tx<-ORFs_tx[[i]]$transcript_id
        ORFs_tx[[i]]$compatible_ORF_id_tr<-ORFs_tx[[i]]$ORF_id_tr
        compats<-elementNROWS(ORFs_tx[[i]]$compatible_with)>1
        comp_txs<-ORFs_tx[[i]]$compatible_with[compats]
        if(length(comp_txs)>0){
            compid<-sapply(comp_txs,function(x){
                txs<-sapply(strsplit(x,split = "_"),function(x){paste(x[-((length(x)-1):length(x))],collapse="_")})
                btps<-Annotation$trann$transcript_biotype[match(txs,Annotation$trann$transcript_id)]
                btps[is.na(btps)]<-"Not_found"
                pcd<-btps=="protein_coding"
                if(sum(pcd,na.rm = T)>0){c(sort(x[pcd])[1],sort(txs[pcd])[1],"protein_coding")}else{c(x[1],txs[1],btps[1])}
            })
            ORFs_tx[[i]]$compatible_ORF_id_tr[compats]<-t(compid)[,1]
            ORFs_tx[[i]]$compatible_tx[compats]<-t(compid)[,2]
            ORFs_tx[[i]]$compatible_biotype[compats]<-t(compid)[,3]
            
        }
        ORFs_tx[[i]]$compatible_with_longest<-ORFs_tx[[i]]$compatible_with
        ORFs_tx[[i]]$compatible_biotype_longest<-ORFs_tx[[i]]$compatible_biotype
        ORFs_tx[[i]]$compatible_tx_longest<-ORFs_tx[[i]]$compatible_tx
        ORFs_tx[[i]]$compatible_ORF_id_tr_longest<-ORFs_tx[[i]]$compatible_ORF_id_tr
        lng<-ORFs_tx[[i]]$longest_ORF
        lng$ORF_id_tr<-paste(seqnames(lng),start(lng),end(lng),sep="_")
        if(start(ORFs_tx[[i]])!=start(lng)){
            x<-from_tx_togen(ORFs = lng,exons = Annotation$exons_txs[ORFs_tx[[i]]$transcript_id],introns = Annotation$introns_txs[[ORFs_tx[[i]]$transcript_id]])[[1]]
            mapp<-mapToTranscripts(x,transcripts = annotated_exons_tx)
            redmapp<-reduce(split(mapp,seqnames(mapp)))
            
            redmapp<-redmapp[which(sum(width(redmapp))==sum(width(x)))]
            redmapp<-redmapp[elementNROWS(redmapp)==1]
            
            comp<-sapply(redmapp,function(x){paste(seqnames(x),start(x),end(x),sep="_")})
            names(comp)<-NULL
            ORFs_tx[[i]]$compatible_with_longest<-CharacterList(unlist(comp))
            ORFs_tx[[i]]$compatible_biotype_longest<-ORFs_tx[[i]]$transcript_biotype
            ORFs_tx[[i]]$compatible_tx_longest<-ORFs_tx[[i]]$transcript_id
            ORFs_tx[[i]]$compatible_ORF_id_tr_longest<-paste(as.character(seqnames(lng)[1]),start(lng),end(lng),sep = "_")
            
            compats<-elementNROWS(ORFs_tx[[i]]$compatible_with)>1
            comp_txs<-ORFs_tx[[i]]$compatible_with_longest[compats]
            if(length(comp_txs)>0){
                compid<-sapply(comp_txs,function(x){
                    txs<-sapply(strsplit(x,split = "_"),function(x){paste(x[-((length(x)-1):length(x))],collapse="_")})
                    btps<-Annotation$trann$transcript_biotype[match(txs,Annotation$trann$transcript_id)]
                    btps[is.na(btps)]<-"Not_found"
                    pcd<-btps=="protein_coding"
                    if(sum(pcd,na.rm = T)>0){c(sort(x[pcd])[1],sort(txs[pcd])[1],"protein_coding")}else{c(x[1],txs[1],btps[1])}
                })
                ORFs_tx[[i]]$compatible_ORF_id_tr_longest[compats]<-t(compid)[,1]
                ORFs_tx[[i]]$compatible_tx_longest[compats]<-t(compid)[,2]
                ORFs_tx[[i]]$compatible_biotype_longest[compats]<-t(compid)[,3]
                
            }        
            
        }
        
    }
    
    for(i in names(ORFs_gen)){
        
        compss<-ORFs_tx[[i]]$compatible_with[[1]]
        compss_txs<-sapply(strsplit(compss,split = "_"),function(x){paste(x[-((length(x)-1):length(x))],collapse="_")})
        ok_id<-compss[compss_txs==ORFs_tx[[i]]$compatible_tx_longest]
        comp_ln<-ORFs_tx[[i]]$compatible_biotype_longest
        comp_prev<-ORFs_tx[[i]]$compatible_biotype
        if(comp_ln=="protein_coding" | (comp_prev!="protein_coding" & comp_ln!="protein_coding") ){
            ORFs_tx[[i]]$compatible_ORF_id_tr<-ok_id
            ORFs_tx[[i]]$compatible_tx<-ORFs_tx[[i]]$compatible_tx_longest
            ORFs_tx[[i]]$compatible_biotype<-ORFs_tx[[i]]$compatible_biotype_longest
        }
        ORFs_tx[[i]]$compatible_tx_longest<-NULL
        ORFs_tx[[i]]$compatible_biotype_longest<-NULL
        
    }
    
    #annotate NMD candidates
    exs<-annotated_exons_tx[unique(unlist(GRangesList(ORFs_tx))$transcript_id)]
    
    strands_exs<-sapply(strand(exs),function(x){x@values[1]})
    exs_pos<-exs[strands_exs=="+"]
    exs_neg<-exs[strands_exs=="-"]
    
    
    last_ex_pos<-which.max(start(exs_pos))
    last_ex_pos<-exs_pos[splitAsList(unname(last_ex_pos), names(last_ex_pos))]
    last_ex_pos<-last_ex_pos[match(names(exs_pos),names(last_ex_pos))]
    txs_pos<-pmapToTranscripts(last_ex_pos,transcripts = exs_pos)
    
    
    last_ex_neg<-which.min(start(exs_neg))
    last_ex_neg<-exs_neg[splitAsList(unname(last_ex_neg), names(last_ex_neg))]
    last_ex_neg<-last_ex_neg[match(names(exs_neg),names(last_ex_neg))]
    txs_neg<-pmapToTranscripts(last_ex_neg,transcripts = exs_neg)
    txs_all<-unlist(c(txs_pos,txs_neg))
    last_exexs<-start(txs_all)[match(unlist(GRangesList(ORFs_tx))$transcript_id,names(txs_all))]
    Distance_EJCs<-last_exexs-end(unlist(GRangesList(ORFs_tx)))
    
    #NMD_compat
    unlORFs<-unlist(GRangesList(ORFs_tx))
    unlORFs_comp<-sapply(strsplit(unlist(unlORFs$compatible_with),"_"),function(x){paste(x[-((length(x)-1):length(x))],collapse="_")})
    exs<-annotated_exons_tx[unlORFs_comp]
    
    strands_exs<-sapply(strand(exs),function(x){x@values[1]})
    exs_pos<-exs[strands_exs=="+"]
    exs_neg<-exs[strands_exs=="-"]
    
    last_ex_pos<-which.max(start(exs_pos))
    last_ex_pos<-exs_pos[splitAsList(unname(last_ex_pos), names(last_ex_pos))]
    last_ex_pos<-last_ex_pos[match(names(exs_pos),names(last_ex_pos))]
    txs_pos<-pmapToTranscripts(last_ex_pos,transcripts = exs_pos)
    
    
    last_ex_neg<-which.min(start(exs_neg))
    last_ex_neg<-exs_neg[splitAsList(unname(last_ex_neg), names(last_ex_neg))]
    last_ex_neg<-last_ex_neg[match(names(exs_neg),names(last_ex_neg))]
    txs_neg<-pmapToTranscripts(last_ex_neg,transcripts = exs_neg)
    txs_all<-unlist(c(txs_pos,txs_neg)[unlORFs_comp])
    
    last_exexs<-start(txs_all)
    
    endcompat<-as.numeric(unlist(lapply(strsplit(unlist(unlORFs$compatible_with),"_"),function(x){x[length(x)]})))
    
    Distance_EJCs_compat<-last_exexs-endcompat
    nrws<-elementNROWS(unlORFs$compatible_with)
    fcttr<-rep(1:length(nrws),nrws)
    Distance_EJCs_compat<-CharacterList(unname(split(Distance_EJCs_compat,fcttr)))
    
    ORFs_tx<-lapply(ORFs_tx,function(x){
        cols<-mcols(x)
        cols[,c("ref_id","ref_id_maxORF","NC_protein_isoform","ORF_category_Tx","ORF_category_Tx_compatible","ORF_category_Gen",
                "NMD_candidate","Distance_to_lastExEx")]<-NA
        cols[,c("NMD_candidate_compatible_txs","Distance_to_lastExEx_compatible_txs")]<-CharacterList("")
        mcols(x)<-cols
        x
    })
    
    for(i in 1:length(ORFs_tx)){
        nmd<-FALSE
        if(Distance_EJCs[i]>0){
            nmd<-TRUE
        }
        ORFs_tx[[i]]$NMD_candidate<-nmd
        ORFs_tx[[i]]$Distance_to_lastExEx<-Distance_EJCs[i]
        
        nmd<-Distance_EJCs_compat[i]>0
        ORFs_tx[[i]]$NMD_candidate_compatible_txs<-nmd
        ORFs_tx[[i]]$Distance_to_lastExEx_compatible_txs<-Distance_EJCs_compat[i]
    }
    
    
    
    annotated_cds_tx2<-annotated_cds_tx
    
    #here bulk of work
    
    for(i in 1:length(ORFs_tx)){
        orf_tx<-ORFs_tx[[i]]
        ORFs_splice_feats[[orf_tx$ORF_id_tr]]<-GRanges()
        ORFs_splice_feats_tomaxORF[[orf_tx$ORF_id_tr]]<-GRanges()
        
        #annotate tx position
        
        annotated_ORF<-Annotation$cds_txs_coords[as.character(seqnames(Annotation$cds_txs_coords))==orf_tx$transcript_id]
        annotated_ORF_compatible<-Annotation$cds_txs_coords[as.character(seqnames(Annotation$cds_txs_coords))==orf_tx$compatible_tx]
        
        if(length(annotated_ORF)==0){
            ORFs_tx[[i]]$ORF_category_Tx<-"novel"
        }
        if(length(annotated_ORF_compatible)==0){
            ORFs_tx[[i]]$ORF_category_Tx_compatible<-"novel"
        }
        
        if(length(annotated_ORF)>0){
            
            ann_sta<-start(annotated_ORF)
            ann_sto<-end(annotated_ORF)-3
            sta<-start(orf_tx)
            sto<-end(orf_tx)
            if(sto==ann_sto){
                
                if(sta==ann_sta){ORFs_tx[[i]]$ORF_category_Tx<-"ORF_annotated"}
                if(sta<ann_sta){ORFs_tx[[i]]$ORF_category_Tx<-"N_extension"}
                if(sta>ann_sta){ORFs_tx[[i]]$ORF_category_Tx<-"N_truncation"}
                
            }
            if(sto!=ann_sto){
                
                if(sta<ann_sta & sto<ann_sto){ORFs_tx[[i]]$ORF_category_Tx<-"overl_uORF"}
                if(sta<ann_sta & sto<ann_sta){ORFs_tx[[i]]$ORF_category_Tx<-"uORF"}
                if(sta<ann_sta & sto>ann_sto){ORFs_tx[[i]]$ORF_category_Tx<-"NC_extension"}
                if(sta>ann_sta & sto>ann_sto){ORFs_tx[[i]]$ORF_category_Tx<-"overl_dORF"}
                if(sta>ann_sto & sto>ann_sto){ORFs_tx[[i]]$ORF_category_Tx<-"dORF"}
                if(sta>ann_sta & sto<ann_sto){ORFs_tx[[i]]$ORF_category_Tx<-"nested_ORF"}
                if(sta==ann_sta & sto<ann_sto){ORFs_tx[[i]]$ORF_category_Tx<-"C_truncation"}
                if(sta==ann_sta & sto>ann_sto){ORFs_tx[[i]]$ORF_category_Tx<-"C_extension"}
                
            }
        }
        
        
        if(length(annotated_ORF_compatible)>0){
            
            ann_sta<-start(annotated_ORF_compatible)
            ann_sto<-end(annotated_ORF_compatible)-3
            #change
            sta<-as.numeric(sapply(strsplit(orf_tx$compatible_ORF_id_tr,"_"),function(x){x[length(x)-1]}))
            sto<-as.numeric(sapply(strsplit(orf_tx$compatible_ORF_id_tr,"_"),function(x){x[length(x)]}))
            if(sto==ann_sto){
                
                if(sta==ann_sta){ORFs_tx[[i]]$ORF_category_Tx_compatible<-"ORF_annotated"}
                if(sta<ann_sta){ORFs_tx[[i]]$ORF_category_Tx_compatible<-"N_extension"}
                if(sta>ann_sta){ORFs_tx[[i]]$ORF_category_Tx_compatible<-"N_truncation"}
                
            }
            if(sto!=ann_sto){
                
                if(sta<ann_sta & sto<ann_sto){ORFs_tx[[i]]$ORF_category_Tx_compatible<-"overl_uORF"}
                if(sta<ann_sta & sto<ann_sta){ORFs_tx[[i]]$ORF_category_Tx_compatible<-"uORF"}
                if(sta<ann_sta & sto>ann_sto){ORFs_tx[[i]]$ORF_category_Tx_compatible<-"NC_extension"}
                if(sta>ann_sta & sto>ann_sto){ORFs_tx[[i]]$ORF_category_Tx_compatible<-"overl_dORF"}
                if(sta>ann_sto & sto>ann_sto){ORFs_tx[[i]]$ORF_category_Tx_compatible<-"dORF"}
                if(sta>ann_sta & sto<ann_sto){ORFs_tx[[i]]$ORF_category_Tx_compatible<-"nested_ORF"}
                if(sta==ann_sta & sto<ann_sto){ORFs_tx[[i]]$ORF_category_Tx_compatible<-"C_truncation"}
                if(sta==ann_sta & sto>ann_sto){ORFs_tx[[i]]$ORF_category_Tx_compatible<-"C_extension"}
                
            }
        }
        
        
        #annotate splice, genomic position and protein termini based on max cds and max pct
        
        orf_gen<-ORFs_gen[[orf_tx$ORF_id_tr]]
        
        #1 of multiple overlapping cds per gene
        
        if(length(annotated_cds)>1){
            moreg<-sapply(annotated_cds,FUN=function(x){
                orf_gen%over%x
            })
            if(length(orf_gen)==1){moreg<-t(data.frame(moreg,stringsAsFactors=F))}
            annotated_cds2<-reduce(unlist(annotated_cds[[colnames(moreg)[which.max(colSums(moreg))]]]))
            
            #if(length(annotated_cds_tx2)>1){annotated_cds_tx<-annotated_cds_tx2[annotated_cds_tx2%over%annotated_cds2]}
            
        }
        
        #otherwise list with 1 gene
        if(length(annotated_cds)==1){annotated_cds2<-reduce(unlist(annotated_cds))}
        if(length(annotated_cds)==0){annotated_cds2<-annotated_cds}
        overl<-orf_gen%over%annotated_cds2
        
        
        if(sum(overl)==0){
            ORFs_tx[[i]]$ORF_category_Gen<-"novel"
            if(length(annotated_cds2)>0){
                nearest_cds<-annotated_cds[[names(unlist(annotated_cds))[nearest(orf_gen,unlist(annotated_cds))[1]]]]
                overl_whole<-orf_gen@ranges%over%IRanges(start=min(start(nearest_cds)),end=max(end(nearest_cds)))
                
                if(as.vector(strand(orf_gen[1]))=="+"){
                    
                    if(sum(overl_whole)==0){
                        if(min(start(orf_gen))<min(start(nearest_cds))){
                            ORFs_tx[[i]]$ORF_category_Gen<-"novel_Upstream"
                        }
                        if(min(start(orf_gen))>max(end(nearest_cds))){
                            ORFs_tx[[i]]$ORF_category_Gen<-"novel_Downstream"
                        }
                    }
                    if(sum(overl_whole)>0){
                        ORFs_tx[[i]]$ORF_category_Gen<-"novel_Internal"
                    }
                    
                }
                if(as.vector(strand(orf_gen[1]))=="-"){
                    if(sum(overl_whole)==0){
                        if(max(end(orf_gen))>max(end(nearest_cds))){
                            ORFs_tx[[i]]$ORF_category_Gen<-"novel_Upstream"
                        }
                        if(min(start(orf_gen))<min(start(nearest_cds))){
                            ORFs_tx[[i]]$ORF_category_Gen<-"novel_Downstream"
                        }
                    }
                    if(sum(overl_whole)>0){
                        ORFs_tx[[i]]$ORF_category_Gen<-"novel_Internal"
                    }                                        
                }
            }
        }
        
        #if overlaps CDS
        spl_ran<-GRanges()
        if(sum(overl)>0){
            
            #annotated NC wrt maxcds
            
            ORFs_tx[[i]]$ORF_category_Gen<-"overlaps_CDS"
            max_cdsok<-annotated_cds_tx[[max_cds[ORFs_tx[[i]]$gene_id]]]
            if(length(max_cdsok)==0){
                max_cdsok<-annotated_cds_tx[[which.max(sum(width(annotated_cds_tx)))]]
            }
            ORFs_tx[[i]]$ref_id<-max_cds[ORFs_tx[[i]]$gene_id]
            max_cds_seq<-unlist(getSeq(x=genome_sequence,max_cdsok))
            segm<-10
            if(length(max_cds_seq)<33 | nchar(orf_tx$Protein)<10){
                segm<-min(c(as.integer(length(max_cds_seq)/3)-1,nchar(orf_tx$Protein)))
            }
            N_pr<-AAString(unlist(orf_tx$Protein))[1:segm]
            C_pr<-AAString(unlist(orf_tx$Protein))[(nchar(unlist(orf_tx$Protein))-(segm-1)):nchar(unlist(orf_tx$Protein))]
            N_ann<-translate(max_cds_seq[1:(segm*3)],genetic.code = genetic_code,if.fuzzy.codon="solve")
            C_ann<-translate(head(tail(max_cds_seq,(segm*3+3)),(segm*3)),genetic.code = genetic_code,if.fuzzy.codon="solve")
            
            
            ORFs_tx[[i]]$NC_protein_isoform<-"N_C"
            if(N_pr==N_ann){
                if(C_pr==C_ann){
                    ORFs_tx[[i]]$NC_protein_isoform<-"same"
                } else {ORFs_tx[[i]]$NC_protein_isoform<-"C"}
            }
            if(C_pr==C_ann & N_pr!=N_ann){
                ORFs_tx[[i]]$NC_protein_isoform<-"N"                                        
            }
            
            #ann genomic
            
            if(as.vector(strand(orf_gen[1]))=="+"){
                
                gen_sta<-min(start(max_cdsok))
                gen_sto<-max(end(max_cdsok))
                
                sta_or<-min(start(orf_gen))
                sto_or<-max(end(orf_gen))+3
                if(sto_or==gen_sto){
                    if(sta_or==gen_sta){ORFs_tx[[i]]$ORF_category_Gen<-"exact_start_stop"}
                    if(sta_or<gen_sta){ORFs_tx[[i]]$ORF_category_Gen<-"Alt5_start"}
                    if(sta_or>gen_sta){ORFs_tx[[i]]$ORF_category_Gen<-"Alt3_start"}
                    
                }
                
                if(sto_or!=gen_sto){
                    if(sta_or<gen_sta & sto_or<gen_sto){ORFs_tx[[i]]$ORF_category_Gen<-"Alt5_start_Alt5_stop"}
                    if(sta_or<gen_sta & sto_or>gen_sto){ORFs_tx[[i]]$ORF_category_Gen<-"Alt5_start_Alt3_stop"}
                    if(sta_or>gen_sta & sto_or>gen_sto){ORFs_tx[[i]]$ORF_category_Gen<-"Alt3_start_Alt3_stop"}
                    if(sta_or>gen_sta & sto_or<gen_sto){ORFs_tx[[i]]$ORF_category_Gen<-"Alt3_start_Alt5_stop"}
                    if(sta_or==gen_sta & sto_or<gen_sto){ORFs_tx[[i]]$ORF_category_Gen<-"Alt5_stop"}
                    if(sta_or==gen_sta & sto_or>gen_sto){ORFs_tx[[i]]$ORF_category_Gen<-"Alt3_stop"}
                    
                }
            }
            if(as.vector(strand(orf_gen[1]))=="-"){
                
                gen_sta<-max(end(max_cdsok))
                gen_sto<-min(start(max_cdsok))
                
                sta_or<-max(end(orf_gen))
                sto_or<-min(start(orf_gen))-3
                
                if(sto_or==gen_sto){
                    if(sta_or==gen_sta){ORFs_tx[[i]]$ORF_category_Gen<-"exact_start_stop"}
                    if(sta_or<gen_sta){ORFs_tx[[i]]$ORF_category_Gen<-"Alt3_start"}
                    if(sta_or>gen_sta){ORFs_tx[[i]]$ORF_category_Gen<-"Alt5_start"}
                    
                }
                
                if(sto_or!=gen_sto){
                    if(sta_or>gen_sta & sto_or>gen_sto){ORFs_tx[[i]]$ORF_category_Gen<-"Alt5_start_Alt5_stop"}
                    
                    if(sta_or>gen_sta & sto_or<gen_sto){ORFs_tx[[i]]$ORF_category_Gen<-"Alt5_start_Alt3_stop"}
                    if(sta_or<gen_sta & sto_or<gen_sto){ORFs_tx[[i]]$ORF_category_Gen<-"Alt3_start_Alt3_stop"}
                    if(sta_or<gen_sta & sto_or>gen_sto){ORFs_tx[[i]]$ORF_category_Gen<-"Alt3_start_Alt5_stop"}
                    if(sta_or==gen_sta & sto_or>gen_sto){ORFs_tx[[i]]$ORF_category_Gen<-"Alt5_stop"}
                    if(sta_or==gen_sta & sto_or<gen_sto){ORFs_tx[[i]]$ORF_category_Gen<-"Alt3_stop"}
                    
                }
                
                
                
            }
        }     
        
        ORFs_splice_feats[[orf_tx$ORF_id_tr]]<-annotate_splicing(orf_gen = orf_gen,ref_cds = max_cdsok)
        #to maxORF
        max_pct<-ORFs_gen[[maxORF_orf[ORFs_tx[[i]]$gene_id]]]
        ORFs_tx[[i]]$ref_id_maxORF<-maxORF_orf[ORFs_tx[[i]]$gene_id]
        ORFs_splice_feats_tomaxORF[[orf_tx$ORF_id_tr]]<-annotate_splicing(orf_gen = orf_gen,ref_cds = max_pct)
        
        
        
    }
    
    results_ORFs$ORFs_tx_position<-ORFs_tx
    list_spl_res<-list(ORFs_splice_feats,ORFs_splice_feats_tomaxORF)
    names(list_spl_res)<-c("annotation_wrt_longest","annotation_wrt_maxORF")
    results_ORFs$ORFs_splice_feats<-list_spl_res
    return(results_ORFs)
}


#' Detection, quantification and annotation of translated ORFs in a genomic region
#'
#' This function detects, quantifies and annotates actively translated ORF in a genomic region
#' @details A set of transcripts, together with genome sequence and Ribo-signal are analyzed to extract translated ORFs
#' @keywords ORFquant
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param region GRanges object with genomic coordinates of the genomic region analyzed
#' @param for_ORFquant "for_ORFquant" Robject containing P_sites positions and junction reads
#' @param genetic_code_region GENETIC_CODE table to use
#' @param orf_find.all_starts \code{get_all_starts} parameter for the \code{detect_translated_orfs} function
#' @param orf_find.nostarts \code{Stop_Stop} parameter for the \code{detect_translated_orfs} function
#' @param orf_find.start_sel_cutoff \code{cutoff} parameter for the \code{detect_translated_orfs} function
#' @param orf_find.start_sel_cutoff_ave \code{cutoff_ave} parameter for the \code{detect_translated_orfs} function
#' @param orf_find.cutoff_fr_ave \code{cutoff} parameter for the \code{detect_translated_orfs} function
#' @param orf_quant.cutoff_cums \code{cutoff_cums} parameter for the \code{select_quantify_ORFs} function
#' @param orf_quant.cutoff_pct \code{cutoff_pct} parameter for the \code{select_quantify_ORFs} function
#' @param orf_quant.cutoff_P_sites \code{cutoff_P_sites} parameter for the \code{select_quantify_ORFs} function
#' @param unique_reads Use only signal from uniquely mapping reads? Defaults to \code{FALSE}.
#' @param orf_quant.scaling \code{scaling} parameter for the \code{select_quantify_ORFs} function. Defaults to total_Psites 
#' @return A list containing transcript coordinates, exonic coordinates and annotation for each ORF.\cr\cr
#' The description for each list object is as follows:\cr\cr
#' \code{ORFs_tx}: transcript coordinates of the detected ORFs.\cr
#' \code{ORFs_gen}: genomic (exon) coordinates of the detected ORFs.\cr
#' \code{ORFs_feat}: list of ORF features together with mapping reads and uniqueness.\cr
#' \code{ORFs_txs_feats}: list of transcript features present in the genomic region, together with mapping reads and uniqueness.\cr
#' \code{ORFs_spl_feat_longest}: splicing annotation for each ORF exon, with respect to the longest annotated coding transcript for each gene.\cr
#' \code{ORFs_spl_feat_maxORF}: splicing annotation for each ORF exon, with respect to the most translated ORF in each gene.\cr
#' \code{selected_txs}: character vector containing the transcript ids of the selected transcripts.\cr
#' \code{ORFs_readthroughs}: (Beta) transcript coordinates of the detected ORFs readthroughs.\cr
#' @seealso \code{\link{select_txs}}, \code{\link{detect_translated_orfs}}, \code{\link{select_quantify_ORFs}}, \code{\link{annotate_ORFs}}, \code{\link{detect_readthrough}}
#' @export

ORFquant<-function(region,for_ORFquant,genetic_code_region,
                 orf_find.all_starts=T,orf_find.nostarts=F,orf_find.start_sel_cutoff = NA,orf_find.start_sel_cutoff_ave = .5,
                 orf_find.cutoff_fr_ave=.5,orf_quant.cutoff_cums = NA,orf_quant.cutoff_pct = 2,orf_quant.cutoff_P_sites=NA,unique_reads=F,orf_quant.scaling="total_Psites"){
    if(!orf_quant.scaling%in%c("total_Psites","average_coverage")){stop(paste("orf_quant.scaling parameter must be either total_Psites (recommended) or average_coverage"),date())}
    
    P_sites_region<-for_ORFquant$P_sites_all[for_ORFquant$P_sites_all%over%region]
    P_sites_uniq_region<-for_ORFquant$P_sites_uniq[for_ORFquant$P_sites_uniq%over%region]
    P_sites_uniq_mm_region<-for_ORFquant$P_sites_uniq_mm[for_ORFquant$P_sites_uniq_mm%over%region]
    
    res_orfs<-list()
    minimum_reads<-length(P_sites_region)>4
    if(unique_reads){length(P_sites_uniq_region)>4}
    
    if(minimum_reads){
        
        selected_transcripts<-select_txs(region = region,P_sites = P_sites_region,P_sites_uniq = P_sites_uniq_region,annotation = GTF_annotation,junction_counts=for_ORFquant$junctions,uniq_signal = unique_reads)
        if(length(selected_transcripts)>0){
            res_orfs<-suppressWarnings(detect_translated_orfs(selected_txs = selected_transcripts,genome_sequence = genome_seq,annotation = GTF_annotation,
                                                              P_sites = P_sites_region,P_sites_uniq = P_sites_uniq_region,P_sites_uniq_mm = P_sites_uniq_mm_region,
                                                              genomic_region=region,genetic_code=genetic_code_region,
                                                              all_starts=orf_find.all_starts,nostarts=,orf_find.nostarts,
                                                              start_sel_cutoff=orf_find.start_sel_cutoff,
                                                              start_sel_cutoff_ave=orf_find.start_sel_cutoff_ave,
                                                              cutoff_fr_ave=orf_find.cutoff_fr_ave,uniq_signal = unique_reads))
        }
        if(length(res_orfs)>0){
            res_orfs<-select_quantify_ORFs(results_ORFs=res_orfs,P_sites = P_sites_region,P_sites_uniq = P_sites_uniq_region,
                                           cutoff_cums = orf_quant.cutoff_cums,cutoff_pct = orf_quant.cutoff_pct,cutoff_P_sites=orf_quant.cutoff_P_sites,uniq_signal = unique_reads,scaling = orf_quant.scaling)
        }
        if(length(res_orfs)>0){
            
            res_orfs<-annotate_ORFs(results_ORFs=res_orfs,Annotation=GTF_annotation,genome_sequence = genome_seq,region=region,genetic_code=genetic_code_region)
            res_orfs[["readthrough"]]<-suppressWarnings(detect_readthrough(results_orf = res_orfs,P_sites = P_sites_region,P_sites_uniq = P_sites_uniq_region,P_sites_uniq_mm = P_sites_uniq_mm_region,
                                                                           genome_sequence=genome_seq, annotation = GTF_annotation,genetic_code_table=genetic_code_region,cutoff_fr_ave=orf_find.cutoff_fr_ave,uniq_signal = unique_reads))
            
            
        }
        res_orfs[["genomic_features"]]<-selected_transcripts
        
    }
    
    return(res_orfs)
}

#' Run the ORFquant pipeline
#'
#' This wrapper function runs the entire ORFquant pipeline
#' @details A set of transcripts, together with genome sequence and Ribo-signal are analyzed to extract translated ORFs
#' @keywords ORFquant
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param annotation_file REQUIRED - path to the *Rannot R file in the annotation directory used in the \code{prepare_annotation_files function}
#' @param for_ORFquant_file REQUIRED - path to the "for_ORFquant" file containing P_sites positions and junction reads
#' @param n_cores REQUIRED - number of cores to use
#' @param prefix prefix to use for the output files. Defaults to same as \code{for_ORFquant_file} (appends to its filename)
#' @param gene_name \code{character} vector of gene names to analyze.
#' @param gene_id \code{character} vector of gene ids to analyze
#' @param genomic_region \code{GRanges} object with genomic regions to analyze
#' @param write_temp_files write temporary files. Defaults to \code{TRUE}
#' @param write_GTF_file write a GTF files with the ORF coordinates. Defaults to \code{TRUE}
#' @param write_protein_fasta write a protein fasta file. Defaults to \code{TRUE}
#' @param interactive should put R object in global environment? Defaults to \code{TRUE}
#' @param stn.orf_find.all_starts \code{orf_find.all_starts} parameter for the \code{ORFquant} function
#' @param stn.orf_find.nostarts \code{orf_find.nostarts} parameter for the \code{ORFquant} function
#' @param stn.orf_find.start_sel_cutoff \code{orf_find.start_sel_cutoff} parameter for the \code{ORFquant} function
#' @param stn.orf_find.start_sel_cutoff_ave \code{orf_find.start_sel_cutoff_ave} parameter for the \code{ORFquant} functio
#' @param stn.orf_find.cutoff_fr_ave \code{orf_find.cutoff_fr_ave} parameter for the \code{ORFquant} function
#' @param stn.orf_quant.cutoff_cums \code{orf_quant.cutoff_cums} parameter for the \code{ORFquant} function
#' @param stn.orf_quant.cutoff_pct \code{orf_quant.cutoff_pct} parameter for the \code{ORFquant} function
#' @param stn.orf_quant.cutoff_P_sites \code{orf_quant.cutoff_P_sites} parameter for the \code{ORFquant} function
#' @param unique_reads_only Use only signal from uniquely mapping reads? Defaults to \code{FALSE}.
#' @param stn.orf_quant.scaling \code{orf_quant.scaling} parameter for the \code{ORFquant} function. Defaults to total_Psites
#' @param canonical_start_only Use only the canonical start codon (no alternative initiation codons)? Defaults to \code{TRUE}.
#' @return A set of output files containing transcript coordinates, exonic coordinates and annotation for each ORF, including optional GTF and protein fasta files.\cr\cr
#' The description for each list object is as follows:\cr\cr
#' \code{tmp_ORFquant_results}: (Optional) RData object file containing the entire set of results for each genomic region.\cr
#' \code{final_ORFquant_results}: RData object file containing the final ORFquant results, see \code{ORFquant}.\cr
#' \code{Protein_sequences.fasta}: (Optional) Fasta file containing the set of translated proteins .\cr
#' \code{Detected_ORFs.gtf}: GTF file containing coordinates of the detected ORFs.\cr\cr
#' In addition, new columns are added in the ORFs_tx file:\cr\cr
#' \code{ORFs_pM}: number of P_sites for each ORF, divided by ORF length and summing up to a million (akin to TPM).\cr
#' @seealso \code{\link{prepare_annotation_files}}, \code{\link{load_annotation}}, \code{\link{ORFquant}}
#' @export

run_ORFquant<-function(for_ORFquant_file,annotation_file,n_cores,prefix=for_ORFquant_file,gene_name=NA,gene_id=NA,genomic_region=NA,
                     write_temp_files=T,write_GTF_file=T,write_protein_fasta=T,interactive=T,
                     stn.orf_find.all_starts=T,stn.orf_find.nostarts=F,stn.orf_find.start_sel_cutoff = NA,
                     stn.orf_find.start_sel_cutoff_ave = .5,stn.orf_find.cutoff_fr_ave=.5,
                     stn.orf_quant.cutoff_cums = NA,stn.orf_quant.cutoff_pct = 2,stn.orf_quant.cutoff_P_sites=NA,unique_reads_only=F,canonical_start_only=T,stn.orf_quant.scaling="total_Psites"){    
    
    if(!stn.orf_quant.scaling%in%c("total_Psites","average_coverage")){stop(paste("stn.orf_quant.scaling parameter must be either total_Psites (recommended) or average_coverage"),date())}
    
    if(n_cores>1){
        registerDoMC(n_cores)
    }
    
    for (f in c(for_ORFquant_file,annotation_file)){
        if(file.access(f, 0)==-1) {
            stop("The following files don't exist:\n",
                 f, "\n")
        }
    }
    
    cat(paste("Loading annotation and Ribo-seq signal ... ",date(),"\n",sep = ""))
    
    load_annotation(annotation_file)
    for_ORFquant_data<-get(load(for_ORFquant_file))
    
    genes_red<-reduce(unlist(GTF_annotation$txs_gene))
    
    if(!is.na(gene_name[1])){
        gnid<-unique(GTF_annotation$trann$gene_id[GTF_annotation$trann$gene_name%in%gene_name])
        genes_red<-genes_red[genes_red%over%GTF_annotation$genes[gnid]]
    }
    
    if(!is.na(gene_id[1])){
        genes_red<-genes_red[genes_red%over%GTF_annotation$genes[gene_id]]
    }
    
    if(!is.na(genomic_region[1])){
        genes_red<-genes_red[genes_red%over%genomic_region]
    }
    
    if(length(genes_red)==0){
        stop(paste("Incorrect gene_id | gene_name | genomic_region! ",date(),sep = ""))
    }
    
    cat(paste("Loading annotation and Ribo-seq signal --- Done! ",date(),"\n",sep = ""))
    
    ovs_genesred<-summarizeOverlaps(genes_red,reads = for_ORFquant_data$P_sites_all)
    
    genes_red<-genes_red[which(assay(ovs_genesred)>4)]
    
    if(length(genes_red)==0){
        stop(paste("Not enough P_sites signal over genomic regions! ",date(),sep = ""))
    }
    
    cat(paste("Summoning ORFquant with ", sum(for_ORFquant_data$P_sites_all%over%genes_red)," P_sites positions over ", length(genes_red), " genomic regions using ",n_cores," processor(s) ... ",date(),"\n",sep = ""))
    
    if(n_cores>1){
        
        ORFs_found<-foreach(g=1:length(genes_red),.packages='GenomicRanges') %dopar%{
            
            gen_region<-genes_red[g]
            genetcd<-GTF_annotation$genetic_codes$genetic_code[rownames(GTF_annotation$genetic_codes)==as.character(seqnames(gen_region))]
            genetcd<-getGeneticCode(genetcd)
            if(canonical_start_only){
                attributes(genetcd)$alt_init_codons<-names(which(genetcd=="M"))
            }
            
            
            ORFquant(region=gen_region,for_ORFquant=for_ORFquant_data,genetic_code_region=genetcd,
                   orf_find.all_starts=stn.orf_find.all_starts,orf_find.nostarts=stn.orf_find.nostarts,
                   orf_find.start_sel_cutoff = stn.orf_find.start_sel_cutoff,orf_find.start_sel_cutoff_ave = stn.orf_find.start_sel_cutoff_ave,
                   orf_find.cutoff_fr_ave=stn.orf_find.cutoff_fr_ave,orf_quant.cutoff_cums = stn.orf_quant.cutoff_cums,
                   orf_quant.cutoff_pct = stn.orf_quant.cutoff_pct,orf_quant.cutoff_P_sites=stn.orf_quant.cutoff_P_sites,unique_reads = unique_reads_only,orf_quant.scaling = stn.orf_quant.scaling)
        }
        
    }
    
    if(n_cores==1){
        ORFs_found<-list()
        for(g in 1:length(genes_red)){
            
            gen_region<-genes_red[g]
            genetcd<-GTF_annotation$genetic_codes$genetic_code[rownames(GTF_annotation$genetic_codes)==as.character(seqnames(gen_region))]
            genetcd<-getGeneticCode(genetcd)
            
            ORFs_found[[g]]<-ORFquant(region=gen_region,for_ORFquant=for_ORFquant_data,genetic_code=genetcd,
                                    orf_find.all_starts=stn.orf_find.all_starts,orf_find.nostarts=stn.orf_find.nostarts,
                                    orf_find.start_sel_cutoff = stn.orf_find.start_sel_cutoff,orf_find.start_sel_cutoff_ave = stn.orf_find.start_sel_cutoff_ave,
                                    orf_find.cutoff_fr_ave=stn.orf_find.cutoff_fr_ave,orf_quant.cutoff_cums = stn.orf_quant.cutoff_cums,
                                    orf_quant.cutoff_pct = stn.orf_quant.cutoff_pct,orf_quant.cutoff_P_sites=stn.orf_quant.cutoff_P_sites,unique_reads = unique_reads_only,orf_quant.scaling = stn.orf_quant.scaling)
            
            
        }
        
    }
    cat(paste("Summoning ORFquant --- Done! ",date(),"\n",sep = ""))
    
    cat(paste("Exporting ORFquant results ... ",date(),"\n",sep = ""))
    
    if(write_temp_files){
        save(ORFs_found,file=paste(prefix,"tmp_ORFquant_results",sep="_"))
    }
    
    lens<-elementNROWS(ORFs_found)
    
    ORFs_found<-ORFs_found[lens>0]
    if(length(ORFs_found)==0){stop(paste("No ORFs found! Please check sub-codon resolution of Ribo-seq reads and ensure the annotation is correct --- ",date(),"\n"))}
    
    ORFs_txs_feats<-unlist(GRangesList(lapply(ORFs_found,function(x){unlist(x$genomic_features)})))
    ORFs_txs_feats<-ORFs_txs_feats[!duplicated(mcols(ORFs_txs_feats)) | !duplicated(ORFs_txs_feats)]
    
    selected_txs<-sort(unique(unlist(ORFs_txs_feats$txs_selected)))
    
    lens<-elementNROWS(ORFs_found)
    
    ORFs_found<-ORFs_found[lens>1]
    if(length(ORFs_found)==0){stop(paste("No ORFs found! Please check sub-codon of Ribo-seq reads or that the annotation is correct --- ",date(),"\n"))}
    ORFs_tx<-unlist(GRangesList(unlist(sapply(ORFs_found,function(x){unlist(x$ORFs_tx_position)}))))
    
    ORFs_feat<-unlist(sapply(ORFs_found,function(x){unlist(x$selected_ORFs_features)}))
    ORFs_feat<-GRangesList(sapply(ORFs_feat,function(x){
        x$X$tx_name<-NULL
        return(x)
    }))
    
    
    #toAdd - use the features unique to ORFs as confidence score
    
    na_ps<-is.na(ORFs_tx$P_sites)
    ORFs_tx$P_sites_pN<-NA
    ORFs_tx$P_sites_pN[!na_ps]<-ORFs_tx$P_sites[!na_ps]/(width(ORFs_tx)[!na_ps])
    ORFs_tx$ORFs_pM<-NA
    ORFs_tx$ORFs_pM[!na_ps]<-ORFs_tx$P_sites_pN[!na_ps]*(1000000/(sum(ORFs_tx$P_sites_pN[!na_ps])))
    ORFs_tx$P_sites_pN<-NULL
    
    
    ORFs_gen<-unlist(GRangesList(sapply(ORFs_found,function(x){unlist(x$ORFs_genomic_position)})))
    
    ORFs_spl_feat_longest<-unlist(GRangesList(sapply(ORFs_found,function(x){unlist(x$ORFs_splice_feats$annotation_wrt_longest)})))
    ORFs_spl_feat_maxORF<-unlist(GRangesList(sapply(ORFs_found,function(x){unlist(x$ORFs_splice_feats$annotation_wrt_maxORF)})))
    ORFs_readthroughs<-unlist(GRangesList(unlist(sapply(ORFs_found,function(x){unlist(x$readthrough)}))))
    if(length(ORFs_readthroughs)>0){
        ORFs_readthroughs<-ORFs_readthroughs[order(ORFs_readthroughs$P_sites_raw,decreasing = T)]
    }
    ORFquant_results<-list(ORFs_tx,ORFs_gen,ORFs_feat,ORFs_spl_feat_longest,ORFs_spl_feat_maxORF,ORFs_readthroughs,ORFs_txs_feats,selected_txs)
    names(ORFquant_results)<-c("ORFs_tx","ORFs_gen","ORFs_feat","ORFs_spl_feat_longest","ORFs_spl_feat_maxORF","ORFs_readthroughs","ORFs_txs_feats","selected_txs")
    
    #added for seqinfo problem
    

    x<-ORFquant_results$ORFs_readthroughs
    seqf<-Seqinfo(seqnames = names(GTF_annotation$exons_txs),seqlengths = sum(width(GTF_annotation$exons_txs)),isCircular = NA,genome = NA)
    x@seqnames<-Rle(factor(as.character(x@seqnames),levels = seqlevels(seqf)))
    x@seqinfo<-seqf
    ORFquant_results$ORFs_readthroughs<-x
    
    x<-ORFquant_results$ORFs_tx
    seqf<-Seqinfo(seqnames = names(GTF_annotation$exons_txs),seqlengths = sum(width(GTF_annotation$exons_txs)),isCircular = NA,genome = NA)
    x@seqnames<-Rle(factor(as.character(x@seqnames),levels = seqlevels(seqf)))
    x@seqinfo<-seqf
    
    x$longest_ORF@seqnames<-Rle(factor(as.character(x$longest_ORF@seqnames),levels = seqlevels(seqf)))
    x$longest_ORF@seqinfo<-seqf
    
    ORFquant_results$ORFs_tx<-x
    
    
    save(ORFquant_results ,file = paste(prefix,"final_ORFquant_results",sep="_"))
    
    if(write_GTF_file){
        map_tx_genes<-mcols(ORFs_tx)[,c("ORF_id_tr","gene_id","gene_biotype","gene_name","transcript_id","transcript_biotype","P_sites","ORF_pct_P_sites","ORF_pct_P_sites_pN","ORFs_pM")]
        
        match_ORF<-match(names(ORFs_gen),map_tx_genes$ORF_id_tr)
        
        ORFs_gen$transcript_id<-map_tx_genes[match_ORF,"transcript_id"]
        
        match_tx<-match(ORFs_gen$transcript_id,map_tx_genes$transcript_id)
        
        ORFs_gen$transcript_biotype<-map_tx_genes[match_tx,"transcript_biotype"]
        ORFs_gen$gene_id<-map_tx_genes[match_tx,"gene_id"]
        ORFs_gen$gene_biotype<-map_tx_genes[match_tx,"gene_biotype"]
        ORFs_gen$gene_name<-map_tx_genes[match_tx,"gene_name"]
        ORFs_gen$ORF_id<-map_tx_genes[match_tx,"ORF_id_tr"]
        
        ORFs_gen$P_sites<-round(map_tx_genes[match_ORF,"P_sites"],digits=4)
        ORFs_gen$ORF_pct_P_sites<-round(map_tx_genes[match_ORF,"ORF_pct_P_sites"],digits=4)
        ORFs_gen$ORF_pct_P_sites_pN<-round(map_tx_genes[match_ORF,"ORF_pct_P_sites_pN"],digits=4)
        ORFs_gen$ORFs_pM<-round(map_tx_genes[match_ORF,"ORFs_pM"],digits=4)
        
        proteins_readthrough<-AAStringSet(ORFs_readthroughs$Protein)
        if(length(proteins_readthrough)>0){
            names(proteins_readthrough)<-paste(ORFs_readthroughs$ORF_id_tr,ORFs_readthroughs$gene_biotype,ORFs_readthroughs$gene_id,"readthrough","readthrough",sep="|")
            proteins_readthrough<-narrow(proteins_readthrough,start = start(proteins_readthrough)[1]+1)
            proteins_readthrough<-AAStringSet(gsub(proteins_readthrough,pattern = "[*]",replacement = "X"))
        }
        
        proteins<-AAStringSet(ORFs_tx$Protein)
        names(proteins)<-paste(ORFs_tx$ORF_id_tr,ORFs_tx$gene_biotype,ORFs_tx$gene_id,ORFs_tx$ORF_category_Gen,ORFs_tx$ORF_category_Tx_compatible,sep="|")
        proteins<-c(proteins,proteins_readthrough)
        if(write_protein_fasta){
            writeXStringSet(proteins,filepath= paste(prefix,"Protein_sequences.fasta",sep="_"))
        }
        
        map_tx_genes<-GTF_annotation$trann
        #ORFs_gen$transcript_id<-names(ORFs_gen)
        ORFs_gen$type="CDS"
        exs_gtf<-unlist(GTF_annotation$exons_txs[selected_txs])
        mcols(exs_gtf)<-NULL
        exs_gtf$transcript_id<-names(exs_gtf)
        exs_gtf$transcript_biotype<-map_tx_genes[match(exs_gtf$transcript_id,map_tx_genes$transcript_id),"transcript_biotype"]
        exs_gtf$gene_id<-map_tx_genes[match(exs_gtf$transcript_id,map_tx_genes$transcript_id),"gene_id"]
        exs_gtf$gene_biotype<-map_tx_genes[match(exs_gtf$transcript_id,map_tx_genes$transcript_id),"gene_biotype"]
        exs_gtf$gene_name<-map_tx_genes[match(exs_gtf$transcript_id,map_tx_genes$transcript_id),"gene_name"]
        mcols(exs_gtf)[,names(mcols(ORFs_gen))[!names(mcols(ORFs_gen))%in%names(mcols(exs_gtf))]]<-NA
        exs_gtf$type<-"exon"
        mcols(ORFs_gen)<-mcols(ORFs_gen)[,names(mcols(exs_gtf))]
        all<-sort(c(exs_gtf,ORFs_gen))
        all$`source`="ORFquant"
        names(all)<-NULL
        suppressWarnings(export.gff2(object=all,con=paste(prefix,"Detected_ORFs.gtf",sep="_")))
    }
    
    if(interactive){
        for_ORFquant<<-for_ORFquant_data
        ORFquant_results<<-ORFquant_results
    }
    cat(paste("Exporting ORFquant results --- Done! ",date(),"\n",sep = ""))
    ORFquant_results$psite_data_file <- for_ORFquant_file
    ORFquant_results   
}


#' Load genomic features and genome sequence
#'
#' This function loads the annotation created by the \code{prepare_annotation_files function}
#' @keywords ORFquant, Ribo-seQC
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param path Full path to the *Rannot R file in the annotation directory used in the \code{prepare_annotation_files function}
#' @return introduces a \code{GTF_annotation} object and a \code{genome_seq} object in the parent environment
#' @seealso \code{\link{prepare_annotation_files}}
#' @export

load_annotation<-function(path){
    GTF_annotation<-get(load(path))
    if(is(GTF_annotation$genome,'FaFile')){
        genome_sequence <- GTF_annotation$genome            
    }else{
        library(GTF_annotation$genome_package,character.only = T)
        genome_sequence<-get(GTF_annotation$genome_package)
    }
    GTF_annotation<<-GTF_annotation
    genome_seq<<-genome_sequence
}




#' Prepare comprehensive sets of annotated genomic features
#'
#' This function processes a gtf file and a twobit file (created using faToTwoBit from ucsc tools: http://hgdownload.soe.ucsc.edu/admin/exe/ ) to create a comprehensive set of genomic regions of interest in genomic and transcriptomic space (e.g. introns, UTRs, start/stop codons).
#'    In addition, by linking genome sequence and annotation, it extracts additional info, such as gene and transcript biotypes, genetic codes for different organelles, or chromosomes and transcripts lengths.
#' @keywords ORFquant, RiboseQC
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param annotation_directory The target directory which will contain the output files
#' @param twobit_file Full path to the genome file in twobit format
#' @param gtf_file Full path to the annotation file in GTF format
#' @param scientific_name A name to give to the organism studied; must be two words separated by a ".", defaults to Homo.sapiens
#' @param annotation_name A name to give to annotation used; defaults to genc25
#' @param export_bed_tables_TxDb Export coordinates and info about different genomic regions in the annotation_directory? It defaults to \code{TRUE}
#' @param forge_BSgenome Forge and install a \code{BSgenome} package? It defaults to \code{TRUE}
#' @param genome_seq Fasta file to use for genome seq if not forging a BSgenome package
#' @param circ_chroms Chromosomes to make circular in the genome sequence - defaults to DEFAULT_CIRC_SEQS
#' @param create_TxDb Create a \code{TxDb} object and a *Rannot object? It defaults to \code{TRUE}
#' @details This function uses the \code{makeTxDbFromGFF} function to  create a TxDb object and extract
#' genomic regions and other info to a *Rannot R file; the \code{mapToTranscripts} and \code{mapFromTranscripts} functions are used to 
#' map features to genomic or transcript-level coordinates. GTF file mist contain "exon" and "CDS" lines,
#' where each line contains "transcript_id" and "gene_id" values. Additional values such as "gene_biotype" or "gene_name" are also extracted.
#' Regarding sequences, the twobit file, together with input scientific and annotation names, is used to forge and install a 
#' BSgenome package using the \code{forgeBSgenomeDataPkg} function.\cr\cr
#' The resulting GTF_annotation object (obtained after runnning \code{load_annotation}) contains:\cr\cr
#' \code{txs}: annotated transcript boundaries.\cr
#' \code{txs_gene}: GRangesList including transcript grouped by gene.\cr
#' \code{seqinfo}: indicating chromosomes and chromosome lengths.\cr
#' \code{start_stop_codons}: the set of annotated start and stop codon, with respective transcript and gene_ids.
#' reprentative_mostcommon,reprentative_boundaries and reprentative_5len represent the most common start/stop codon,
#' the most upstream/downstream start/stop codons and the start/stop codons residing on transcripts with the longest 5'UTRs\cr
#' \code{cds_txs}: GRangesList including CDS grouped by transcript.\cr
#' \code{introns_txs}: GRangesList including introns grouped by transcript.\cr
#' \code{cds_genes}: GRangesList including CDS grouped by gene.\cr
#' \code{exons_txs}: GRangesList including exons grouped by transcript.\cr
#' \code{exons_bins}: the list of exonic bins with associated transcripts and genes.\cr
#' \code{junctions}: the list of annotated splice junctions, with associated transcripts and genes.\cr
#' \code{genes}: annotated genes coordinates.\cr
#' \code{threeutrs}: collapsed set of 3'UTR regions, with correspinding gene_ids. This set does not overlap CDS region.\cr
#' \code{fiveutrs}: collapsed set of 5'UTR regions, with correspinding gene_ids. This set does not overlap CDS region.\cr
#' \code{ncIsof}: collapsed set of exonic regions of protein_coding genes, with correspinding gene_ids. This set does not overlap CDS region.\cr
#' \code{ncRNAs}: collapsed set of exonic regions of non_coding genes, with correspinding gene_ids. This set does not overlap CDS region.\cr
#' \code{introns}: collapsed set of intronic regions, with correspinding gene_ids. This set does not overlap exonic region.\cr
#' \code{intergenicRegions}: set of intergenic regions, defined as regions with no annotated genes on either strand.\cr
#' \code{trann}: DataFrame object including (when available) the mapping between gene_id, gene_name, gene_biotypes, transcript_id and transcript_biotypes.\cr
#' \code{cds_txs_coords}: transcript-level coordinates of ORF boundaries, for each annotated coding transcript. Additional columns are the same as as for the \code{start_stop_codons} object.\cr
#' \code{genetic_codes}: an object containing the list of genetic code ids used for each chromosome/organelle. see GENETIC_CODE_TABLE for more info.\cr
#' \code{genome_package}: the name of the forged BSgenome package. Loaded with \code{load_annotation} function.\cr
#' \code{stop_in_gtf}: stop codon, as defined in the annotation.\cr
#' @return a TxDb file and a *Rannot files are created in the specified \code{annotation_directory}. 
#' In addition, a BSgenome object is forged, installed, and linked to the *Rannot object
#' @seealso \code{\link{load_annotation}}, \code{\link{forgeBSgenomeDataPkg}}, \code{\link{makeTxDbFromGFF}}, \code{\link{run_ORFquant}}.
#' @export

prepare_annotation_files<-function(annotation_directory,twobit_file=NULL,gtf_file,scientific_name="Homo.sapiens",annotation_name="genc25",export_bed_tables_TxDb=TRUE,forge_BSgenome=TRUE,genome_seq=NULL,circ_chroms=DEFAULT_CIRC_SEQS,create_TxDb=TRUE){


  
    DEFAULT_CIRC_SEQS <- unique(c("chrM","MT","MtDNA","mit","Mito","mitochondrion",
                                  "dmel_mitochondrion_genome","Pltd","ChrC","Pt","chloroplast",
                                  "Chloro","2micron","2-micron","2uM",
                                  "Mt", "NC_001879.2", "NC_006581.1","ChrM","mitochondrion_genome"))

    #adjust variable names (some chars not permitted)
    annotation_name<-gsub(annotation_name,pattern = "_",replacement = "")
    annotation_name<-gsub(annotation_name,pattern = "-",replacement = "")
    if(!dir.exists(annotation_directory)){dir.create(path = annotation_directory,recursive = TRUE)}
    annotation_directory<-normalizePath(annotation_directory)
    gtf_file<-normalizePath(gtf_file)

    filestotest <- c(gtf_file)
    if(forge_BSgenome) filestotest <- c(filestotest,twobit_file)
    for (f in filestotest){
        if(file.access(f, 0)==-1) {
            stop("
                 The following files don't exist:\n",
                 f, "\n")
        }
        }
    }
    
   
    if(forge_BSgenome){

    scientific_name_spl<-strsplit(scientific_name,"[.]")[[1]]
    ok<-length(scientific_name_spl)==2
    if(!ok){stop("\"scientific_name\" must be two words separated by a \".\", like \"Homo.sapiens\"")}
    
    #get circular sequences
    
    seqinfotwob<-seqinfo(TwoBitFile(twobit_file))
    circss<-seqnames(seqinfotwob)[which(seqnames(seqinfotwob)%in%DEFAULT_CIRC_SEQS)]
    seqinfotwob@is_circular[which(seqnames(seqinfotwob)%in%DEFAULT_CIRC_SEQS)]<-TRUE
    pkgnm<-paste("BSgenome",scientific_name,annotation_name,sep=".")
    
    circseed<-circss
    if(length(circseed)==0){circseed<-NULL}
    
    #Forge a BSGenome package
    
    if(forge_BSgenome){
        cat(paste("Creating the BSgenome package ... ",date(),"\n",sep = ""))
        seed_text<-paste("Package: BSgenome.",scientific_name,".",annotation_name,"\n",
                         "Title: Full genome sequences for ",scientific_name,", ",annotation_name,"\n",
                         "Description: Full genome sequences for ",scientific_name,", ",annotation_name,"\n",
                         "Version: 1.0","\n",
                         "organism: ",scientific_name,"\n",
                         "common_name: ",scientific_name,"\n",
                         "provider: NA","\n",
                         "provider_version: ",annotation_name,"\n",
                         "release_date: NA","\n",
                         "release_name: NA","\n",
                         "source_url: NA","\n",
                         "organism_biocview: ", scientific_name,"\n",
                         "BSgenomeObjname: ",scientific_name,"\n",
                         "seqs_srcdir: ",dirname(twobit_file),"\n",
                         "seqfile_name: ",basename(twobit_file),sep="")
        
        
        seed_dest<-paste(annotation_directory,"/",basename(twobit_file),"_",scientific_name,"_seed",sep = "")
        
        if(length(circseed)==0){
            writeLines(text = seed_text,con = seed_dest)
        }
        if(length(circseed)==1){
            seed_text<-paste(seed_text,"\n",
                             "circ_seqs: \"",circseed,"\"",sep="")
            writeLines(text = seed_text,con = seed_dest)
        }
        
        if(length(circseed)>1){
            circseed<-paste('c("',paste(circseed,collapse=","),'")',sep="")
            circseed<-gsub(circseed,pattern = ",",replacement='","')
            
            cat(seed_text,"\n","circ_seqs: ",circseed,"\n",sep="",file = seed_dest)
        }
        
        
        
        unlink(paste(annotation_directory,pkgnm,sep="/"),recursive=T)
        
        forgeBSgenomeDataPkg(x=seed_dest,destdir=annotation_directory,seqs_srcdir=dirname(twobit_file))
        cat(paste("Creating the BSgenome package --- Done! ",date(),"\n",sep = ""))
        
        cat(paste("Installing the BSgenome package ... ",date(),"\n",sep = ""))
        
        install(paste(annotation_directory,pkgnm,sep="/"),upgrade = F)
        cat(paste("Installing the BSgenome package --- Done! ",date(),"\n",sep = ""))

        
        }else{
            if(!is(genome_seq,'FaFile')){
                genome_seq <- Rsamtools::FaFile(genome_seq)
            }
            if(!is(genome_seq,'FaFile_Circ')){
                genome_seq <- FaFile_Circ(genome_seq,circularRanges=circ_chroms)
            }
            seqinfo(genome_seq)@is_circular[which(seqnames(seqinfo_genome)%in%circ_chroms)]<-TRUE
            seqinfotwob<-seqinfo(genome_seq)
        }
        
    }
    
    #Create the TxDb from GTF and BSGenome info
    
    if(create_TxDb){
        cat(paste("Creating the TxDb object ... ",date(),"\n",sep = ""))
        
        annotation<-makeTxDbFromGFF(file=gtf_file,format="gtf",chrominfo = seqinfotwob)
        
        saveDb(annotation, file=paste(annotation_directory,"/",basename(gtf_file),"_TxDb",sep=""))
        cat(paste("Creating the TxDb object --- Done! ",date(),"\n",sep = ""))
        cat(paste("Extracting genomic regions ... ",date(),"\n",sep = ""))
        
        genes<-genes(annotation)
        exons_ge<-exonsBy(annotation,by="gene")
        exons_ge<-reduce(exons_ge)
        
        cds_gen<-cdsBy(annotation,"gene")
        cds_ge<-reduce(cds_gen)
        
        
        #define regions not overlapping CDS ( or exons when defining introns)
        
        threeutrs<-reduce(GenomicRanges::setdiff(unlist(threeUTRsByTranscript(annotation)),unlist(cds_ge),ignore.strand=FALSE))
        
        fiveutrs<-reduce(GenomicRanges::setdiff(unlist(fiveUTRsByTranscript(annotation)),unlist(cds_ge),ignore.strand=FALSE))
        
        introns<-reduce(GenomicRanges::setdiff(unlist(intronsByTranscript(annotation)),unlist(exons_ge),ignore.strand=FALSE))
        
        nc_exons<-reduce(GenomicRanges::setdiff(unlist(exons_ge),reduce(c(unlist(cds_ge),fiveutrs,threeutrs)),ignore.strand=FALSE))
        
        #assign gene ids (mutiple when overlapping multiple genes)
        ov<-findOverlaps(threeutrs,genes)
        ov<-split(subjectHits(ov),queryHits(ov))
        threeutrs$gene_id<-CharacterList(lapply(ov,FUN = function(x){names(genes)[x]}))
        ov<-findOverlaps(fiveutrs,genes)
        ov<-split(subjectHits(ov),queryHits(ov))
        fiveutrs$gene_id<-CharacterList(lapply(ov,FUN = function(x){names(genes)[x]}))
        ov<-findOverlaps(introns,genes)
        ov<-split(subjectHits(ov),queryHits(ov))
        introns$gene_id<-CharacterList(lapply(ov,FUN = function(x){names(genes)[x]}))
        ov<-findOverlaps(nc_exons,genes)
        ov<-split(subjectHits(ov),queryHits(ov))
        nc_exons$gene_id<-CharacterList(lapply(ov,FUN = function(x){names(genes)[x]}))
        
        intergenicRegions<-genes
        strand(intergenicRegions)<-"*"
        intergenicRegions <- gaps(reduce(intergenicRegions))
        intergenicRegions<-intergenicRegions[strand(intergenicRegions)=="*"]
        
        cds_tx<-cdsBy(annotation,"tx",use.names=T)
        txs_gene<-transcriptsBy(annotation,by="gene")
        genes_red<-reduce(sort(genes(annotation)))
        
        exons_tx<-exonsBy(annotation,"tx",use.names=T)
        
        transcripts_db<-transcripts(annotation)
        intron_names_tx<-intronsByTranscript(annotation,use.names=T)
        
        
        #define exonic bins, including regions overlapping multiple genes
        nsns<-disjointExons(annotation,aggregateGenes=T)
        
        
        
        #define tx_coordinates of ORF boundaries
        
        exsss_cds<-exons_tx[names(cds_tx)]
        chunks<-seq(1,length(cds_tx),by = 20000)
        if(chunks[length(chunks)]<length(cds_tx)){chunks<-c(chunks,length(cds_tx))}
        mapp<-GRangesList()
        for(i in 1:(length(chunks)-1)){
            if(i!=(length(chunks)-1)){
                mapp<-suppressWarnings(c(mapp,pmapToTranscripts(cds_tx[chunks[i]:(chunks[i+1]-1)],transcripts = exsss_cds[chunks[i]:(chunks[i+1]-1)])))
            }
            if(i==(length(chunks)-1)){
                mapp<-suppressWarnings(c(mapp,pmapToTranscripts(cds_tx[chunks[i]:(chunks[i+1])],transcripts = exsss_cds[chunks[i]:(chunks[i+1])])))
            }
        }
        cds_txscoords<-unlist(mapp)
        
        
        #extract biotypes and ids
        
        cat(paste("Extracting ids and biotypes ... ",date(),"\n",sep = ""))
        
        trann<-unique(mcols(import.gff2(gtf_file,colnames=c("gene_id","gene_biotype","gene_type","gene_name","gene_symbol","transcript_id","transcript_biotype","transcript_type"))))
        trann<-trann[!is.na(trann$transcript_id),]
        trann<-data.frame(unique(trann),stringsAsFactors=F)
        
        if(sum(!is.na(trann$transcript_biotype))==0 & sum(!is.na(trann$transcript_type))==0 ){
            trann$transcript_biotype<-"no_type"
        }
        if(sum(!is.na(trann$transcript_biotype))==0){trann$transcript_biotype<-NULL}
        if(sum(!is.na(trann$transcript_type))==0){trann$transcript_type<-NULL}
        
        
        if(sum(!is.na(trann$gene_biotype))==0 & sum(!is.na(trann$gene_type))==0 ){
            
            trann$gene_type<-"no_type"
            
        }
        if(sum(!is.na(trann$gene_name))==0 & sum(!is.na(trann$gene_symbol))==0 ){
            
            trann$gene_name<-"no_name"
            
        }
        if(sum(!is.na(trann$gene_biotype))==0){trann$gene_biotype<-NULL}
        if(sum(!is.na(trann$gene_type))==0){trann$gene_type<-NULL}
        if(sum(!is.na(trann$gene_name))==0){trann$gene_name<-NULL}
        if(sum(!is.na(trann$gene_symbol))==0){trann$gene_symbol<-NULL}
        colnames(trann)<-c("gene_id","gene_biotype","gene_name","transcript_id","transcript_biotype")
        
        trann<-DataFrame(trann)
        
        
        
        #introns and transcript_ids/gene_ids
        unq_intr<-sort(unique(unlist(intron_names_tx)))
        names(unq_intr)<-NULL
        all_intr<-unlist(intron_names_tx)
        
        ov<-findOverlaps(unq_intr,all_intr,type="equal")
        ov<-split(subjectHits(ov),queryHits(ov))
        a_nam<-CharacterList(lapply(ov,FUN = function(x){unique(names(all_intr)[x])}))
        
        unq_intr$type="J"
        unq_intr$tx_name<-a_nam
        
        
        mat_genes<-match(unq_intr$tx_name,trann$transcript_id)
        g<-unlist(apply(cbind(1:length(mat_genes),Y = elementNROWS(mat_genes)),FUN =function(x) rep(x[1],x[2]),MARGIN = 1))
        g2<-split(trann[unlist(mat_genes),"gene_id"],g)
        unq_intr$gene_id<-CharacterList(lapply(g2,unique))
        
        
        #filter ncRNA and ncIsof regions
        ncrnas<-nc_exons[!nc_exons%over%genes[trann$gene_id[trann$gene_biotype=="protein_coding"]]]
        ncisof<-nc_exons[nc_exons%over%genes[trann$gene_id[trann$gene_biotype=="protein_coding"]]]
        
        
        # define genetic codes to use
        # IMPORTANT : modify if needed (e.g. different organelles or species) check ids of GENETIC_CODE_TABLE for more info
        
        ifs<-seqinfo(annotation)
        translations<-as.data.frame(ifs)
        translations$genetic_code<-"1"
        
        #insert new codes for chromosome name
        
        #Mammalian mito
        translations$genetic_code[rownames(translations)%in%c("chrM","MT","MtDNA","mit","mitochondrion")]<-"2"
        
        #Yeast mito
        translations$genetic_code[rownames(translations)%in%c("Mito")]<-"3"
        
        #Drosophila mito
        translations$genetic_code[rownames(translations)%in%c("dmel_mitochondrion_genome")]<-"5"
        
        circs<-ifs@seqnames[which(ifs@is_circular)]
        
        
        #define start and stop codons (genome space)
        
        suppressPackageStartupMessages(library(pkgnm,character.only=TRUE))
        genome<-get(pkgnm)
        tocheck<-as.character(runValue(seqnames(cds_tx)))
        tocheck<-cds_tx[!tocheck%in%circs]
        seqcds<-extractTranscriptSeqs(genome,transcripts = tocheck)
        cd<-unique(translations$genetic_code[!rownames(translations)%in%circs])
        trsl<-suppressWarnings(translate(seqcds,genetic.code = getGeneticCode(cd),if.fuzzy.codon = "solve"))
        trslend<-as.character(narrow(trsl,end = width(trsl),width = 1))
        stop_inannot<-NA
        if(names(sort(table(trslend),decreasing = T)[1])=="*"){stop_inannot<-"*"}
        
        cds_txscoords$gene_id<-trann$gene_id[match(as.vector(seqnames(cds_txscoords)),trann$transcript_id)]
        cds_cc<-cds_txscoords
        strand(cds_cc)<-"*"
        sta_cc<-resize(cds_cc,width = 1,"start")
        sta_cc<-unlist(pmapFromTranscripts(sta_cc,exons_tx[seqnames(sta_cc)],ignore.strand=F))
        sta_cc$gene_id<-trann$gene_id[match(names(sta_cc),trann$transcript_id)]
        sta_cc<-sta_cc[sta_cc$hit]
        strand(sta_cc)<-structure(as.vector(strand(transcripts_db)),names=transcripts_db$tx_name)[names(sta_cc)]
        sta_cc$type<-"start_codon"
        mcols(sta_cc)<-mcols(sta_cc)[,c("exon_rank","type","gene_id")]
        
        sto_cc<-resize(cds_cc,width = 1,"end")
        #stop codon is the 1st nt, e.g. U of the UAA
        #To-do: update with regards to different organelles, and different annotations
        sto_cc<-shift(sto_cc,-2)
        if(is.na(stop_inannot)){sto_cc<-resize(trim(shift(sto_cc,3)),width = 1,fix = "end")}
        
        sto_cc<-unlist(pmapFromTranscripts(sto_cc,exons_tx[seqnames(sto_cc)],ignore.strand=F))
        sto_cc<-sto_cc[sto_cc$hit]
        sto_cc$gene_id<-trann$gene_id[match(names(sto_cc),trann$transcript_id)]
        strand(sto_cc)<-structure(as.vector(strand(transcripts_db)),names=transcripts_db$tx_name)[names(sto_cc)]
        sto_cc$type<-"stop_codon"
        mcols(sto_cc)<-mcols(sto_cc)[,c("exon_rank","type","gene_id")]
        
        
        #define most common, most upstream/downstream
        
        cat(paste("Defining most common start/stop codons ... ",date(),"\n",sep = ""))
        
        start_stop_cc<-sort(c(sta_cc,sto_cc))
        start_stop_cc$transcript_id<-names(start_stop_cc)
        start_stop_cc$most_up_downstream<-FALSE
        start_stop_cc$most_frequent<-FALSE
        
        df<-cbind.DataFrame(start(start_stop_cc),start_stop_cc$type,start_stop_cc$gene_id)
        colnames(df)<-c("start_pos","type","gene_id")
        upst<-by(df$start_pos,INDICES = df$gene_id,function(x){x==min(x) | x==max(x)})
        start_stop_cc$most_up_downstream<-unlist(upst[unique(df$gene_id)])
        
        mostfr<-by(df[,c("start_pos","type")],INDICES = df$gene_id,function(x){
            mfreq<-table(x)
            x$start_pos%in%as.numeric(names(which(mfreq[,1]==max(mfreq[,1])))) | x$start_pos%in%as.numeric(names(which(mfreq[,2]==max(mfreq[,2]))))
        })
        
        start_stop_cc$most_frequent<-unlist(mostfr[unique(df$gene_id)])
        
        names(start_stop_cc)<-NULL
        
        
        
        #define transcripts as containing frequent start/stop codons or most upstream ones, in relation with 5'UTR length
        
        mostupstr_tx<-sum(LogicalList(split(start_stop_cc$most_up_downstream,start_stop_cc$transcript_id)))[as.character(seqnames(cds_txscoords))]
        cds_txscoords$upstr_stasto<-mostupstr_tx
        mostfreq_tx<-sum(LogicalList(split(start_stop_cc$most_frequent,start_stop_cc$transcript_id)))[as.character(seqnames(cds_txscoords))]
        cds_txscoords$mostfreq_stasto<-mostfreq_tx
        cds_txscoords$lentx<-sum(width(exons_tx[as.character(seqnames(cds_txscoords))]))
        df<-cbind.DataFrame(as.character(seqnames(cds_txscoords)),width(cds_txscoords),start(cds_txscoords),cds_txscoords$mostfreq_stasto,cds_txscoords$gene_id)
        colnames(df)<-c("txid","cdslen","utr5len","var","gene_id")
        repres_freq<-by(df[,c("txid","cdslen","utr5len","var")],df$gene_id,function(x){
            x<-x[order(x$var,x$utr5len,x$cdslen,decreasing = T),]
            x<-x[x$var==max(x$var),]
            ok<-x$txid[which(x$cdslen==max(x$cdslen) & x$utr5len==max(x$utr5len) & x$var==max(x$var))][1]
            if(length(ok)==0 | is.na(ok[1])){ok<-x$txid[1]}
            ok
        })
        
        df<-cbind.DataFrame(as.character(seqnames(cds_txscoords)),width(cds_txscoords),start(cds_txscoords),cds_txscoords$upstr_stasto,cds_txscoords$gene_id)
        colnames(df)<-c("txid","cdslen","utr5len","var","gene_id")
        repres_upstr<-by(df[,c("txid","cdslen","utr5len","var")],df$gene_id,function(x){
            x<-x[order(x$var,x$utr5len,x$utr5len,decreasing = T),]
            x<-x[x$var==max(x$var),]
            ok<-x$txid[which(x$cdslen==max(x$cdslen) & x$utr5len==max(x$utr5len) & x$var==max(x$var))][1]
            if(length(ok)==0 | is.na(ok[1])){ok<-x$txid[1]}
            ok
        })
        df<-cbind.DataFrame(as.character(seqnames(cds_txscoords)),width(cds_txscoords),start(cds_txscoords),cds_txscoords$upstr_stasto,cds_txscoords$gene_id)
        colnames(df)<-c("txid","cdslen","utr5len","var","gene_id")
        repres_len5<-by(df[,c("txid","cdslen","utr5len","var")],df$gene_id,function(x){
            x<-x[order(x$utr5len,x$var,x$cdslen,decreasing = T),]
            ok<-x$txid[which(x$utr5len==max(x$utr5len) & x$var==max(x$var))][1]
            if(length(ok)==0 | is.na(ok[1])){ok<-x$txid[1]}
            ok
        })
        
        cds_txscoords$reprentative_mostcommon<-as.character(seqnames(cds_txscoords))%in%unlist(repres_freq)
        cds_txscoords$reprentative_boundaries<-as.character(seqnames(cds_txscoords))%in%unlist(repres_upstr)
        cds_txscoords$reprentative_5len<-as.character(seqnames(cds_txscoords))%in%unlist(repres_len5)
        unq_stst<-start_stop_cc
        mcols(unq_stst)<-NULL
        unq_stst<-sort(unique(unq_stst))
        ov<-findOverlaps(unq_stst,start_stop_cc,type="equal")
        ov<-split(subjectHits(ov),queryHits(ov))
        unq_stst$type<-CharacterList(lapply(ov,FUN = function(x){unique(start_stop_cc$type[x])}))
        unq_stst$transcript_id<-CharacterList(lapply(ov,FUN = function(x){start_stop_cc$transcript_id[x]}))
        unq_stst$gene_id<-CharacterList(lapply(ov,FUN = function(x){unique(start_stop_cc$gene_id[x])}))
        
        unq_stst$reprentative_mostcommon<-sum(!is.na(match(unq_stst$transcript_id,unlist(as(repres_freq,"CharacterList")))))>0
        unq_stst$reprentative_boundaries<-sum(!is.na(match(unq_stst$transcript_id,unlist(as(repres_upstr,"CharacterList")))))>0
        unq_stst$reprentative_5len<-sum(!is.na(match(unq_stst$transcript_id,unlist(as(repres_len5,"CharacterList")))))>0
        
        
        #put in a list
        GTF_annotation<-list(transcripts_db,txs_gene,ifs,unq_stst,cds_tx,intron_names_tx,cds_gen,exons_tx,nsns,unq_intr,genes,threeutrs,fiveutrs,ncisof,ncrnas,introns,intergenicRegions,trann,cds_txscoords,translations,pkgnm,stop_inannot)
        names(GTF_annotation)<-c("txs","txs_gene","seqinfo","start_stop_codons","cds_txs","introns_txs","cds_genes","exons_txs","exons_bins","junctions","genes","threeutrs","fiveutrs","ncIsof","ncRNAs","introns","intergenicRegions","trann","cds_txs_coords","genetic_codes","genome_package","stop_in_gtf")
        
        txs_all<-unique(GTF_annotation$trann$transcript_id)
        txs_exss<-unique(names(GTF_annotation$exons_txs))
        
        txs_notok<-txs_all[!txs_all%in%txs_exss]
        if(length(txs_notok)>0){
            set.seed(666)
            cat(paste("Warning: ",length(txs_notok)," txs with incorrect/unspecified exon boundaries - e.g. trans-splicing events, examples: "
                      ,paste(txs_notok[sample(1:length(txs_notok),size = min(3,length(txs_notok)),replace = F)],collapse=", ")," - ",date(),"\n",sep = ""))
        }
        
        
        #Save as a RData object
        save(GTF_annotation,file=paste(annotation_directory,"/",basename(gtf_file),"_Rannot",sep=""))
        cat(paste("Rannot object created!   ",date(),"\n",sep = ""))
        
        
        #create tables and bed files (with colnames, so with header)
        if(export_bed_tables_TxDb==T){
            cat(paste("Exporting annotation tables ... ",date(),"\n",sep = ""))
            for(bed_file in c("fiveutrs","threeutrs","ncIsof","ncRNAs","introns","cds_txs_coords")){
                bf<-GTF_annotation[[bed_file]]
                bf_t<-bf
                if(length(bf)>0){
                    bf_t<-data.frame(chromosome=seqnames(bf),start=start(bf),end=end(bf),name=".",score=width(bf),strand=strand(bf))
                    meccole<-mcols(bf)
                    for(mecc in names(meccole)){
                        if(is(meccole[,mecc],"CharacterList") | is(meccole[,mecc],"NumericList") | is(meccole[,mecc],"IntegerList")){
                            meccole[,mecc]<-paste(meccole[,mecc],collapse=";")
                        }
                    }
                    bf_t<-cbind.data.frame(bf_t,meccole)
                }
                write.table(bf_t,file = paste(annotation_directory,"/",bed_file,"_similbed.bed",sep=""),sep="\t",quote = FALSE,row.names = FALSE,col.names = F)
                
            }
            
            write.table(GTF_annotation$trann,file = paste(annotation_directory,"/table_gene_tx_IDs",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
            seqi<-as.data.frame(GTF_annotation$seqinfo)
            seqi$chromosome<-rownames(seqi)
            write.table(seqi,file = paste(annotation_directory,"/seqinfo",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
            
            gen_cod<-as.data.frame(GTF_annotation$genetic_codes)
            gen_cod$chromosome<-rownames(gen_cod)
            write.table(gen_cod,file = paste(annotation_directory,"/genetic_codes",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
            cat(paste("Exporting annotation tables --- Done! ",date(),"\n",sep = ""))
            
        }
        
    }
    return(annot_file)
}



#' Offset spliced reads on plus strand
#'
#' This function calculates P-sites positions for spliced reads on the plus strand
#' @keywords Ribo-seQC, ORFquant
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param x a \code{GAlignments} object with a cigar string
#' @param cutoff number representing the offset value
#' @return a \code{GRanges} object with offset reads
#' @seealso \code{\link{prepare_for_ORFquant}}
#' @export

get_ps_fromspliceplus<-function(x,cutoff){
    rang<-cigarRangesAlongReferenceSpace(cigar(x), pos=start(x),ops="M")
    cs<-lapply(rang,function(x){cumsum(x@width)})
    rangok<-lapply(which(IntegerList(cs)>cutoff),"[[",1)
    rangok<-unlist(rangok)
    gr<-as(x,"GRanges")
    
    ones<-rangok==1
    mores<-rangok>1
    psmores<-GRanges()
    psones<-GRanges()
    
    if(sum(ones)>0){
        psones<-shift(resize(gr[ones],width=1,fix="start"),shift=cutoff)
    }
    
    if(sum(mores)>0){
        rangmore<-rang[mores]
        rangok<-rangok[mores]
        cms<-cumsum(width(rangmore))
        shft<-c()
        for(i in 1:length(rangok)){shft<-c(shft,cutoff-cms[[i]][rangok[i]-1])}
        stt<-start(rangmore)
        stok<-c()
        for(i in 1:length(shft)){stok<-c(stok,stt[[i]][rangok[i]]+shft[i])}
        psmores<-GRanges(IRanges(start=stok,width = 1),seqnames = seqnames(x[mores]),strand=strand(x[mores]),seqlengths=seqlengths(x[mores]))
        
    }
    ps<-sort(c(psones,psmores))
    return(ps)
}

#' Offset spliced reads on minus strand
#'
#' This function calculates P-sites positions for spliced reads on the minus strand
#' @keywords Ribo-seQC, ORFquant
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param x a \code{GAlignments} object with a cigar string
#' @param cutoff number representing the offset value
#' @return a \code{GRanges} object with offset reads
#' @seealso \code{\link{prepare_for_ORFquant}}
#' @export

get_ps_fromsplicemin<-function(x,cutoff){
    rang<-cigarRangesAlongReferenceSpace(cigar(x), pos=start(x),ops="M")
    rang<-endoapply(rang,rev)
    cs<-lapply(rang,function(x){cumsum(x@width)})
    rangok<-lapply(which(IntegerList(cs)>cutoff),"[[",1)
    rangok<-unlist(rangok)
    gr<-as(x,"GRanges")
    
    ones<-rangok==1
    mores<-rangok>1
    psmores<-GRanges()
    psones<-GRanges()
    
    
    if(sum(ones)>0){
        psones<-shift(resize(gr[ones],width=1,fix="start"),shift=-cutoff)
    }
    
    
    if(sum(mores)>0){
        rangmore<-rang[mores]
        rangok<-rangok[mores]
        cms<-cumsum(width(rangmore))
        shft<-c()
        for(i in 1:length(rangok)){shft<-c(shft,cutoff-cms[[i]][rangok[i]-1])}
        #start?
        stt<-end(rangmore)
        stok<-c()
        for(i in 1:length(shft)){stok<-c(stok,stt[[i]][rangok[i]]-shft[i])}
        psmores<-GRanges(IRanges(start=stok,width = 1),seqnames = seqnames(x[mores]),strand=strand(x[mores]),seqlengths=seqlengths(x[mores]))
        
    }
    ps<-sort(c(psones,psmores))
    return(ps)
    
    
    if(rangok!=1){
        
        ps<-shift(resize(GRanges(rang[rangok],seqnames=seqnames(x),strand=strand(x),seqlengths=seqlengths(x)),width=1,fix="start"),shift=-(cutoff-sum(rang[1:(rangok-1)]@width)))
    }
    return(ps)
}

#' Prepare the "for_ORFquant" file
#'
#' 
#' @details This function uses a list of pre-determined read lengths, cutoffs and compartments to calculate P_sites positions.\cr
#' Alternatively, bigwig files containing P_sites position for each strand can be specified. Optional bigwig files for uniquely mapping P_sites position (with and without mismatches)
#' can be specified to obtain more statistics on the ORFquant-identified ORFs
#' @keywords ORFquant
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param annotation_file Full path to the annotation file (*Rannot)
#' @param bam_file Full path to the bam file
#' @param chunk_size the number of alignments to read at each iteration, defaults to 5000000, increase when more RAM is available
#' @param path_to_rl_cutoff_file path to the rl_cutoff_file file specifying in 3 columns the read lengths, cutoffs and compartments ("nucl" for standard chromosomes)
#' @param path_to_P_sites_plus_bw path to a bigwig file containing P_sites positions on the plus strand
#' @param path_to_P_sites_minus_bw path to a bigwig file containing P_sites positions on the minus strand
#' @param path_to_P_sites_uniq_plus_bw (Optional) path to a bigwig file containing uniquely mapping P_sites positions on the plus strand
#' @param path_to_P_sites_uniq_minus_bw (Optional) path to a bigwig file containing uniquely mapping P_sites positions on the minus strand
#' @param path_to_P_sites_uniq_mm_plus_bw (Optional) path to a bigwig file containing uniquely mapping (with mismatches) P_sites positions on the plus strand
#' @param path_to_P_sites_uniq_mm_minus_bw (Optional) path to a bigwig file containing uniquely mapping (with mismatches) P_sites positions on the minus strand
#' @param dest_name prefix to use for the output files. Defaults to same as \code{bam_file} (appends "for_ORFquant" to its filename)
#' @seealso \code{\link{run_ORFquant}}
#' @export

prepare_for_ORFquant<-function(annotation_file,bam_file,path_to_rl_cutoff_file=NA,chunk_size=5000000,path_to_P_sites_plus_bw=NA,
                             path_to_P_sites_minus_bw=NA,path_to_P_sites_uniq_plus_bw=NA,path_to_P_sites_uniq_minus_bw=NA,
                             path_to_P_sites_uniq_mm_plus_bw=NA,path_to_P_sites_uniq_mm_minus_bw=NA,
                             dest_name=NA){
    
    load_annotation(annotation_file)
    
    if(is.na(dest_name)){dest_name=bam_file}
    
    if(is.na(path_to_rl_cutoff_file) & is.na(path_to_P_sites_plus_bw) & is.na(path_to_P_sites_minus_bw)){
        stop(paste("Please input either the paths to the P_sites bw files, or the path a suitable rl_cutoff table! ", date(),sep=""))
    }
    
    if(!is.na(path_to_rl_cutoff_file) & !is.na(path_to_P_sites_plus_bw) & !is.na(path_to_P_sites_minus_bw)){
        stop(paste("Please input either the paths to the P_sites bw files, or the path a suitable rl_cutoff table! ", date(),sep=""))
    }
    
    if(!is.na(path_to_rl_cutoff_file)){
        rl_cutoff<-read.table(path_to_rl_cutoff_file,header = T,sep = "\t",stringsAsFactors = F)
        
        if(dim(rl_cutoff)[2]!=3){stop(
            paste("Error: please format the rl_cutoff file correctly, using 3 tab-separated columns with 'rl', 'cutoff' and 'compartment' as column names! ",date(),sep="")
        )}
        
        rl_cutoffs_comp<-split(rl_cutoff,rl_cutoff[,3])
        compnms<-names(rl_cutoffs_comp)
        
        for(compar in compnms){
            cat(paste("Using ",paste(rl_cutoffs_comp[[compar]][,1],collapse=","), " nt long footprints with ",paste(rl_cutoffs_comp[[compar]][,2],collapse=",")," as cutoffs, '", compar,"' compartment ... ","\n",sep=""))
        }
    }
    
    opts <- BamFile(file=bam_file, yieldSize=chunk_size) 
    circs_seq<-seqnames(GTF_annotation$seqinfo)[which(isCircular(GTF_annotation$seqinfo))]
    seqs <- seqinfo(opts)
    circs <- seqs@seqnames[which(seqs@seqnames%in%circs_seq)]
    
    param <- ScanBamParam(flag=scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE),what=c("mapq"),tag = "MD")
    seqllll<-seqlevels(GTF_annotation$seqinfo)
    seqleee<-seqlengths(GTF_annotation$seqinfo)
    
    input_P_sites<-GRanges()
    seqlevels(input_P_sites)<-seqllll
    seqlengths(input_P_sites)<-seqleee
    input_P_sites_uniq<-input_P_sites
    input_P_sites_uniq_mm<-input_P_sites
    
    if(!is.na(path_to_P_sites_plus_bw)){
        input_P_sites_mn<-GRanges()
        input_P_sites_pl<-import(path_to_P_sites_plus_bw)
        strand(input_P_sites_pl)<-"+"
        if(!is.na(path_to_P_sites_minus_bw)){
            input_P_sites_mn<-import(path_to_P_sites_minus_bw)
            strand(input_P_sites_mn)<-"-"
        }
        suppressWarnings(input_P_sites<-sort(c(input_P_sites_pl,input_P_sites_mn)))
        seqlevels(input_P_sites,pruning.mode="coarse")<-seqllll
        seqlengths(input_P_sites)<-seqleee
        
    }
    
    if(!is.na(path_to_P_sites_uniq_plus_bw)){
        input_P_sites_uniq_mn<-GRanges()
        input_P_sites_uniq_pl<-import(path_to_P_sites_uniq_plus_bw)
        strand(input_P_sites_uniq_pl)<-"+"
        if(!is.na(path_to_P_sites_uniq_minus_bw)){
            input_P_sites_uniq_mn<-import(path_to_P_sites_uniq_minus_bw)
            strand(input_P_sites_uniq_mn)<-"-"
        }
        suppressWarnings(input_P_sites_uniq<-sort(c(input_P_sites_uniq_pl,input_P_sites_uniq_mn)))
        seqlevels(input_P_sites_uniq,pruning.mode="coarse")<-seqllll
        seqlengths(input_P_sites_uniq)<-seqleee
        
    }
    
    if(!is.na(path_to_P_sites_uniq_mm_plus_bw)){
        input_P_sites_uniq_mm_mn<-GRanges()
        input_P_sites_uniq_mm_pl<-import(path_to_P_sites_uniq_mm_plus_bw)
        strand(input_P_sites_uniq_mm_pl)<-"+"
        if(!is.na(path_to_P_sites_uniq_mm_minus_bw)){
            input_P_sites_uniq_mm_mn<-import(path_to_P_sites_uniq_mm_minus_bw)
            strand(input_P_sites_uniq_mm_mn)<-"-"
        }
        suppressWarnings(input_P_sites_uniq<-sort(c(input_P_sites_uniq_mm_pl,input_P_sites_uniq_mm_mn)))
        seqlevels(input_P_sites_uniq,pruning.mode="coarse")<-seqllll
        seqlengths(input_P_sites_uniq)<-seqleee
        
    }
    
    reduc<-function(x,y){
        all_ps<-x[["P_sites_all"]]
        uniq_ps<-x[["P_sites_uniq"]]
        uniq_mm_ps<-x[["P_sites_uniq_mm"]]
        
        if(!is.na(path_to_rl_cutoff_file)){
            all_ps<-GRangesList()
            rls<-unique(c(names(x[["P_sites_all"]]),names(y[["P_sites_all"]])))
            
            for(rl in rls){
                reads_x<-GRanges()
                reads_y<-GRanges()
                
                seqlevels(reads_x)<-seqllll
                seqlevels(reads_y)<-seqllll
                
                seqlengths(reads_x)<-seqleee
                seqlengths(reads_y)<-seqleee
                
                if(sum(rl%in%names(x[["P_sites_all"]]))>0){reads_x<-x[["P_sites_all"]][[rl]]}
                if(sum(rl%in%names(y[["P_sites_all"]]))>0){reads_y<-y[["P_sites_all"]][[rl]]}
                
                plx<-reads_x[strand(reads_x)=="+"]
                mnx<-reads_x[strand(reads_x)=="-"]
                ply<-reads_y[strand(reads_y)=="+"]
                mny<-reads_y[strand(reads_y)=="-"]
                if(length(plx)>0){covv_pl<-coverage(plx,weight = plx$score)}else{covv_pl<-coverage(plx)}
                if(length(ply)>0){covv_pl<-covv_pl+coverage(ply,weight = ply$score)}
                
                covv_pl<-GRanges(covv_pl)
                covv_pl<-covv_pl[covv_pl$score>0]
                
                if(length(mnx)>0){covv_min<-coverage(mnx,weight = mnx$score)}else{covv_min<-coverage(mnx)}
                if(length(mny)>0){covv_min<-covv_min+coverage(mny,weight = mny$score)}
                
                covv_min<-GRanges(covv_min)
                covv_min<-covv_min[covv_min$score>0]
                
                strand(covv_pl)<-"+"
                strand(covv_min)<-"-"
                
                all_ps[[rl]]<-sort(c(covv_pl,covv_min))
                
            }
            
            uniq_ps<-GRangesList()
            rls<-unique(c(names(x[["P_sites_uniq"]]),names(y[["P_sites_uniq"]])))
            
            for(rl in rls){
                reads_x<-GRanges()
                reads_y<-GRanges()
                
                seqlevels(reads_x)<-seqllll
                seqlevels(reads_y)<-seqllll
                
                seqlengths(reads_x)<-seqleee
                seqlengths(reads_y)<-seqleee
                if(sum(rl%in%names(x[["P_sites_uniq"]]))>0){reads_x<-x[["P_sites_uniq"]][[rl]]}
                if(sum(rl%in%names(y[["P_sites_uniq"]]))>0){reads_y<-y[["P_sites_uniq"]][[rl]]}
                
                plx<-reads_x[strand(reads_x)=="+"]
                mnx<-reads_x[strand(reads_x)=="-"]
                ply<-reads_y[strand(reads_y)=="+"]
                mny<-reads_y[strand(reads_y)=="-"]
                if(length(plx)>0){covv_pl<-coverage(plx,weight = plx$score)}else{covv_pl<-coverage(plx)}
                if(length(ply)>0){covv_pl<-covv_pl+coverage(ply,weight = ply$score)}
                
                covv_pl<-GRanges(covv_pl)
                covv_pl<-covv_pl[covv_pl$score>0]
                
                if(length(mnx)>0){covv_min<-coverage(mnx,weight = mnx$score)}else{covv_min<-coverage(mnx)}
                if(length(mny)>0){covv_min<-covv_min+coverage(mny,weight = mny$score)}
                
                covv_min<-GRanges(covv_min)
                covv_min<-covv_min[covv_min$score>0]
                
                strand(covv_pl)<-"+"
                strand(covv_min)<-"-"
                
                uniq_ps[[rl]]<-sort(c(covv_pl,covv_min))
                
                
            }
            
            uniq_mm_ps<-GRangesList()
            rls<-unique(c(names(x[["P_sites_uniq_mm"]]),names(y[["P_sites_uniq_mm"]])))
            
            for(rl in rls){
                reads_x<-GRanges()
                reads_y<-GRanges()
                
                seqlevels(reads_x)<-seqllll
                seqlevels(reads_y)<-seqllll
                
                seqlengths(reads_x)<-seqleee
                seqlengths(reads_y)<-seqleee
                if(sum(rl%in%names(x[["P_sites_uniq_mm"]]))>0){reads_x<-x[["P_sites_uniq_mm"]][[rl]]}
                if(sum(rl%in%names(y[["P_sites_uniq_mm"]]))>0){reads_y<-y[["P_sites_uniq_mm"]][[rl]]}
                
                plx<-reads_x[strand(reads_x)=="+"]
                mnx<-reads_x[strand(reads_x)=="-"]
                ply<-reads_y[strand(reads_y)=="+"]
                mny<-reads_y[strand(reads_y)=="-"]
                if(length(plx)>0){covv_pl<-coverage(plx,weight = plx$score)}else{covv_pl<-coverage(plx)}
                if(length(ply)>0){covv_pl<-covv_pl+coverage(ply,weight = ply$score)}
                
                covv_pl<-GRanges(covv_pl)
                covv_pl<-covv_pl[covv_pl$score>0]
                
                if(length(mnx)>0){covv_min<-coverage(mnx,weight = mnx$score)}else{covv_min<-coverage(mnx)}
                if(length(mny)>0){covv_min<-covv_min+coverage(mny,weight = mny$score)}
                
                covv_min<-GRanges(covv_min)
                covv_min<-covv_min[covv_min$score>0]
                
                strand(covv_pl)<-"+"
                strand(covv_min)<-"-"
                
                uniq_mm_ps[[rl]]<-sort(c(covv_pl,covv_min))
                
            }
        }
        
        rang_jun<-x$junctions
        rang_jun$reads<-rang_jun$reads+y$junctions$reads
        rang_jun$unique_reads<-rang_jun$unique_reads+y$junctions$unique_reads
        
        list_res<-list(all_ps,uniq_ps,uniq_mm_ps,rang_jun)
        names(list_res)<-c("P_sites_all","P_sites_uniq","P_sites_uniq_mm","junctions")
        
        
        return(list_res)
    }
    
    #what to do with each chunk (read as alignment file)
    
    yiel<-function(x){
        readGAlignments(x,param = param)
    }
    
    #operations on the chunk (here count reads and whatnot)
    
    mapp<-function(x){
        
        mcols(x)$MD[which(is.na(mcols(x)$MD))]<-"NO"
        
        x<-x[seqnames(x)%in%seqnames(GTF_annotation$seqinfo)]
        seqlevels(x)<-seqlevels(GTF_annotation$seqinfo)
        x_I<-x[grep("I",cigar(x))]
        
        if(length(x_I)>0){
            x<-x[grep("I",cigar(x),invert=T)]
            
        }
        x_D<-x[grep("D",cigar(x))]
        if(length(x_D)>0){
            x<-x[grep("D",cigar(x),invert=T)]
            
        }
        
        # softclipping
        
        clipp <- width(cigarRangesAlongQuerySpace(x@cigar, ops="S"))
        clipp[elementNROWS(clipp)==0] <- 0
        len_adj <- qwidth(x)-sum(clipp)
        mcols(x)$len_adj <- len_adj
        
        # Remove S from Cigar (read positions/length are already adjusted)
        # it helps calculating P-sites positions for spliced reads
        
        cigg<-cigar(x)
        cigg_s<-grep(cigg,pattern = "S")
        if(length(cigg_s)>0){
            cigs<-cigg[cigg_s]
            cigs<-gsub(cigs,pattern = "^[0-9]+S",replacement = "")
            cigs<-gsub(cigs,pattern = "[0-9]+S$",replacement = "")
            cigg[cigg_s]<-cigs
            x@cigar<-cigg
        }
        mcols(x)$cigar_str<-x@cigar
        x_uniq<-x[x@elementMetadata$mapq>50]
        
        pos<-x[strand(x)=="+"]
        neg<-x[strand(x)=="-"]
        
        uniq_pos<-x_uniq[strand(x_uniq)=="+"]
        uniq_neg<-x_uniq[strand(x_uniq)=="-"]
        
        
        
        # junctions
        
        juns<-summarizeJunctions(x)
        juns_pos<-juns
        juns_neg<-juns
        mcols(juns_pos)<-NULL
        mcols(juns_neg)<-NULL
        juns_pos$reads<-juns$plus_score
        juns_neg$reads<-juns$minus_score
        strand(juns_pos)<-"+"
        strand(juns_neg)<-"-"
        juns<-sort(c(juns_pos,juns_neg))
        juns<-juns[juns$reads>0]
        
        uniq_juns<-summarizeJunctions(x_uniq)
        uniq_juns_pos<-uniq_juns
        uniq_juns_neg<-uniq_juns
        mcols(uniq_juns_pos)<-NULL
        mcols(uniq_juns_neg)<-NULL
        uniq_juns_pos$reads<-uniq_juns$plus_score
        uniq_juns_neg$reads<-uniq_juns$minus_score
        strand(uniq_juns_pos)<-"+"
        strand(uniq_juns_neg)<-"-"
        uniq_juns<-sort(c(uniq_juns_pos,uniq_juns_neg))
        uniq_juns<-uniq_juns[uniq_juns$reads>0]
        if(length(juns)>0){
            juns$unique_reads<-0
            mat<-match(uniq_juns,juns)
            juns$unique_reads[mat]<-uniq_juns$reads
        }
        rang_jun<-GTF_annotation$junctions
        rang_jun$reads<-0
        rang_jun$unique_reads<-0
        if(length(juns)>0){
            mat<-match(juns,rang_jun)
            juns<-juns[!is.na(mat)]
            mat<-mat[!is.na(mat)]
            rang_jun$reads[mat]<-juns$reads
            rang_jun$unique_reads[mat]<-juns$unique_reads
        }
        
        
        # P-sites calculation
        all_ps_comps<-GRanges()
        seqlevels(all_ps_comps)<-seqllll
        seqlengths(all_ps_comps)<-seqleee
        
        uniq_ps_comps<-GRanges()
        seqlevels(uniq_ps_comps)<-seqllll
        seqlengths(uniq_ps_comps)<-seqleee
        
        uniq_mm_ps_comps<-GRanges()
        seqlevels(uniq_mm_ps_comps)<-seqllll
        seqlengths(uniq_mm_ps_comps)<-seqleee
        
        if(!is.na(path_to_rl_cutoff_file)){
            list_pss<-list()
            for(comp in names(rl_cutoffs_comp)){
                all_rl_ps<-GRangesList()
                uniq_rl_ps<-GRangesList()
                uniq_rl_mm_ps<-GRangesList()
                
                seqlevels(all_rl_ps)<-seqllll
                seqlevels(uniq_rl_ps)<-seqllll
                seqlevels(uniq_rl_mm_ps)<-seqllll
                
                seqlengths(all_rl_ps)<-seqleee
                seqlengths(uniq_rl_ps)<-seqleee
                seqlengths(uniq_rl_mm_ps)<-seqleee
                
                
                chroms<-comp
                
                if(comp=="nucl"){chroms=seqlevels(x)[!seqlevels(x)%in%circs]}
                resul<-rl_cutoffs_comp[[comp]]
                
                for(i in seq_along(resul$read_length)){
                    
                    all_ps<-GRangesList()
                    uniq_ps<-GRangesList()
                    uniq_mm_ps<-GRangesList()
                    seqlevels(all_ps)<-seqllll
                    seqlevels(uniq_ps)<-seqllll
                    seqlevels(uniq_mm_ps)<-seqllll
                    
                    seqlengths(all_ps)<-seqleee
                    seqlengths(uniq_ps)<-seqleee
                    seqlengths(uniq_mm_ps)<-seqleee
                    
                    rl<-as.numeric(resul$read_length[i])
                    ct<-as.numeric(resul$cutoff[i])
                    ok_reads<-pos[mcols(pos)$len_adj%in%rl]
                    ok_reads<-ok_reads[as.vector(seqnames(ok_reads))%in%chroms]
                    
                    ps_plus<-GRanges()
                    seqlevels(ps_plus)<-seqllll
                    seqlengths(ps_plus)<-seqleee
                    ps_plus_uniq<-ps_plus
                    ps_plus_uniq_mm<-ps_plus
                    
                    if(length(ok_reads)>0){
                        unspl<-ok_reads[grep(pattern="N",x=cigar(ok_reads),invert=T)]
                        
                        ps_unspl<-shift(resize(GRanges(unspl),width=1,fix="start"),shift=ct)
                        
                        spl<-ok_reads[grep(pattern="N",x=cigar(ok_reads))]
                        firstb<-as.numeric(sapply(strsplit(cigar(spl),"M"),"[[",1))
                        lastb<-as.numeric(sapply(strsplit(cigar(spl),"M"),function(x){gsub(x[length(x)],pattern="^[^_]*N",replacement="")}))
                        firstok<-spl[firstb>ct]
                        firstok<-shift(resize(GRanges(firstok),width=1,fix="start"),shift=ct)
                        
                        lastok<-spl[lastb>=rl-ct]
                        lastok<-shift(resize(GRanges(lastok),width=1,fix="end"),shift=-(rl-ct-1))
                        
                        
                        multi<-spl[firstb<=ct & lastb<rl-ct]
                        
                        
                        ps_spl<-GRanges()
                        seqlevels(ps_spl)<-seqllll
                        seqlengths(ps_spl)<-seqleee
                        
                        if(length(multi)>0){
                            ps_spl<-get_ps_fromspliceplus(multi,cutoff=ct)
                            
                        }
                        mcols(ps_spl)<-mcols(multi)
                        
                        seqlevels(firstok)<-seqllll
                        seqlevels(lastok)<-seqllll
                        seqlevels(ps_unspl)<-seqllll
                        seqlevels(ps_spl)<-seqllll
                        
                        seqlengths(firstok)<-seqleee
                        seqlengths(lastok)<-seqleee
                        seqlengths(ps_unspl)<-seqleee
                        seqlengths(ps_spl)<-seqleee
                        
                        
                        ps_plus<-c(ps_unspl,firstok,lastok,ps_spl)
                        ps_plus_uniq<-ps_plus[mcols(ps_plus)$mapq>50]
                        ps_plus_uniq_mm<-ps_plus[mcols(ps_plus)$mapq>50 & nchar(mcols(ps_plus)$MD)>3]
                        
                        mcols(ps_plus)<-NULL
                        mcols(ps_plus_uniq)<-NULL
                        mcols(ps_plus_uniq_mm)<-NULL
                        
                    }
                    ok_reads<-neg[mcols(neg)$len_adj%in%rl]
                    ok_reads<-ok_reads[as.vector(seqnames(ok_reads))%in%chroms]
                    
                    ps_neg<-GRanges()
                    seqlevels(ps_neg)<-seqllll
                    seqlengths(ps_neg)<-seqleee
                    ps_neg_uniq<-ps_neg
                    ps_neg_uniq_mm<-ps_neg
                    
                    if(length(ok_reads)>0){
                        unspl<-ok_reads[grep(pattern="N",x=cigar(ok_reads),invert=T)]
                        
                        ps_unspl<-shift(resize(GRanges(unspl),width=1,fix="start"),shift=-ct)
                        
                        spl<-ok_reads[grep(pattern="N",x=cigar(ok_reads))]
                        
                        firstb<-as.numeric(sapply(strsplit(cigar(spl),"M"),"[[",1))
                        lastb<-as.numeric(sapply(strsplit(cigar(spl),"M"),function(x){gsub(x[length(x)],pattern="^[^_]*N",replacement="")}))
                        lastok<-spl[lastb>ct]
                        lastok<-shift(resize(GRanges(lastok),width=1,fix="start"),shift=-ct)
                        
                        firstok<-spl[firstb>=rl-ct]
                        firstok<-shift(resize(GRanges(firstok),width=1,fix="end"),shift=(rl-ct-1))
                        
                        multi<-spl[firstb<rl-ct & lastb<=ct]
                        
                        
                        ps_spl<-GRanges()
                        seqlevels(ps_spl)<-seqllll
                        seqlengths(ps_spl)<-seqleee
                        
                        
                        if(length(multi)>0){
                            ps_spl<-get_ps_fromsplicemin(multi,cutoff=ct)
                        }
                        mcols(ps_spl)<-mcols(multi)
                        
                        seqlevels(firstok)<-seqllll
                        seqlevels(lastok)<-seqllll
                        seqlevels(ps_unspl)<-seqllll
                        seqlevels(ps_spl)<-seqllll
                        
                        seqlengths(firstok)<-seqleee
                        seqlengths(lastok)<-seqleee
                        seqlengths(ps_unspl)<-seqleee
                        seqlengths(ps_spl)<-seqleee
                        
                        ps_neg<-c(ps_unspl,firstok,lastok,ps_spl)
                        ps_neg_uniq<-ps_neg[mcols(ps_neg)$mapq>50]
                        ps_neg_uniq_mm<-ps_neg[mcols(ps_neg)$mapq>50 & nchar(mcols(ps_neg)$MD)>3]
                        
                        mcols(ps_neg)<-NULL
                        mcols(ps_neg_uniq)<-NULL
                        mcols(ps_neg_uniq_mm)<-NULL
                    }
                    
                    all_ps<-sort(c(ps_plus,ps_neg))
                    uniq_ps<-sort(c(ps_plus_uniq,ps_neg_uniq))
                    uniq_mm_ps<-sort(c(ps_plus_uniq_mm,ps_neg_uniq_mm))
                    if(length(all_ps)>0){
                        ps_res<-unique(all_ps)
                        ps_res$score<-countOverlaps(ps_res,all_ps,type="equal")
                        all_ps<-ps_res
                    }
                    if(length(uniq_ps)>0){
                        ps_res<-unique(uniq_ps)
                        ps_res$score<-countOverlaps(ps_res,uniq_ps,type="equal")
                        uniq_ps<-ps_res
                        
                    }
                    if(length(uniq_mm_ps)>0){
                        ps_res<-unique(uniq_mm_ps)
                        ps_res$score<-countOverlaps(ps_res,uniq_mm_ps,type="equal")
                        uniq_mm_ps<-ps_res
                        
                    }
                    all_rl_ps[[as.character(rl)]]<-all_ps
                    uniq_rl_ps[[as.character(rl)]]<-uniq_ps
                    uniq_rl_mm_ps[[as.character(rl)]]<-uniq_mm_ps
                    
                    
                }
                #here comps
                list_rlct<-list(all_rl_ps,uniq_rl_ps,uniq_rl_mm_ps)
                names(list_rlct)<-c("P_sites_all","P_sites_uniq","P_sites_uniq_mm")
                list_pss[[comp]]<-list_rlct
            }
            #for rl, merge psites
            
            
            all_ps_comps<-GRangesList()
            seqlevels(all_ps_comps)<-seqllll
            seqlengths(all_ps_comps)<-seqleee
            rls_comps<-unique(unlist(lapply(list_pss,FUN=function(x) names(x[["P_sites_all"]]) )))
            for(rl in rls_comps){
                reads_rl_comp<-GRanges()
                seqlevels(reads_rl_comp)<-seqllll
                seqlengths(reads_rl_comp)<-seqleee
                for(comp in names(list_pss)){
                    if(sum(rl%in%names(list_pss[[comp]][["P_sites_all"]]))>0){
                        oth<-list_pss[[comp]][["P_sites_all"]][[rl]]
                        
                        if(!is.null(oth)){
                            seqlevels(oth)<-seqllll
                            seqlengths(oth)<-seqleee
                            reads_rl_comp<-c(reads_rl_comp,oth)
                        }
                    }          
                    all_ps_comps[[rl]]<-reads_rl_comp
                }
                
            }
            
            uniq_ps_comps<-GRangesList()
            seqlevels(uniq_ps_comps)<-seqllll
            seqlengths(uniq_ps_comps)<-seqleee
            rls_comps<-unique(unlist(lapply(list_pss,FUN=function(x) names(x[["P_sites_uniq"]]) )))
            for(rl in rls_comps){
                reads_rl_comp<-GRanges()
                seqlevels(reads_rl_comp)<-seqllll
                seqlengths(reads_rl_comp)<-seqleee
                for(comp in names(list_pss)){
                    if(sum(rl%in%names(list_pss[[comp]][["P_sites_uniq"]]))>0){
                        oth<-list_pss[[comp]][["P_sites_uniq"]][[rl]]
                        
                        if(!is.null(oth)){
                            seqlevels(oth)<-seqllll
                            seqlengths(oth)<-seqleee
                            reads_rl_comp<-c(reads_rl_comp,oth)
                        }
                    }          
                    uniq_ps_comps[[rl]]<-reads_rl_comp
                }
                
            }
            
            uniq_mm_ps_comps<-GRangesList()
            seqlevels(uniq_mm_ps_comps)<-seqllll
            seqlengths(uniq_mm_ps_comps)<-seqleee
            rls_comps<-unique(unlist(lapply(list_pss,FUN=function(x) names(x[["P_sites_uniq_mm"]]) )))
            for(rl in rls_comps){
                reads_rl_comp<-GRanges()
                seqlevels(reads_rl_comp)<-seqllll
                seqlengths(reads_rl_comp)<-seqleee
                for(comp in names(list_pss)){
                    if(sum(rl%in%names(list_pss[[comp]][["P_sites_uniq_mm"]]))>0){
                        oth<-list_pss[[comp]][["P_sites_uniq_mm"]][[rl]]
                        if(!is.null(oth)){
                            seqlevels(oth)<-seqllll
                            seqlengths(oth)<-seqleee
                            reads_rl_comp<-c(reads_rl_comp,oth)
                        }
                    }          
                    uniq_mm_ps_comps[[rl]]<-reads_rl_comp
                }
                
            }
            
        }
        
        list_res<-list(all_ps_comps,uniq_ps_comps,uniq_mm_ps_comps,rang_jun)
        names(list_res)<-c("P_sites_all","P_sites_uniq","P_sites_uniq_mm","junctions")
        
        return(list_res)
    }
    
    cat(paste("Calculating P-sites positions and junctions ...", date(),"\n"))
    
    for_ORFquant<-reduceByYield(X=opts,YIELD=yiel,MAP=mapp,REDUCE=reduc)
    
    if(length(for_ORFquant$P_sites_all)>0){
        merged_all_ps<-unlist(for_ORFquant$P_sites_all)
        
        if(length(merged_all_ps)>0){
            covv_pl<-coverage(merged_all_ps[strand(merged_all_ps)=="+"],weight = merged_all_ps[strand(merged_all_ps)=="+"]$score)
            covv_pl<-GRanges(covv_pl)
            covv_pl<-covv_pl[covv_pl$score>0]
            covv_min<-coverage(merged_all_ps[strand(merged_all_ps)=="-"],weight = merged_all_ps[strand(merged_all_ps)=="-"]$score)
            covv_min<-GRanges(covv_min)
            covv_min<-covv_min[covv_min$score>0]
            strand(covv_pl)<-"+"
            strand(covv_min)<-"-"
            
            merged_all_ps<-sort(c(covv_pl,covv_min))
            
        }
        for_ORFquant$P_sites_all<-merged_all_ps
    }
    
    
    if(length(for_ORFquant$P_sites_uniq)>0){
        merged_uniq_ps<-unlist(for_ORFquant$P_sites_uniq)
        if(length(merged_uniq_ps)>0){
            
            covv_pl<-coverage(merged_uniq_ps[strand(merged_uniq_ps)=="+"],weight = merged_uniq_ps[strand(merged_uniq_ps)=="+"]$score)
            covv_pl<-GRanges(covv_pl)
            covv_pl<-covv_pl[covv_pl$score>0]
            
            covv_min<-coverage(merged_uniq_ps[strand(merged_uniq_ps)=="-"],weight = merged_uniq_ps[strand(merged_uniq_ps)=="-"]$score)
            covv_min<-GRanges(covv_min)
            covv_min<-covv_min[covv_min$score>0]
            strand(covv_pl)<-"+"
            strand(covv_min)<-"-"
            merged_uniq_ps<-sort(c(covv_pl,covv_min))
            
        }
        for_ORFquant$P_sites_uniq<-merged_uniq_ps
    }
    
    if(length(for_ORFquant$P_sites_uniq_mm)>0){
        merged_uniq_mm_ps<-unlist(for_ORFquant$P_sites_uniq_mm)
        if(length(merged_uniq_mm_ps)>0){
            
            covv_pl<-coverage(merged_uniq_mm_ps[strand(merged_uniq_mm_ps)=="+"],weight = merged_uniq_mm_ps[strand(merged_uniq_mm_ps)=="+"]$score)
            covv_pl<-GRanges(covv_pl)
            covv_pl<-covv_pl[covv_pl$score>0]
            
            covv_min<-coverage(merged_uniq_mm_ps[strand(merged_uniq_mm_ps)=="-"],weight = merged_uniq_mm_ps[strand(merged_uniq_mm_ps)=="-"]$score)
            covv_min<-GRanges(covv_min)
            covv_min<-covv_min[covv_min$score>0]
            strand(covv_pl)<-"+"
            strand(covv_min)<-"-"
            merged_uniq_mm_ps<-sort(c(covv_pl,covv_min))
        }
        for_ORFquant$P_sites_uniq_mm<-merged_uniq_mm_ps
    }
    
    
    if(!is.na(path_to_P_sites_plus_bw) | !is.na(path_to_P_sites_minus_bw) ){
        for_ORFquant$P_sites_all<-input_P_sites
    }
    
    if(!is.na(path_to_P_sites_uniq_plus_bw) | !is.na(path_to_P_sites_uniq_minus_bw) ){
        for_ORFquant$P_sites_uniq<-input_P_sites_uniq
    }
    
    if(!is.na(path_to_P_sites_uniq_mm_plus_bw) | !is.na(path_to_P_sites_uniq_mm_minus_bw) ){
        for_ORFquant$P_sites_uniq_mm<-input_P_sites_uniq_mm
    }
    
    cat(paste("Calculating P-sites positions and junctions --- Done!", date(),"\n"))
    
    save(for_ORFquant,file = paste(dest_name,"for_ORFquant",sep = "_"))
    paste(dest_name,"for_ORFquant",sep = "_")
}




#' Plot general statistics about ORFquant results
#'
#' This function produces a series of plots and statistics about the set ORFs called by ORFquant compared to the annotation.
#' IMPORTANT: Use only on transcriptome-wide ORFquant results. See \code{run_ORFquant}
#' @keywords ORFquant, Ribo-seQC
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param for_ORFquant_file path to the "for_ORFquant" file containing P_sites positions and junction reads
#' @param ORFquant_output_file Full path to the "_final_ORFquant_results" RData object output by ORFquant. See \code{run_ORFquant}
#' @param annotation_file Full path to the *Rannot R file in the annotation directory used in the \code{prepare_annotation_files function}
#' @param coverage_file_plus Full path to a Ribo-seq coverage (no P-sites but read coverage) bigwig file (plus strand), as the ones created by \code{RiboseQC}
#' @param coverage_file_minus Full path to a Ribo-seq coverage (no P-sites but read coverage) bigwig file (minus strand), as the ones created by \code{RiboseQC}
#' @param output_plots_path Full path to the directory where plots in .pdf format are stored.
#' @param prefix prefix appended to output filenames
#' @return the function exports a RData object (*ORFquant_plots_RData) containing data to produce all plots, and produces different QC plots in .pdf format.
#' The plots created are as follows:\cr\cr
#' \code{ORFs_found}: Number of ORF categories detected per gene biotype.\cr
#' \code{ORFs_found_pct_tr}: Distribution of ORF_pct_P_sites (% of gene translation) for different ORF categories and gene biotypes.\cr
#' \code{ORFs_found_ORFs_pM}: Distribution of ORFs_pM (ORFs per Million, similar to TPM) for different ORF categories and gene biotypes.\cr
#' \code{ORFs_found_len}: Distribution of ORF length for different ORF categories and gene biotypes.\cr
#' \code{ORFs_genes}: Number of detected ORFs per gene.\cr
#' \code{ORFs_genes_tpm}: Gene level TPM values, plotted by number of ORFs detected.\cr
#' \code{ORFs_maxiso}:  Number of genes plotted against the percentages of gene translation of their most translated ORF.\cr
#' \code{ORFs_maxiso_tpm}: Gene level TPM values, plotted against the percentages of gene translation of their most translated ORF.\cr
#' \code{Sel_txs_genes}: Number of genes plotted against the number of selected transcripts.\cr
#' \code{Sel_txs_genes_tpm}: Gene level TPM values, plotted against the number of selected transcripts.\cr
#' \code{Sel_txs_genes_pct}: Percentages of annotated trascripts per gene, plotted against the number of selected transcripts.\cr
#' \code{Sel_txs_bins_juns}: Percentages of covered exonic bins or junctions, using all annotated transcripts, coding transcripts only, or the set of selected transcripts.\cr
#' \code{Meta_splicing_coverage}: Aggregate signal of Ribo-seq coverage and normalized ORF coverage across different splice sites combinations, with different mixtures of translated overlapping ORFs.
#' @seealso \code{\link{run_ORFquant}}
#' @export



plot_ORFquant_results<-function(for_ORFquant_file,ORFquant_output_file,annotation_file,coverage_file_plus=NA,coverage_file_minus=NA,output_plots_path=NA,prefix=NA){
    
    
    cat(paste("Plotting ORFquant results for ",ORFquant_output_file," ... ",date(),"\n",sep = ""))
    
    if(is.na(prefix)){prefix<-gsub(gsub(basename(ORFquant_output_file),pattern = "_final_ORFquant_results",replacement = ""),pattern = "final_ORFquant_results",replacement = "")}
    
    for (f in c(for_ORFquant_file,ORFquant_output_file,annotation_file)){
        if(file.access(f, 0)==-1) {
            stop("The following files don't exist:\n",
                 f, "\n")
        }
    }
    
    if(is.na(output_plots_path)){output_plots_path<-paste0(ORFquant_output_file, "_plots")}
    dir.create(output_plots_path,recursive=TRUE, showWarnings=FALSE)
    
    load(annotation_file)
    load(ORFquant_output_file)
    ORFs_tx<-ORFquant_results$ORFs_tx
    ORFs_gen<-ORFquant_results$ORFs_gen
    selected_txs<-ORFquant_results$selected_txs
    ORFs_txs_feats<-ORFquant_results$ORFs_txs_feats
    ORFs_feat<-ORFquant_results$ORFs_feat
    
    for_ORFquant<-get(load(for_ORFquant_file))
    
    list_ORFquant_plots<-list()
    
    
    gids<-ORFs_tx$gene_id
    cat_biot<-ORFs_tx$gene_biotype
    cat_biot[grep(cat_biot,pattern = "pseudo")]<-"pseudogene"
    cat_biot[!cat_biot%in%c("protein_coding","pseudogene")]<-"non-coding RNA"
    tbid<-table(gids)
    tbid[tbid>3]<-">3"
    #boxplot(log(tpms)~tbid)
    tb_bio<-cat_biot[match(names(tbid),gids)]
    
    
    cat_tx<-ORFs_tx$ORF_category_Tx_compatible
    cat_tx[cat_tx=="novel"]<-"not_annotated"
    cat_tx[cat_tx=="overl_uORF"]<-"uORF"
    cat_tx[cat_tx=="overl_dORF"]<-"other"
    cat_tx[cat_tx=="dORF"]<-"other"
    
    cat_gen<-ORFs_tx$ORF_category_Gen
    cat_biot<-ORFs_tx$compatible_biotype
    cat_biot_gene<-ORFs_tx$gene_biotype
    cat_biot[cat_biot!="protein_coding" & cat_biot_gene=="protein_coding"]<-"non-coding isoform"
    cat_biot[grep(cat_biot,pattern = "pseudo")]<-"pseudogene"
    cat_biot[!cat_biot%in%c("protein_coding","pseudogene","non-coding isoform")]<-"non-coding RNA"
    
    #cat_tx[intersect(grep(ORFs_tx$gene_biotype,pattern = "pseudogene"),which(cat_tx=="novel"))]<-"novel_pseudogene"
    cat_tx[cat_tx%in%c("C_extension","C_truncation","NC_extension","N_extension","nested_ORF","overl_dORF")]<-"other"
    
    cat_gen[grep(cat_gen,pattern = "novel")]<-"non-overlapping\nCDS regions"
    cat_gen[grep(cat_gen,pattern = "non-overlapping\nCDS regions",invert = T)]<-"overlapping\nCDS regions"
    levvs<-c("ORF_annotated","N_truncation","uORF","other","not_annotated")
    
    df<-melt(table(cat_tx,cat_biot))
    
    df$value[df$value==0]<-NA
    df$cat_tx<-factor(df$cat_tx,levels=levvs)
    df$cat_biot<-factor(df$cat_biot,levels=c("protein_coding","non-coding isoform","pseudogene","non-coding RNA"))
    
    
    colli<-c("red","orange","dark blue","cornflowerblue")
    
    
    dfiso<-data.frame(cat_tx,cat_biot,ORFs_tx$ORF_pct_P_sites)
    dfiso$cat_tx<-factor(dfiso$cat_tx,levels=levvs)
    dfiso$cat_biot<-factor(dfiso$cat_biot,levels=c("protein_coding","non-coding isoform","pseudogene","non-coding RNA"))
    
    dfpnpm<-data.frame(cat_tx,cat_biot,ORFs_tx$ORFs_pM)
    #dfpnpm$value[dfpnpm$value==0]<-NA
    dfpnpm$cat_tx<-factor(dfpnpm$cat_tx,levels=levvs)
    dfpnpm$cat_biot<-factor(dfpnpm$cat_biot,levels=c("protein_coding","non-coding isoform","pseudogene","non-coding RNA"))
    
    dfwid<-data.frame(cat_tx,cat_biot,width(ORFs_tx))
    dfwid$cat_tx<-factor(dfwid$cat_tx,levels=levvs)
    dfwid$cat_biot<-factor(dfwid$cat_biot,levels=c("protein_coding","non-coding isoform","pseudogene","non-coding RNA"))
    
    df<-df[!is.na(df$value),]
    
    
    a<-ggplot(df,aes(x=cat_biot,y=value,fill=cat_biot))
    a<-a + geom_bar(stat="identity",position = "dodge",colour="black")
    a<-a + facet_grid(. ~ cat_tx ,drop = T,scales = "free_x")
    a<-a + ylim(0,max(df$value,na.rm=T)*1.1)
    a<-a + geom_text(aes(x=cat_biot, y=value, hjust=.5,vjust=-.2,label=value),check_overlap = TRUE,colour="black",position = position_dodge(width = .9),size=3.5)
    a<-a + theme_bw()
    a<-a + ylab("n of ORFs")
    a<-a + xlab("")
    a<-a + scale_fill_manual(values = colli,"biotype")
    a<-a + theme(axis.title.x = element_blank(),axis.text.x  = element_blank())
    a<-a + theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=45, vjust=0.5, size=15))
    orfs_found_n<-a + theme(strip.text.x = element_text(size=14, face="bold"),strip.text.y = element_text(size=9),strip.background = element_rect(colour="black", fill=c(rep("darkkhaki",8),rep("red",8))))
    
    list_ORFquant_plots[["ORFs_found"]]<-orfs_found_n
    list_ORFquant_plots[["ORFs_found"]][["pars"]]<-c(14,4)
    
    b<-ggplot(dfiso,aes(x=cat_biot,y=ORFs_tx.ORF_pct_P_sites,fill=cat_biot))
    b<-b + geom_violin(scale="width",draw_quantiles=.5,adjust = 1)
    b<-b + theme_bw()
    b<-b + ylab("ORF_pct_P-sites")
    b<-b + xlab("")
    b<-b + facet_grid(. ~  cat_tx ,drop = T,scales = "free_x")
    #b<-b + theme(legend.position="none")
    b<-b + scale_fill_manual(values = colli,"biotype")
    b<-b + theme(axis.title.x = element_blank(),axis.text.x  = element_blank())
    b<-b + theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=45, vjust=0.5, size=15))
    orfs_found_pct<-b + theme(strip.text.x = element_text(size=14, face="bold"),strip.text.y = element_text(size=9),strip.background = element_rect(colour="black", fill=c(rep("darkkhaki",8),rep("red",8))))
    
    list_ORFquant_plots[["ORFs_found_pct_tr"]]<-orfs_found_pct
    list_ORFquant_plots[["ORFs_found_pct_tr"]][["pars"]]<-c(14,4)
    
    
    c<-ggplot(dfpnpm,aes(x=cat_biot,y=ORFs_tx.ORFs_pM+1,fill=cat_biot))
    c<-c + geom_violin(scale="width",draw_quantiles=.5,adjust = 1)
    #c<-c + geom_jitter(aes(x=cat_tx,y=ORFs_tx.ORFs_pM),position = position_jitterdodge(jitter.width = .5), alpha = 0.2)
    c<-c + scale_y_log10(breaks=c(1,11,101,1001),limits=c(1,max(dfpnpm$ORFs_tx.ORFs_pM)*1.1),labels=c(1,11,101,1001)-1)
    #c<-c + geom_text(aes(x=cat_tx, y=ORFs_tx.ORFs_pM, hjust="top",label=ORFs_tx.ORFs_pM),colour="black",position = position_dodge(width = 1),size=5)
    c<-c + facet_grid(. ~ cat_tx ,drop = T,scales = "free_x")
    c<-c + theme_bw()
    c<-c + xlab("")
    c<-c + ylab("P-sites_pNpM")
    c<-c + theme(axis.title.x = element_blank(),axis.text.x = element_blank())
    #c<-c + theme(legend.position="none")
    c<-c + scale_fill_manual(values = colli,"biotype")
    c<-c + theme(axis.title.x = element_blank(),axis.text.x  = element_blank())
    c<-c + theme(axis.title.y = element_text(size=18),axis.text.y  = element_text(angle=45, vjust=0.5, size=16))
    c<-c + theme(strip.text.x = element_text(size=16, face="bold"),strip.text.y = element_text(size=20),strip.background = element_rect(colour="black", fill=c("darkkhaki")))
    orfs_found_pspn<-c + theme(strip.text.x = element_text(size=14, face="bold"),strip.text.y = element_text(size=9),strip.background = element_rect(colour="black", fill=c(rep("darkkhaki",8),rep("red",8))))
    
    list_ORFquant_plots[["ORFs_found_ORFs_pM"]]<-orfs_found_pspn
    list_ORFquant_plots[["ORFs_found_ORFs_pM"]][["pars"]]<-c(14,4)
    
    d<-ggplot(dfwid,aes(x=cat_biot,y=width.ORFs_tx.,fill=cat_biot))
    d<-d + geom_violin(scale="width",draw_quantiles=.5,adjust = 1)
    #c<-c + geom_jitter(aes(x=cat_tx,y=ORFs_tx.ORFs_pM),position = position_jitterdodge(jitter.width = .5), alpha = 0.2)
    d<-d + scale_y_log10(breaks=c(10,100,1000,10000),labels=c(10,100,1000,10000),limits=c(min(dfwid$width.ORFs_tx.),max(dfwid$width.ORFs_tx.)*1.1))
    #c<-c + geom_text(aes(x=cat_tx, y=ORFs_tx.ORFs_pM, hjust="top",label=ORFs_tx.ORFs_pM),colour="black",position = position_dodge(width = 1),size=5)
    d<-d + theme_bw()
    d<-d + facet_grid(. ~ cat_tx ,drop = T,scales = "free_x")
    d<-d + xlab("")
    d<-d + ylab("ORF length (nt)")
    d<-d + theme(axis.title.x = element_blank(),axis.text.x = element_blank())
    #c<-c + theme(legend.position="none")
    d<-d + scale_fill_manual(values = colli,"biotype")
    d<-d + theme(axis.title.x = element_blank(),axis.text.x  = element_blank())
    d<-d + theme(axis.title.y = element_text(size=20),axis.text.y  = element_text(angle=45, vjust=0.5, size=16))
    orfs_found_len<-d + theme(strip.text.x = element_text(size=14, face="bold"),strip.text.y = element_text(size=9),strip.background = element_rect(colour="black", fill=c(rep("darkkhaki",8),rep("red",8))))
    
    list_ORFquant_plots[["ORFs_found_len"]]<-orfs_found_len
    list_ORFquant_plots[["ORFs_found_len"]][["pars"]]<-c(14,4)
    
    
    red_ex <- unlist(GTF_annotation$exons_txs)
    red_ex$gene_id<-GTF_annotation$trann$gene_id[match(names(red_ex),GTF_annotation$trann$transcript_id)]
    red_ex<-reduce(split(red_ex,red_ex$gene_id))
    
    ovv<-findOverlaps(red_ex,for_ORFquant$P_sites_all)
    aggo<-aggregate(for_ORFquant$P_sites_all$score[ovv@to],by=list(ovv@from),sum)
    cnts<-aggo[,2]
    names(cnts)<-names(red_ex)[aggo[,1]]
    cnts0<-rep(0,sum(!names(red_ex)%in%names(cnts)))
    names(cnts0)<-names(red_ex)[!names(red_ex)%in%names(cnts)]
    cnts<-c(cnts,cnts0)[names(red_ex)]
    
    cnts<-DataFrame(gene_id=names(cnts),P_sites=as.numeric(cnts))
    cnts<-cnts[match(cnts$gene_id,names(red_ex)),]
    cnts$RPKM<-cnts$P_sites/(sum(cnts$P_sites)/1e06)
    cnts$RPKM<-round(cnts$RPKM/(sum(width(red_ex))/1000),digits = 4)
    cnts$TPM<-cnts$P_sites/(sum(width(red_ex))/1000)
    cnts$TPM<-round(cnts$TPM/(sum(cnts$TPM)/1e06),digits = 4)
    
    
    tpms<-cnts[match(names(tbid),cnts$gene_id),"TPM"]
    gids<-ORFs_tx$gene_id
    cat_biot<-ORFs_tx$gene_biotype
    cat_biot[grep(cat_biot,pattern = "pseudo")]<-"pseudogene"
    cat_biot[!cat_biot%in%c("protein_coding","pseudogene")]<-"non-coding RNA"
    tbid<-table(gids)
    tbid[tbid>3]<-">3"
    #boxplot(log(tpms)~tbid)
    tb_bio<-cat_biot[match(names(tbid),gids)]
    max_ORF<-aggregate(ORFs_tx$ORF_pct_P_sites,by=list(ORFs_tx$gene_id),max)
    colnames(max_ORF)<-c("gene_id","ORF_pct_P_sites")
    
    multg<-names(which(tbid!="1"))
    #mult_max_ORF<-max_ORF[max_ORF$gene_id%in%multg,]
    mult_max_ORF<-max_ORF
    
    maxiso<-cut(mult_max_ORF$ORF_pct_P_sites,breaks = seq(0,100,by = 10),include.lowest = T,right = T)
    maxiso[is.na(maxiso)]<-"(90,100]"
    maxiso<-gsub(maxiso,pattern = ",",replacement = "-")
    maxiso<-gsub(maxiso,pattern = "\\[",replacement = "")
    maxiso<-gsub(maxiso,pattern = "]",replacement = "")
    maxiso<-gsub(maxiso,pattern = "\\(",replacement = "")
    tpms_mult<-tpms[mult_max_ORF$gene_id]
    #boxplot(log(tpms_mult+1)~maxiso)
    df<-data.frame(tbid,tb_bio,tpms)
    
    df$Freq<-factor(df$Freq,levels=c("1","2","3",">3"))
    
    df$tb_bio<-factor(df$tb_bio,levels=c("protein_coding","pseudogene","non-coding RNA"))
    df<-table(df$Freq)
    df<-melt(df)
    
    df$Var1<-factor(df$Var1,levels=c("1","2","3",">3"))
    
    colli<-c("red","dark blue","cornflowerblue")
    a<-ggplot(df,aes(x=Var1,y=value,fill="dark grey"))
    a<-a + geom_bar(stat="identity",position = "stack",colour="black")
    #a<-a + scale_y_log10(breaks=c(1,10,100,1000,10000),limits=c(1,15000))
    a<-a + geom_text(aes(x=Var1, y=value, hjust=.5,vjust=-.8,label=value),check_overlap = TRUE,colour="black",size=4)
    a<-a + theme_bw()
    a<-a + theme(legend.position="none")
    a<-a  + ylim(0,max(df$value)*1.1)
    a<-a + ylab("n of genes")
    a<-a + xlab("n of ORFs")
    a<-a + scale_fill_manual(values = "dark grey","biotype")
    a<-a + theme(axis.title.y = element_text(size=20),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))
    orfsn_genes<-a + theme(axis.title.x = element_text(size=20),axis.text.x  = element_text(angle=45, vjust=0.5, size=18))
    
    list_ORFquant_plots[["ORFs_genes"]]<-orfsn_genes
    list_ORFquant_plots[["ORFs_genes"]][["pars"]]<-c(7,6)
    
    
    df<-data.frame(tbid,tb_bio,tpms)
    df$tb_bio<-factor(df$tb_bio,levels=names(sort(table(df$tb_bio),decreasing = T)))
    df$Freq<-factor(df$Freq,levels=c("1","2","3",">3"))
    #df<-df[df$Freq!="1",]
    
    b<-ggplot(df,aes(x=Freq,y=tpms+1,fill="dark grey"))
    b<-b + geom_violin(scale="width",draw_quantiles=.5,adjust = 1)
    b<-b + scale_y_log10(breaks=c(1,11,101,1001),limits=c(1,max(df$tpms+1)*1.1),labels=c(0,10,100,1000))
    #b<-b + geom_text(aes(x=cat_tx, y=value, hjust="top",label=value),colour="black",position = position_dodge(width = 1),size=5)
    b<-b + theme_bw()
    b<-b + xlab("n of ORFs")
    b<-b + theme(legend.position="none")
    
    b<-b + ylab("TPM (gene level)")
    b<-b + scale_fill_manual(values ="dark grey")
    b<-b + theme(axis.title.y = element_text(size=20),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))
    orfsn_genes_tpm<-b + theme(axis.title.x = element_text(size=20),axis.text.x  = element_text(angle=45, vjust=0.5, size=20))
    
    list_ORFquant_plots[["ORFs_genes_tpm"]]<-orfsn_genes_tpm
    list_ORFquant_plots[["ORFs_genes_tpm"]][["pars"]]<-c(7,6)
    
    
    df<-melt(table(maxiso))
    df$maxiso<-factor(df$maxiso,levels = rev(df$maxiso))
    c<-ggplot(df,aes(x=maxiso,y=value,fill="dark grey"))
    c<-c + geom_bar(stat="identity",position = "stack",colour="black")
    #a<-a + scale_y_log10(breaks=c(1,10,100,1000,10000),limits=c(1,15000))
    c<-c + geom_text(aes(x=maxiso, y=value, hjust=.5,vjust=-.8,label=value),check_overlap = TRUE,colour="black",size=4)
    c<-c + theme_bw()
    c<-c + theme(legend.position="none")
    c<-c + ylim(0,max(df$value)*1.1)
    c<-c + ylab("n of genes")
    c<-c + xlab("ORF_pct_P_sites (max ORF)")
    c<-c + scale_fill_manual(values = "dark grey")
    c<-c + theme(axis.title.y = element_text(size=20),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))
    orfsn_genes_maxiso<-c + theme(axis.title.x = element_text(size=20),axis.text.x  = element_text(angle=45, vjust=0.5, size=20))
    
    list_ORFquant_plots[["ORFs_maxiso"]]<-orfsn_genes_maxiso
    list_ORFquant_plots[["ORFs_maxiso"]][["pars"]]<-c(7,6)
    
    
    df<-data.frame(tpms_mult,maxiso)
    df$maxiso<-factor(df$maxiso,levels = rev(levels(df$maxiso)))
    
    d<-ggplot(df,aes(x=maxiso,y=tpms_mult+1,fill="dark grey"))
    d<-d + geom_violin(scale="width",draw_quantiles=.5,adjust = 1)
    d<-d + scale_y_log10(breaks=c(1,11,101,1001),limits=c(1,max(tpms_mult+1)*1.1),labels=c(0,10,100,1000))
    #b<-b + geom_text(aes(x=cat_tx, y=value, hjust="top",label=value),colour="black",position = position_dodge(width = 1),size=5)
    d<-d + theme_bw()
    d<-d + theme(legend.position="none")
    
    d<-d + xlab("ORF_pct_P_sites (max ORF)")
    d<-d + ylab("TPM (gene level)")
    d<-d + theme(axis.title.x = element_blank(),axis.text.x = element_blank())
    d<-d + scale_fill_manual(values ="dark grey")
    
    d<-d + theme(axis.title.x = element_text(size=20),axis.text.x  = element_text(angle=45, vjust=0.5, size=20))
    orfsn_genes_maxiso_tpm<-d + theme(axis.title.y = element_text(size=20),axis.text.y  = element_text(angle=45, vjust=0.5, size=18))
    
    list_ORFquant_plots[["ORFs_maxiso_tpm"]]<-orfsn_genes_maxiso_tpm
    list_ORFquant_plots[["ORFs_maxiso_tpm"]][["pars"]]<-c(7,6)
    
    
    
    ch_txs_sel<-CharacterList(split(selected_txs,GTF_annotation$trann$gene_id[match(selected_txs,GTF_annotation$trann$transcript_id)]))
    tot_n_tx<-table(GTF_annotation$trann$gene_id)[names(ch_txs_sel)]
    genes_m<-names(which(tot_n_tx>1))
    genes_m<-genes_m[genes_m%in%ORFs_tx$gene_id]
    #genes_m<-unique(ORFs_tx$gene_id)
    tot_n_tx<-tot_n_tx[genes_m]
    ch_txs_sel<-ch_txs_sel[genes_m]
    
    tpms<-cnts[match(names(ch_txs_sel),cnts$gene_id),"TPM"]
    pct_sel<-(elementNROWS(ch_txs_sel))/tot_n_tx*100
    nsels<-elementNROWS(ch_txs_sel)
    qnt<-cut(nsels,breaks = c(0,3,6,9,max(nsels)),include.lowest = T)
    qnt<-gsub(qnt,pattern = ",",replacement = "-")
    qnt<-gsub(qnt,pattern = "\\[",replacement = "")
    qnt<-gsub(qnt,pattern = "]",replacement = "")
    qnt<-gsub(qnt,pattern = "\\(",replacement = "")
    qnt[qnt=="0-3"]<-"1-3"
    qnt[qnt=="3-6"]<-"4-6"
    qnt[qnt=="6-9"]<-"7-9"
    qnt[qnt=="9-49"]<-"10-49"
    qnt<-factor(qnt,levels = names(sort(table(qnt),decreasing = T)))
    df<-melt(table(qnt))
    
    a<-ggplot(df,aes(x=qnt,y=value,fill="dark grey"))
    a<-a + geom_bar(stat="identity",position = "stack",colour="black")
    #a<-a + scale_y_log10(breaks=c(1,10,100,1000,10000),limits=c(1,15000))
    a<-a + geom_text(aes(x=qnt, y=value, hjust=.5,vjust=-.8,label=value),check_overlap = TRUE,colour="black",size=4)
    a<-a + theme_bw()
    a<- a + ylim(0,max(df$value)*1.1)
    a<-a + theme(legend.position="none")
    a<-a + ylab("n of genes")
    a<-a + xlab("selected Txs per gene")
    a<-a + scale_fill_manual(values = "dark grey","biotype")
    a<-a + theme(axis.title.y = element_text(size=16),axis.text.y  = element_text(angle=45, vjust=0.5, size=16))
    sel_txs_genes<-a + theme(axis.title.x = element_text(size=16),axis.text.x  = element_text(angle=45, vjust=0.5, size=16))
    
    list_ORFquant_plots[["Sel_txs_genes"]]<-sel_txs_genes
    list_ORFquant_plots[["Sel_txs_genes"]][["pars"]]<-c(7,6)
    
    df<-data.frame(tpm=tpms[names(ch_txs_sel)],elementNROWS(ch_txs_sel))
    df$qnt<-qnt
    b<-ggplot(df,aes(x=qnt,y=tpms+1,fill="dark grey"))
    b<-b + geom_violin(scale="width",draw_quantiles=.5,adjust = 1)
    b<-b + scale_y_log10(breaks=c(1,11,101,1001,10001),limits=c(1,max(df$tpm+1)*1.1),labels=c(0,10,100,1000,10000))
    #b<-b + geom_text(aes(x=cat_tx, y=value, hjust="top",label=value),colour="black",position = position_dodge(width = 1),size=5)
    b<-b + theme_bw()
    b<-b + xlab("selected Txs per gene")
    b<-b + theme(legend.position="none")
    b<-b + ylab("TPM (gene)")
    b<-b + scale_fill_manual(values ="dark grey")
    b<-b + theme(axis.title.y = element_text(size=16),axis.text.y  = element_text(angle=45, vjust=0.5, size=16))
    sel_txs_genes_tpm<-b + theme(axis.title.x = element_text(size=16),axis.text.x  = element_text(angle=45, vjust=0.5, size=16))
    
    list_ORFquant_plots[["Sel_txs_genes_tpm"]]<-sel_txs_genes_tpm
    list_ORFquant_plots[["Sel_txs_genes_tpm"]][["pars"]]<-c(7,6)
    
    
    df<-data.frame(pct_sel=as.numeric(pct_sel),elementNROWS(ch_txs_sel))
    df$qnt<-qnt
    c<-ggplot(df,aes(x=qnt,y=pct_sel,fill="dark grey"))
    c<-c + geom_violin(scale="width",draw_quantiles=.5,adjust = 1)
    #b<-b + scale_y_log10(breaks=c(1,11,101,1001,10001),limits=c(1,10001),labels=c(0,10,100,1000,10000))
    #b<-b + geom_text(aes(x=cat_tx, y=value, hjust="top",label=value),colour="black",position = position_dodge(width = 1),size=5)
    c<-c + theme_bw()
    c<-c + xlab("selected Txs per gene")
    c<-c + theme(legend.position="none")
    c<-c + ylab("% annotated Txs")
    c<-c + scale_fill_manual(values ="dark grey")
    c<-c + theme(axis.title.y = element_text(size=16),axis.text.y  = element_text(angle=45, vjust=0.5, size=16))
    sel_txs_genes_pct<-c + theme(axis.title.x = element_text(size=18),axis.text.x  = element_text(angle=45, vjust=0.5, size=16))
    
    list_ORFquant_plots[["Sel_txs_genes_pct"]]<-sel_txs_genes_pct
    list_ORFquant_plots[["Sel_txs_genes_pct"]][["pars"]]<-c(7,6)
    
    gens_sel<-unique(GTF_annotation$trann$gene_id[GTF_annotation$trann$transcript_id%in%selected_txs])
    tx_all<-GTF_annotation$trann$transcript_id[GTF_annotation$trann$gene_id%in%gens_sel]
    tx_sel<-selected_txs
    if(sum(!tx_sel%in%names(GTF_annotation$exons_txs))>0){warning(paste("some selected Txs not present in annotation!, like:",head(tx_sel[!tx_sel%in%names(GTF_annotation$exons_txs)],1)))}
    tx_sel<-tx_sel[tx_sel%in%names(GTF_annotation$exons_txs)]
    
    tx_cds<-names(GTF_annotation$cds_txs)
    #NOCDS, USE ALL
    b<-disjoin(unlist(reduce(GTF_annotation$exons_txs[GTF_annotation$trann$transcript_id[GTF_annotation$trann$gene_id%in%gens_sel]])),with.revmap=T)
    b$ps<-assay(summarizeOverlaps(b,reads = for_ORFquant$P_sites_all,ignore.strand=F,mode="Union",inter.feature=FALSE))
    
    #b<-b[b%over%GTF_annotation$exons_txs[tx_sel]]
    
    b<-b[as.vector(b$ps>0)]
    
    ov<-findOverlaps(b,GTF_annotation$exons_txs[tx_sel])
    ov<-split(subjectHits(ov),queryHits(ov))
    b$tx_sel<-0
    b$tx_sel[as.numeric(names(ov))]<-elementNROWS(ov)
    
    ov<-findOverlaps(b,GTF_annotation$exons_txs[tx_all])
    ov<-split(subjectHits(ov),queryHits(ov))
    b$tx_all<-0
    b$tx_all[as.numeric(names(ov))]<-elementNROWS(ov)
    
    ov<-findOverlaps(b,GTF_annotation$exons_txs[tx_cds])
    ov<-split(subjectHits(ov),queryHits(ov))
    b$tx_cds<-0
    b$tx_cds[as.numeric(names(ov))]<-elementNROWS(ov)
    
    ov<-findOverlaps(b,GTF_annotation$genes[gens_sel])
    ov<-split(subjectHits(ov),queryHits(ov))
    b$gene<-ov
    nms<-names(GTF_annotation$genes)
    b$gene_id<-CharacterList(lapply(b$gene,function(x){nms[x]}))
    nts_all<-c()
    nts_cds<-c()
    nts_sel<-c()
    
    
    for(i in c(0,1,2,3,4,5,"more")){
        if(i!="more"){
            nts_all<-c(nts_all,length(b[(b$tx_all)==i]))
            nts_cds<-c(nts_cds,length(b[(b$tx_cds)==i]))
            nts_sel<-c(nts_sel,length(b[(b$tx_sel)==i]))
        }
        if(i=="more"){
            nts_all<-c(nts_all,length(b[(b$tx_all)>5]))
            nts_cds<-c(nts_cds,length(b[(b$tx_cds)>5]))
            nts_sel<-c(nts_sel,length(b[(b$tx_sel)>5]))
        }
    }
    
    juns<-for_ORFquant$junctions
    juns<-juns[juns$reads>0]
    juns$tx_cds<-match(juns$tx_name,tx_cds)
    juns$tx_sel<-match(juns$tx_name,tx_sel)
    juns$tx_cds<-CharacterList(lapply(juns$tx_cds,function(x){x[!is.na(x)]}))
    juns$tx_sel<-CharacterList(lapply(juns$tx_sel,function(x){x[!is.na(x)]}))
    #juns_ok<-juns[elementNROWS(juns$tx_sel)>0 & elementNROWS(juns$tx_cds)>0]
    juns_ok<-juns[sapply(juns$gene_id,function(x){sum(x%in%gens_sel)>0})]
    
    juns_all<-c()
    juns_cds<-c()
    juns_sel<-c()
    for(i in c(0,1,2,3,4,5,"more")){
        if(i!="more"){
            juns_all<-c(juns_all,sum((elementNROWS(juns_ok$tx_name)==i)))
            juns_cds<-c(juns_cds,sum((elementNROWS(juns_ok$tx_cds)==i)))
            juns_sel<-c(juns_sel,sum((elementNROWS(juns_ok$tx_sel)==i)))
        }
        if(i=="more"){
            juns_all<-c(juns_all,sum((elementNROWS(juns_ok$tx_name)>5)))
            juns_cds<-c(juns_cds,sum((elementNROWS(juns_ok$tx_cds)>5)))
            juns_sel<-c(juns_sel,sum((elementNROWS(juns_ok$tx_sel)>5)))     
        }
    }
    
    juns_all<-round(juns_all/sum(juns_all)*100,digits = 2)
    juns_cds<-round(juns_cds/sum(juns_cds)*100,digits = 2)
    juns_sel<-round(juns_sel/sum(juns_sel)*100,digits = 2)
    
    nts_all<-round(nts_all/sum(nts_all)*100,digits = 2)
    nts_cds<-round(nts_cds/sum(nts_cds)*100,digits = 2)
    nts_sel<-round(nts_sel/sum(nts_sel)*100,digits = 2)
    
    
    
    juns_nt_mapping<-data.frame(t(rbind(juns_all,juns_cds,juns_sel)))
    juns_nt_mapping<-suppressMessages(melt(juns_nt_mapping))
    juns_nt_mapping$type="covered junctions"
    juns_nt_mapping$n_tx<-rep(c(0:5,">5"),3)
    juns_nt_mapping$ref<-c(rep("All Txs",7),rep("Coding Txs",7),rep("Selected Txs",7))
    cds_nt_mapping<-data.frame(t(rbind(nts_all,nts_cds,nts_sel)))
    cds_nt_mapping<-suppressMessages(melt(cds_nt_mapping))
    cds_nt_mapping$type="covered exon bins"
    cds_nt_mapping$n_tx<-rep(c(0:5,">5"),3)
    cds_nt_mapping$ref<-c(rep("All Txs",7),rep("Coding Txs",7),rep("Selected Txs",7))
    
    
    
    df<-rbind(cds_nt_mapping,juns_nt_mapping)
    df$n_tx<-factor(df$n_tx,levels=c("0",">5","5","4","3","2","1"))
    colss<-rev(c(colorRampPalette(colors = c("black","white"))(6),"red"))
    a<-ggplot(df,aes(x=ref,y=value,fill=n_tx))
    a<-a + geom_bar(stat="identity",position = "stack",colour="black")
    #a<-a + scale_y_log10(breaks=c(1,10,100,1000,10000),limits=c(1,15000))
    a<-a + theme_bw()
    a<-a + ylab("percentage")
    a<-a + xlab("")
    a<-a + facet_grid(type~.,scale="free",)
    a<-a + scale_fill_manual(values = colss,"mapping Txs")
    a<-a + theme(axis.title.y = element_text(size=16),axis.text.y  = element_text(angle=45, vjust=0.5, size=16),strip.text.y =  element_text(size=12))
    selection_bins_juns<-a + theme(axis.title.x = element_text(size=12),axis.text.x  = element_text(angle=45, vjust=0.5, size=16))
    
    list_ORFquant_plots[["Sel_txs_bins_juns"]]<-selection_bins_juns
    list_ORFquant_plots[["Sel_txs_bins_juns"]][["pars"]]<-c(5,4.5)
    
    
    
    if(!is.na(coverage_file_plus) & !is.na(coverage_file_minus)){
        
        cat(paste("Plotting alternative splice sites profiles ... ",date(),"\n",sep = ""))
        
        aggregate_regions<-function(range,coverage_plus,coverage_minus,norm_x=NA,norm_y=NA,nozero=F){
            pl<-range[strand(range)=="+"]
            min<-range[strand(range)=="-"]
            mat_pl<-t(sapply(coverage_plus[pl],as.vector))
            if(nozero){
                mat_pl<-mat_pl[rowSums(mat_pl)>0,]
            }
            if(!is.na(norm_x)){
                mat_pl<-mat_pl/(rowSums(mat_pl)/norm_x)
            }
            if(!is.na(norm_y)){
                mat_pl<-mat_pl/(colSums(mat_pl)/norm_y)
            }
            mat_min<-t(sapply(coverage_minus[min],as.vector))
            if(nozero){
                mat_min<-mat_min[rowSums(mat_min)>0,]
            }
            if(!is.na(norm_x)){
                mat_min<-mat_min/(rowSums(mat_min)/norm_x)
            }
            if(!is.na(norm_y)){
                mat_min<-mat_min/(colSums(mat_min)/norm_y)
            }
            
            mat<-rbind(mat_pl,mat_min[,dim(mat_min)[2]:1])
            rownames(mat)<-c(names(pl),names(min))
            mat<-mat[match(names(range),rownames(mat)),]
            mat
        }
        
        ORFs_gen<-ORFquant_results$ORFs_gen
        ORFs_spl_feat_maxORF<-ORFquant_results$ORFs_spl_feat_maxORF
        
        seqlevels(ORFs_gen)<-seqlevels(for_ORFquant$P_sites_all)
        seqlengths(ORFs_gen)<-seqlengths(for_ORFquant$P_sites_all)
        #nope, use coverage
        
        isovals<-cbind.data.frame(ORFs_tx$ORF_id_tr,ORFs_tx$ORF_pct_P_sites_pN)
        ORFs_gen$ORF_pct_P_sites_pN<-isovals[match(names(ORFs_gen),rownames(isovals)),2]
        ORFs_gen$ORF_pct_P_sites_pN[which(is.na(ORFs_gen$ORF_pct_P_sites_pN))]<-0
        ORFs_gen$ORF_pct_P_sites_pN<-round(ORFs_gen$ORF_pct_P_sites_pN,digits = 4)
        covisopl<-coverage(ORFs_gen[strand(ORFs_gen)=="+"],weight = ORFs_gen[strand(ORFs_gen)=="+"]$ORF_pct_P_sites_pN)
        covisomn<-coverage(ORFs_gen[strand(ORFs_gen)=="-"],weight = ORFs_gen[strand(ORFs_gen)=="-"]$ORF_pct_P_sites_pN)
        
        
        covpspl<-import(coverage_file_plus)
        strand(covpspl)<-"+"
        covpspl<-coverage(x = covpspl,weight = as.numeric(covpspl$score))
        covpsmn<-import(coverage_file_minus)
        strand(covpsmn)<-"-"
        covpsmn<-coverage(x = covpsmn,weight = as.numeric(covpsmn$score))
        
        
        #fix the previous annotation, then re-run
        listvals<-list()
        listiso<-list()
        listrang<-list()
        #CHECK DISTS
        types<-c("up_5ss","down_5ss","same_5ss","up_3ss","down_3ss","same_3ss","CDS_spanning")
        types_alt<-types[grep(types,pattern = "same",invert = T)]
        orfs_alt<-ORFs_spl_feat_maxORF[sapply(strsplit(ORFs_spl_feat_maxORF$spl_type,";"),function(x){sum(x%in%types_alt)>0})]
        ORFs_ok<-names(orfs_alt)
        genes_ok<-unique(ORFs_tx$gene_id[ORFs_tx$ORF_id_tr%in%ORFs_ok])
        ORFs_ok<-ORFs_tx$ORF_id_tr[ORFs_tx$gene_id%in%genes_ok]
        ORFs_ok<-ORFs_spl_feat_maxORF[names(ORFs_spl_feat_maxORF)%in%ORFs_ok]
        
        
        for(i in types){
            region<-ORFs_ok[grep(i,ORFs_ok$spl_type)]
            region<-region[!duplicated(GRanges(region))]
            if(i%in%c("up_3ss","up_lastCDS")){
                region<-region[elementNROWS(region$ref)==1]
                set.seed(666)
                if(length(region)>1000){region<-region[sample(x = 1:length(region),size = 1000,replace = F)]}
                region$dist<-end(region$ref)-end(region)
                region$dist[as.vector(strand(region)=="-")]<-start(region[as.vector(strand(region)=="-")])-start(region[as.vector(strand(region)=="-")]$ref)
                region$dist<-as.numeric(region$dist)
                dst<-region$dist
                nms<-names(region)
                region1<-promoters(resize(region,width = 1,fix = "end"),upstream = 25,downstream = 0)
                mcols(region1)<-mcols(region)
                
                region2<-promoters(resize(region,width = 1,fix = "end"),upstream = 0,downstream = 25)
                mcols(region2)<-mcols(region)
                
                region<-promoters(resize(region,width = 1,fix = "end"),upstream = 25,downstream = 25)
                mcols(region)<-mcols(region1)
                
                region$dist<-dst
                names(region)<-nms
            }
            
            if(i%in%c("down_3ss","down_lastCDS","same_3ss","same_lastCDS")){
                region<-region[elementNROWS(region$ref)==1]
                set.seed(666)
                if(length(region)>1000){region<-region[sample(x = 1:length(region),size = 1000,replace = F)]}
                region$dist<-end(region)-end(region$ref)
                region$dist[as.vector(strand(region)=="-")]<-start(region[as.vector(strand(region)=="-")]$ref)-start(region[as.vector(strand(region)=="-")])
                region$dist<-as.numeric(region$dist)
                dst<-region$dist
                nms<-names(region)
                region1<-promoters(resize(unlist(region$ref),width = 1,fix = "end"),upstream = 25,downstream = 0)
                mcols(region1)<-mcols(region)
                
                region2<-promoters(resize(unlist(region$ref),width = 1,fix = "end"),upstream = 0,downstream = 25)
                mcols(region2)<-mcols(region)
                
                region<-promoters(resize(unlist(region$ref),width = 1,fix = "end"),upstream = 25,downstream = 25)
                mcols(region)<-mcols(region2)
                
                region$dist<-dst
                names(region)<-nms
            }
            
            
            if(i%in%c("down_5ss","down_firstCDS")){
                region<-region[elementNROWS(region$ref)==1]
                set.seed(666)
                if(length(region)>1000){region<-region[sample(x = 1:length(region),size = 1000,replace = F)]}
                region$dist<-start(region)-start(region$ref)
                region$dist[as.vector(strand(region)=="-")]<-end(region[as.vector(strand(region)=="-")]$ref)-end(region[as.vector(strand(region)=="-")])
                region$dist<-as.numeric(region$dist)
                #invert numbers for donwstr-upstr
                dst<-region$dist
                nms<-names(region)
                region2<-promoters(resize(region,width = 1,fix = "start"),upstream = 25,downstream = 0)
                mcols(region2)<-mcols(region)
                
                region1<-promoters(resize(region,width = 1,fix = "start"),upstream = 0,downstream = 25)
                mcols(region1)<-mcols(region)
                
                region<-promoters(resize(region,width = 1,fix = "start"),upstream = 25,downstream = 25)
                mcols(region)<-mcols(region1)
                
                region$dist<-dst
                names(region)<-nms
            }
            
            if(i%in%c("up_5ss","up_firstCDS","same_5ss","same_firstCDS")){
                region<-region[elementNROWS(region$ref)==1]
                set.seed(666)
                if(length(region)>1000){region<-region[sample(x = 1:length(region),size = 1000,replace = F)]}
                region$dist<-start(region$ref)-start(region)
                region$dist[as.vector(strand(region)=="-")]<-end(region[as.vector(strand(region)=="-")])-end(region[as.vector(strand(region)=="-")]$ref)
                region$dist<-as.numeric(region$dist)
                dst<-region$dist
                nms<-names(region)
                region2<-promoters(resize(unlist(region$ref),width = 1,fix = "start"),upstream = 25,downstream = 0)
                mcols(region2)<-mcols(region)
                region1<-promoters(resize(unlist(region$ref),width = 1,fix = "start"),upstream = 0,downstream = 25)
                mcols(region1)<-mcols(region)
                region<-promoters(resize(unlist(region$ref),width = 1,fix = "start"),upstream = 25,downstream = 25)
                mcols(region)<-mcols(region1)
                region$dist<-dst
                names(region)<-nms
            }
            
            if(!i%in%c("CDS_spanning")){
                region1$iso<-0
                region1$iso[as.vector(strand(region1)=="+")]<-mean(covisopl[region1[as.vector(strand(region1)=="+")]])
                region1$iso[as.vector(strand(region1)=="-")]<-mean(covisomn[region1[as.vector(strand(region1)=="-")]])
                region2$iso<-0
                region2$iso[as.vector(strand(region2)=="+")]<-mean(covisopl[region2[as.vector(strand(region2)=="+")]])
                region2$iso[as.vector(strand(region2)=="-")]<-mean(covisomn[region2[as.vector(strand(region2)=="-")]])
                region$iso_1<-region1$iso
                region$iso_2<-region2$iso
                region$delta_iso<-region$iso_1-region$iso_2
                #check for those depending on region type
                #region<-region[region$iso_1<101 & region$iso_2<101]
                #region<-region[region$delta_iso<1]
                listrang[[i]]<-region
                listvals[[i]]<-aggregate_regions(range=region,coverage_plus=covpspl,coverage_minus=covpsmn)
                listiso[[i]]<-aggregate_regions(range=region,coverage_plus=covisopl,coverage_minus=covisomn)
                
                
            }
            
            if(i%in%c("CDS_spanning")){
                region<-region[elementNROWS(region$ref)==2]
                ress<-unlist(region$ref)
                start(ress[seq(1,length(region)*2,by = 2)])<-end(ress[seq(1,length(region)*2,by = 2)])
                
                region1_a<-promoters(ress[seq(1,length(region)*2,by = 2)],upstream = 25,downstream = 0)
                region2_a<-promoters(ress[seq(1,length(region)*2,by = 2)],upstream = 0,downstream = 25)
                region2_b<-promoters(ress[seq(2,length(region)*2,by = 2)],upstream = 25,downstream = 0)
                region1_b<-promoters(ress[seq(2,length(region)*2,by = 2)],upstream = 0,downstream = 25)
                region_a<-promoters(ress[seq(1,length(region)*2,by = 2)],upstream = 25,downstream = 25)
                names(region_a)<-names(region)
                region_b<-promoters(ress[seq(2,length(region)*2,by = 2)],upstream = 25,downstream = 25)
                names(region_b)<-names(region)
                region_a$dist<-sapply(region$ref,function(x){min(width(gaps(x)))})
                region_b$dist<-sapply(region$ref,function(x){min(width(gaps(x)))})
                
                #check for upstr and dowstr, or maybe main and secondary
                
                region1_a$iso<-0
                region1_a$iso[as.vector(strand(region1_a)=="+")]<-mean(covisopl[region1_a[as.vector(strand(region1_a)=="+")]])
                region1_a$iso[as.vector(strand(region1_a)=="-")]<-mean(covisomn[region1_a[as.vector(strand(region1_a)=="-")]])
                region1_b$iso<-0
                region1_b$iso[as.vector(strand(region1_b)=="+")]<-mean(covisopl[region1_b[as.vector(strand(region1_b)=="+")]])
                region1_b$iso[as.vector(strand(region1_b)=="-")]<-mean(covisomn[region1_b[as.vector(strand(region1_b)=="-")]])
                region2_a$iso<-0
                region2_a$iso[as.vector(strand(region2_a)=="+")]<-mean(covisopl[region2_a[as.vector(strand(region2_a)=="+")]])
                region2_a$iso[as.vector(strand(region2_a)=="-")]<-mean(covisomn[region2_a[as.vector(strand(region2_a)=="-")]])
                region2_b$iso<-0
                region2_b$iso[as.vector(strand(region2_b)=="+")]<-mean(covisopl[region2_b[as.vector(strand(region2_b)=="+")]])
                region2_b$iso[as.vector(strand(region2_b)=="-")]<-mean(covisomn[region2_b[as.vector(strand(region2_b)=="-")]])
                
                
                
                region_a$iso_1<-region1_a$iso
                region_a$iso_2<-region2_a$iso
                region_b$iso_1<-region1_b$iso
                region_b$iso_2<-region2_b$iso
                
                region_a$delta_iso<-region_a$iso_1-region_a$iso_2
                region_b$delta_iso<-region_b$iso_1-region_b$iso_2
                names(region_a)<-paste(names(region),"left",sep = "_")
                names(region_b)<-paste(names(region),"right",sep = "_")
                region<-c(region_a,region_b)
                
                listrang[[i]]<-region
                listvals[[i]]<-aggregate_regions(range=region,coverage_plus=covpspl,coverage_minus=covpsmn)
                listiso[[i]]<-aggregate_regions(range=region,coverage_plus=covisopl,coverage_minus=covisomn)
                
            }
        }
        
        same_all<-ORFs_ok[grep("same_5ss;same_3ss",ORFs_ok$spl_type)]
        same_all<-same_all[width(same_all)>110]
        same_all$dist<-50
        same_all$iso_1<-50
        same_all$iso_2<-50
        same_all$delta_iso<-50
        set.seed(666)
        if(length(same_all)>1000){same_all<-same_all[sample(x = 1:length(same_all),size = 1000,replace = F)]}
        same_all_5<-promoters(resize(same_all,width = 1,fix = "center"),upstream = 50,downstream = 0)
        same_all_3<-promoters(resize(same_all,width = 1,fix = "center"),upstream = 0,downstream = 50)
        
        listrang[["sameall_5"]]<-same_all_5
        listvals[["sameall_5"]]<-aggregate_regions(range=same_all_5,coverage_plus=covpspl,coverage_minus=covpsmn)
        listiso[["sameall_5"]]<-aggregate_regions(range=same_all_5,coverage_plus=covisopl,coverage_minus=covisomn)
        
        listrang[["sameall_3"]]<-same_all_3
        listvals[["sameall_3"]]<-aggregate_regions(range=same_all_3,coverage_plus=covpspl,coverage_minus=covpsmn)
        listiso[["sameall_3"]]<-aggregate_regions(range=same_all_3,coverage_plus=covisopl,coverage_minus=covisomn)
        
        
        topl_ok<-c("up_5ss","down_3ss","down_5ss","up_3ss","CDS_spanning")
        topl_cons<-c("same_5ss","same_3ss","sameall_5","sameall_3","CDS_spanning")
        
        plots_ribo<-list()
        plots_iso<-list()
        plots_sketch<-list()
        
        for(i in 1:length(topl_ok)){
            rans<-listrang[[topl_ok[i]]]
            vals<-listvals[[topl_ok[i]]]
            vals_contr<-listvals[[topl_cons[i]]]
            isos_contr<-listiso[[topl_cons[i]]]
            rans_contr<-listrang[[topl_cons[i]]]
            if(topl_ok[i]=="CDS_spanning"){
                vals_contr<-rbind(listvals[["same_3ss"]],listvals[["same_5ss"]])
                isos_contr<-rbind(listiso[["same_3ss"]],listiso[["same_5ss"]])
                rans_contr<-c(listrang[["same_3ss"]],listrang[["same_5ss"]])
                names(rans_contr)<-paste(names(rans_contr),c(rep("left",length(listrang[["same_3ss"]])),rep("right",length(listrang[["same_3ss"]]))),sep = "_")
                
            }
            
            isos<-listiso[[topl_ok[i]]]
            rans$id<-paste(names(rans),GRanges(rans),sep = "_")
            #ok<-intersect(intersect(which(abs(rans$dist)>25),which(rans$iso_1<100 & rans$iso_2<100)),which(rans$delta_iso<0))
            ok<-intersect(which(abs(rans$dist)>25),which(rans$iso_1<101 & rans$iso_2<101))
            ok<-intersect(ok,which(rans$delta_iso>0))
            ok<-intersect(ok,which(rowSums(vals)>25))
            
            rans_contr<-rans_contr[rowSums(vals_contr)>25]
            isos_contr<-isos_contr[rowSums(vals_contr)>25,]
            vals_contr<-vals_contr[rowSums(vals_contr)>25,]
            vals<-vals[ok,]
            isos<-isos[ok,]
            rans<-rans[ok]
            rownames(vals)<-rans$id
            rownames(isos)<-rans$id
            rans$delta_iso<-round(as.numeric(abs(rans$delta_iso)),digits=2)
            rans$delta_iso2<-abs(rans$iso_1-rans$delta_iso)
            if(topl_ok[i]%in%c("up_3ss","down_5ss","down_firstCDS","up_lastCDS")){
                rans$delta_iso2<-abs(rans$iso_1-rans$iso_2)
            }
            rans$delta_iso<-rans$delta_iso2
            qnt<-unique(quantile(rans$delta_iso,probs=seq(0,1,length.out = 4),include.lowest = T))
            
            rans$group<-as.character(cut(rans$delta_iso,breaks = qnt,include.lowest = T))
            
            
            
            tb<-aggregate(rans$delta_iso,by=list(rans$group),mean)
            rans$group<-round(tb$x[match(rans$group,tb$Group.1)],digits = 2)
            ok<-which(!is.na(rans$group))
            vals<-vals[ok,]
            isos<-isos[ok,]
            rans<-rans[ok]
            rans_contr$id<-paste(names(rans_contr),GRanges(rans_contr),sep = "_")
            rans_contr$group<-"No_mixture"
            rownames(isos_contr)<-rans_contr$id
            rownames(vals_contr)<-rans_contr$id
            levss<-sort(unique(rans$group),decreasing = F)
            levss<-c("No_mixture",levss)
            if(topl_ok[i]=="CDS_spanning"){
                rans$spanning<-"downstream"
                rans$spanning[grep(names(rans),pattern = "left")]<-"upstream"
                rans_contr$spanning<-"downstream"
                rans_contr$spanning[grep(names(rans_contr),pattern = "left")]<-"upstream"
            }
            notcol<-names(mcols(rans))[!names(mcols(rans))%in%names(mcols(rans_contr))]
            rans$ref<-NULL
            mcols(rans_contr)[,notcol]<-""
            mcols(rans_contr)<-mcols(rans_contr)[,names(mcols(rans))]
            rans<-c(rans,rans_contr)
            rans$group<-factor(rans$group,levels=levss)
            
            
            #HERE IT WAS THE PROBLEM
            #levss<-levels(rans$group)[order(as.numeric(sapply(strsplit(unique(rans$group),"-"),"[[",1)),decreasing = F)]
            #levels(rans$group)<-levss
            
            
            #vals<-t(apply(vals,1,scale))
            #vals<-vals/apply(vals,1,max)
            
            vals<-rbind(vals,vals_contr)
            isos<-rbind(isos,isos_contr)
            
            vals<-vals/rowSums(vals)
            #vals<-scale(vals)
            colnames(vals)<--25:24
            vls<-DataFrame(melt(vals))
            vls$delta_iso<-rans$delta_iso[match(vls$Var1,rans$id)]
            vls$group<-rans$group[match(vls$Var1,rans$id)]
            #isos<-t(apply(isos,1,scale))
            #isos<-isos/apply(isos,1,max)
            
            isos<-isos/rowSums(isos)
            #isos<-scale(isos)
            colnames(isos)<--25:24
            
            vls2<-DataFrame(melt(isos))
            vls2$delta_iso<-rans$delta_iso[match(vls2$Var1,rans$id)]
            
            vls2$group<-rans$group[match(vls2$Var1,rans$id)]
            vls$type<-"Ribo-seq_coverage"
            vls2$type<-"Iso_values"
            if(topl_ok[i]=="CDS_spanning"){
                vls2$spanning<-rans$spanning[match(vls2$Var1,rans$id)]
                vls$spanning<-rans$spanning[match(vls$Var1,rans$id)]
            }
            
            vlsall<-rbind(vls,vls2)
            vlsall[,"Position_from_splice_site"]<-vlsall$Var2
            vlsall$type<-factor(vlsall$type)
            
            
            dat<-as.data.frame(vlsall)
            cols<- colorRampPalette(c("dark blue","blue","red"))(8)[c(1,3,6,7)]
            cols<-cols[(4-length(qnt)+1):4]
            
            means<-aggregate(dat$value,by=list(dat$group,dat$Position_from_splice_site,dat$type),mean)
            colnames(means)<-c("group","Position_from_splice_site","type","value")
            
            if(topl_ok[i]=="CDS_spanning"){
                means<-aggregate(dat$value,by=list(dat$group,dat$Position_from_splice_site,dat$type,dat$spanning),mean)
                colnames(means)<-c("group","Position_from_splice_site","type","spanning","value")
                ok<-ok[1:(length(ok)/2)]
                means$spanning<-factor(means$spanning,levels=c("upstream","downstream"))
            }
            
            
            means$groups_all<-paste(means$group,means$type)
            means$shape<-means$type
            
            means$color<-cols[as.numeric(means$group)]
            means_ribo<-means[means$type!="Iso_values",]
            means_Iso<-means[means$type=="Iso_values",]
            #ss<-aes(aes(ymax = value + sd, ymin = value - sd,group = groups_all,color=groups_all,shape=type,linetype=type))
            nm_titl<-topl_ok[i]
            
            a<-ggplot(means_ribo, aes(x = Position_from_splice_site,y =  value, group = group,color=group)) +
                geom_line(size=1.8) +
                geom_vline(xintercept = 0,col="black",lty=2) +
                theme_classic() +
                #geom_ribbon(aes(ymin=value-sds, ymax=value+sds, x = Position_from_splice_site, fill = group), alpha = 0.01) +
                ylab("Normalized coverage \nRibo-seq") +
                theme(axis.title.x = element_text(size=20),axis.text.x  = element_text(angle=45, vjust=0.5, size=15)) +
                theme(axis.title.y = element_text(size=20),axis.text.y  = element_text(angle=45, vjust=0.5, size=15))+
                scale_color_manual(values=cols,name="Avg.\ncontribution\nadditional ORF(s)") +
                ggtitle(paste("Event =",nm_titl,"\nn =",length(ok)))
            
            b<-ggplot(means_Iso, aes(x = Position_from_splice_site,y =  value, group = group,color=group)) +
                geom_line(size=1.8) +
                geom_vline(xintercept = 0,col="black",lty=2) +
                theme_classic() +
                #geom_ribbon(aes(ymin=value-sds, ymax=value+sds, x = Position_from_splice_site, fill = group), alpha = 0.01) +
                theme(axis.title.x = element_text(size=20),axis.text.x  = element_text(angle=45, vjust=0.5, size=15)) +
                theme(axis.title.y = element_text(size=20),axis.text.y  = element_text(angle=45, vjust=0.5, size=15))+
                scale_color_manual(values=cols,name="Avg.\ncontribution\nadditional ORF(s)") +
                ggtitle(paste("Event =",nm_titl,"\nn =",length(ok)))
            if(topl_ok[i]=="CDS_spanning"){
                a<-a+facet_grid(facets = .~spanning) + xlab("Position from Splice Sites") +  ylab("Normalized coverage\nRibo-seq")
                b<-b+facet_grid(facets = .~spanning) + xlab("Position from Splice Sites") +  ylab("Normalized coverage\nORFquant ORFs")
            }
            if(topl_ok[i]%in%c("up_5ss")){
                a<-a + xlab("Position from Splice Site") +  ylab("Normalized coverage\nRibo-seq")
                b<-b + xlab("Position from Splice Site") +  ylab("Normalized coverage\nORFquant ORFs")
            }
            if(topl_ok[i]%in%c("up_firstCDS")){
                a<-a + xlab("Position from Start codon") +  ylab("Normalized coverage\nRibo-seq")
                b<-b + xlab("Position from Start codon") +  ylab("Normalized coverage\nORFquant ORFs")
            }
            
            if(topl_ok[i]%in%c("up_lastCDS")){
                a<-a + xlab("Position from Stop codon") +  ylab("")
                b<-b + xlab("Position from Stop codon") +  ylab("")
            }
            
            
            if(!topl_ok[i]%in%c("up_lastCDS","up_firstCDS","up_5ss","CDS_spanning")){
                a<-a + xlab("") +  ylab("")
                b<-b + xlab("") +  ylab("")
            }
            
            
            if(topl_ok[i]=="down_5ss"){
                blocks=data.frame(x1=c(.5,0,0,.14,.28,.42), x2=c(1,1,.09,.23,.37,.5), y1=c(.5,.05,.675,.675,.675,.675), y2=c(.85,.4,.685,.685,.685,.685),coll=c("red","blue","red","red","red","red"),stringsAsFactors = F)
            }
            
            if(topl_ok[i]=="down_firstCDS"){
                blocks=data.frame(x1=c(.5,0,0), x2=c(1,1,1), y1=c(.5,.05,.575), y2=c(.85,.4,.785),coll=c("red","blue","red"),stringsAsFactors = F)
            }
            
            if(topl_ok[i]=="up_5ss"){
                blocks=data.frame(x1=c(0,.5,0,.14,.28,.42), x2=c(1,1,.09,.23,.37,.5), y1=c(.5,.05,.225,.225,.225,.225), y2=c(.85,.4,.235,.235,.235,.235),coll=c("red","blue","blue","blue","blue","blue"),stringsAsFactors = F)
            }
            
            if(topl_ok[i]=="up_firstCDS"){
                blocks=data.frame(x1=c(0,.5,0), x2=c(1,1,1), y1=c(.5,.05,.125), y2=c(.85,.4,.330),coll=c("red","blue","blue"),stringsAsFactors = F)
            }
            
            if(topl_ok[i]=="down_3ss"){
                blocks=data.frame(x1=c(0,0,.5,.64,.78,.92), x2=c(1,.5,.59,.73,.87,1), y1=c(.5,.05,.225,.225,.225,.225), y2=c(.85,.4,.235,.235,.235,.235),coll=c("red","blue","blue","blue","blue","blue"),stringsAsFactors = F)
            }
            
            if(topl_ok[i]=="down_lastCDS"){
                blocks=data.frame(x1=c(0,0,.5), x2=c(1,.5,1), y1=c(.5,.05,.125), y2=c(.85,.4,.330),coll=c("red","blue","blue"),stringsAsFactors = F)
            }
            
            
            if(topl_ok[i]=="up_3ss"){
                blocks=data.frame(x1=c(0,0,.5,.64,.78,.92), x2=c(.5,1,.59,.73,.87,1), y1=c(.5,.05,.675,.675,.675,.675), y2=c(.85,.4,.685,.685,.685,.685),coll=c("red","blue","red","red","red","red"),stringsAsFactors = F)
            }
            
            if(topl_ok[i]=="up_lastCDS"){
                blocks=data.frame(x1=c(0,0,.5), x2=c(.5,1,1), y1=c(.5,.05,.575), y2=c(.85,.4,.785),coll=c("red","blue","red"),stringsAsFactors = F)
            }
            
            
            if(topl_ok[i]=="CDS_spanning"){
                blocks=data.frame(x1=c(0,0,0.25,0.39,0.53,0.67,.75), x2=c(1,.25,0.34,0.48,0.62,0.75,1), y1=c(.5,.05,.225,.225,.225,.225,.05), y2=c(.85,.4,.235,.235,.235,.235,.4),coll=c("red","blue","blue","blue","blue","blue","blue"),stringsAsFactors = F)
            }
            
            
            colss<-blocks$coll
            blocks$coll<-factor(blocks$coll,levels=c("red","blue","white"))
            df<-blocks
            bl<-ggplot() + geom_rect(data=df, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2,fill=coll,colour=coll))+
                scale_fill_manual(values = colss) +
                scale_color_manual(values = colss) + 
                theme_nothing() + ylim(c(0,1))
            
            plots_ribo[[i]]<-a
            plots_iso[[i]]<-b
            plots_sketch[[i]]<-bl
            
        }
        
        x1s<-c(0,.25,.5,.75,.3)
        x1s_iso<-x1s
        x1s_sk<-x1s
        
        wids<-c(rep(.25,4),.5)
        wids_iso<-wids
        wids_sk<-wids
        
        y2s<-c(rep(.8,4),.3)
        y2s_iso<-y2s-.2
        y2s_sk<-round(as.numeric(y2s_iso-.1),digits = 3)
        
        heis<-c(rep(.2,4),.2)
        heis_iso<-heis
        heis_sk<-rep(.1,10)
        
        all<-ggdraw()
        for(i in 1:length(plots_ribo)){
            all<-all+draw_plot(plots_ribo[[i]], x = x1s[i], y = y2s[i], width = wids[i], height = heis[i])  +
                draw_plot(plots_iso[[i]],  x = x1s_iso[i], y = y2s_iso[i], width = wids_iso[i], height = heis_iso[i])  +
                draw_plot(plots_sketch[[i]],  x = x1s_sk[i], y = y2s_sk[i], width = wids_sk[i], height = heis_sk[i])
            
        }
        
        all_metapl<-all
        list_ORFquant_plots[["Meta_splicing_coverage"]]<-all_metapl
        list_ORFquant_plots[["Meta_splicing_coverage"]][["pars"]]<-c(22,16)
    }
    
    save(list_ORFquant_plots,file = paste0(output_plots_path,"/",prefix,"_","ORFquant_plots_RData"))
    for(i in names(list_ORFquant_plots)){
        pdf(file = paste0(output_plots_path,"/",prefix,"_",i,".pdf"),width = list_ORFquant_plots[[i]][["pars"]][1] ,height = list_ORFquant_plots[[i]][["pars"]][2])
        print(list_ORFquant_plots[[i]])
        dev.off()
    }
    cat(paste("Plotting ORFquant results for ",ORFquant_output_file,"  --- Done! ",date(),"\n",sep = ""))
    
}



#' Create an html report summarizing ORFquant results
#'
#' This function creates an html report showing summary statistics for ORFquant-detected ORFs.
#' 
#' @keywords ORFquant
#' @author Lorenzo Calviello, \email{calviello.bio@@gmail.com}
#' 
#' @param input_files Character vector with full paths to plot files (*ORFquant_plots_RData) 
#' generated with \code{plot_ORFquant_results}. 
#' Must be of same length as \code{input_sample_names}.
#' 
#' @param input_sample_names Character vector containing input names.
#' Must be of same length as \code{input_files}.
#' 
#' @param output_file String; full path to html report file.
#' @return The function saves the html report file with the file path \code{output_file}.
#' @details This function creates the html report visualizing final ORFquant results. \cr \cr
#' Input are two lists of the same length: \cr \cr
#' a) \code{input_files}: list of full paths to one or multiple input files 
#' (*ORFquant_plots_RData files generated with \code{plot_ORFquant_results}) and \cr \cr
#' b) \code{input_sample_names}: list of corresponding names describing the file content (these are used as names in the report). \cr \cr
#' For the report, a RMarkdown file is rendered as html document, saved as \code{output_file}. \cr \cr
#' @seealso \code{\link{plot_ORFquant_results}}, \code{\link{run_ORFquant}}
#' @export

create_ORFquant_html_report <- function(input_files, input_sample_names, output_file){
    
    # get input and output file paths
    input_files <- paste(normalizePath(dirname(input_files)),basename(input_files),sep="/")
    output_file <- paste(normalizePath(dirname(output_file)),basename(output_file),sep="/")
    
    # get path to RMarkdown file (to be rendered)
    rmd_path <- paste(system.file(package="ORFquant"),"/rmd/ORFquant_template.Rmd",sep="")
    
    sink(file = paste(output_file,"_ORFquant_report_output.txt",sep = ""))
    # render RMarkdown file > html report
    suppressWarnings(render(rmd_path, 
                            params = list(input_files = input_files,
                                          input_sample_names = input_sample_names),
                            output_file = output_file))
    sink()
}


#' Create a plot of the ORFquant results at a locus
#'
#' Create a plot of the ORFquant results at a locus. Uses the info form the orfquant results about where the psite data is located.
#' 
#' @keywords ORFquant
#' @author Lorenzo Calviello, \email{calviello.bio@@gmail.com}
#' 
#' @param input_files Character vector with full paths to plot files (*ORFquant_plots_RData) 
#' generated with \code{plot_ORFquant_results}. 
#' Must be of same length as \code{input_sample_names}.
#' 
#' @param input_sample_names Character vector containing input names.
#' Must be of same length as \code{input_files}.
#' 
#' @param locus String; a gene name, must be present in names(orfquant_results$ORFs_gen)
#' @param orfquant_results A list containing processed output from ORFquant
#' @param bam_files bam files, (or pre-processed bam data from RiboseQC) to be plotted
#' @param plotfile the file into which the plot will be saved as a pdf
#' @return returns the value of plotfile if successfull.
#' @export
plot_orfquant_locus<-function(locus,orfquant_results, bam_files=NULL,plotfile='locusplot.pdf', col ='green' ){
    
    if(!is.null(orfquant_results$psite_data_file)){stop("this object looks like it's form an old version of ORFquant, it doesn't list the psite data file")}
    if(length(orfquant_results$psite_data_file)>1){stop("locus plots aren't supported for multiple psite tracks - either unify the psite tracks or modify the input object to have only one track")}
    riboseqcoutput<-get(load(orfquant_results$psite_data_file))
    selgene <- locus
    stopifnot(length(locus)==1)
    stopifnot(locus %in% names(orfquant_results$ORFs_gen))
    orfs_quantified_gen <-  orfquant_results$ORFs_gen%>%subset(.,str_detect(names(.),selgene))
    orfquantgrscores = orfs_quantified_gen$TrP_pNpM[match(names(orfquantgr),orfs_quantified_gen$ORF_id_tr)]
    mcols(orfs_quantified_gen) <- mcols(orfquant_results$ORFs_tx)[match(names(orfs_quantified_gen),orfquant_results$ORFs_tx$ORF_id_tr),]
    seltxs <- orfs_quantified_gen$transcript_id%>%unique
    orfs_quantified_gen$feature <- 'CDS'
    orfs_quantified_gen$transcript=orfs_quantified_gen$transcript_id
    orfs_quantified_tr <- anno$exons_tx
    orfs_quantified_tr <- orfs_quantified_tr[unique(orfs_quantified_gen$transcript_id)]
    #get the orfs, then get the negative coverage
    #add transcript info to the ORFs_tx object
    seqinf <- Seqinfo(names(anno$exons_tx),anno$exons_tx%>%width%>%sum)
    #Now get the negatives for each ORF
    utrs <- orfquant_results$ORFs_tx%>%subset(gene_id==selgene)%>%keepSeqlevels(seltxs)%>%{seqinfo(.)<-seqinf[seltxs];.}%>%
        coverage%>%
        as('GRanges')%>%subset(score==0)%>%mapFromTranscripts(anno$exons_tx)%>%
        {.$transcript <- names(anno$exons_tx)[.$transcriptsHits];.}%>%
        {.$feature='utr';.}

    orfquantgr <- c(
        orfs_quantified_gen[,c('feature','transcript')]
        ,utrs[,c('feature','transcript')]
    )
    orfquantgr$feature[orfquantgr$feature=='CDS'] <- names(orfquantgr)[orfquantgr$feature=='CDS']
    #get correct col name
    if('ORF_pM' %in% metacols) quantcol = 'ORF_pM' else quantcol = 'TrP_pNpM'
    orfscores<- orfquantgr$feature%>%unique%>%setNames(match(.,orfs_quantified_gen$ORF_id_tr)%>%orfs_quantified_gen[[quantcol]][.],.)
    orfcols <- orfscores%>%{./max(na.omit(.))}%>%
        c(0,.)%>%
        map_chr(~possibly(rgb,'white')(0,.,0))%>%setNames(c('0',orfquantgr$feature%>%unique))
    ###Define non selected
    disctxs<-anno$txs_gene[selgene]%>%unlist%>%.$tx_name%>%unique%>%setdiff(seltxs)
    disc_orfquantgr <- anno$cds_txs_coords%>%
        keepSeqlevels(disctxs,'coarse')%>%
        {seqinfo(.)<-seqinf[disctxs];.}%>%  coverage%>%
        as('GRanges')%>%
        subset(.$score==0)%>%
        {   
            txgr = .
            out = mapFromTranscripts(txgr,anno$exons_tx)
            out$score = txgr$score[out$xHits]
            out
        }%>%
        {.$transcript <- names(anno$exons_tx)[.$transcriptsHits];.}%>%
        {.$feature=ifelse(.$score==0,'utr','CDS');.}
    disc_orfquantgr <- disc_orfquantgr%>% c(.,anno$cds_txs[disctxs]%>%unlist%>%{.$feature=rep('CDS',length(.));.$transcript=names(.);.})
    discORFnames<-paste0(disctxs,'_',start(anno$cds_txs_coords[disctxs]),'_',end(anno$cds_txs_coords[disctxs]))%>%setNames(disctxs)
    disc_orfquantgr$symbol = discORFnames[disc_orfquantgr$transcript]
    fakejreads <- riboseqcoutput$junctions%>%subset(any(gene_id==selgene))%>%resize(width(.)+2,'center')%>%
        {.$cigar <- paste0('1M',width(.)-2,'N','1M');.}
    fakejreads <- fakejreads[map2(seq_along(fakejreads[]),fakejreads$reads,rep)%>%unlist]
    ncols <- 2
    nrows <- 1
    orfcols <- orfcols[order(-orfscores[names(orfcols)])]
    orfquantgr_sorted <- orfquantgr[order(orfscores[names(orfquantgr)])]
    fix_utrs <- function(orfquantgr_sorted){
        orfquantgr_sorted$symbol = names(orfquantgr_sorted)
        #add utrs for each ORF
        #for each selected ORF
        orftrpairs<-orfquantgr_sorted%>%subset(feature!='utr')%>%mcols%>%as.data.frame%>%distinct
        orfutrs <- orfquantgr_sorted%>%subset(feature=='utr')
        orfutrs<-lapply(1:nrow(orftrpairs),function(i){
            orfutrs <- orfutrs%>%subset(transcript==orftrpairs$transcript[i])
            orfutrs$symbol = orftrpairs$feature[i]
            names(orfutrs) = orfutrs$symbol 
            orfutrs
        })%>%GRangesList%>%unlist
        orfquantgr_sorted<-orfquantgr_sorted%>%subset(feature!='utr')%>%c(.,orfutrs)
        orfquantgr_sorted
    }
    orfquantgr_sorted<-fix_utrs(orfquantgr_sorted)
    # disc_orfquantgrfix<-fix_utrs(disc_orfquantgr)
    disc_orfquantgrfix<-(disc_orfquantgr)
    #dimenions, extent of the plot
    plotstart = start(selgenerange) - (0.2 * (end(selgenerange)-start(selgenerange)))
    plotend = end(selgenerange) + (0 * (end(selgenerange)-start(selgenerange)))
    library(Gviz)
    legendwidth=1/10
    plottitle <- paste0('ORFquant: ',selgene)
    #write to pdf
    pdf(plotfile,width=14+2,h=7)
    #code for arranging legend next to the locus plot
    grid.newpage()
    vp1 <- viewport(x = 0, y = 0, w = 1-legendwidth*1.5, h = 1,
    just = c("left", "bottom"), name = "vp1")
    vp2 <- viewport(x = 1-legendwidth*1.5, y = 0, w = legendwidth*1.5, h = 1,
    just = c("left", "bottom"))
    vp3 <- viewport(x = 1-legendwidth*1.5, y = 0, w = legendwidth*1.5, h = 1/7,
    just = c("left", "bottom"))
    pushViewport(vp1)
    #finally plot the locus
    plotTracks(main=plottitle,cex.main=2,legend=TRUE,add=TRUE,
    from=plotstart,to=plotend,#zoomed in on the orf in question
    sizes=c(1,1,1,1,1,1,1),rot.title=0,cex.title=1,title.width=2.5,
    c(
        GenomeAxisTrack(range=selgenerange),
        # rnaseqtrack, # plot the riboseq signal
        # txs_discarded_Track,
        # txs_selected_track,
        # DataTrack(riboseqcoutput$P_sites_all%>%subsetByOverlaps(selgenerange),type='hist'),
        GeneRegionTrack(name='discarded\ntranscripts',anno$exons_tx[disctxs]%>%unlist%>%{.$transcript=names(.);.$feature=rep('exon',length(.));.},fill='#F7CAC9',
                transcriptAnnotation='transcript'),
        GeneRegionTrack(exon='forestgreen',name='selected\ntranscripts',anno$exons_tx[seltxs]%>%unlist%>%{.$transcript=names(.);.$feature=rep('exon',length(.));.},fill='#F7CAC9',
                transcriptAnnotation='transcript'),
        DataTrack(legend=TRUE,name='\t\t P-Sites',col.histogram='forestgreen',riboseqcoutput$P_sites_all%>%subsetByOverlaps(selgenerange),type='hist'),
        AlignmentsTrack(name='\nJunction Reads\n\n\n',col.sashimi='forestgreen',fakejreads[,'cigar'],type='sashimi',sashimiNumbers=TRUE),
        # GeneRegionTrack(discarded_orfs_gen),
        GeneRegionTrack(name='Discarded\nORFs',disc_orfquantgrfix,
            transcriptAnnotation='symbol',collapse=FALSE,thinBoxFeature='utr',CDS='blue',utr='white'),
        GeneRegionTrack(name='Selected\nORFs',
                range=orfquantgr_sorted,collapse=FALSE,
            thinBoxFeature='utr',CDS='red',utr='white',
            transcriptAnnotation='symbol'
        )%>%
        # identity
        {displayPars(.)[names(orfcols)]<-orfcols;.}
    ),
    col.labels='black',
    chr=seqnames(selgenerange)
    )
    #create barchart of intensities
    popViewport(1)
    pushViewport(vp2)
    cols = I(c(orfcols[which.min(orfscores)],orfcols[which.max(orfscores)]))
    grid.draw(g_legend(qplot(x=1:2,y=1:2,color=range(orfscores,na.rm=T))+
    scale_color_gradient(name='Normalized ORF Expr\n(ORFs_pM)',
        breaks = setNames(sort(na.omit(orfscores)),floor(na.omit(sort(orfscores)))%>% format(big.mark=",",scientific=FALSE) ),
        low=cols[1],high=cols[2])+theme(text=element_text(size=14),legend.key.size=unit(.5,'inches'))
    ))
    popViewport(1)
    pushViewport(vp3)
    dev.off()
    normalizePath(plotfile)
    #return file name
    return(plotfile)
} 
}