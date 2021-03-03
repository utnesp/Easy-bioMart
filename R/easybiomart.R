#' @title get cDNA sequence
#' @param input input
#' @export
#' @import biomaRt
#' @import jsonlite
#' @import httr
#' @import xml2
#' @examples getSequence_cDNA()
getSequence_cDNA <- function(input, which.mart = mart, which.type = "external_gene_name", which.seq.type = "cdna", fasta.output = "") {
    df.temp <- biomaRt::getSequence(id = input, type = which.type, seqType = which.seq.type, mart = which.mart)
    return(df.temp)
    print(paste("Retrieved", nrow(df.temp), "sequence(s) from getSequence. Input had", length(input), "entries."))

    if (fasta.output != "") {
        print("Moving any existing file to trash folder using trash command")
        system(paste("trash", fasta.output))
        exportFASTA(df.temp, fasta.output)
        fasta.temp <- gsub(".fasta", ".temp.fasta", fasta.output)
        print(fasta.out)
        print(fasta.output)
        command <- capture.output(cat("awk '", '/^>/{print s? s"\\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}', "' ", fasta.output, " > ", fasta.temp, sep = ""))
        system(command)
        system(paste("mv",  fasta.temp, fasta.output))
        system(paste("trash", fasta.temp))
    }
}

#' @title get peptide sequence
#' @param input input
#' @export
#' @import biomaRt
#' @examples getSequence_peptide()
getSequence_peptide <- function(input, which.mart = mart, which.type = "external_gene_name", which.seq.type = "peptide", fasta.output = "") {
    df.temp <- biomaRt::getSequence(id = input, type = which.type, seqType = which.seq.type, mart = which.mart)
    return(df.temp)
    print(paste("Retrieved", nrow(df.temp), "sequence(s) from getSequence. Input had", length(input), "entries."))

    if (fasta.output != "") {
        print("Moving any existing file to trash folder using trash command")
        system(paste("trash", fasta.output))
        exportFASTA(df.temp, fasta.output)
        fasta.temp <- gsub(".fasta", ".temp.fasta", fasta.output)
        print(fasta.out)
        print(fasta.output)
        command <- capture.output(cat("awk '", '/^>/{print s? s"\\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}', "' ", fasta.output, " > ", fasta.temp, sep = ""))
        system(command)
        system(paste("mv",  fasta.temp, fasta.output))
        system(paste("trash", fasta.temp))
    }
}

#' @title get 3' UTR sequence
#' @param ensg input
#' @export
#' @import biomaRt
#' @examples getSequence_3UTR()
getSequence_3UTR <- function(ensg, which.mart = mart, which.type = "ensembl_gene_id", which.seq.type = "3utr", fasta.output = "") {
    df.temp <- biomaRt::getSequence(id = ensg, type = which.type, seqType = which.seq.type, mart = which.mart)
    return(df.temp)
    print(paste("Retrieved", nrow(df.temp), "sequence(s) from getSequence. Input had", length(input), "entries."))

    if (fasta.output != "") {
        print("Moving any existing file to trash folder using trash command")
        system(paste("trash", fasta.output))
        exportFASTA(df.temp, fasta.output)
        fasta.temp <- gsub(".fasta", ".temp.fasta", fasta.output)
        print(fasta.out)
        print(fasta.output)
        command <- capture.output(cat("awk '", '/^>/{print s? s"\\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}', "' ", fasta.output, " > ", fasta.temp, sep = ""))
        system(command)
        system(paste("mv",  fasta.temp, fasta.output))
        system(paste("trash", fasta.temp))
    }
}

#' @title get 5' UTR sequence
#' @param ensg input
#' @export
#' @import biomaRt
#' @examples getSequence_5UTR()
getSequence_5UTR <- function(ensg, which.mart = mart, which.type = "ensembl_gene_id", which.seq.type = "5utr", fasta.output = "") {
    df.temp <- biomaRt::getSequence(id = ensg, type = which.type, seqType = which.seq.type, mart = which.mart)
    return(df.temp)
    print(paste("Retrieved", nrow(df.temp), "sequence(s) from getSequence. Input had", length(input), "entries."))

    if (fasta.output != "") {
        print("Moving any existing file to trash folder using trash command")
        system(paste("trash", fasta.output))
        exportFASTA(df.temp, fasta.output)
        fasta.temp <- gsub(".fasta", ".temp.fasta", fasta.output)
        print(fasta.out)
        print(fasta.output)
        command <- capture.output(cat("awk '", '/^>/{print s? s"\\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}', "' ", fasta.output, " > ", fasta.temp, sep = ""))
        system(command)
        system(paste("mv",  fasta.temp, fasta.output))
        system(paste("trash", fasta.temp))
    }
}

#' @title get upstream sequence
#' @param input input
#' @export
#' @import biomaRt
#' @examples getUpstream()
getUpstream <- function(input, bp_upstream = 2000, which.mart = mart, which.type = "external_gene_name", which.seq.type = "gene_flank", fasta.output = "") {
    df.temp <- biomaRt::getSequence(id = input, type = which.type, seqType = which.seq.type, upstream = bp_upstream, mart = which.mart)
    print(paste("Retrieved", nrow(df.temp), "sequence(s) from getSequence. Input had", length(input), "entries."))
    df.temp <- data.frame(df.temp)
    print(paste(length(grep("N", df.temp$gene_flank)), "entries contain one or more N's"))
    df.temp2 <- gsub('N', '', df.temp$gene_flank)
    df.temp3 <- data.frame(cbind(df.temp2, df.temp$external_gene_name, deparse.level = F))
    df.temp3 <- df.temp3[nchar(as.character(df.temp3$X1)) > (bp_upstream * 0.8), ]
    print(paste("Left with", nrow(df.temp3), "when removing sequences shorter than", (bp_upstream * 0.8), "basepairs"))
    df.temp4 <- df.temp3[!duplicated(df.temp3$X2), ]
    print(paste("Removed", nrow(df.temp3) - nrow(df.temp4), "duplicate entries"))
    df.temp5 <- df.temp4[df.temp4$X2 %in% input, ]
    # df.temp3 <- merge(df.temp3, df, by.x = "X2", by.y = "external_gene_name", all.y = T)
    print(paste("Final fasta contains", nrow(df.temp5), "entries compared to input:", length(input)))
    # note before exporting: exportFASTA appends to existing files

    if (fasta.output != "") {
        print("Moving any existing file to trash folder using trash command")
        system(paste("trash", fasta.output))
        exportFASTA(df.temp5, fasta.output)
        fasta.temp <- gsub(".fasta", ".temp.fasta", fasta.output)
        print(fasta.out)
        print(fasta.output)
        command <- capture.output(cat("awk '", '/^>/{print s? s"\\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}', "' ", fasta.output, " > ", fasta.temp, sep = ""))
        system(command)
        system(paste("mv",  fasta.temp, fasta.output))
        system(paste("trash", fasta.temp))
    }

    return(df.temp5)
    print("Note: exportFASTA appends to any existing files!")
    print("To delete newline, see awk one-liner in the function code")
    # # fix newline: awk '/^>/{print s? s"\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}'
}

#' @title Convert hugo names to ensembl gene identifiers
#' @param hgnc input
#' @export
#' @import biomaRt
#' @examples hugo2ensg()
hugo2ensg <- function(hgnc, biomart = mart, combine = F, df2 = "", by.x = "hgnc_symbol", by.y = "hgnc_symbol", all = F){
    # if(is.object(hgnc)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F & class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame" ) hgnc <- row.names(hgnc)
    
    df <- res <- getBM(
        filters= "hgnc_symbol",
        attributes= c("hgnc_symbol", "ensembl_gene_id"),
        values= hgnc,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}

#' @title Biomart conversion
#' @param ext_name input
#' @export
#' @import biomaRt
#' @examples ext_name2ensg()
ext_name2ensg <- function(ext_name, biomart = mart, combine = F, df2 = "", by.x = "external_gene_name", by.y = "external_gene_name", all = F){
    # if(is.object(ext_name)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F &  class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame"  ) ext_name <- row.names(ext_name)

    df <- getBM(
        filters= "external_gene_name",
        attributes= c("external_gene_name", "ensembl_gene_id"),
        values= ext_name,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}

#' @title Biomart conversion
#' @param ext_name input
#' @export
#' @import biomaRt
#' @examples ext_name2enst()
ext_name2enst <- function(ext_name, biomart = mart, combine = F, df2 = "", by.x = "external_gene_name", by.y = "external_gene_name", all = F){
    # if(is.object(ext_name)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F &  class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame"  ) ext_name <- row.names(ext_name)

    df <- getBM(
        filters= "external_gene_name",
        attributes= c("external_gene_name", "ensembl_transcript_id"),
        values= ext_name,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}

#' @title Biomart conversion
#' @param ext_name input
#' @export
#' @import biomaRt
#' @examples ext_name2chr_band()
ext_name2chr_band <- function(ext_name, biomart = mart, combine = F, df2 = "", by.x = "external_gene_name", by.y = "external_gene_name", all = F){
    # if(is.object(ext_name)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F &  class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame"  ) ext_name <- row.names(ext_name)

    df <- getBM(
        filters= "external_gene_name",
        attributes= c("external_gene_name", "chromosome_name", "band"),
        values= ext_name,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}

#' @title Biomart conversion
#' @param ext_name input
#' @export
#' @import biomaRt
#' @examples ext_name2chr_band()
ext_name2chr_TSS_strand <- function(ext_name, biomart = mart, combine = F, df2 = "", by.x = "external_gene_name", by.y = "external_gene_name", all = F){

    df <- getBM(
        filters= "external_gene_name",
        attributes= c("external_gene_name", "chromosome_name", "transcription_start_site","strand"),
        values= ext_name,
        mart = biomart)


    
    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}

#' @title Biomart conversion
#' @param ext_name input
#' @export
#' @import biomaRt
#' @examples ext_name2chr_band()
ext_name2chr_TSS_geneStartEnd_strand <- function(ext_name, biomart = mart, combine = F, df2 = "", by.x = "external_gene_name", by.y = "external_gene_name", all = F){

    df <- getBM(
        filters= "external_gene_name",
        attributes= c("external_gene_name", "chromosome_name", "transcription_start_site", "cdna_coding_start", "cdna_coding_end", "strand"),
        values= ext_name,
        mart = biomart)


    
    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}


#' @title Biomart conversion
#' @param df input
#' @export
#' @import biomaRt
#' @examples ext_name2chr_band_df()
ext_name2chr_band_df <- function(df, biomart = mart, combine = F, df2 = "", by.x = "external_gene_name", by.y = "external_gene_name", all = F){
    df.band <- getBM(
        filters= "external_gene_name",
        attributes= c("external_gene_name", "chromosome_name", "band"),
        values= df$external_gene_name,
        mart = biomart)

    df.band$band <- paste(df.band$chromosome_name, df.band$band, sep="")
    df.band_2 <- df.band[duplicated(df.band$band), ]
    df.band_sub <- df.band[df.band$band %in% df.band_2$band, ]
    df.band_sub3 <- df.band_sub[order(df.band_sub$band),]
    df.band_sub3 <- merge(df.band_sub3, df, by = "external_gene_name")
    df.band_sub3 <- df.band_sub3[order(df.band_sub3$band),]
    df.band_sub3$chromosome_name <- NULL

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}

#' @title Biomart conversion
#' @param goid input
#' @export
#' @import biomaRt
#' @examples goid2goname()
goid2goname <- function(goid, biomart = mart, combine = F, df2 = "", by.x = "go_id", by.y = "go_id", all = F){
    # if(is.object(goid)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F &  class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame"  ) goid <- row.names(goid)

    df <- getBM(
        filters= "go",
        attributes= c("go_id", "name_1006"),
        values= goid,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}

#' @title Biomart conversion
#' @param goid input
#' @export
#' @import biomaRt
#' @examples goid2ensg_extName_biotype()
goid2ensg_extName_biotype <- function(goid, biomart = mart, combine = F, df2 = "", by.x = "go_id", by.y = "go_id", all = F){
    # if(is.object(goid)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F &  class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame"  ) goid <- row.names(goid)

    df <- getBM(
        filters= "go",
        attributes= c("go_id", "ensembl_gene_id", "external_gene_name", "gene_biotype"),
        values= goid,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}

#' @title Biomart conversion
#' @param goid input
#' @export
#' @import biomaRt
#' @examples goid2ensg_extName_biotype_name_evidence()
goid2ensg_extName_biotype_name_evidence <- function(goid, biomart = mart, combine = F, df2 = "", by.x = "go_id", by.y = "go_id", all = F){
    # if(is.object(goid)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F &  class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame"  ) goid <- row.names(goid)

    df <- getBM(filters= "go",
                attributes= c("go_id", "ensembl_gene_id", "external_gene_name", "gene_biotype", "go_linkage_type"),
                values= goid,
                mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}

#' @title Biomart conversion
#' @param ensg input
#' @export
#' @import biomaRt
#' @examples enst2ensg()
enst2ensg <- function(ensg, biomart = mart, combine = F, df2 = "", by.x = "ensembl_transcript_id", by.y = "ensembl_transcript_id", all = F){
    # # if(is.object(ensg)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F &  class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame"  ) ensg <- row.names(ensg)

    df <- getBM(
        filters= "ensembl_transcript_id",
        attributes= c("ensembl_gene_id", "ensembl_transcript_id"),
        values= ensg,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}

#' @title Biomart conversion
#' @param ensg input
#' @export
#' @import biomaRt
#' @examples enst2ensg()
enst2utr <- function(ensg, biomart = mart, combine = F, df2 = "", by.x = "ensembl_transcript_id", by.y = "ensembl_transcript_id", all = F){
    # # if(is.object(ensg)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F &  class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame"  ) ensg <- row.names(ensg)

    df <- getBM(
        filters= "ensembl_transcript_id",
        attributes= c("ensembl_transcript_id", "3_utr_start", "3_utr_end", "5_utr_start", "5_utr_end"),
        values= ensg,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}

#' @title Biomart conversion
#' @param enst input
#' @export
#' @import biomaRt
#' @examples enst2ensg()
enst2ensg.extName <- function(enst = NULL, biomart = mart, combine = F, df2 = "", by.x = "ensembl_transcript_id", by.y = "ensembl_transcript_id", all = F){
    # if(is.object(enst)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F &  class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame"  ) ensg <- row.names(ensg)

    if (is.null(enst)) {
          df <- getBM(
            attributes= c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"),
            mart = biomart)
    } else {
        df <- getBM(
            filters= "ensembl_transcript_id",
            attributes= c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"),
            values= enst,
            mart = biomart)
    }
    
    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}

#' @title Biomart conversion
#' @param enst input
#' @export
#' @import biomaRt
#' @examples enst2ensg()
enst2ensg.extName.with.version <- function(enst = NULL, biomart = mart, combine = F, df2 = "", by.x = "ensembl_transcript_id", by.y = "ensembl_transcript_id", all = F){
    # if(is.object(enst)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F &  class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame"  ) ensg <- row.names(ensg)

    if (is.null(enst)) {
          df <- getBM(
            attributes= c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name", "version", "transcript_version"),
            mart = biomart)
    } else {
        df <- getBM(
            filters= "ensembl_transcript_id",
            attributes= c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name", "version", "transcript_version"),
            values= enst,
            mart = biomart)
    }
    
    df$ensembl_transcript_id <- paste(df$ensembl_transcript_id, df$transcript_version, sep = ".")
    df$ensembl_gene_id <- paste(df$ensembl_gene_id, df$version, sep = ".")
    df$version <- NULL; df$transcript_version <- NULL

    
    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}

#' @title Biomart conversion
#' @param ensg input
#' @export
#' @import biomaRt
#' @examples ensg2hugo()
ensg2hugo <- function(ensg, biomart = mart, combine = F, df2 = "", by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all = F){
    # # if(is.object(ensg)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F &  class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame"  ) ensg <- row.names(ensg)

    df <- getBM(
        filters= "ensembl_gene_id",
        attributes= c("ensembl_gene_id", "hgnc_symbol"),
        values= ensg,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }

            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}

#' @title Biomart conversion
#' @param ensg input
#' @export
#' @import biomaRt
#' @examples ensg2hugo()
ext_name2wiggleplotr <- function(ext_name, biomart = mart, combine = F, df2 = "", by.x = "gene_name", by.y = "external_gene_name", all = F){
    # # if(is.object(ensg)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F &  class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame"  ) ensg <- row.names(ensg)

    df <- getBM(
        filters= "external_gene_name",
        attributes= c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name", "strand"),
        values= ext_name,
        mart = biomart)

    colnames(df) <- c("transcript_id", "gene_id", "gene_name", "strand")
    
    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }

            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}

#' @title Biomart conversion
#' @param ensg input
#' @export
#' @import biomaRt
#' @examples ensg2source_status()
ensg2source_status <- function(ensg, biomart = mart, combine = F, df2 = "", by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all = T) {
    # # if(is.object(ensg)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F &  class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame"  ) ensg <- row.names(ensg)

    df <- getBM(
        filters= "ensembl_gene_id",
        attributes= c("ensembl_gene_id", "source", "status"),
        values= ensg,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}

#' @title Biomart conversion
#' @param ensg input
#' @export
#' @import biomaRt
#' @examples ensg2reactome()
ensg2reactome <- function(ensg, biomart = mart, combine = F, df2 = "", by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all = F){
    # # if(is.object(ensg)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F & class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame" ) ensg <- row.names(ensg)

    df <- getBM(
        filters= "ensembl_gene_id",
        attributes= c("ensembl_gene_id", "external_gene_name", "reactome", "reactome_gene"),
        values= ensg,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}

#' @title Biomart conversion
#' @param ensg input
#' @export
#' @import biomaRt
#' @examples ensg2chr_start_end()
ensg2chr_start_end <- function(ensg, biomart = mart, combine = F, df2 = "", by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all = F, as.IGV = F){
    # # if(is.object(ensg)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F & class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame" ) ensg <- row.names(ensg)

    df <- getBM(
        filters= "ensembl_gene_id",
        attributes= c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "band"),
        values= ensg,
        mart = biomart)

    if (as.IGV == T) {
        df$chromosome_name <- paste(df$chromosome_name, ":", df$start_position, "-", df$end_position, sep = "")
        colnames(df)[colnames(df) == "chromosome_name"] <- "position"
        df$start_position <- NULL
        df$end_position <- NULL
    }

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}

#' @title Biomart conversion
#' @param ensg input
#' @export
#' @import biomaRt
#' @examples ext_name2previous_ext_name()
ext_name2previous_ext_name <- function(ext_name, hg38.to.hg19 = T, from.mart = NULL, to.mart = NULL, combine = F, df2 = "", by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all = F, remove.bracket.info = T) {
    # # # if(is.object(ensg)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F & class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame" ) ensg <- row.names(ensg)

    if(hg38.to.hg19 == T) {
        if(is.null(from.mart)) from.mart = useMart("ENSEMBL_MART_ENSEMBL", dataset='hsapiens_gene_ensembl')
        if(is.null(to.mart)) to.mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")   
    } else {
        if(is.null(to.mart)) to.mart = useMart("ENSEMBL_MART_ENSEMBL", dataset='hsapiens_gene_ensembl')
        if(is.null(from.mart)) from.mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")   
    }
    
    gene_id <- getBM(
        filters= "external_gene_name",
        attributes= c("external_gene_name", "ensembl_gene_id"),
        values= ext_name,
        mart = from.mart)
    if(hg38.to.hg19 == T) colnames(gene_id)[1] <- "external_gene_name_hg38"
    if(hg38.to.hg19 != T) colnames(gene_id)[1] <- "external_gene_name_hg19"
    
    df <- getBM(
        filters= "ensembl_gene_id",
        attributes= c("ensembl_gene_id", "external_gene_name"),
        values= gene_id$ensembl_gene_id,
        mart = to.mart)
    if(hg38.to.hg19 == T) colnames(df)[2] <- "external_gene_name_hg19"
    if(hg38.to.hg19 != T) colnames(df)[2] <- "external_gene_name_hg38"
    

    df <- merge(df, gene_id, by = "ensembl_gene_id", all = T)
    
    
    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}


#' @title Biomart conversion
#' @param ensg input
#' @export
#' @import biomaRt
#' @examples ensg2description()
ensg2description <- function(ensg, biomart = mart, combine = F, df2 = "", by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all = F, remove.bracket.info = T) {
    # # # if(is.object(ensg)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F & class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame" ) ensg <- row.names(ensg)

    df <- getBM(
        filters= "ensembl_gene_id",
        attributes= c("ensembl_gene_id", "description"),
        values= ensg,
        mart = biomart)

    if(remove.bracket.info == T) df$description <- gsub(" \\[.*","", df$description)
    
    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}

#' @title Biomart conversion
#' @param ensg input
#' @export
#' @import biomaRt
#' @examples ensg2uniprot()
ensg2uniprot <- function(ensg, biomart = mart, combine = F, df2 = "", by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all = F) {
   # # if(is.object(ensg)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F & class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame" ) ensg <- row.names(ensg)

   df <- getBM(
        filters= "ensembl_gene_id",
        attributes= c("ensembl_gene_id", "uniprotsptrembl", "uniprotswissprot"),
        values= ensg,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}

#' @title Biomart conversion
#' @param uniprotswissprot input
#' @export
#' @import biomaRt
#' @examples uniprot2ensg()
uniprot2ensg <- function(uniprotswissprot, useSwissprot = T, biomart = mart, combine = F, df2 = "", by.x = "uniprotswissprot", by.y = "uniprotswissprot", all = F) {
    # if(is.object(uniprotswissprot)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F & class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame" ) uniprotswissprot <- row.names(uniprotswissprot)

    if (useSwissprot == T) {
        print("Using Swiss-Prot, a manually annotated and reviewed database")
        df <- getBM(
            filters= "uniprotswissprot",
            attributes= c("ensembl_gene_id", "uniprotswissprot", "uniprot_gn", "external_gene_name"),
            values= uniprotswissprot,
            mart = biomart)

        if (combine == T) {
            if (df2 == "") {
                if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
                return(df.anno)
            } else {
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                return(df.anno)
            }
        } else {
            return(df)
        }


    } else if (useSwissprot == F) {
        by.x = "uniprotsptrembl"
        print("Using TrEMBL, an automatically annotated and not reviewed database")
        df <- getBM(
            filters= "uniprotsptrembl",
            attributes= c("ensembl_gene_id", "uniprotsptrembl", "uniprot_gn", "external_gene_name"),
            values= uniprotswissprot,
            mart = biomart)

        if (combine == T) {
            if (df2 == "") {
                if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
                return(df.anno)
            } else {
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                return(df.anno)
            }
        } else {
            return(df)
        }

    } else if (useSwissprot == "both") {

        print("Using both Swiss-Prot and TrEMBL")

        df <- getBM(
            filters= "uniprotswissprot",
            attributes= c("ensembl_gene_id", "uniprotswissprot", "uniprot_gn", "external_gene_name"),
            values= uniprotswissprot,
            mart = biomart)

        colnames(df) <- c("ensembl_gene_id", "uniprotswissport_sptrembl", "uniprot_gn", "external_gene_name")

        df3 <- getBM(
            filters= "uniprotsptrembl",
            attributes= c("ensembl_gene_id", "uniprotsptrembl", "uniprot_gn", "external_gene_name"),
            values= uniprotswissprot,
            mart = biomart)

        colnames(df3) <- c("ensembl_gene_id", "uniprotswissport_sptrembl", "uniprot_gn", "external_gene_name")

        df <- rbind(df, df3)

        by.x = "uniprotswissport_sptrembl"

        if (combine == T) {
            if (df2 == "") {
                if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
                return(df.anno)
            } else {
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                return(df.anno)
            }
        } else {
            return(df)
        }
    }
}

#' @title Biomart conversion
#' @param probe_id input
#' @export
#' @import biomaRt
#' @examples agilent44k2ensg_chr()
agilent44k2ensg_chr <- function(probe_id, biomart = mart, combine = F, df2 = "", by.x = "efg_agilent_wholegenome_4x44k_v1", by.y = "efg_agilent_wholegenome_4x44k_v1", all = F) {
    # if(is.object(probe_id)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F & class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame" ) probe_id <- row.names(probe_id)

    df <- getBM(
        filters= "efg_agilent_wholegenome_4x44k_v1",
        attributes= c("ensembl_gene_id", "external_gene_id", ),
        values= probe_id,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}

#' @title Biomart conversion
#' @param ensg input
#' @export
#' @import biomaRt
#' @examples ensg2GO()
ensg2GO <- function(ensg, biomart = mart, combine = F, df2 = "", by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all = F){
    # # if(is.object(ensg)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F & class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame" ) ensg <- row.names(ensg)

    df <- getBM(
        filters= "ensembl_gene_id",
        attributes= c("ensembl_gene_id", "hgnc_symbol", "go_id", "name_1006", "namespace_1003"),
        values= ensg,
        mart = biomart)

    if(is.na(df$hgnc_symbol == TRUE)) {
        df <- getBM(
            filters= "ensembl_gene_id",
            attributes= c("ensembl_gene_id", "external_gene_name",  "go_id", "name_1006", "namespace_1003"),
            values= ensg,
            mart = biomart)
    }

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}

#' @title Biomart conversion
#' @param ext_name input
#' @export
#' @import biomaRt
#' @examples hugo2GO()
hugo2GO <- function(ext_name, biomart = mart, combine = F, df2 = "", by.x = "external_gene_name", by.y = "external_gene_name", all = F){
    # if(is.object(ext_name)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F & class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame" ) ext_name <- row.names(ext_name)

    df <- getBM(
        filters= "external_gene_name",
        attributes= c("external_gene_name", "name_1006", "namespace_1003"),
        values= ext_name,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}

#' @title Biomart conversion
#' @param df input
#' @export
#' @import biomaRt
#' @examples ensg2entrez_fc()
ensg2entrez_fc <- function(df, ensg_col = row.names(df), fc_col = "log2FC", biomart = mart,
                           combine = F, df2 = "", by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all = F) {
    fc_col <- match(fc_col, colnames(MEG3))
    if (ensg_col != row.names(df)) {
        ensg_col <- match(ensg_col, colnames(MEG3))
        t <- data.frame(ensembl_gene_id = df[ensg_col], fc = df[fc_col])
    } else {
        t <- data.frame(ensembl_gene_id = row.names(df), fc = df[fc_col])
    }
    z <- ensg2entrez(t$ensembl_gene_id, biomart)
    x <- merge(t, z, by = "ensembl_gene_id")
    row.names(x) <- x$entrezgene_id
    x$entrezgene_id <- NULL
    x$ensembl_gene_id <- NULL
    names(x) <- ""
    y <- t(x[1])
    return(y)
}

#' @title Biomart conversion
#' @param ense_id input
#' @export
#' @import biomaRt
#' @examples exon2ensg.is_constitutive()
exon2ensg.is_constitutive <- function(ense_id, biomart = mart, combine = F, df2 = "", by.x = "ensembl_exon_id", by.y = "ensembl_exon_id", all = F){
    # if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F & class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame" ) ense_id <- row.names(ense_id)

    df <- getBM(
        filters= "ensembl_exon_id",
        attributes= c("ensembl_exon_id", "ensembl_gene_id", "is_constitutive"),
        values= ense_id,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}

#' @title Biomart conversion
#' @param ense_id input
#' @export
#' @import biomaRt
#' @examples exon2ensg.is_constitutive.rank()
exon2ensg.is_constitutive.rank <- function(ense_id, biomart = mart, combine = F, df2 = "", by.x = "ensembl_exon_id", by.y = "ensembl_exon_id", all = F){
    # if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F & class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame" ) ense_id <- row.names(ense_id)

    df <- getBM(
        filters= "ensembl_exon_id",
        attributes= c("ensembl_exon_id", "ensembl_gene_id", "is_constitutive", "rank"),
        values= ense_id,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}



#' @title Biomart conversion
#' @param ense_id input
#' @export
#' @import biomaRt
#' @examples exon2gencode.basic()
exon2gencode.basic <- function(ense_id, biomart = mart, combine = F, df2 = "", by.x = "ensembl_exon_id", by.y = "ensembl_exon_id", all = F){
    # if(is.object(ense_id)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F & class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame" ) ense_id <- row.names(ense_id)

    df <- getBM(
        filters= "ensembl_exon_id",
        attributes= c("ensembl_exon_id", "ensembl_gene_id", "transcript_gencode_basic"),
        values= ense_id,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }

            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}

#' @title Biomart conversion
#' @param ense_id input
#' @export
#' @import biomaRt
#' @examples exon2appris()
exon2appris <- function(ense_id, biomart = mart, combine = F, df2 = "", by.x = "ensembl_exon_id", by.y = "ensembl_exon_id", all = F){
   # if(is.object(ense_id)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F & class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame" ) ense_id <- row.names(ense_id)

   df <- getBM(
        filters= "ensembl_exon_id",
        attributes= c("ensembl_exon_id", "ensembl_gene_id", "transcript_appris"),
        values= ense_id,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}

#' @title Biomart conversion
#' @param ense_id input
#' @export
#' @import biomaRt
#' @examples exon2ensg()
exon2ensg <- function(ense_id, biomart = mart, combine = F, df2 = "", by.x = "ensembl_exon_id", by.y = "ensembl_exon_id", all = F){
    # if(is.object(ense_id)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F & class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame" ) ense_id <- row.names(ense_id)
    df <- getBM(
        filters= "ensembl_exon_id",
        attributes= c("ensembl_exon_id", "ensembl_gene_id"),
        values= ense_id,
        mart = biomart)


    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}

#' @title Biomart conversion
#' @param ensg input
#' @export
#' @import biomaRt
#' @examples ensg2gene_biotype()
ensg2gene_biotype <- function(ensg, biomart = mart, combine = F, df2 = "", by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all = F){
    # # if(is.object(ensg)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F & class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame" ) ensg <- row.names(ensg)
    df <- getBM(
        filters= "ensembl_gene_id",
        attributes= c("ensembl_gene_id", "gene_biotype"),
        values= ensg,
        mart = biomart)


    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}

#' @title Biomart conversion
#' @param ensg input
#' @export
#' @import biomaRt
#' @examples ensg2ext_name()
ensg2ext_name <- function(ensg, biomart = mart, combine = F, df2 = "", by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all = F){
    # # if(is.object(ensg)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F & class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame" ) ensg <- row.names(ensg)

    df <- getBM(
        filters= "ensembl_gene_id",
        attributes= c("ensembl_gene_id", "external_gene_name"),
        values= ensg,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }

            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}

#' @title Biomart conversion
#' @param ensg input
#' @export
#' @import biomaRt
#' @examples ensg2ext_name_biotype()
ensg2ext_name_biotype <- function(ensg, biomart = mart, combine = F, df2 = "", by.x = "ensembl_gene_id", by.y = "row.names", all = F, title = "", ...) {
    # # if(is.object(ensg)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F & class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame" ) ensg <- row.names(ensg)

    df <- getBM(
        filters= "ensembl_gene_id",
        attributes= c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
        values= ensg,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {
                ## by.y is called from function, and by default is row.names
                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}

#' @title Biomart conversion
#' @param ensg input
#' @export
#' @import biomaRt
#' @examples ensg2ext_name_biotype()
ensg2hugo_biotype <- function(ensg, biomart = mart, combine = F, df2 = "", by.x = "ensembl_gene_id", by.y = "row.names", all = F, title = "", ...) {
    # # if(is.object(ensg)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F & class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame" ) ensg <- row.names(ensg)

    df <- getBM(
        filters= "ensembl_gene_id",
        attributes= c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
        values= ensg,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {
                ## by.y is called from function, and by default is row.names
                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}

#' @title Biomart conversion
#' @param ensg input
#' @export
#' @import biomaRt
#' @examples ensg2ext_name_biotype()
ensg2ext_name_biotype_gencode <- function(ensg, biomart = mart, combine = F, df2 = "", by.x = "ensembl_gene_id", by.y = "row.names", all = F, title = "", ...) {
    # # if(is.object(ensg)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F & class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame" ) ensg <- row.names(ensg)

    df <- getBM(
        filters= "ensembl_gene_id",
        attributes= c("ensembl_gene_id", "external_gene_name", "gene_biotype", "transcript_gencode_basic"),
        values= ensg,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {
                ## by.y is called from function, and by default is row.names
                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}

#' @title Biomart conversion
#' @param ensg input
#' @export
#' @import biomaRt
#' @examples enst2ext_name_biotype()
enst2ext_name_biotype <- function(ensg, biomart = mart, combine = F, df2 = "", by.x = "ensembl_transcript_id", by.y = "ensembl_transcript_id", all = F){
    # # if(is.object(ensg)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F & class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame" ) ensg <- row.names(ensg)

    ensg <- gsub("\\..*", "", ensg) # remove version numbers if they exist

    
    df <- getBM(
        filters= "ensembl_transcript_id",
        attributes= c("ensembl_transcript_id", "transcript_biotype", "external_transcript_name", "external_gene_name", "gene_biotype"),
        values= ensg,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])
    
                df2[, which(colnames(df2) == by.y) ] <- ensg

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}

#' @title Biomart conversion
#' @param ensg input
#' @export
#' @import biomaRt
#' @examples enst2ext_name_biotype()
enst2ensg2_ext_name_biotype <- function(ensg, biomart = mart, combine = F, df2 = "", by.x = "ensembl_transcript_id", by.y = "ensembl_transcript_id", all = F){
    # # if(is.object(ensg)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F & class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame" ) ensg <- row.names(ensg)

    df <- getBM(
        filters= "ensembl_transcript_id",
        attributes= c("ensembl_transcript_id", "ensembl_gene_id", "transcript_biotype", "external_transcript_name", "external_gene_name", "gene_biotype"),
        values= ensg,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}
    
#' @title Biomart conversion
#' @param entrez_ids input
#' @export
#' @import biomaRt
#' @examples entrez2ext_name()
entrez2ext_name <- function(entrez_ids, biomart = mart, combine = F, df2 = "", by.x = "entrezgene_id", by.y = "entrezgene_id", all = F)
{
    # if(is.object(entrez_ids)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F & class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame" ) entrez_ids <- row.names(entrez_ids)
    df <- getBM(
        filters= "entrezgene_id",
        attributes= c("entrezgene_id", "external_gene_name"),
        values= entrez_ids,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }

            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}

#' @title Biomart conversion
#' @param entrez_ids input
#' @export
#' @import biomaRt
#' @examples entrez2ensg()
entrez2ensg <- function(entrez_ids, biomart = mart, combine = F, df2 = "", by.x = "entrezgene_id", by.y = "entrezgene_id", all = F)
{
    # if(is.object(entrez_ids)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F & class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame" ) entrez_ids <- row.names(entrez_ids)
    df <- getBM(
        filters= "entrezgene_id",
        attributes= c("entrezgene_id", "ensembl_gene_id"),
        values= entrez_ids,
        mart = biomart)


    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}

#' @title Biomart conversion
#' @param hgnc input
#' @export
#' @import biomaRt
#' @examples hugo2entrez()
hugo2entrez <- function(hgnc, biomart = mart, combine = F, df2 = "", by.x = "hgnc_symbol", by.y = "hgnc_symbol", all = F){
    # if(is.object(hgnc)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F & class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame" ) hgnc <- row.names(hgnc)
    df <- getBM(
        filters= "hgnc_symbol",
        attributes= c("hgnc_symbol", "entrezgene_id"),
        values= hgnc,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}

#' @title Biomart conversion
#' @param ENSG_v82 input
#' @export
#' @import biomaRt
#' @examples ensg2entrez()
ensg2entrez <- function(ENSG_v82, biomart = mart, combine = F, df2 = "", by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all = F){
   #  # if(is.object(ENSG_v82)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F & class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame" ) ENSG_v82 <- row.names(ENSG_v82)

    df <- getBM(
        filters= "ensembl_gene_id",
        attributes= c("ensembl_gene_id", "entrezgene_id"),
        values= ENSG_v82,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }

}

#' @title Biomart conversion
#' @param ENSG_v82 input
#' @export
#' @import biomaRt
#' @examples ensg2entrez()
enst2entrez <- function(ENSG_v82, biomart = mart, combine = F, df2 = "", by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all = F){
    # if(is.object(ENSG_v82)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F & class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame" ) ENSG_v82 <- row.names(ENSG_v82)

    df <- getBM(
        filters= "ensembl_transcript_id",
        attributes= c("ensembl_transcript_id", "entrezgene_id"),
        values= ENSG_v82,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }

}

#' @title Biomart conversion
#' @param external_gene_name input
#' @export
#' @import biomaRt
#' @examples ext_name2entrez()
ext_name2reactome <- function(external_gene_name, biomart = mart, combine = F, df2 = "", by.x = "external_gene_name", by.y = "external_gene_name", all = F){
    # if(is.object(external_gene_name)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F & class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame" ) external_gene_name <- row.names(external_gene_name)

    df <- getBM(
        filters= "external_gene_name",
        attributes= c("reactome", "external_gene_name"),
        values= external_gene_name,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }

}



#' @title Biomart conversion
#' @param external_gene_name input
#' @export
#' @import biomaRt
#' @examples ext_name2entrez()
ext_name2source <- function(external_gene_name, biomart = mart, combine = F, df2 = "", by.x = "external_gene_name", by.y = "external_gene_name", all = F){

    df <- getBM(
        filters= "external_gene_name",
        attributes= c("external_gene_name", "external_gene_source"),
        values= external_gene_name,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }

}

#' @title Biomart conversion
#' @param external_gene_name input
#' @export
#' @import biomaRt
#' @examples ext_name2entrez()
ext_name2entrez <- function(external_gene_name, biomart = mart, combine = F, df2 = "", by.x = "external_gene_name", by.y = "external_gene_name", all = F){
    # if(is.object(external_gene_name)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F & class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame" ) external_gene_name <- row.names(external_gene_name)

    df <- getBM(
        filters= "external_gene_name",
        attributes= c("external_gene_name", "entrezgene_id"),
        values= external_gene_name,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }

}

#' @title Biomart conversion
#' @param entrez_ids input
#' @export
#' @import biomaRt
entrez2ext_name <- function(entrez_ids, biomart = mart, combine = F, df2 = "", by.x = "entrezgene_id", by.y = "entrezgene_id", all = F){
    # if(is.object(entrez_ids)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F & class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame" ) entrez_ids <- row.names(entrez_ids)

    df <- getBM(
        filters= "entrezgene_id",
        attributes= c("entrezgene_id", "external_gene_name"),
        values= entrez_ids,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }

}

#' @title Biomart conversion
#' @param refseq_ids input
#' @export
#' @import biomaRt
refseq2ext_name <- function(refseq_ids, biomart = mart, combine = F, df2 = "", by.x = "refseq", by.y = "refseq", all = F){
    # if(is.object(entrez_ids)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F & class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame" ) entrez_ids <- row.names(entrez_ids)

    df <- getBM(
        filters= "refseq_mrna",
        attributes= c("refseq_mrna", "external_gene_name"),
        values= refseq_ids,
        mart = biomart)
    
    df2 <- getBM(
        filters= "refseq_ncrna",
        attributes= c("refseq_ncrna", "external_gene_name"),
        values= refseq_ids,
        mart = biomart)

    colnames(df) <- c("refseq", "external_gene_name")
    colnames(df2) <- colnames(df)
    
    df <- rbind(df, df2)
    
    if (combine == T) {
        
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

                return(df.anno)
            } else {
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                
                return(df.anno)
        }
    } else {
        return(df)
    }

}

ensg2refseq_ext_name <- function(ensg_ids, biomart = mart, combine = F, df2 = "", by.x = "refseq", by.y = "refseq", all = F) {
    # if(is.object(entrez_ids)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F & class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame" ) entrez_ids <- row.names(entrez_ids)
    
    df <- getBM(
        filters= "ensembl_gene_id",
        attributes= c("ensembl_gene_id", "external_gene_name", "refseq_mrna", "refseq_ncrna"),
        values= ensg_ids,
        mart = biomart)
    
      df2 <- getBM(
        filters= "ensembl_gene_id",
        attributes= c("ensembl_gene_id", "external_gene_name", "refseq_mrna_predicted", "refseq_ncrna_predicted"),
        values= ensg_ids,
        mart = biomart)
      
      
        df3 <- getBM(
        filters= "ensembl_gene_id",
        attributes= c("ensembl_gene_id", "external_gene_name", "refseq_peptide", "refseq_peptide_predicted"),
        values= ensg_ids,
        mart = biomart)
    
    
    colnames(df2) <- colnames(df)
    colnames(df3) <- colnames(df)
    df$type <- "canonical"
    df2$type <- "predicted"
    df3$type <- "peptide"
    df <- rbind(df,df2,df3)
        
    if (combine == T) {
        
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

                return(df.anno)
            } else {
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                
                return(df.anno)
        }
    } else {
        return(df)
    }

}

#' @title Biomart conversion
#' @param refseq_ids input
#' @export
#' @import biomaRt
refseq2ensg_ext_name <- function(refseq_ids, biomart = mart, combine = F, df2 = "", by.x = "refseq", by.y = "refseq", all = F){
    # if(is.object(entrez_ids)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F & class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame" ) entrez_ids <- row.names(entrez_ids)
    
    cat("Searchin refseq_mrna \n")
    df <- getBM(
        filters= "refseq_mrna",
        attributes= c("refseq_mrna", "ensembl_gene_id", "external_gene_name"),
        values= refseq_ids,
        mart = biomart)
    
    cat("Searchin refseq_ncrna \n")
    df2 <- getBM(
        filters= "refseq_ncrna",
        attributes= c("refseq_ncrna", "ensembl_gene_id", "external_gene_name"),
        values= refseq_ids,
        mart = biomart)
    
    cat("Searchin refseq_mrna_predicted \n")
    df3 <- getBM(
        filters= "refseq_mrna_predicted",
        attributes= c("refseq_ncrna", "ensembl_gene_id", "external_gene_name"),
        values= refseq_ids,
        mart = biomart)
    
    cat("Searchin refseq_ncrna_predicted \n")
    df4 <- getBM(
        filters= "refseq_ncrna_predicted",
        attributes= c("refseq_ncrna", "ensembl_gene_id", "external_gene_name"),
        values= refseq_ids,
        mart = biomart)
    
    cat("Searchin refseq_peptide \n")
    df5 <- getBM(
        filters= "refseq_peptide",
        attributes= c("refseq_ncrna", "ensembl_gene_id", "external_gene_name"),
        values= refseq_ids,
        mart = biomart)
    
    cat("Searchin refseq_peptide_predicted \n")
    df6 <- getBM(
        filters= "refseq_peptide_predicted",
        attributes= c("refseq_ncrna", "ensembl_gene_id", "external_gene_name"),
        values= refseq_ids,
        mart = biomart)

    colnames(df) <- c("refseq", "ensembl_gene_id", "external_gene_name")
    colnames(df2) <- colnames(df)
    colnames(df3) <- colnames(df)
    colnames(df4) <- colnames(df)
    colnames(df5) <- colnames(df)
    colnames(df6) <- colnames(df)
    
    df <- rbind(df, df2, df3, df4, df5, df6)
    
    if (combine == T) {
        
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

                return(df.anno)
            } else {
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                
                return(df.anno)
        }
    } else {
        return(df)
    }

}

#' @title Biomart conversion
#' @param array_ids input
#' @export
#' @import biomaRt
array2ensg <- function(array_ids, array = NULL, biomart = mart, combine = F, df2 = "", by.x = array, by.y = array, all = F){
    # if(is.object(ext_name)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F &  class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame"  ) ext_name <- row.names(ext_name)
    
    if(is.null(array)) {
        t <- listAttributes(mart)
        print(t[grep("probe", t$description), ])
    } else {
    
    df <- getBM(
        filters= array,
        attributes= c(array, "ensembl_gene_id"),
        values= array_ids,
        mart = biomart)
    

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
    }
}

#' @title Biomart conversion
#' @param reactome_id input
#' @export
#' @import biomaRt
#' @examples reactome2ensg()
reactome2ensg <- function(reactome_id, biomart = mart, combine = F, df2 = "", by.x = "reactome", by.y = "reactome", all = F){
    # if(is.object(ext_name)) if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F &  class(get(strsplit(as.character(match.call()), "=")[[2]])) == "data.frame"  ) ext_name <- row.names(ext_name)

    df <- getBM(
        filters= "reactome",
        attributes= c("reactome", "ensembl_gene_id", "external_gene_name"),
        values= reactome_id,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                class <- class(df2)
                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
                class(df.anno) <- class

            } else {

                df2 <- get(strsplit(as.character(match.call()), "=")[[2]]  )
                by.y = "row.names"

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            }
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            return(df.anno)
        }
    } else {
        return(df)
    }
}


#' @title Biomart conversion
#' @param input either a data.frama or a GRanges object
#' @export
#' @import biomaRt
#' @import jsonlite
#' @import httr
#' @import xml2
#' @import doParallel
#' @import foreach
#' @import GenomicRanges
#' @examples chrStartEnd2seq()
chrStartEndStrand2seq <- function(input, cols = c("chromosome_name", "start_position", "end_position", "strand"), server = "http://rest.ensembl.org/", mask = "-", coord_system = "GRCh38"){

cat("See optional arguments at http://rest.ensembl.org/documentation/info/sequence_region\n")
    
    if (coord_system == "GRCh38") cat("coord_system_version is set to GRCh38\ncoord_system_version may also be set to GRCh37\n")
    if(mask == "-") cat('\nSequence(s) are not masked for repeat sequences.\nMask may be set to "soft" or "hard".\nHard will mask all repeats as Ns.\nSoft will mask repeats as lower cased characters.\n')
    
    if(class(input) == "GRanges") {
        for (i in 1:length(input)) {
            pos <- as.data.frame(cbind(as.data.frame(input)$seqnames[i], as.data.frame(ranges(input)[i])[,1:2], as.data.frame(input)$strand[i]))
            colnames(pos) <- cols
            pos$chromosome_name <- gsub("chr", "", pos$chromosome_name) 
            if(!is.numeric(pos$strand) & pos$strand == "+") pos$strand = 1
            if(!is.numeric(pos$strand) & pos$strand == "-") pos$strand = -1
            ext <- paste("/sequence/region/human/", pos$chromosome_name, ":", pos$start, "..", pos$end, ":", pos$strand, "??mask=", mask, "?coord_system_version==", coord_system, sep = "")
            r <- GET(paste(server, ext, sep = ""), content_type("text/plain"))
            mcols(input)$seq[i] <- httr::content(r)
            cat(paste("Completed", i, "of", length(input), "sequences.\n"))
        }
    }
    
    if("data.frame" %in% class(input)) {
        for (i in 1:nrow(input)) {
            pos <- input[, cols]
            pos$chromosome_name <- gsub("chr", "", pos$chromosome_name)
            if(!is.numeric(pos$strand) & pos$strand == "+") pos$strand = 1
            if(!is.numeric(pos$strand) & pos$strand == "-") pos$strand = -1
            ext <- paste("/sequence/region/human/", pos$chromosome_name, ":", pos$start, "..", pos$end, ":", pos$strand, "??mask=", mask, "?coord_system_version==", coord_system, sep = "")
            r <- GET(paste(server, ext, sep = ""), content_type("text/plain"))
            input$seq[i] <- httr::content(r)
            cat(paste("Completed", i, "of", nrow(input), "sequences.\n"))
        }
    }
    
    return(input)    
}

#' @title Biomart conversion
#' @param input either a data.frama or a GRanges object
#' @export
#' @import biomaRt
#' @import jsonlite
#' @import httr
#' @import xml2
#' @import doParallel
#' @import foreach
#' @import GenomicRanges
#' @examples chrStartEnd2seq()
chrStartEndStrand2seq.par <- function(input, cols = c("chromosome_name", "start_position", "end_position", "strand"), server = "http://rest.ensembl.org/", mask = "-", coord_system = "GRCh38", cores = NULL){
    if(is.null(cores)) cores=detectCores()
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    
    cat("See optional arguments at http://rest.ensembl.org/documentation/info/sequence_region\n")
    if (coord_system == "GRCh38") cat("coord_system_version is set to GRCh38\ncoord_system_version may also be set to GRCh37\n")
    if(mask == "-") cat('\nSequence(s) are not masked for repeat sequences.\nMask may be set to "soft" or "hard".\nHard will mask all repeats as Ns.\nSoft will mask repeats as lower cased characters.\n')
    
    # Convert to GRanges if is data.frame
    if("data.frame" %in% class(input)) input <- GRanges(Rle(input$seqnames), IRanges(input$start_position, input$end_position), strand = input$strand)

    if(class(input) == "GRanges") {
        input$seq <- foreach(i=icount(length(input)), .combine =c, .inorder = T, .packages="httr") %dopar% {
            pos <- as.data.frame(cbind(as.data.frame(input)$seqnames[i], as.data.frame(ranges(input)[i])[,1:2], as.data.frame(input)$strand[i]))
            colnames(pos) <- c("chromosome_name", "start_position", "end_position", "strand")
            pos$chromosome_name <- gsub("chr", "", pos$chromosome_name) 
            if(!is.numeric(pos$strand) & pos$strand == "+") pos$strand = 1
            if(!is.numeric(pos$strand) & pos$strand == "-") pos$strand = -1
            ext <- paste("/sequence/region/human/", pos$chromosome_name, ":", pos$start, "..", pos$end, ":", pos$strand, "??mask=", mask, "?coord_system_version==", coord_system, sep = "")
            r <- httr::GET(paste(server, ext, sep = ""), content_type("text/plain"))
            httr::content(r)
        }
    } 
    
    if("data.frame" %in% class(input)) {
        input$seq <- foreach (i=icount(nrow(input)), .combine =c, .packages=c("httr", "data.table")) %dopar% {
            if("data.table" %in% class(input))  pos <- input[, cols, with = F]
            if(!"data.table" %in% class(input)) pos <- input[, cols]
            colnames(pos) <- c("chromosome_name", "start_position", "end_position", "strand")
            pos$chromosome_name <- gsub("chr", "", pos$chromosome_name)
            if(!is.numeric(pos$strand) & pos$strand == "+") pos$strand = 1
            if(!is.numeric(pos$strand) & pos$strand == "-") pos$strand = -1
            ext <- paste("/sequence/region/human/", pos$chromosome_name, ":", pos$start, "..", pos$end, ":", pos$strand, "??mask=", mask, "?coord_system_version==", coord_system, sep = "")
            r <- httr::GET(paste(server, ext, sep = ""), content_type("text/plain"))
            httr::content(r)
        }
    }
    
    stopImplicitCluster()
    stopCluster(cl)
    
    return(input)    
}
