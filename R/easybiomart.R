## Default mart:
# if ( exists("mart") == "FALSE") {
#    mart = useMart("ENSEMBL_MART_ENSEMBL", dataset='hsapiens_gene_ensembl')
# }
#
## GRCh38.p3
# mart =  useMart("ENSEMBL_MART_ENSEMBL", dataset='hsapiens_gene_ensembl', host="jul2015.archive.ensembl.org")
## Sometimes biomart is down for maintenance, then we can switch to:
# mart =  useMart("ENSEMBL_MART_ENSEMBL", dataset='hsapiens_gene_ensembl', host="useast.ensembl.org")
# mart =  useMart("ENSEMBL_MART_ENSEMBL", dataset='hsapiens_gene_ensembl', host="uswest.ensembl.org")
## Using GRCh37
# mart =  useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
# GRCh37.p12
# mart =  useMart(host="sep2013.archive.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
#
#
#' @title get cDNA sequence
#' @param input input
#' @export
#' @import biomaRt
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
#' @param input input
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
#' @param input input
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
#' @param input input
#' @export
#' @import biomaRt
#' @examples hugo2ensg()
hugo2ensg <- function(hgnc, biomart = mart, combine = F, df2 = "", by.x = "hgnc_symbol", by.y = "hgnc_symbol", all = F){
    if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F) hgnc <- row.names(hgnc)

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

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)

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
#' @param input input
#' @export
#' @import biomaRt
#' @examples ext_name2ensg()
ext_name2ensg <- function(ext_name, biomart = mart, combine = F, df2 = "", by.x = "external_gene_name", by.y = "external_gene_name", all = F){
    if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F) ext_name <- row.names(ext_name)

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

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)

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
#' @param input input
#' @export
#' @import biomaRt
#' @examples ext_name2chr_band()
ext_name2chr_band <- function(ext_name, biomart = mart, combine = F, df2 = "", by.x = "external_gene_name", by.y = "external_gene_name", all = F){
    if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F) ext_name <- row.names(ext_name)

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

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)

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
#' @param input input
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

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)

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
#' @param input input
#' @export
#' @import biomaRt
#' @examples goid2goname()
goid2goname <- function(goid, biomart = mart, combine = F, df2 = "", by.x = "go_id", by.y = "go_id", all = F){
    if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F) goid <- row.names(goid)

    df <- getBM(
        filters= "go_id",
        attributes= c("go_id", "name_1006"),
        values= goid,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)

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
#' @param input input
#' @export
#' @import biomaRt
#' @examples goid2ensg_extName_biotype()
goid2ensg_extName_biotype <- function(goid, biomart = mart, combine = F, df2 = "", by.x = "go_id", by.y = "go_id", all = F){
    if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F) goid <- row.names(goid)

    df <- getBM(
        filters= "go_id",
        attributes= c("go_id", "ensembl_gene_id", "external_gene_name", "gene_biotype"),
        values= goid,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)

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
#' @param input input
#' @export
#' @import biomaRt
#' @examples goid2ensg_extName_biotype_name_evidence()
goid2ensg_extName_biotype_name_evidence <- function(goid, biomart = mart, combine = F, df2 = "", by.x = "go_id", by.y = "go_id", all = F){
    if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F) goid <- row.names(goid)

    df <- getBM(filters= "go_id",
                attributes= c("go_id", "ensembl_gene_id", "external_gene_name", "gene_biotype", "go_linkage_type"),
                values= goid,
                mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)

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
#' @param input input
#' @export
#' @import biomaRt
#' @examples enst2ensg()
enst2ensg <- function(ensg, biomart = mart, combine = F, df2 = "", by.x = "ensembl_transcript_id", by.y = "ensembl_transcript_id", all = F){
    if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F) ensg <- row.names(ensg)

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

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)

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
#' @param input input
#' @export
#' @import biomaRt
#' @examples ensg2hugo()
ensg2hugo <- function(ensg, biomart = mart, combine = F, df2 = "", by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all = F){
    if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F) ensg <- row.names(ensg)

    df <- getBM(
        filters= "ensembl_gene_id",
        attributes= c("ensembl_gene_id", "hgnc_symbol"),
        values= ensg,
        mart = biomart)

    if(is.na(df$hgnc_symbol == TRUE)) {
        df.na <- getBM(
            filters= "ensembl_gene_id",
            attributes= c("ensembl_gene_id", "external_gene_name"),
            values= ensg,
            mart = biomart)
        df <- df.na
    }

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)

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
#' @param input input
#' @export
#' @import biomaRt
#' @examples ensg2source_status()
ensg2source_status <- function(ensg, biomart = mart, combine = F, df2 = "", by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all = T) {
    if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F) ensg <- row.names(ensg)

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

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)

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
#' @param input input
#' @export
#' @import biomaRt
#' @examples ensg2reactome()
ensg2reactome <- function(ensg, biomart = mart, combine = F, df2 = "", by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all = F){
    if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F) ensg <- row.names(ensg)

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

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)

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
#' @param input input
#' @export
#' @import biomaRt
#' @examples ensg2chr_start_end()
ensg2chr_start_end <- function(ensg, biomart = mart, combine = F, df2 = "", by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all = F, as.IGV = F){
    if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F) ensg <- row.names(ensg)

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

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)

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
#' @param input input
#' @export
#' @import biomaRt
#' @examples ensg2description()
ensg2description <- function(ensg, biomart = mart, combine = F, df2 = "", by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all = F) {
    if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F) ensg <- row.names(ensg)

    df <- getBM(
        filters= "ensembl_gene_id",
        attributes= c("ensembl_gene_id", "description"),
        values= ensg,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)

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
#' @param input input
#' @export
#' @import biomaRt
#' @examples ensg2uniprot()
ensg2uniprot <- function(ensg, biomart = mart, combine = F, df2 = "", by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all = F) {
   if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F) ensg <- row.names(ensg)

   df <- getBM(
        filters= "ensembl_gene_id",
        attributes= c("ensembl_gene_id", "uniprot_sptrembl", "uniprot_swissprot"),
        values= ensg,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)

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
#' @param input input
#' @export
#' @import biomaRt
#' @examples uniprot2ensg()
uniprot2ensg <- function(uniprot_swissprot, useSwissprot = T, biomart = mart, combine = F, df2 = "", by.x = "uniprot_swissprot", by.y = "uniprot_swissprot", all = F) {
    if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F) uniprot_swissprot <- row.names(uniprot_swissprot)

    if (useSwissprot == T) {
        print("Using Swiss-Prot, a manually annotated and reviewed database")
        df <- getBM(
            filters= "uniprot_swissprot",
            attributes= c("ensembl_gene_id", "uniprot_swissprot", "uniprot_genename", "uniprot_genename", "external_gene_name"),
            values= uniprot_swissprot,
            mart = biomart)

        if (combine == T) {
            if (df2 == "") {
                if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)

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
        by.x = "uniprot_sptrembl"
        print("Using TrEMBL, an automatically annotated and not reviewed database")
        df <- getBM(
            filters= "uniprot_sptrembl",
            attributes= c("ensembl_gene_id", "uniprot_sptrembl", "uniprot_genename", "external_gene_name"),
            values= uniprot_swissprot,
            mart = biomart)

        if (combine == T) {
            if (df2 == "") {
                if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)

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
            filters= "uniprot_swissprot",
            attributes= c("ensembl_gene_id", "uniprot_swissprot", "uniprot_genename", "external_gene_name"),
            values= uniprot_swissprot,
            mart = biomart)

        colnames(df) <- c("ensembl_gene_id", "uniprot_swissport_sptrembl", "uniprot_genename", "external_gene_name")

        df3 <- getBM(
            filters= "uniprot_sptrembl",
            attributes= c("ensembl_gene_id", "uniprot_sptrembl", "uniprot_genename", "external_gene_name"),
            values= uniprot_swissprot,
            mart = biomart)

        colnames(df3) <- c("ensembl_gene_id", "uniprot_swissport_sptrembl", "external_gene_name", "uniprot_genename")

        df <- rbind(df, df3)

        by.x = "uniprot_swissport_sptrembl"

        if (combine == T) {
            if (df2 == "") {
                if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)

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
#' @param input input
#' @export
#' @import biomaRt
#' @examples agilent44k2ensg_chr()
agilent44k2ensg_chr <- function(probe_id, biomart = mart, combine = F, df2 = "", by.x = "efg_agilent_wholegenome_4x44k_v1", by.y = "efg_agilent_wholegenome_4x44k_v1", all = F) {
    if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F) probe_id <- row.names(probe_id)

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

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)

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
#' @param input input
#' @export
#' @import biomaRt
#' @examples ensg2GO()
ensg2GO <- function(ensg, biomart = mart, combine = F, df2 = "", by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all = F){
    if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F) ensg <- row.names(ensg)

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

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)

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
#' @param input input
#' @export
#' @import biomaRt
#' @examples hugo2GO()
hugo2GO <- function(ext_name, biomart = mart, combine = F, df2 = "", by.x = "external_gene_name", by.y = "external_gene_name", all = F){
    if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F) ext_name <- row.names(ext_name)

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

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)

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
#' @param input input
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
    row.names(x) <- x$entrezgene
    x$entrezgene <- NULL
    x$ensembl_gene_id <- NULL
    names(x) <- ""
    y <- t(x[1])
    return(y)
}

#' @title Biomart conversion
#' @param input input
#' @export
#' @import biomaRt
#' @examples exon2ensg.is_constitutive()
exon2ensg.is_constitutive <- function(ense_id, biomart = mart, combine = F, df2 = "", by.x = "ensembl_exon_id", by.y = "ensembl_exon_id", all = F){
    if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F) ense_id <- row.names(ense_id)

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

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)

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
#' @param input input
#' @export
#' @import biomaRt
#' @examples exon2gencode.basic()
exon2gencode.basic <- function(ense_id, biomart = mart, combine = F, df2 = "", by.x = "ensembl_exon_id", by.y = "ensembl_exon_id", all = F){
    if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F) ense_id <- row.names(ense_id)

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

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)

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
#' @param input input
#' @export
#' @import biomaRt
#' @examples exon2appris()
exon2appris <- function(ense_id, biomart = mart, combine = F, df2 = "", by.x = "ensembl_exon_id", by.y = "ensembl_exon_id", all = F){
   if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F) ense_id <- row.names(ense_id)

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

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)

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
#' @param input input
#' @export
#' @import biomaRt
#' @examples exon2ensg()
exon2ensg <- function(ense_id, biomart = mart, combine = F, df2 = "", by.x = "ensembl_exon_id", by.y = "ensembl_exon_id", all = F){
    if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F) ense_id <- row.names(ense_id)
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

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)

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
#' @param input input
#' @export
#' @import biomaRt
#' @examples ensg2gene_biotype()
ensg2gene_biotype <- function(ensg, biomart = mart, combine = F, df2 = "", by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all = F){
    if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F) ensg <- row.names(ensg)
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

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)

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
#' @param input input
#' @export
#' @import biomaRt
#' @examples ensg2ext_name()
ensg2ext_name <- function(ensg, biomart = mart, combine = F, df2 = "", by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all = F){
    if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F) ensg <- row.names(ensg)

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

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)

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
#' @param input input
#' @export
#' @import biomaRt
#' @examples ensg2ext_name_biotype()
ensg2ext_name_biotype <- function(ensg, biomart = mart, combine = F, df2 = "", by.x = "ensembl_gene_id", by.y = "row.names", all = F, title = "", ...) {
    if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F) ensg <- row.names(ensg)

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

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)

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
#' @param input input
#' @export
#' @import biomaRt
#' @examples enst2ext_name_biotype()
enst2ext_name_biotype <- function(ensg, biomart = mart, combine = F, df2 = "", by.x = "ensembl_transcript_id", by.y = "ensembl_transcript_id", all = F){
    if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F) ensg <- row.names(ensg)

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

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)

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
#' @param input input
#' @export
#' @import biomaRt
#' @examples reactome.ensg()
reactome.ensg <- function(ensg.diff.exp, biomart = mart, pvalueCutoff=0.05, readable = T) {
try(entrez.ids <- ensg2entrez(ensg.diff.exp, biomart))
    try(print(paste("Returned ", nrow(entrez.ids), " entrez.ids from ", length(ensg.diff.exp), " ensembl_gene_ids", sep = "")))
    try(entrez.ids <- entrez.ids[!is.na(entrez.ids$entrezgene), ])
    try(t <- enrichPathway(gene=entrez.ids$entrezgene, pvalueCutoff=pvalueCutoff, readable=readable))
    try(x <- summary(t))
    try(x$Description <- gsub(" ", "_", x$Description))
    try(print(x))
    return(t)
}

#' @title Biomart conversion
#' @param input input
#' @export
#' @import biomaRt
#' @import limma
#' @examples goana.ensg()
goana.ensg <- function(ensg.diff.exp, ensg.universe = F, biomart = mart, pval = "", ont = "",
                               combine = F, df2 = "", by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all = F
                               ) {

    entrez.ids <- ensg2entrez(ensg.diff.exp, biomart)

    if (ensg.universe == F) {
        ensg.universe <- universe
    } else {
        ensg.universe <- ensg2entrez(ensg.diff.exp, biomart)
    }

    df <- goana(entrez.ids$entrezgene, ensg.universe$entrezgene, species = "Hs")

    df <- df[order(df$P.DE), ]

    if (pval != "") {
        df <- df[df$P.DE < pval, ]
    }

    if (ont != "") {
        df <- df[df$Ont %in% ont, ]
    }


    if (combine == T) {

        k <- ifelse(10 > nrow(df), nrow(df), 10)

        for (i in 1:k) {

            if (i == 1) {
                df.anno <- goid2ensg_extName_biotype(row.names(df)[i])
            } else {
                m <- goid2ensg_extName_biotype(row.names(df)[i])
                df.anno <- rbind(df.anno, m)
            }
        }

        df.anno <- merge(df.anno, df, by.x = "go_id", by.y = "row.names")

        if (df2 == "") {
            # insted of this, we automate this expr to the line below df2 <- get(gsub("\\$.*","", deparse(substitute(formalArgs(exon2ensg.is_constitutive)[1]))))
            # df2 is the first argument passed to this function:
            df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
            # by.y is the column selected in the first argument sent to this function
            by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

            df.anno <- merge(df.anno, df2, by = "ensembl_gene_id")
            df.anno <- df.anno[order(df.anno$P.DE), ]
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            df.anno <- df.anno[order(df.anno$P.DE), ]
            return(df.anno)
        }
    } else {
        return(df)
    }
}

#' @title Biomart conversion
#' @param input input
#' @export
#' @import biomaRt
#' @import limma
#' @examples kegga.ensg()
kegga.ensg <- function(ensg.diff.exp, ensg.universe = F, biomart = mart) {
    entrez.ids <- ensg2entrez(ensg.diff.exp, biomart)
    if (ensg.universe == F) {
        ensg.universe <- universe
    } else {
        ensg.universe <- ensg2entrez(ensg.diff.exp, biomart)
    }
    return(kegga(entrez.ids$entrezgene, ensg.universe$entrezgene))
}

#' @title Biomart conversion
#' @param input input
#' @export
#' @import biomaRt
#' @import ReactomePA
#' @examples reactome.ext_name()
reactome.ext_name <- function(ext_name.diff.exp, return.annotated.counts = F, ranks = 0, biomart = mart, pvalueCutoff=0.05, readable = T,
                              df2 = "", by.x = "external_gene_name", by.y = "external_gene_name", all = F, org = "human", universe = "") {
    try(entrez.ids <- ext_name2entrez(ext_name.diff.exp))
    try(print(paste("Returned ", nrow(entrez.ids), " entrez.ids from ", length(ext_name.diff.exp), " external_gene_names", sep = "")))
    try(entrez.ids <- entrez.ids[!is.na(entrez.ids$entrezgene), ])
    if (universe == "") {
        universe <- as.character(universe.entrez$entrezgene[!is.na(universe.entrez$entrezgene)])
    }
    print(universe)
    try(t <- enrichPathway(gene=entrez.ids$entrezgene, pvalueCutoff=pvalueCutoff, readable=readable, organism = org, universe = universe))
    try(x <- summary(t))
    try(x$Description <- gsub(" ", "_", x$Description))

    if (ranks != 0) {
        try(x <- x[1:ranks,])
    } else {
        try(print(x))
    }

    if (return.annotated.counts == T) {
        for (i in 1:nrow(x)) {
            genes <- unlist(strsplit(x[i, "geneID"], "/"))

            if (i == 1) {
                df <- data.frame(external_gene_name = genes, reactome = rep(x$Description[i], length(genes)), reactome.qval = rep(x$qvalue[i], length(genes)), rank = rep(i, length(genes)))
            } else {
                df <- rbind(df, data.frame(external_gene_name = genes, reactome = rep(x$Description[i], length(genes)), reactome.qval = rep(x$qvalue[i], length(genes)), rank = rep(i, length(genes))))
            }
        }

        if (df2 == "") {
            df2 <- get(gsub("$external_gene_name", "", deparse(substitute(ext_name.diff.exp)), fixed = T))
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            df.anno <- df.anno[order(df.anno$rank), ]
            return(df.anno)
        } else {
            df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)
            df.anno <- df.anno[order(df.anno$rank), ]
            return(df.anno)
        }

    } else {
        return(x)
    }
}

#' @title Biomart conversion
#' @param input input
#' @export
#' @import biomaRt
#' @import limma
#' @examples goana.ext_name()
goana.ext_name<- function(ext_name.diff.exp, ensg.universe = F, biomart = mart) {
    entrez.ids <- ext_name2entrez(ext_name.diff.exp, biomart)
    if (ensg.universe == F) {
        ensg.universe <- universe
    } else {
        ensg.universe <- ext_name2entrez(ext_name.diff.exp, biomart)
    }
    return(goana(entrez.ids$entrezgene, ensg.universe$entrezgene, species = "Hs"))
}

#' @title Biomart conversion
#' @param input input
#' @export
#' @import biomaRt
#' @import limma
#' @examples kegga.ext_name()
kegga.ext_name <- function(ext_name.diff.exp, ensg.universe = F, biomart = mart) {
    entrez.ids <- ext_name2entrez(ext_name.diff.exp, biomart)
    if (ensg.universe == F) {
        ensg.universe <- universe
    } else {
        ensg.universe <- ext_name2entrez(ext_name.diff.exp, biomart)
    }
    return(kegga(entrez.ids$entrezgene, ensg.universe$entrezgene))
}

#' @title Biomart conversion
#' @param input input
#' @export
#' @import biomaRt
#' @examples entrez2ext_name()
entrez2ext_name <- function(entrez_ids, biomart = mart, combine = F, df2 = "", by.x = "entrezgene", by.y = "entrezgene", all = F)
{
    if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F) entrez_ids <- row.names(entrez_ids)
    df <- getBM(
        filters= "entrezgene",
        attributes= c("entrezgene", "external_gene_name"),
        values= entrez_ids,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)

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
#' @param input input
#' @export
#' @import biomaRt
#' @examples entrez2ensg()
entrez2ensg <- function(entrez_ids, biomart = mart, combine = F, df2 = "", by.x = "entrezgene", by.y = "entrezgene", all = F)
{
    if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F) entrez_ids <- row.names(entrez_ids)
    df <- getBM(
        filters= "entrezgene",
        attributes= c("entrezgene", "ensembl_gene_id"),
        values= entrez_ids,
        mart = biomart)


    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)

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
#' @param input input
#' @export
#' @import biomaRt
#' @examples hugo2entrez()
hugo2entrez <- function(hgnc, biomart = mart, combine = F, df2 = "", by.x = "hgnc_symbol", by.y = "hgnc_symbol", all = F){
    if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F) hgnc <- row.names(hgnc)
    df <- getBM(
        filters= "hgnc_symbol",
        attributes= c("hgnc_symbol", "entrezgene"),
        values= hgnc,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)

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
#' @param input input
#' @export
#' @import biomaRt
#' @examples ensg2entrez()
ensg2entrez <- function(ENSG_v82, biomart = mart, combine = F, df2 = "", by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all = F){
    if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F) ENSG_v82 <- row.names(ENSG_v82)

    df <- getBM(
        filters= "ensembl_gene_id",
        attributes= c("ensembl_gene_id", "entrezgene"),
        values= ENSG_v82,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)

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
#' @param input input
#' @export
#' @import biomaRt
#' @examples ext_name2entrez()
ext_name2entrez <- function(external_gene_name, biomart = mart, combine = F, df2 = "", by.x = "external_gene_name", by.y = "external_gene_name", all = F){
    if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F) external_gene_name <- row.names(external_gene_name)

    df <- getBM(
        filters= "external_gene_name",
        attributes= c("external_gene_name", "entrezgene"),
        values= external_gene_name,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)

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
#' @param input input
#' @export
#' @import biomaRt
entrez2ext_name <- function(entrez_ids, biomart = mart, combine = F, df2 = "", by.x = "entrezgene", by.y = "entrezgene", all = F){
    if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == F) entrez_ids <- row.names(entrez_ids)

    df <- getBM(
        filters= "entrezgene",
        attributes= c("entrezgene", "external_gene_name"),
        values= entrez_ids,
        mart = biomart)

    if (combine == T) {
        if (df2 == "") {
            if(grepl("$", strsplit(as.character(match.call()), "=")[[2]], fixed = T) == T) {
                df2 <- get(gsub("\\$.*","", strsplit(as.character(match.call()), "=")[[2]]  ))
                # by.y is the column selected in the first argument sent to this function
                by.y = gsub(".*\\$","", strsplit(as.character(match.call()), "=")[[2]])

                df.anno <- merge(df, df2, by.x = by.x, by.y = by.y, all = all)

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