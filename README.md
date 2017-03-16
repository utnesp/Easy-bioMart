# Easy-bioMart
Provide some functions to make biomart easier to handle

# Installation
```R
library(devtools)
devtools::install_github("utnesp/Easy-bioMart")
library(easybiomart)
```

# User guide

## Init mart(s)

The first step would be to initialize a mart:

```R
# Default mart:
if ( exists("mart") == "FALSE") {
    mart = useMart("ENSEMBL_MART_ENSEMBL", dataset='hsapiens_gene_ensembl')
}
```

Alternative mirrors:
```R
## GRCh38.p3
mart =  useMart("ENSEMBL_MART_ENSEMBL", dataset='hsapiens_gene_ensembl', host="jul2015.archive.ensembl.org")
## Sometimes biomart is down for maintenance, then we can switch to:
mart =  useMart("ENSEMBL_MART_ENSEMBL", dataset='hsapiens_gene_ensembl', host="useast.ensembl.org")
mart =  useMart("ENSEMBL_MART_ENSEMBL", dataset='hsapiens_gene_ensembl', host="uswest.ensembl.org")
## Using GRCh37
mart =  useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
# GRCh37.p12
mart =  useMart(host="sep2013.archive.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl") 
```
You can connect to different marts, by changing host="" in the above code snippet with some of the marts found [here](http://www.ensembl.org/info/website/archives/index.html).




## Using the functions

```R
> ensg2ext_name_biotype("ENSG00000136997", biomart = mart)
  ensembl_gene_id external_gene_name   gene_biotype
1 ENSG00000136997                MYC protein_coding
```

You do not need to specify the mart (default = mart):
```R
> ensg2ext_name_biotype("ENSG00000136997")
  ensembl_gene_id external_gene_name   gene_biotype
1 ENSG00000136997                MYC protein_coding
```

If you have a data.frame, you can figure out what the genes are like this:
```R
> head(test)
  ensembl_gene_id expression
1 ENSG00000067601          6
2 ENSG00000073905          4
3 ENSG00000078319          8
4 ENSG00000080947          7
5 ENSG00000088340          4
6 ENSG00000099251          4

> ensg2ext_name_biotype(test$ensembl_gene_id)
  ensembl_gene_id external_gene_name                       gene_biotype
1 ENSG00000067601             PMS2P4 transcribed_unprocessed_pseudogene
2 ENSG00000073905            VDAC1P1               processed_pseudogene
3 ENSG00000078319             PMS2P1             unprocessed_pseudogene
4 ENSG00000080947            CROCCP3 transcribed_unprocessed_pseudogene
5 ENSG00000088340             FER1L4                 unitary_pseudogene
6 ENSG00000099251          HSD17B7P2 transcribed_unprocessed_pseudogene

# You can combine you intial data.frame with the results passing combine = T:
> ensg2ext_name_biotype(test$ensembl_gene_id, combine = T)
  ensembl_gene_id external_gene_name                       gene_biotype expression
1 ENSG00000067601             PMS2P4 transcribed_unprocessed_pseudogene          6
2 ENSG00000073905            VDAC1P1               processed_pseudogene          4
3 ENSG00000078319             PMS2P1             unprocessed_pseudogene          8
4 ENSG00000080947            CROCCP3 transcribed_unprocessed_pseudogene          7
5 ENSG00000088340             FER1L4                 unitary_pseudogene          4
6 ENSG00000099251          HSD17B7P2 transcribed_unprocessed_pseudogene          4
```

If for instance the ensembl_gene_id's are set as row.names, you can just type: 
```R
> ensg2ext_name_biotype(test)
  ensembl_gene_id external_gene_name                       gene_biotype
1 ENSG00000067601             PMS2P4 transcribed_unprocessed_pseudogene
2 ENSG00000073905            VDAC1P1               processed_pseudogene
3 ENSG00000078319             PMS2P1             unprocessed_pseudogene
4 ENSG00000080947            CROCCP3 transcribed_unprocessed_pseudogene
5 ENSG00000088340             FER1L4                 unitary_pseudogene
6 ENSG00000099251          HSD17B7P2 transcribed_unprocessed_pseudogene

## Whereas typing the below code would not work 
> ensg2ext_name_biotype(row.names(test))
```

The functions should be self-explanatory. If you have any problems, please do not hesitate to ask.
