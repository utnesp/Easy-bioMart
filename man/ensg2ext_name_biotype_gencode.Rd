% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/easybiomart.R
\name{ensg2ext_name_biotype_gencode}
\alias{ensg2ext_name_biotype_gencode}
\title{Biomart conversion}
\usage{
ensg2ext_name_biotype_gencode(
  ensg,
  biomart = mart,
  combine = F,
  df2 = "",
  by.x = "ensembl_gene_id",
  by.y = "row.names",
  all = F,
  title = "",
  ...
)
}
\arguments{
\item{ensg}{input}
}
\description{
Biomart conversion
}
\examples{
ensg2ext_name_biotype()
}
