% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download_ala.R
\name{download_occurrences}
\alias{download_occurrences}
\title{Download occurrences from the ALA}
\usage{
download_occurrences(taxon, wkt, fields, qa, method = "offline", email,
  download_reason_id, dwcHeaders = TRUE, dst = NULL, read = TRUE,
  dry_run = FALSE, verbose = TRUE)
}
\arguments{
\item{taxon}{(string) Can be a free text search, or taxonomic rank can be used to focus search.}

\item{wkt}{(WKT formatted string) Bounding box to filter geographic locations of occurrences to (optional).}

\item{fields}{(CSV) Fields to return. Defaults to subset of fields.}

\item{qa}{(CSV) Quality assertion fields. Currently set to 'includeall.'}

\item{method}{(string) One of 'offline' or indexed. 'offline' only method available currently.}

\item{email}{(string) Email address to receive path to .zip archive.}

\item{dwcHeaders}{(bool) Return column names as Darwin Core formatted (where available). See website. Default = TRUE.}

\item{dst}{(string) File path to a directory to save downloads in. If NULL, defaults to tempdir().}

\item{read}{(bool) Read file into memory once downloaded? Only relevant if method = indexed, and dst is supplied. Default = TRUE.}

\item{dry_run}{(bool) Return URL to be requested without sending the request. Default: FALSE}

\item{verbose}{(bool) Switch print statements and status checking from the ALA servers on. D
Default: TRUE.}

\item{download_reason}{(int) Valid ALA download reason. See website. Accepts numeric codes only (0:12).}
}
\value{
if verbose, query status
}
\description{
Simplified and customised ALA occurrence download function. 
Largely a hack of the ALA4R function 'occurrences'
}
\note{
strips away much of the functionality of the ALA4R function
}
\seealso{
ALA4R
}
