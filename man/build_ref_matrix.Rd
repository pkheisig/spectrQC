\name{build_ref_matrix}
\alias{build_ref_matrix}
\title{Build Reference Matrix from Single Stains}
\usage{build_ref_matrix(files, detectors = NULL)}
\arguments{
\item{files}{Named list of FCS file paths}
\item{detectors}{Vector of detector column names}
}
\value{Matrix (Fluorophores x Detectors), normalized to max=1}
