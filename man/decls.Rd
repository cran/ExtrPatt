% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TPDM.R
\name{decls}
\alias{decls}
\title{Declustering}
\usage{
decls(x, th, k)
}
\arguments{
\item{x}{Real vector}

\item{th}{Threshold}

\item{k}{Cluster length}
}
\value{
numeric vector of declustered threshold exceedances
}
\description{
Declustering routine, which will can be applied on radial component r in estimation of the TPDM. Subroutine of \link{est.tpdm}.
}
\references{
Jiang & Cooley (2020) <doi:10.1175/JCLI-D-19-0413.1>
}
\seealso{
\link{est.tpdm}
}
\author{
Yuing Jiang, Dan Cooley
}
