\name{CreateNetSciDataCompanionObject}
\alias{CreateNetSciDataCompanionObject}

\title{
Constructor for the \code{NetSciDataCompanion} object.
}
\description{
This function is used to construct a \code{NetSciDataCompanion} object. The member functions of this object are the functions of this package.
}
\usage{
CreateNetSciDataCompanionObject(clinical_patient_file, project_name)
}

\arguments{
  \item{clinical_patient_file}{Path to a comma-separated file containing clinical data for the samples of interest.}
  \item{project_name}{A character string that identifies the project.}
}

\value{A \code{NetSciDataCompanion} object.
}

\author{
Panagiotis Mandros (pmandros@hsph.harvard.edu)

}
