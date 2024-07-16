\name{CreateNetSciDataCompanionObject}
\alias{CreateNetZooDataCompanionObject}

\title{
Constructor for the \code{NetZooDataCompanion} object.
}
\description{
This function is used to construct a \code{NetZooDataCompanion} object. The member functions of this object are the functions of this package.
}
\usage{
CreateNetZooDataCompanionObject(clinical_patient_file, project_name)
}

\arguments{
  \item{clinical_patient_file}{Path to a comma-separated file containing clinical data for the samples of interest.}
  \item{project_name}{A character string that identifies the project.}
}

\value{A \code{NetZooDataCompanion} object.
}

\author{
Panagiotis Mandros (pmandros@hsph.harvard.edu)

}