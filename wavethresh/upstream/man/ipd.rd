\name{ipd}
\docType{data}
\alias{ipd}
\title{Inductance plethysmography data.}
\description{
Inductance plethysmography trace. 
}
\usage{data(ipd)}
\source{
This data set contains 4096 observations of inductance plethsymography data sampled at 50Hz starting at 1229.98 seconds. This is a regular time series object. 

I am grateful to David Moshal and Andrew Black of the Department of Anaesthesia, University of Bristol for permission to include this data set. 

This data set was used in Nason, 1996 to illustrate noise reduction with wavelet shrinkage and using cross-validation for choosing the threshold. 

A plethysmograph is an apparatus for measuring variations in the size of parts of the body. In this experiment the inductance plethysmograph consists of a coil of wire encapsulated in a belt. A radio-frequency carrier signal is passed through the wire and size variations change the inductance of the coil that can be detected as a change in voltage. When properly calibrated the output voltage of the inductance plethysmograph is proportional to the change in volume of the part of the body under examination. 

It is of both clinical and scientific interest to discover how anaesthetics or analgesics may alter normal breathing patterns post-operatively. Sensors exist that measure blood oxygen saturation but by the time they indicate critically low levels the patient is often apnoeic (cease breathing) and in considerable danger. It is possible for a nurse to continually observe a patient but this is expensive, prone to error and requires training. In this examples the plethysmograph is arranged around the chest and abdomen of a set of patients and is used to measure the flow of air during breathing. The recordings below were made by the Department of Anaesthesia at the Bristol Royal Infirmary after the patients had undergone surgery under general anaesthetic. The data set (shown below) shows a section of plethysmograph recording lasting approximately 80 seconds. The two main sets of regular oscillations correspond to normal breathing. The disturbed behaviour in the centre of the plot where the normal breathing pattern disappears corresponds to the patient vomiting. 
}
\examples{
#
data(ipd)
\dontrun{ts.plot(ipd)}
}
\keyword{datasets}
\author{G P Nason}
