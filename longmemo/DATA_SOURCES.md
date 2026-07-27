# Dataset provenance

The CSV files in `data/` were converted from datasets distributed with the R
package `longmemo` 1.1-4. The original package credits Jan Beran, Martin
Maechler, and Brandon Whitcher and is distributed under GPL-2.0-or-later. This
file preserves the dataset-specific provenance and references supplied in the
original package documentation.

## `NBSdiff1kg.csv`

National Bureau of Standards weight measurements: deviations from one kilogram
in micrograms, listed chronologically. The original documentation identifies
the source as Jan Beran and Brandon Whitcher by email in fall 1995 and cites:

- H. P. Graf, F. R. Hampel, and J. Tacier (1984), “The problem of unsuspected
  serial correlations,” in *Robust and Nonlinear Time Series Analysis*, Lecture
  Notes in Statistics 26, pp. 127–145, Springer.
- M. Pollak, C. Croakin, and C. Hagwood (1993), *Surveillance schemes with
  applications to mass calibration*, NIST Report 5158.

## `NhemiTemp.csv`

Monthly Northern Hemisphere temperatures for 1854–1989 from the Climate
Research Unit, University of East Anglia. Values are temperature differences
from monthly averages over 1950–1979. The original documentation identifies
the source as Jan Beran and Brandon Whitcher by email in fall 1995 and cites:

- P. D. Jones and K. R. Briffa (1992), “Global surface air temperature
  variations during the twentieth century, part 1,” *The Holocene* 2,
  pp. 165–179.
- Jan Beran (1994), dataset 5, pp. 29–31.

## `NileMin.csv`

Yearly minimum Nile water levels measured at the Roda gauge near Cairo for
years 622–1284. The original documentation says the data supplied by Beran
contained 500 observations and that remaining observations were transcribed
from the cited sources. It cites:

- O. Tousson (1925), *Mémoire sur l'Histoire du Nil*, Mémoire de l'Institut
  d'Egypte.
- Jan Beran (1994), dataset 1, pp. 20–22.

## `ethernetTraffic.csv`

Ethernet traffic from a Bellcore LAN in Morristown, listed chronologically.
The original documentation identifies Leland et al. (1993), Leland and Wilson
(1991), and Jan Beran and Brandon Whitcher by email in fall 1995 as its
provenance.

## `videoVBR.csv`

Variable-bit-rate coded information per frame for a video sequence. The
original documentation cites:

- H. Heeke (1991), “Statistical multiplexing gain for variable bit rate codecs
  in ATM networks,” *International Journal of Digital and Analog Communication
  Systems* 4, pp. 261–268.
- D. Heyman, A. Tabatabai, and T. V. Lakshman (1991), “Statistical analysis and
  simulation of video teleconferencing in ATM networks,” *IEEE Transactions on
  Circuits and Systems for Video Technology* 2, pp. 49–59.
- Jan Beran (1994), dataset 2, pp. 22–23.

The source documentation for all five datasets remains available in the
original `longmemo` package. Conversion to CSV changes representation only and
does not change the applicable package license or original attribution.

