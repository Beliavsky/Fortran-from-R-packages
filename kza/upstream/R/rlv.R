##Authors: Mingzeng Sun	 Igor G Zurbenko
##Department of Epidemiology and Biostatistics, State University of New York at Albany, New York, USA
##Email: msun@albany.edu; igorg.zurbenko@gmail.com

# Rolling local variance atlas.
#
# The 4.2.0 implementation computes every window's sample variance from
# cumulative-sum tables -- n/(n-1) * (E[x^2] - E[x]^2) with box sums and
# counts built per axis -- instead of calling sd() once per cell, which
# made rlv() O(cells * krnl^d) interpreted R (a 256x256 image at krnl = 9
# was ~65k sd() calls). Two boundary policies:
#
#   pad = "clamp" (default): windows are clamped at the data edge and the
#     variance uses the real count of cells present -- every value in the
#     atlas is computed from data only.
#   pad = "zero": the 4.1.x behavior, exactly -- conceptually the data is
#     zero-padded by the window radius, so windows near the border mix
#     zeros into the variance, biasing it toward the (0 - mean)^2 spread.
#     Kept for reproducing old results.
#
# krnl == 1 selects the legacy one-sided mode (the max over the 2/4/8
# corner windows that include the center), unchanged from 4.1.x.

# ---- legacy 4.1.x implementation ------------------------------------------
# Kept verbatim: it is the krnl == 1 code path, the fallback for data
# containing NA, and the behavioral oracle for pad = "zero" in the tests.
.rlv_legacy <- function(inpt, krnl){
inpt_chk=ncol(rbind(dim(inpt)))
if (length(inpt_chk)==0) {
# 1D vector
#note: make the working space with 2*krnl bigger than the data space at all directions
l=length(inpt)+2*krnl
imgn0= rep(0,l)
imgn1=imgn0
imgn1[(krnl+1):(l-krnl)]=inpt

rolv=imgn0
for (i in (krnl+1):(length(imgn0)-krnl)) {

ssp0 <- imgn1[(i-(0.5*(krnl-1))):(i+(0.5*(krnl-1)))]
ssp0a <- imgn1[(i-1):i] ## left, including self
ssp0b <- imgn1[i:(i+1)] ## right, including self
if (krnl>1) {vrn0=sd(ssp0)*sd(ssp0)
}
if (krnl==1) {vrn0a=sd(ssp0a)*sd(ssp0a)
vrn0b=sd(ssp0b)*sd(ssp0b)
vrn0=max(vrn0a, vrn0b)
}
rolv[i]=vrn0
}
rolVariance=rolv[(krnl+1):(length(rolv)-krnl)]
} else if (length(inpt_chk)==1 & inpt_chk==2) {
# 2D matrix
#note: make the working space with 2*krnl bigger than the data space at both directions
l=nrow(inpt)+2*krnl
w=ncol(inpt)+2*krnl
imgn0= matrix(rep(0,l*w),nrow=l)
imgn1=imgn0
imgn1[(krnl+1):(l-krnl), (krnl+1):(w-krnl)]=inpt

rolv=imgn0
# starting loop
for (i in (krnl+1):(nrow(imgn0)-krnl)) {
for (j in (krnl+1):(ncol(imgn0)-krnl)) {
ssp0 <- imgn1[(i-(0.5*(krnl-1))):(i+(0.5*(krnl-1))), (j-(0.5*(krnl-1))):(j+(0.5*(krnl-1)))]
ssp0a <- imgn1[(i-1):i, j:(j+1)] ## right-above corner, including self
ssp0b <- imgn1[(i-1):i, (j-1):j] ## left-above corner, including self
ssp0c <- imgn1[i:(i+1), (j-1):j] ## left-bottom corner, including self
ssp0d <- imgn1[i:(i+1), j:(j+1)] ## right-bottom corner, including self

if (krnl>1) {vrn0=sd(ssp0)*sd(ssp0)
}
if (krnl==1) {vrn0a=sd(ssp0a)*sd(ssp0a)
vrn0b=sd(ssp0b)*sd(ssp0b)
vrn0c=sd(ssp0c)*sd(ssp0c)
vrn0d=sd(ssp0d)*sd(ssp0d)
vrn0=max(vrn0a, vrn0b, vrn0c, vrn0d)
}
rolv[i,j]=vrn0
}
}
rolVariance=rolv[(krnl+1):(nrow(rolv)-krnl), (krnl+1):(ncol(rolv)-krnl)]

} else if (length(inpt_chk)==1 & inpt_chk==3) {
# 3D array
##it is required that krnl is odd number
#note: make the working space with 2*krnl bigger than the data space at both directions
l=nrow(inpt)+2*krnl
w=ncol(inpt)+2*krnl
h=2*krnl+length(inpt)/(nrow(inpt)*ncol(inpt))
imgn0= array(rep(0,l*w*h),dim=c(l, w, h))
imgn1=imgn0
imgn1[(krnl+1):(l-krnl), (krnl+1):(w-krnl), (krnl+1):(h-krnl)]=inpt

rolv=imgn0
# starting loop 1
for (i in (krnl+1):(l-krnl)) {
for (j in (krnl+1):(w-krnl)) {
for (k in (krnl+1):(h-krnl)) {
ssp0 <- imgn1[(i-(0.5*(krnl-1))):(i+(0.5*(krnl-1))), (j-(0.5*(krnl-1))):(j+(0.5*(krnl-1))), (k-(0.5*(krnl-1))):(k+(0.5*(krnl-1)))]
ssp0a <- imgn1[(i-1):i, (j-1):j, (k-1):k] ## one of the 8 corner, including self
ssp0b <- imgn1[(i-1):i, (j-1):j, k:(k+1)] ## one of the 8 corner, including self
ssp0c <- imgn1[(i-1):i, j:(j+1), (k-1):k] ## one of the 8 corner, including self
ssp0d <- imgn1[(i-1):i, j:(j+1), k:(k+1)] ## one of the 8 corner, including self
ssp0e <- imgn1[i:(i+1), (j-1):j, (k-1):k] ## one of the 8 corner, including self
ssp0f <- imgn1[i:(i+1), (j-1):j, k:(k+1)] ## one of the 8 corner, including self
ssp0g <- imgn1[i:(i+1), j:(j+1), (k-1):k] ## one of the 8 corner, including self
ssp0h <- imgn1[i:(i+1), j:(j+1), k:(k+1)] ## one of the 8 corner, including self

if (krnl>1) {vrn0=sd(ssp0)*sd(ssp0)
}
if (krnl==1) {vrn0a=sd(ssp0a)*sd(ssp0a)
vrn0b=sd(ssp0b)*sd(ssp0b)
vrn0c=sd(ssp0c)*sd(ssp0c)
vrn0d=sd(ssp0d)*sd(ssp0d)
vrn0e=sd(ssp0e)*sd(ssp0e)
vrn0f=sd(ssp0f)*sd(ssp0f)
vrn0g=sd(ssp0g)*sd(ssp0g)
vrn0h=sd(ssp0h)*sd(ssp0h)
vrn0=max(vrn0a, vrn0b, vrn0c, vrn0d, vrn0e, vrn0f, vrn0g, vrn0h)
}
rolv[i,j,k]=vrn0
}
}
}
rolVariance=rolv[(krnl+1):(l-krnl), (krnl+1):(w-krnl), (krnl+1):(h-krnl)]
}
return(rolVariance)
}

# ---- cumulative-sum machinery ---------------------------------------------

# edge-clamped box sum of a vector over [i - r, i + r]
.box1 <- function(v, r) {
    n <- length(v)
    C <- c(0, cumsum(v))
    C[pmin(seq_len(n) + r, n) + 1] - C[pmax(seq_len(n) - r - 1, 0) + 1]
}

# per-position count of cells in the clamped 1D window
.cnt1 <- function(n, r) {
    pmin(seq_len(n) + r, n) - pmax(seq_len(n) - r, 1) + 1
}

.box2 <- function(m, r) {
    tmp <- apply(m, 2, .box1, r = r)      # along rows
    t(apply(tmp, 1, .box1, r = r))        # along columns
}

.box3 <- function(a, r) {
    a <- apply(a, c(2, 3), .box1, r = r)                    # along dim 1
    a <- aperm(apply(a, c(1, 3), .box1, r = r), c(2, 1, 3)) # along dim 2
    aperm(apply(a, c(1, 2), .box1, r = r), c(2, 3, 1))      # along dim 3
}

rlv <- function(inpt, krnl, pad = c("clamp", "zero")) {
    pad <- match.arg(pad)
    if (krnl == 1) return(.rlv_legacy(inpt, krnl))
    if (krnl %% 2 == 0) stop("krnl must be an odd number (the window is centered on each cell); got ", krnl)
    if (anyNA(inpt)) {
        if (pad == "clamp") stop("rlv with pad = \"clamp\" requires complete data; use pad = \"zero\" for the legacy NA behavior")
        return(.rlv_legacy(inpt, krnl))
    }

    r <- (krnl - 1) / 2
    d <- dim(inpt)
    ndim <- length(d)
    a <- inpt
    storage.mode(a) <- "double"

    if (ndim <= 1) {
        a <- as.vector(a)
        S1 <- .box1(a, r); S2 <- .box1(a * a, r)
        n <- if (pad == "zero") krnl else .cnt1(length(a), r)
    } else if (ndim == 2) {
        S1 <- .box2(a, r); S2 <- .box2(a * a, r)
        n <- if (pad == "zero") krnl^2 else outer(.cnt1(d[1], r), .cnt1(d[2], r))
    } else if (ndim == 3) {
        S1 <- .box3(a, r); S2 <- .box3(a * a, r)
        n <- if (pad == "zero") krnl^3 else outer(outer(.cnt1(d[1], r), .cnt1(d[2], r)), .cnt1(d[3], r))
    } else stop("rlv supports 1-, 2-, or 3-dimensional input")

    # sample variance from the sums; clip the tiny negative values float
    # cancellation can produce on flat regions
    pmax((S2 - S1 * S1 / n) / (n - 1), 0)
}
