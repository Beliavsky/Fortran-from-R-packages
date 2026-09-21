

qn <- function (x, corrFact)
{
	if (missing (corrFact))	
		corrFact = 1 / (sqrt(2) * qnorm(5/8))

	ret.C <- .C (C_qn,
			as.integer (length (x)),
			as.double (corrFact),
			parOutD = double (1),
			x = as.double (x))
	ret.C$parOutD
}
