.restr2_eigenv <- function (autovalues, ni.ini, restr.fact, zero.tol)
{
    ev <- autovalues

###### function parameters:
###### ev: matrix containin eigenvalues									#n proper naming - changed autovalue to ev (eigenvalue)
###### ni.ini: current sample size of the clusters						#n proper naming - changed ni.ini to cSize (cluster size)
###### factor: the factor parameter in tclust program
###### init: 1 if we are applying restrictions during the smart inicialization
######       0 if we are applying restrictions during the C steps execution

###### some inicializations

	if (!is.matrix (ev))												#n	checking for malformed ev - argument (e.g. a vector of determinants instead of a matrix)
		if (is.atomic (ev))												#n
			ev <- t (as.numeric (ev))									#n
		else															#n
			ev <- as.matrix (ev)										#n

	stopifnot (ncol (ev) == length (ni.ini))							#n	check wether the matrix autovalues and ni.ini have the right dimension.

	d <- t (ev)

    p <- nrow (ev)
	K <- ncol (ev)	

	n <- sum(ni.ini)

	nis <- matrix(data=ni.ini, nrow=K, ncol=p)
																		#m	MOVED: this block has been moved up a bit. see "old position of "block A" for it's old occurrence
	idx.nis.gr.0 <- nis > zero.tol										#n	as nis is in R we have to be carefull when checking against 0
	used.ev <- ni.ini > zero.tol										#n
	ev.nz <- ev[,used.ev]												#n	non-zero eigenvalues
																		#m
#o	if ((max (d[nis > 0]) <= zero.tol))									#m
#i	if ((max (d[idx.nis.gr.0]) <= zero.tol))							#n	"idx.nis.gr.0" instead of (nis > 0)
	if ((max (ev.nz) <= zero.tol))										#n	simplify syntax
		return (matrix (0, nrow = p, ncol = K))							#m
																		#m
	###### we check if the  eigenvalues verify the restrictions			#m
																		#m
#o	if (max (d[nis > 0]) / min (d[nis > 0]) <= restr.fact)				#m
#i	if (max (d[idx.nis.gr.0]) / min (d[idx.nis.gr.0]) <= restr.fact)	#n	"idx.nis.gr.0" instead of (nis > 0)
	if (max (ev.nz) / min (ev.nz) <= restr.fact)						#n	simplify syntax

	{																	#m
#o		d[!idx.nis.gr.0] <- mean (d[idx.nis.gr.0])						#n	"idx.nis.gr.0" instead of (nis > 0)
		ev[,!used.ev] <- mean (ev.nz)									#n	simplify syntax
		return (ev)														#m
	}																	#m

	###### d_ is the ordered set of values in which the restriction objective function change the definition
	###### points in d_ correspond to  the frontiers for the intervals in which this objective function has the same definition
	###### ed is a set with the middle points of these intervals

#o	d_ <- sort (c (d, d / restr.fact))
	d_ <- sort (c (ev, ev / restr.fact))								#n	using ev instead of d
	dim <- length (d_)													##2do: continue here cleaning up .restr2_eigenv
	d_1 <- d_
	d_1[dim+1] <- d_[dim] * 2
	d_2 <- c (0, d_)
	ed <- (d_1 + d_2) / 2
	dim <- dim + 1;

##o
##o	old position of "block A"
##o

	###### the only relevant eigenvalues are those belong to a clusters with sample size greater than 0.
	###### eigenvalues corresponding to a clusters with 0 individuals has no influence in the objective function
	###### if all the eigenvalues are 0 during the smart initialization we assign to all the eigenvalues the value 1

	###### we build the sol array
	###### sol[1],sol[2],.... this array contains the critical values of the interval functions which defines the m objective function
	###### we use the centers of the interval to get a definition for the function in each interval
	###### this set with the critical values (in the array sol) contains the optimum m value

	t <- s <- r <- array(0, c(K, dim))
	sol <- sal <- array(0, c(dim))

	for (mp_ in 1:dim)
	{
		for (i in 1:K)
		{
			r[i,mp_] <- sum ((d[i,] < ed[mp_])) + sum((d[i,] > ed[mp_]*restr.fact))
			s[i,mp_] <- sum (d[i,]*(d[i,] < ed[mp_]))
			t[i,mp_] <- sum (d[i,]*(d[i,] > ed[mp_] * restr.fact))
		}

		sol[mp_] <- sum (ni.ini / n * (s[,mp_] + t[,mp_] / restr.fact)) / (sum(ni.ini /n * (r[, mp_])))

		e <-	sol[mp_] * (d < sol[mp_]) +
				d * (d >= sol[mp_]) * (d <= restr.fact * sol[mp_]) +
				(restr.fact*sol[mp_]) * (d > restr.fact * sol[mp_])
		o <- -1/2 * nis / n * (log(e) + d / e)

		sal[mp_] <- sum(o)
	}

	###### m is the optimum value for the eigenvalues procedure
#o	eo <- which.max (c (sal))						## remove c ()
#	m <- sol[eo]
	m <- sol[which.max (sal)]						#n
	###### based on the m value we get the restricted eigenvalues

	t (m * (d < m) + d * (d >= m) * (d <= restr.fact * m) + (restr.fact * m) * (d > restr.fact * m))	##	the return value
}

.restr2_deter_old <- function(autovalues, ni.ini, restr.fact, zero.tol=1e-16)
{
    ###### function parameters:
    ###### autovalues: matrix containing eigenvalues
    ###### ni.ini: current sample size of the clusters
    ###### factor: the factor parameter in tclust program
    ###### some initializations

	p = nrow (autovalues)
	K = ncol (autovalues)

	if(p == 1)
		return (.restr2_eigenv (autovalues, ni.ini, restr.fact, zero.tol))

	es = apply (autovalues, 2, prod)

	idx.ni.ini.gr.0 <- ni.ini > zero.tol										#n	as ni.ini is in R we have to be carefull when checking against 0

				######	we check if all the determinants in no empty populations are 0
#o	if (max(es[ni.ini > 0]) <= zero.tol)	##	all eigenvalues are somehow zero.
	if (max(es[idx.ni.ini.gr.0]) <= zero.tol)	#n	"idx.ni.ini.gr.0" instead of (ni.ini > 0)
		return (matrix(0, p, K))										 		##		-> return zero mat

				######  we put in d the determinants of the populations (converting es into a matrix of dim 1 x K)
	d = t(es)	#### --> dim (d) = 1 x K (has once been "d <- matrix (es, nrow = 1)")
																				#n	put this block into a function (for improved readability)
	autovalues_det <- .HandleSmallEv (autovalues, zero.tol)						##	handling close to zero eigenvalues here

print(autovalues)
print(autovalues_det)

#cat ("\n1d^(1/p):\t", d^(1/p), "\n")

					######	we check if all the determinants verify the restrictions
#o	if (max (d[ni.ini > 0]) / min (d[ni.ini > 0]) <= restr.fact)
	if (max (d[idx.ni.ini.gr.0]) / min (d[idx.ni.ini.gr.0]) <= restr.fact)		#n	"idx.ni.ini.gr.0" instead of (ni.ini > 0)
	{
#o		d [ni.ini == 0] <-  mean (d[ni.ini > 0])								## and get the mean - determinants for all clusters without observations.
		d [!idx.ni.ini.gr.0] <-  mean (d[idx.ni.ini.gr.0])						#n	"idx.ni.ini.gr.0" instead of (ni.ini > 0)

		dfin <- d^(1/p)
	}
	else {
		dfin <- .restr2_eigenv(d^(1/p), ni.ini, restr.fact^(1/p), zero.tol)
        cat("\ndfin coming from restr2_eigenv()...\n")
        print(dfin)
    }
					######  we apply the restriction to the determinants by using the .restr2_eigenv function
					######  In order to apply this function is neccessary to transform d and factor with the power (1/p)

#cat ("\nfin:\t", dfin, "\n")


	return(autovalues_det %*% diag(as.numeric(dfin))) 							## autovalues_det %*% diag (dfin)
}

.HandleSmallEv <- function (autovalues, zero.tol)							    #n	a part of .restr2_deter_, which handles almost zero eigenvalues
{	

    ##	handling close to zero eigenvalues here
    ######  populations with one eigenvalue close to 0 are very close to be contained in a hyperplane
    ######  autovalues2 is equal to autovalues except for the columns corresponding to populations close to singular
    ######  for these populations we put only one eigenvalue close to 0 and the rest far from 0   	
	
	p <- nrow (autovalues)													#n
	autovalues[autovalues <= zero.tol] <- zero.tol							#n	"<= zero.tol" for checking for zero

	mi <- apply(autovalues, 2, min)											##	the minimum eigenvalue of each cluster
	ma <- apply(autovalues, 2, max)											##	the maximum eigenvalue of each cluster

	idx.iter <- which(mi/ma <= zero.tol)									#n	making more obvious for what to check!

    ##  cat("\nmi/ma\n")
    ##  print(mi/ma)
    ##  cat("\nmi/zero.tol\n")
    ##  print(mi/zero.tol)
    ##  cat("\nidx.iter\n")
    ##  print(idx.iter)
    
    ##	for each of these clusters. set all "normal" eigenvalues to a high value
	for(i in idx.iter)	{
        ##  cat("\nTo modify index: ", which(autovalues[,i] > mi[i] / zero.tol), "\n")    
        autovalues[autovalues[,i] > mi[i] / zero.tol, i] <- mi[i] / zero.tol
    }
    
    ##  cat("\nModified autovalues\n")
    ##  print(autovalues)

	det = apply(autovalues, 2, prod)											#n	the new determinants

	######	autovalues_det contains the autovalues corrected by the determinant
	######	the product of the eigenvalues of each column in autovalues_det is equal to 1

	autovalues_det <- autovalues %*% diag(det^(-1/p))  					
	return (autovalues_det)
}

.restr2_deter_ <- function(autovalues, ni.ini, restr.factor, cshape=1e10, zero.tol = 1e-16) {
    ## autovalues: matrix containing eigenvalues
    ## ni.ini: current sample size of the clusters
    ## implicitly we are using two restriction parameters: 
    ##  - restr.factor  - constraint level for the determinant constraints
    ##  - cshape        - constraint level for the eigenvalues constraints, default is 1e10


    p <- nrow(autovalues)
    K <- ncol(autovalues)

    autovalues[autovalues < 1e-16] <- 0
    autovalues_ <- autovalues

    if(p == 1)
        return(.restr2_eigenv(autovalues, ni.ini, restr.factor, zero.tol))

    for(k in 1:K)                                
        autovalues_[,k] <- .restr2_eigenv(autovalues=cbind(autovalues[,k]), ni.ini=1, restr.fact=cshape, zero.tol)

    ##  cat("\nautovalues_\n"); print(autovalues_)
    
    es <- apply(autovalues_, 2, prod)
    es[es==0] <- 1
    
    ##  cat("\nes\n"); print(es)
    
    ##  cat("\nesmat\n"); print(matrix((es^(1/p)), ncol=K, nrow=p, byrow=TRUE))

    gm <- autovalues_/matrix((es^(1/p)), ncol=K, nrow=p, byrow=TRUE)

    ##  cat("\ngm\n"); print(gm)

    d <- rbind(apply(autovalues/gm, 2, sum)/p)
    d[is.nan(d)] <- 0

    dfin <- .restr2_eigenv(d, ni.ini, restr.factor^(1/p), zero.tol)
    dout <- matrix(dfin, ncol=K, nrow=p, byrow=TRUE) * (gm * (gm>0) + 1 * (gm==0))  
    
    return (dout)
}
    
