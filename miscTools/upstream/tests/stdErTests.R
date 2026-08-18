## testing stdEr() methods
library( "miscTools" )
data( "cars" )

# lm()
lmRes <- lm(dist ~ speed, data=cars)
summary( lmRes )
round( stdEr( lmRes ), 4 )

# nls()
nlsRes <- nls( dist ~ b0 + b1 * speed^b2, start = c( b0=0, b1=1, b2=1 ),
   data = cars )
tmp <- summary( nlsRes )
tmp$convInfo$finTol <- round( tmp$convInfo$finTol, 4 )
print( tmp )
round( stdEr( nlsRes ), 3 )
