summarizeDF <- function( dat, printValues = TRUE, maxLevel = 20,
      file = NULL, ... ) {
   if( !inherits( dat, "data.frame" ) ) {
      stop( "argument 'dat' must be a data.frame" )
   }
   if( !is.null( file ) ) {
      sink( file = file, ... )
   }
   cat( "Summary of data.frame\n" )
   cat( "number of observations:", nrow(dat), "\n" )
   cat( "number of variables:", ncol(dat), "\n" )
   cat( "MD5:", digest(dat), "\n\n" )
   
   for( i in 1:length( dat ) ) {
      cat( "variable:", names( dat )[i], "\n" )
      cat( "MD5:", digest(dat[[i]]), "\n" )
      if( isTRUE( printValues ) ) {
         if( is.numeric( dat[[i]] ) ) {
            print( cbind( summary( dat[[i]] ) ) )
            if( length( unique( dat[[i]] ) ) <= maxLevel ) {
               print( cbind( table( dat[[i]], useNA = "ifany" ) ) )
            }
         } else if( is.character( dat[[i]] ) & 
               length( unique( dat[[i]] ) ) <= maxLevel ) {
            print( cbind( table( dat[[i]], useNA = "ifany" ) ) )
         } else if( is.factor( dat[[i]] ) ) {
            if( length( levels( dat[[i]] )) <= maxLevel ) {
               print( table( dat[[i]], useNA = "ifany" ) )
            }
         } else if( is.logical( dat[[i]] ) ) {
            print( table( dat[[i]], useNA = "ifany" ) )
         }
      } else if( printValues == "mean+sd" ) {
         if( is.numeric( dat[[i]] ) ) {
            cat( "mean:", mean( dat[[i]], na.rm = TRUE ), "\n")
            cat( "sd:", sd( dat[[i]], na.rm = TRUE ), "\n" )
            if( sum( is.na( dat[[i]] ) ) > 0 ) {
               cat( "NAs:", sum( is.na( dat[[i]] ) ), "\n" )
            }
         }
      }
      cat( "\n" )
   }
   if( !is.null( file ) ) {
      sink( )
   }
}
