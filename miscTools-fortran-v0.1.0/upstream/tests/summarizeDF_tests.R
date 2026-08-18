## testing summarizeDF()
library( "miscTools" )
data( "cars" )

summarizeDF( cars )

summarizeDF( cars, maxLevel = 10 )

summarizeDF( cars, printValues = FALSE )

tmpFile <- file()
summarizeDF( cars, file = tmpFile )
tmpLines <- readLines( tmpFile )
close( tmpFile )
print( matrix( tmpLines, ncol = 1 ) )


data( "iris" )

summarizeDF( iris )

summarizeDF( iris, printValues = FALSE )

tmpFile <- file()
summarizeDF( iris, file = tmpFile )
tmpLines <- readLines( tmpFile )
close( tmpFile )
print( matrix( tmpLines, ncol = 1 ) )


data( "swiss" )

summarizeDF( swiss )

summarizeDF( swiss, maxLevel = 10 )

summarizeDF( swiss, printValues = FALSE )

tmpFile <- file()
summarizeDF( swiss, file = tmpFile )
tmpLines <- readLines( tmpFile )
close( tmpFile )
print( matrix( tmpLines, ncol = 1 ) )
