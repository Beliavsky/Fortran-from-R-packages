#   Tests for rotation functions with random starts (GPFRSorth, GPFRSoblq).
#   Tests compare against values computed previously with this code,
#   to ensure there was no accidental change.
#
#   All rotation tests use set.seed() for reproducibility.
#   fuzz    = 1e-5 standard tolerance
#   fuzz_lo = 1e-4 looser tolerance for BLAS-sensitive tests

require("stats")
require("GPArotation")

fuzz    <- 1e-5
fuzz_lo <- 1e-4
all.ok  <- TRUE


# Helper: sort factors by descending SS loadings and fix signs
sortFac <- function(x) {
  vx   <- order(colSums(x$loadings^2), decreasing = TRUE)
  Dsgn <- diag(sign(colSums(x$loadings^3)))[, vx]
  x$Th      <- x$Th %*% Dsgn
  x$loadings <- x$loadings %*% Dsgn
  if ("Phi" %in% names(x))
    x$Phi <- diag(1 / diag(Dsgn)) %*% x$Phi %*% Dsgn
  x
}


# --- Data ---

data(ability.cov)
L <- loadings(factanal(factors = 2, covmat = ability.cov))


# --- Test 1: Varimax vs stats::varimax, normalize=FALSE ---

set.seed(100)
tst <- sortFac(Varimax(L, normalize = FALSE, randomStarts = 100))$loadings
if (0.001 < max(abs(varimax(L, normalize = FALSE)$loadings - tst))) {
  cat("Test 1 Varimax normalize=FALSE: difference exceeds tolerance\n")
  print(unclass(varimax(L, normalize = FALSE)$loadings) - tst, digits = 18)
  all.ok <- FALSE
}


# --- Test 2: Varimax vs stats::varimax, normalize=TRUE ---

set.seed(100)
tst <- sortFac(Varimax(L, normalize = TRUE, randomStarts = 100))$loadings
if (0.01 < max(abs(varimax(L, normalize = TRUE)$loadings - tst))) {
  cat("Test 2 Varimax normalize=TRUE: difference exceeds tolerance\n")
  print(unclass(varimax(L, normalize = TRUE)$loadings) - tst, digits = 18)
  all.ok <- FALSE
}


# --- Test 3: oblimin with random starts ---

set.seed(99)
v <- sortFac(oblimin(L, eps = 1e-8, randomStarts = 100))$loadings
tst <- t(matrix(c(
         0.3863615904740822504,  0.4745127741495974161,
        -0.0110059418769087539,  0.6458720769633764514,
        -0.0262926272350604423,  0.8961141105684561348,
        -0.0180200526810754824,  0.4882928281695405048,
         0.9900944939102318543, -0.0370718282544326011,
         0.7905657274265397438,  0.0526109550054999417
    ), 2, 6))

if (fuzz < max(abs(v - tst))) {
  cat("Test 3 oblimin randomStarts: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 4: quartimin with random starts ---

set.seed(98)
v <- sortFac(quartimin(L, eps = 1e-8, randomStarts = 100))$loadings
tst <- t(matrix(c(
         0.3863615904740822504,  0.4745127741495974161,
        -0.0110059418769087539,  0.6458720769633764514,
        -0.0262926272350604423,  0.8961141105684561348,
        -0.0180200526810754824,  0.4882928281695405048,
         0.9900944939102318543, -0.0370718282544326011,
         0.7905657274265397438,  0.0526109550054999417
    ), 2, 6))

if (fuzz < max(abs(v - tst))) {
  cat("Test 4 quartimin randomStarts: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 5: entropy with random starts ---

set.seed(99)
v <- sortFac(entropy(L, maxit = 3000, eps = 1e-5, randomStarts = 100))$loadings
tst <- t(matrix(c(
        0.528292107548243184, 0.515443945340967824,
        0.189686511729033253, 0.612116304198454975,
        0.252311894464850861, 0.847442931117894815,
        0.133843268148035738, 0.461156452364903380,
        0.964740133927989407, 0.129750551769587635,
        0.795847094000000532, 0.181751199795689433
    ), 2, 6))

if (0.01 < max(abs(v - tst))) {
  cat("Test 5 entropy randomStarts: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 6: quartimax with random starts ---

set.seed(110)
v <- sortFac(quartimax(L, eps = 1e-7, randomStarts = 100))$loadings
tst <- t(matrix(c(
        0.534714740804540178, 0.508778102568043678,
        0.197348140750149392, 0.609689309353509956,
        0.262919828098457153, 0.844212045390758559,
        0.139616102327241837, 0.459441658926639795,
        0.966291466215733252, 0.117641548844535412,
        0.798063848020893585, 0.171756193883937508
    ), 2, 6))

if (fuzz_lo < max(abs(v - tst))) {
  cat("Test 6 quartimax randomStarts: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}

# Same via GPFRSorth
set.seed(110)
v <- sortFac(GPFRSorth(L, eps = 1e-7, method = "quartimax", randomStarts = 100))$loadings
if (fuzz_lo < max(abs(v - tst))) {
  cat("Test 6b GPFRSorth quartimax randomStarts: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 7: Varimax with random starts ---

set.seed(90)
v <- sortFac(Varimax(L, eps = 1e-8, randomStarts = 100))$loadings
tst <- t(matrix(c(
        0.515866523962843160, 0.527879475961036904,
        0.175054634278874244, 0.616460231981747930,
        0.232057748479543163, 0.853211588623112749,
        0.122822468397975171, 0.464213243286899446,
        0.961376376417989453, 0.152689863976982837,
        0.791292800869773050, 0.200653429940987366
    ), 2, 6))

if (fuzz_lo < max(abs(v - tst))) {
  cat("Test 7 Varimax randomStarts: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}

# Same via GPFRSorth
set.seed(90)
v <- sortFac(GPFRSorth(L, eps = 1e-7, method = "varimax", randomStarts = 100))$loadings
if (fuzz_lo < max(abs(v - tst))) {
  cat("Test 7b GPFRSorth varimax randomStarts: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 8: simplimax with random starts ---

set.seed(89)
v <- sortFac(simplimax(L, eps = 1e-5, randomStarts = 100))$loadings
tst <- t(matrix(c(
        0.3384175759313114429, 0.508414890494446547464,
       -0.0654601124161610648, 0.670992229004664153535,
       -0.1016231721735353366, 0.930535379393095940515,
       -0.0589933707274080121, 0.506904360351960181497,
        0.9733094402675376289, 0.000234046050254643859,
        0.7702037184085044341, 0.085651123319384916965
    ), 2, 6))

if (fuzz_lo < max(abs(v - tst[, 2:1]))) {
  cat("Test 8 simplimax randomStarts: value not same as test value. Difference:\n")
  print(unclass(v - tst[, 2:1]), digits = 18)
  all.ok <- FALSE
}

# Same via GPFRSoblq
set.seed(89)
v <- sortFac(GPFRSoblq(L, eps = 1e-5, method = "simplimax", randomStarts = 100))$loadings
if (fuzz_lo < max(abs(v - tst[, 2:1]))) {
  cat("Test 8b GPFRSoblq simplimax randomStarts: value not same as test value. Difference:\n")
  print(unclass(v - tst[, 2:1]), digits = 18)
  all.ok <- FALSE
}


# --- Test 9: bentlerT with random starts ---

set.seed(88)
v <- sortFac(bentlerT(L, eps = 1e-8, randomStarts = 100))$loadings
tst <- t(matrix(c(
        0.523583611303327312, 0.520226117818945788,
        0.184113022124463677, 0.613815719643687197,
        0.244596116053327067, 0.849702038129718673,
        0.129644684715025493, 0.462354355134084738,
        0.963520501269179652, 0.138517057902201340,
        0.794161628656258278, 0.188979901644201559
    ), 2, 6))

if (fuzz_lo < max(abs(v - tst))) {
  cat("Test 9 bentlerT randomStarts: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}

# Same via GPFRSorth
set.seed(88)
v <- sortFac(GPFRSorth(L, eps = 1e-7, method = "bentler", randomStarts = 100))$loadings
if (fuzz_lo < max(abs(v - tst))) {
  cat("Test 9b GPFRSorth bentler randomStarts: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 10: bentlerQ with random starts ---

set.seed(97)
v <- sortFac(bentlerQ(L, eps = 1e-8, randomStarts = 100))$loadings
tst <- t(matrix(c(
         0.3801726240258240241,  0.4741208368044214638,
        -0.0223632969057368826,  0.6514196922540864687,
        -0.0421105927111659756,  0.9039359851665277334,
        -0.0266594447192576613,  0.4925968005718689424,
         0.9961524457620027917, -0.0485973498906049697,
         0.7939648477384558811,  0.0440983921679098251
    ), 2, 6))

if (fuzz_lo < max(abs(v - tst))) {
  cat("Test 10 bentlerQ randomStarts: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 11: tandemI with random starts ---

set.seed(85)
v <- sortFac(tandemI(L, eps = 1e-6, maxit = 1000, randomStarts = 100))$loadings
tst <- t(matrix(c(
       0.615424480780047745,  0.4074649925368262759,
       0.300894306348887419,  0.5658002819054848143,
       0.406455233467338028,  0.7852483408305571677,
       0.217785179074990981,  0.4279590047675180808,
       0.971977129465111611, -0.0530960591067626969,
       0.815800376450207976,  0.0295946184147908228
    ), 2, 6))

if (fuzz < max(abs(v - tst))) {
  cat("Test 11 tandemI randomStarts: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 12: tandemII with random starts ---

set.seed(84)
v <- sortFac(tandemII(L, eps = 1e-5, randomStarts = 100))$loadings
tst <- t(matrix(c(
       0.512160139332842212, 0.531476249107136312,
       0.170736763115044710, 0.617670057812827134,
       0.226081850628144149, 0.854814488884392154,
       0.119571200821562001, 0.465061309851099225,
       0.960284416460420398, 0.159413208985883820,
       0.789869387186175276, 0.206185467095899383
    ), 2, 6))

if (fuzz < max(abs(v - tst))) {
  cat("Test 12 tandemII randomStarts: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 13: geominT with random starts ---
# Note: algorithm = "legacy" used for reproducibility of reference values

set.seed(81)
v <- sortFac(geominT(L, eps = 1e-5, randomStarts = 100,
                     algorithm = "legacy", fwindow = 1))$loadings
tst <- t(matrix(c(
        0.572197044101002361, 0.4662247895688098054,
        0.243573415560656120, 0.5927388411683653935,
        0.326956608263186954, 0.8215352639437966120,
        0.174476792179181994, 0.4473668997335142894,
        0.972471249855535680, 0.0431091626026945812,
        0.808894688433769660, 0.1099794466209375043
    ), 2, 6))

if (fuzz < max(abs(v - tst))) {
  cat("Test 13 geominT randomStarts: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 14: geominQ with random starts ---

set.seed(80)
v <- sortFac(geominQ(L, eps = 1e-5, randomStarts = 100))$loadings
tst <- t(matrix(c(
         0.39672053553904490508,  0.4713295988080449250,
         0.00424452688463150020,  0.6389466007374070555,
        -0.00510976786312981532,  0.8864521406378518265,
        -0.00646959173137159373,  0.4830101828530461994,
         0.98709860078485589518, -0.0318959930081098297,
         0.79011178369962709045,  0.0558689642678330683
    ), 2, 6))

if (fuzz < max(abs(v - tst))) {
  cat("Test 14 geominQ randomStarts: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 15: cfT with random starts ---

set.seed(70)
v <- sortFac(cfT(L, eps = 1e-8, randomStarts = 100))$loadings
tst <- t(matrix(c(
        0.534721263659975854, 0.508771247100584523,
        0.197355957387199576, 0.609686779159006154,
        0.262930651479430233, 0.844208674501022327,
        0.139621992686633722, 0.459439868910532512,
        0.966292974385164483, 0.117629160286744874,
        0.798066049992627313, 0.171745962120156664
    ), 2, 6))

if (fuzz < max(abs(v - tst))) {
  cat("Test 15 cfT randomStarts: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 16: cfQ with random starts ---

set.seed(60)
v <- sortFac(cfQ(L, eps = 1e-8, randomStarts = 100))$loadings
tst <- t(matrix(c(
         0.3863615904740822504,  0.4745127741495974161,
        -0.0110059418769087539,  0.6458720769633764514,
        -0.0262926272350604423,  0.8961141105684561348,
        -0.0180200526810754824,  0.4882928281695405048,
         0.9900944939102318543, -0.0370718282544326011,
         0.7905657274265397438,  0.0526109550054999417
    ), 2, 6))

if (fuzz < max(abs(v - tst))) {
  cat("Test 16 cfQ randomStarts: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 17: infomaxT with random starts ---

set.seed(55)
v <- sortFac(infomaxT(L, eps = 1e-5, randomStarts = 100))$loadings
tst <- t(matrix(c(
        0.495330443338021176, 0.547195361446864537,
        0.151384273205308784, 0.622695868320644275,
        0.199304253086364791, 0.861451466010626055,
        0.105004533733904976, 0.468565194910632365,
        0.954843809781045660, 0.189293503899924942,
        0.783052579543945471, 0.230726576980168713
    ), 2, 6))

if (fuzz < max(abs(v - tst))) {
  cat("Test 17 infomaxT randomStarts: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 18: infomaxQ with random starts ---

set.seed(50)
v <- sortFac(infomaxQ(L, eps = 1e-5, randomStarts = 100))$loadings
tst <- t(matrix(c(
        0.39327554287862442894,  0.4693137508305071925,
       -0.00319802321222481794,  0.6422985517185823001,
       -0.01549245038490981718,  0.8912279460026399924,
       -0.01214605901641467763,  0.4856544522916727002,
        0.99260028929193111491, -0.0433225495465055510,
        0.79356458059567791530,  0.0471559021503157039
    ), 2, 6))

if (fuzz < max(abs(v - tst))) {
  cat("Test 18 infomaxQ randomStarts: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 19: mccammon with random starts ---

set.seed(40)
v <- sortFac(mccammon(L, eps = 1e-5, randomStarts = 100,
                      algorithm = "legacy", fwindow = 1))$loadings
tst <- t(matrix(c(
       0.4293472299617892007, 0.600363196582340275,
       0.0790140496845253004, 0.635943490060206229,
       0.0992523811009183854, 0.878618107277518656,
       0.0506062164774049028, 0.477512622702450096,
       0.9268544198491108776, 0.297488850382792269,
       0.7514463663627769519, 0.318958389348199534
    ), 2, 6))
tst <- tst[, 2:1]

if (fuzz < max(abs(v - tst))) {
  cat("Test 19 mccammon randomStarts: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 20: oblimax with random starts ---

data(Harman, package = "GPArotation")
set.seed(35)
v <- sortFac(oblimax(Harman8, eps = 1e-5, randomStarts = 100,
                     algorithm = "legacy", fwindow = 1))$loadings
tst <- t(matrix(c(
       0.93395421734409445058, -0.0302013026726007383,
       0.99243032312927881300, -0.1121899246869615951,
       0.96509469978483286567, -0.1322258547171115683,
       0.91647702431117861188, -0.0502569243958834178,
       0.08441855308346873921,  0.8875309317276611765,
       0.04427084251510177149,  0.7907585046311147448,
       0.00332736511424391868,  0.7399752420126202157,
       0.14133359391312094733,  0.6483050831171799366
   ), 2, 8))

if (fuzz < max(abs(v - tst))) {
  cat("Test 20 oblimax randomStarts: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


cat("rotationsRS tests completed.\n")
if (!all.ok) stop("some tests FAILED")