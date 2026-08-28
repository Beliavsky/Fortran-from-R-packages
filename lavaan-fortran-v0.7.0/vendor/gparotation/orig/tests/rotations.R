#   Tests here only compare against values computed previously with this code,
#   to ensure there was no accidental change. It would be better to have
#   comparisons with known correct values.
#
#   Test for oblimax is commented out as it appears to be unstable.
#   All rotation tests use algorithm = "legacy", fwindow = 1 for reproducibility.
#   Reference values were computed with standard BLAS on Linux x86_64.

require("stats")
require("GPArotation")

fuzz    <- 1e-6
fuzz_lo <- 1e-4   # looser tolerance for BLAS-sensitive tests
all.ok  <- TRUE


# --- Varimax vs stats::varimax ---

data(ability.cov)
L <- loadings(factanal(factors = 2, covmat = ability.cov))

if (0.001 < max(abs(varimax(L, normalize = FALSE)$loadings -
                    Varimax(L, normalize = FALSE)$loadings))) {
  cat("Test varimax-normalize=FALSE: difference exceeds tolerance\n")
  print(unclass(varimax(L, normalize = FALSE)$loadings -
        Varimax(L, normalize = FALSE)$loadings), digits = 18)
  all.ok <- FALSE
}

if (0.01 < max(abs(varimax(L, normalize = TRUE)$loadings -
                   Varimax(L, normalize = TRUE, eps = 1e-5)$loadings))) {
  cat("Test varimax-normalize=TRUE: difference exceeds tolerance\n")
  print(unclass(varimax(L, normalize = TRUE)$loadings -
        Varimax(L, normalize = TRUE, eps = 1e-5)$loadings), digits = 18)
  all.ok <- FALSE
}


# --- Test 1: oblimin gam=0 (quartimin) ---

v <- oblimin(L, eps = 1e-8, algorithm = "legacy", fwindow = 1)$loadings
tst <- t(matrix(c(
         0.3863615904740822504,  0.4745127741495974161,
        -0.0110059418769087539,  0.6458720769633764514,
        -0.0262926272350604423,  0.8961141105684561348,
        -0.0180200526810754824,  0.4882928281695405048,
         0.9900944939102318543, -0.0370718282544326011,
         0.7905657274265397438,  0.0526109550054999417
    ), 2, 6))

if (fuzz < max(abs(v - tst))) {
  cat("Test 1 oblimin gam=0: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 2: oblimin gam=1 ---

v <- oblimin(L, gam = 1, eps = 1e-8, algorithm = "legacy", fwindow = 1)$loadings
tst <- t(matrix(c(
       0.2160827,  0.5403732,
      -0.4787025,  1.0224006,
      -0.6800410,  1.4244177,
      -0.3758706,  0.7781390,
       1.4517362, -0.5873946,
       1.1002585, -0.3396290
    ), 2, 6))

if (fuzz < max(abs(v - tst))) {
  cat("Test 2 oblimin gam=1: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 3: oblimin gam=0.1 ---

v <- oblimin(L, gam = .1, eps = 1e-8, algorithm = "legacy", fwindow = 1)$loadings
tst <- t(matrix(c(
       0.37893531,  0.47408606,
      -0.02465543,  0.65257986,
      -0.04530330,  0.90557045,
      -0.02840333,  0.49349573,
       0.99740452, -0.05088932,
       0.79467442,  0.04241284
    ), 2, 6))

if (fuzz < max(abs(v - tst))) {
  cat("Test 3 oblimin gam=0.1: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 4: quartimin ---

v <- quartimin(L, eps = 1e-8, algorithm = "legacy", fwindow = 1)$loadings
tst <- t(matrix(c(
         0.3863615904740822504,  0.4745127741495974161,
        -0.0110059418769087539,  0.6458720769633764514,
        -0.0262926272350604423,  0.8961141105684561348,
        -0.0180200526810754824,  0.4882928281695405048,
         0.9900944939102318543, -0.0370718282544326011,
         0.7905657274265397438,  0.0526109550054999417
    ), 2, 6))

if (fuzz < max(abs(v - tst))) {
  cat("Test 4 quartimin: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 5: targetT ---

v <- targetT(L, Target = matrix(c(rep(1, 3), rep(0, 6), rep(1, 3)), 6, 2),
             eps = 1e-5, algorithm = "legacy", fwindow = 1)$loadings
tst <- t(matrix(c(
        0.551529228817982942, 0.4905002767031292898,
        0.217748645523411000, 0.6027046291262584399,
        0.291173432863349457, 0.8348885228488550636,
        0.154994397662456290, 0.4544843569140373241,
        0.969702339393929247, 0.0850652965070581996,
        0.803390575440818822, 0.1448091121037717866
    ), 2, 6))

if (fuzz < max(abs(v - tst))) {
  cat("Test 5 targetT: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}

# Test deprecated L= argument
v <- targetT(L = L, Target = matrix(c(rep(1, 3), rep(0, 6), rep(1, 3)), 6, 2),
             eps = 1e-5, algorithm = "legacy", fwindow = 1)$loadings
if (fuzz < max(abs(v - tst))) {
  cat("Test 5L targetT (deprecated L=): value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 6: targetQ ---

v <- targetQ(L, Target = matrix(c(rep(1, 3), rep(0, 6), rep(1, 3)), 6, 2),
             eps = 1e-5, algorithm = "legacy", fwindow = 1)$loadings
tst <- t(matrix(c(
        0.735795682866631218, 0.565351705145453853,
        0.433590223819374398, 0.664644550038417159,
        0.589924557708411568, 0.920006940799857786,
        0.317543426981046928, 0.500590650032113116,
        1.021758247914384077, 0.155121528590726393,
        0.872521244896209747, 0.208735706420634437
    ), 2, 6))

if (fuzz < max(abs(v - tst))) {
  cat("Test 6 targetQ: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}

v <- targetQ(L = L, Target = matrix(c(rep(1, 3), rep(0, 6), rep(1, 3)), 6, 2),
             eps = 1e-5, algorithm = "legacy", fwindow = 1)$loadings
if (fuzz < max(abs(v - tst))) {
  cat("Test 6L targetQ (deprecated L=): value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 7: pstT (does not converge but loadings stable) ---

v <- pstT(L, W = matrix(c(rep(.4, 6), rep(.6, 6)), 6, 2),
          Target = matrix(c(rep(1, 3), rep(0, 6), rep(1, 3)), 6, 2),
          maxit = 15000, eps = 1e-5, algorithm = "legacy", fwindow = 1)$loadings
tst <- t(matrix(c(
        0.37067889993474656407, 0.638257130653133720,
        0.01855112570739854416, 0.640564749523800270,
        0.01576132191496706567, 0.884065831441111172,
        0.00524531003824213384, 0.480158078874985073,
        0.89458633399812259590, 0.383762977265515448,
        0.71793428958051475064, 0.388556883222951677
    ), 2, 6))

if (fuzz < max(abs(v - tst))) {
  cat("Test 7 pstT: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}

v <- pstT(L = L, W = matrix(c(rep(.4, 6), rep(.6, 6)), 6, 2),
          Target = matrix(c(rep(1, 3), rep(0, 6), rep(1, 3)), 6, 2),
          maxit = 2000, eps = 1e-5, algorithm = "legacy", fwindow = 1)$loadings
if (fuzz < max(abs(v - tst))) {
  cat("Test 7L pstT (deprecated L=): value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 8: pstQ (does not converge but loadings stable) ---

v <- pstQ(L, W = matrix(c(rep(.4, 6), rep(.6, 6)), 6, 2),
          Target = matrix(c(rep(1, 3), rep(0, 6), rep(1, 3)), 6, 2),
          maxit = 2000, eps = 1e-5, algorithm = "legacy", fwindow = 1)$loadings
tst <- t(matrix(c(
        0.573125161748393785, 0.700868331877288475,
        0.214899397066479453, 0.681727425525818886,
        0.286558275327103040, 0.940272379393286339,
        0.152257795885557295, 0.510481967637567036,
        1.029289798076480578, 0.462598702071116141,
        0.850691132520651205, 0.456859727346562328
    ), 2, 6))

if (fuzz < max(abs(v - tst))) {
  cat("Test 8 pstQ: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}

v <- pstQ(L = L, W = matrix(c(rep(.4, 6), rep(.6, 6)), 6, 2),
          Target = matrix(c(rep(1, 3), rep(0, 6), rep(1, 3)), 6, 2),
          maxit = 2000, eps = 1e-5, algorithm = "legacy", fwindow = 1)$loadings
if (fuzz < max(abs(v - tst))) {
  cat("Test 8L pstQ (deprecated L=): value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 9: entropy ---

v <- entropy(L, maxit = 3000, eps = 1e-5, algorithm = "legacy", fwindow = 1)$loadings
tst <- t(matrix(c(
        0.528292107548243184, 0.515443945340967824,
        0.189686511729033253, 0.612116304198454975,
        0.252311894464850861, 0.847442931117894815,
        0.133843268148035738, 0.461156452364903380,
        0.964740133927989407, 0.129750551769587635,
        0.795847094000000532, 0.181751199795689433
    ), 2, 6))

if (0.01 < max(abs(v - tst))) {
  cat("Test 9 entropy: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 10: quartimax ---

v <- quartimax(L, eps = 1e-5, algorithm = "legacy", fwindow = 1)$loadings
tst <- t(matrix(c(
        0.534714740804540178, 0.508778102568043678,
        0.197348140750149392, 0.609689309353509956,
        0.262919828098457153, 0.844212045390758559,
        0.139616102327241837, 0.459441658926639795,
        0.966291466215733252, 0.117641548844535412,
        0.798063848020893585, 0.171756193883937508
    ), 2, 6))

if (fuzz < max(abs(v - tst))) {
  cat("Test 10 quartimax: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 11: Varimax ---

v <- Varimax(L, eps = 1e-8, algorithm = "legacy", fwindow = 1)$loadings
tst <- t(matrix(c(
        0.515866523962843160, 0.527879475961036904,
        0.175054634278874244, 0.616460231981747930,
        0.232057748479543163, 0.853211588623112749,
        0.122822468397975171, 0.464213243286899446,
        0.961376376417989453, 0.152689863976982837,
        0.791292800869773050, 0.200653429940987366
    ), 2, 6))

if (fuzz < max(abs(v - tst))) {
  cat("Test 11 Varimax: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 12: simplimax ---

v <- simplimax(L, eps = 1e-5, algorithm = "legacy")$loadings
tst <- t(matrix(c(
        0.3384175759313114429, 0.508414890494446547464,
       -0.0654601124161610648, 0.670992229004664153535,
       -0.1016231721735353366, 0.930535379393095940515,
       -0.0589933707274080121, 0.506904360351960181497,
        0.9733094402675376289, 0.000234046050254643859,
        0.7702037184085044341, 0.085651123319384916965
    ), 2, 6))

if (fuzz < max(abs(v - tst))) {
  cat("Test 12 simplimax: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 13: bentlerT ---

v <- bentlerT(L, eps = 1e-8, algorithm = "legacy", fwindow = 1)$loadings
tst <- t(matrix(c(
        0.523583611303327312, 0.520226117818945788,
        0.184113022124463677, 0.613815719643687197,
        0.244596116053327067, 0.849702038129718673,
        0.129644684715025493, 0.462354355134084738,
        0.963520501269179652, 0.138517057902201340,
        0.794161628656258278, 0.188979901644201559
    ), 2, 6))

if (fuzz < max(abs(v - tst))) {
  cat("Test 13 bentlerT: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 14: bentlerQ ---

v <- bentlerQ(L, eps = 1e-8, algorithm = "legacy", fwindow = 1)$loadings
tst <- t(matrix(c(
         0.3801726240258240241,  0.4741208368044214638,
        -0.0223632969057368826,  0.6514196922540864687,
        -0.0421105927111659756,  0.9039359851665277334,
        -0.0266594447192576613,  0.4925968005718689424,
         0.9961524457620027917, -0.0485973498906049697,
         0.7939648477384558811,  0.0440983921679098251
    ), 2, 6))

if (fuzz < max(abs(v - tst))) {
  cat("Test 14 bentlerQ: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 15: tandemI ---

v <- tandemI(L, eps = 1e-5, algorithm = "legacy", fwindow = 1)$loadings
tst <- t(matrix(c(
       0.615424480780047745,  0.4074649925368262759,
       0.300894306348887419,  0.5658002819054848143,
       0.406455233467338028,  0.7852483408305571677,
       0.217785179074990981,  0.4279590047675180808,
       0.971977129465111611, -0.0530960591067626969,
       0.815800376450207976,  0.0295946184147908228
    ), 2, 6))

if (fuzz < max(abs(v - tst))) {
  cat("Test 15 tandemI: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 16: tandemII ---

v <- tandemII(L, eps = 1e-5, algorithm = "legacy", fwindow = 1)$loadings
tst <- t(matrix(c(
       0.512160139332842212, 0.531476249107136312,
       0.170736763115044710, 0.617670057812827134,
       0.226081850628144149, 0.854814488884392154,
       0.119571200821562001, 0.465061309851099225,
       0.960284416460420398, 0.159413208985883820,
       0.789869387186175276, 0.206185467095899383
    ), 2, 6))

if (fuzz < max(abs(v - tst))) {
  cat("Test 16 tandemII: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 17: geominT ---

v <- geominT(L, eps = 1e-5, algorithm = "legacy", fwindow = 1)$loadings
tst <- t(matrix(c(
        0.572197044101002361, 0.4662247895688098054,
        0.243573415560656120, 0.5927388411683653935,
        0.326956608263186954, 0.8215352639437966120,
        0.174476792179181994, 0.4473668997335142894,
        0.972471249855535680, 0.0431091626026945812,
        0.808894688433769660, 0.1099794466209375043
    ), 2, 6))

if (fuzz < max(abs(v - tst))) {
  cat("Test 17 geominT: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 18: geominQ ---

v <- geominQ(L, eps = 1e-5, algorithm = "legacy", fwindow = 1)$loadings
tst <- t(matrix(c(
         0.39672053553904490508,  0.4713295988080449250,
         0.00424452688463150020,  0.6389466007374070555,
        -0.00510976786312981532,  0.8864521406378518265,
        -0.00646959173137159373,  0.4830101828530461994,
         0.98709860078485589518, -0.0318959930081098297,
         0.79011178369962709045,  0.0558689642678330683
    ), 2, 6))

if (fuzz < max(abs(v - tst))) {
  cat("Test 18 geominQ: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 19: bigeominT ---

v <- bigeominT(L, eps = 1e-5, delta = 0.01, algorithm = "legacy", fwindow = 1)$loadings
tst <- t(matrix(c(
       0.735864675930537948,  0.0572554836159610558,
       0.537832587849186305,  0.3484299786828288781,
       0.736693622718023078,  0.4889819217894930681,
       0.398234350199199783,  0.2683070932862005598,
       0.823835376232448180, -0.5185113350380095021,
       0.727481839163037991, -0.3703731487890759011
    ), 2, 6))

if (fuzz < max(abs(v - tst))) {
  cat("Test 19 bigeominT: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 20: bigeominQ ---

v <- bigeominQ(L, eps = 1e-5, delta = 0.01, algorithm = "legacy", fwindow = 1)$loadings
tst <- t(matrix(c(
       0.735864459785110725,  0.0572566968425438083,
       0.537831272508815683,  0.3484308654124375071,
       0.736691776787216535,  0.4889831363831567135,
       0.398233337326754422,  0.2683077498589007681,
       0.823837333630433988, -0.5185099767738823306,
       0.727483237333638733, -0.3703719493837349663
    ), 2, 6))

if (fuzz < max(abs(v - tst))) {
  cat("Test 20 bigeominQ: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 21: cfT (kappa=0) = quartimax ---

v <- cfT(L, eps = 1e-8, algorithm = "legacy", fwindow = 1)$loadings
tst <- t(matrix(c(
        0.534721263659975854, 0.508771247100584523,
        0.197355957387199576, 0.609686779159006154,
        0.262930651479430233, 0.844208674501022327,
        0.139621992686633722, 0.459439868910532512,
        0.966292974385164483, 0.117629160286744874,
        0.798066049992627313, 0.171745962120156664
    ), 2, 6))

if (fuzz < max(abs(v - tst))) {
  cat("Test 21 cfT: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 22: cfQ (kappa=0) = quartimin ---

v <- cfQ(L, eps = 1e-8, algorithm = "legacy", fwindow = 1)$loadings
tst <- t(matrix(c(
         0.3863615904740822504,  0.4745127741495974161,
        -0.0110059418769087539,  0.6458720769633764514,
        -0.0262926272350604423,  0.8961141105684561348,
        -0.0180200526810754824,  0.4882928281695405048,
         0.9900944939102318543, -0.0370718282544326011,
         0.7905657274265397438,  0.0526109550054999417
    ), 2, 6))

if (fuzz < max(abs(v - tst))) {
  cat("Test 22 cfQ: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 23: infomaxT ---

v <- infomaxT(L, eps = 1e-5, algorithm = "legacy", fwindow = 1)$loadings
tst <- t(matrix(c(
        0.495330443338021176, 0.547195361446864537,
        0.151384273205308784, 0.622695868320644275,
        0.199304253086364791, 0.861451466010626055,
        0.105004533733904976, 0.468565194910632365,
        0.954843809781045660, 0.189293503899924942,
        0.783052579543945471, 0.230726576980168713
    ), 2, 6))

if (fuzz < max(abs(v - tst))) {
  cat("Test 23 infomaxT: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 24: infomaxQ ---

v <- infomaxQ(L, eps = 1e-5, algorithm = "legacy", fwindow = 1)$loadings
tst <- t(matrix(c(
        0.39327554287862442894,  0.4693137508305071925,
       -0.00319802321222481794,  0.6422985517185823001,
       -0.01549245038490981718,  0.8912279460026399924,
       -0.01214605901641467763,  0.4856544522916727002,
        0.99260028929193111491, -0.0433225495465055510,
        0.79356458059567791530,  0.0471559021503157039
    ), 2, 6))

if (fuzz < max(abs(v - tst))) {
  cat("Test 24 infomaxQ: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 25: mccammon ---

v <- mccammon(L, eps = 1e-5, algorithm = "legacy", fwindow = 1)$loadings
tst <- t(matrix(c(
       0.4293472299617892007, 0.600363196582340275,
       0.0790140496845253004, 0.635943490060206229,
       0.0992523811009183854, 0.878618107277518656,
       0.0506062164774049028, 0.477512622702450096,
       0.9268544198491108776, 0.297488850382792269,
       0.7514463663627769519, 0.318958389348199534
    ), 2, 6))

if (fuzz < max(abs(v - tst))) {
  cat("Test 25 mccammon: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Tests with box26 (3-factor, fuzz_lo = 1e-4 for BLAS sensitivity) ---

data(Thurstone)

# --- Test 26: equamax ---

v <- equamax(box26, eps = 1e-5, maxit = 2000,
             algorithm = "legacy", fwindow = 1)$loadings
tst <- t(matrix(c(
       0.511813618717971597,  0.1252460667724786814,  0.835031881099661200,
       0.211275278125612587,  0.9469860693462274215,  0.024701038786419674,
       0.923671387190205140,  0.1861505968810791833, -0.278366886980007111,
       0.414270797796799317,  0.7243752493532077397,  0.530526346393166759,
       0.927099794400001564,  0.1710560637343615797,  0.314400690653154735,
       0.685509679739711331,  0.6873945075387188908, -0.212674093365320949,
       0.500975325417812756,  0.4985944480056956341,  0.693100497576226382,
       0.350251174602310256,  0.8631423492204841619,  0.303299191676876356,
       0.809196181501955492,  0.1468111894018074293,  0.540855816747015439,
       1.051940508364259674,  0.2023337382785123650,  0.126016765617061266,
       0.528246625368315792,  0.8145581663496035407, -0.154555803579673606,
       0.791784749686200273,  0.5353191515116044741, -0.254010464723911089,
       0.283760830721282831, -0.7132278971933163625,  0.633221728633476699,
      -0.283760830721282831,  0.7132278971933163625, -0.633221728633476699,
      -0.351981708826951678,  0.0145585781278812498,  0.920862598031950474,
       0.351981708826951678, -0.0145585781278812498, -0.920862598031950474,
      -0.641238077659381234,  0.7340358583767647715,  0.211813801195267382,
       0.641238077659381234, -0.7340358583767647715, -0.211813801195267382,
       0.370916272566192251,  0.7781992933002486179,  0.457012011497068604,
       0.943267697340363864,  0.1458935486092693412,  0.269085717994103968,
       0.683769139477491628,  0.6932804480935084168, -0.193612975261152009,
       0.375506314902942506,  0.7683789003013250518,  0.444462454027040654,
       0.921697465732450816,  0.1542330203892136042,  0.244709944799956780,
       0.664806997738585315,  0.6918110118942031317, -0.165931249557543792,
       0.748952844572093657,  0.5985308972371030656,  0.239842451746804020,
       0.716556890444816297,  0.6343221919993241587,  0.139425892477791219
   ), 3, 26))

if (fuzz_lo < max(abs(v - tst))) {
  cat("Test 26 equamax: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}

# equamax = cfT with kappa = m / (2 * p)
v <- cfT(box26, kappa = (3 / (2 * 26)), eps = 1e-5, maxit = 2000,
         algorithm = "legacy", fwindow = 1)$loadings
if (fuzz_lo < max(abs(v - tst))) {
  cat("Test 26b cfT(kappa=equamax): value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 27: parsimax ---

v <- parsimax(box26, eps = 1e-5, maxit = 2000,
              algorithm = "legacy", fwindow = 1)$loadings
tst <- t(matrix(c(
       0.7201835790622810318, -0.2820790149262949464,  0.6137467244615277817,
      -0.0679423851913938670,  0.6010788795762025405,  0.7590243822315081434,
       0.6707172136894012926,  0.7174085874409354968, -0.0277909684381195121,
       0.3564975652920873705,  0.2169149780725644350,  0.8964682806594965747,
       0.8905375652375422391,  0.2961288407436278303,  0.3269136806262873951,
       0.3419238671745732927,  0.8395902605853544642,  0.4072361273101923196,
       0.5551801796495600128,  0.0181442676741417341,  0.8193941381745735164,
       0.1791815177720131602,  0.4232257028261393605,  0.8651329309165379788,
       0.8735154515073240145,  0.0707778916688399651,  0.4481498031114935499,
       0.9254209439635597834,  0.5018907196720771013,  0.2347334274887991901,
       0.1871266535649374341,  0.7975561080494565358,  0.5434380093797205324,
       0.4642991413806251688,  0.8329607560592182658,  0.2619421428072832847,
       0.6793187635833454197, -0.7070938474525871875, -0.1694942722875707464,
      -0.6793187635833454197,  0.7070938474525871875,  0.1694942722875707464,
       0.0126413340577520572, -0.7959999181796318934,  0.5816488003350992475,
      -0.0126413340577520572,  0.7959999181796318934, -0.5816488003350992475,
      -0.7005089039174136056, -0.0349340740508546910,  0.7091733821870598309,
       0.7005089039174136056,  0.0349340740508546910, -0.7091733821870598309,
       0.2764675266718980007,  0.2782887989292237019,  0.8933946782281918519,
       0.8957285854821546156,  0.3212930506383097073,  0.2790825626255922787,
       0.3455615085575033385,  0.8287155760723180498,  0.4236528505493507568,
       0.2788032457489866833,  0.2837361452757249936,  0.8779069142135194070,
       0.8654448252950356357,  0.3331277242202841382,  0.2706070467039994321,
       0.3390107603224090660,  0.7999236919074760310,  0.4396564471388281214,
       0.5855440811683275681,  0.5029973433114328651,  0.6171108503586539840,
       0.5106729412991620753,  0.5782224406112924653,  0.5832066153589001711
   ), 3, 26))

if (fuzz_lo < max(abs(v - tst))) {
  cat("Test 27 parsimax: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}

# parsimax = cfT with kappa = (m-1) / (p+m-2)
v <- cfT(box26, kappa = ((3 - 1) / (26 + 3 - 2)), eps = 1e-5, maxit = 2000,
         algorithm = "legacy", fwindow = 1)$loadings
if (fuzz_lo < max(abs(v - tst))) {
  cat("Test 27b cfT(kappa=parsimax): value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 28: varimin ---

data(Harman, package = "GPArotation")
v <- varimin(Harman8, eps = 1e-5, algorithm = "legacy", fwindow = 1)$loadings
tst <- t(matrix(c(
       0.800626657046876855, -0.452452158825595752,
       0.783606930490612252, -0.524447498313301397,
       0.742635936060292656, -0.522609669324872517,
       0.768357486963803682, -0.455227165519225097,
       0.818696625686402668,  0.444445536696790211,
       0.702064973637186673,  0.410429985249392060,
       0.623283524595303340,  0.401857745935120247,
       0.668480210595655877,  0.287458184858228272
   ), 2, 8))

if (fuzz_lo < max(abs(v - tst))) {
  cat("Test 28 varimin: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 29: oblimax ---

v <- oblimax(Harman8, eps = 1e-5, algorithm = "legacy", fwindow = 1)$loadings
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

if (fuzz_lo < max(abs(v - tst))) {
  cat("Test 29 oblimax: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 30: echelon ---

data(WansbeekMeijer)
fa.unrotated <- factanal(factors = 2, covmat = NetherlandsTV, rotation = "none")
v   <- echelon(fa.unrotated$loadings)$loadings
tst <- matrix(c(
       0.7910866, 0.000000000,
       0.8356693, 0.086562877,
       0.7732175, 0.003599155,
       0.4712891, 0.520430204,
       0.4716313, 0.596050663,
       0.4340904, 0.693182608,
       0.3934293, 0.610972389
   ), ncol = 2, nrow = 7, byrow = TRUE)

if (fuzz_lo < max(abs(v - tst))) {
  cat("Test 30 echelon: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 31: eiv ---

v   <- eiv(fa.unrotated$loadings)$loadings
tst <- matrix(c(
        1.0000000, 0.0000000,
        0.0000000, 1.0000000,
        0.9334902, 0.0415785,
       -5.7552381, 6.0121639,
       -6.6776277, 6.8857539,
       -7.9104168, 8.0078509,
       -6.9585766, 7.0581340
   ), ncol = 2, nrow = 7, byrow = TRUE)

if (fuzz_lo < max(abs(v - tst))) {
  cat("Test 31 eiv: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 32: bifactorT ---

data(WansbeekMeijer)
fa.unrotated <- factanal(factors = 3, covmat = NetherlandsTV, rotation = "none")
v <- bifactorT(fa.unrotated$loadings, algorithm = "legacy", fwindow = 1)$loadings
tst <- matrix(c(
       0.605259846207760410,  0.50926282550276680272, -0.0326441696928007064,
       0.674848533008883589,  0.49892997353743967492,  0.0360559565502543422,
       0.585307624074672961,  0.50452427555859225006, -0.0116655363009391944,
       0.569143997010692737,  0.07907134469847305891,  0.4496106087464453172,
       0.639591951170676909, -0.00317621525167442742,  0.3903956610664964244,
       0.644285080202307459, -0.06194060671151470354,  0.4889273316625122323,
       0.894466584777320994, -0.44135797428894490979, -0.0115183403900823399
   ), ncol = 3, nrow = 7, byrow = TRUE)

if (fuzz_lo < max(abs(v - tst))) {
  cat("Test 32 bifactorT: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


# --- Test 33: bifactorQ ---
# Note: uses fuzz_lo due to BLAS-dependent numerical sensitivity

v <- bifactorQ(fa.unrotated$loadings, algorithm = "legacy", fwindow = 1)$loadings
tst <- matrix(c(
       0.639628269247016989,  0.4617887681341789063, -0.03408526722627130967,
       0.708202536330854282,  0.4531948277009248405,  0.03509478386512845938,
       0.619356995328652293,  0.4604997427039814184, -0.01303680319918767982,
       0.572539420889668138,  0.0819716834028682007,  0.45225086629811012129,
       0.637128483915923693, -0.0107264127434737802,  0.39334646255480443244,
       0.637511978242973232, -0.0601392817604798208,  0.49259101012935263553,
       0.861232867683719872, -0.5044079449050349329, -0.00751942773799910234
   ), ncol = 3, nrow = 7, byrow = TRUE)

if (fuzz_lo < max(abs(v - tst))) {
  cat("Test 33 bifactorQ: value not same as test value. Difference:\n")
  print(unclass(v - tst), digits = 18)
  all.ok <- FALSE
}


cat("rotations tests completed.\n")
if (!all.ok) stop("some tests FAILED")