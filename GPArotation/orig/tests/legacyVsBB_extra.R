# Additional legacy vs BB convergence tests using package datasets.
# Verifies that BB+fw10 finds same or better solution than legacy
# across multiple datasets and criteria using fixed, hardcoded Tmat matrices.

require("GPArotation")

all.ok <- TRUE
fuzz   <- 1e-4   # loose tolerance --- same solution within rounding

data("CCAI",           package = "GPArotation")
data("WansbeekMeijer", package = "GPArotation")
data("GriffithMulaik", package = "GPArotation")
data("Harman",         package = "GPArotation")
data("Thurstone",      package = "GPArotation")

# Helper: check BB finds solution at least as good as legacy
check_bb <- function(label, leg, bb) {
  # Extract final objective value 'f' from the last row of Table
  f.leg <- leg$Table[nrow(leg$Table), "f"]
  f.bb  <- bb$Table[nrow(bb$Table), "f"]
  
  # BB should find same or better (lower) criterion value
  if (f.bb > f.leg + fuzz) {
    cat(label, ": BB found worse solution than legacy\n")
    cat("  legacy f =", f.leg, "  BB f =", f.bb, "\n")
    all.ok <<- FALSE
  }
  
  # Both should converge
  if (!bb$convergence) {
    cat(label, ": BB did not converge\n")
    all.ok <<- FALSE
  }
}


# ==============================================================================
# --- Explicit Hardcoded Starting Matrices ---
# ==============================================================================

T3_ccai <- structure(c(0.431414381712226, -0.894744332740752, 0.115387218876984, 
                       -0.0776174383091876, -0.164239860314741, -0.983361989073258, 
                       -0.898808747466433, -0.415280444172491, 0.140303200839187), dim = c(3L, 3L))

T2_nl <- structure(c(0.434315353073653, -0.900760886187066, 0.900760886187066, 
                     0.434315353073653), dim = c(2L, 2L))

T6_gm <- structure(c(0.181373108294974, -0.376164002957957, 0.0485105259174951, 
                     -0.0295264563740687, -0.232257082268142, -0.876604816009788, 
                     -0.746940052505326, -0.522756224813139, 0.052758861370963, -0.361585308596306, 
                     -0.142312065982704, 0.122581933466333, -0.150353031847131, -0.207811107487956, 
                     -0.860627239445551, 0.408803078608172, 0.156233992577518, -0.0447241168398343, 
                     -0.0432441982549581, 0.287431949041387, -0.215875719270493, 0.0405111134904342, 
                     -0.925918263678788, 0.0997231492023105, -0.580746432095367, 0.375678001257268, 
                     0.311043225538829, 0.574862048564166, 0.0640150636997622, -0.300478667406478, 
                     0.217811478785953, -0.564171986138982, 0.332917064480368, 0.607582727458707, 
                     -0.199930747462431, 0.338090643845403), dim = c(6L, 6L))

T3_box <- structure(c(0.8122394, -0.4721832,  0.3418291,
                      0.5630312,  0.7911029, -0.2381044,
                     -0.1580211,  0.3878192,  0.9080702), dim = c(3L, 3L))
					 
T2_harman <- T2_nl

# ==============================================================================
# --- CCAI 3-factor ---
# ==============================================================================

fa.ccai <- factanal(factors = 3, covmat = CCAI_R,
                    n.obs = 461, rotation = "none")

check_bb("CCAI oblimin",
  oblimin(fa.ccai, Tmat = T3_ccai, algorithm = "legacy", fwindow = 1),
  oblimin(fa.ccai, Tmat = T3_ccai, algorithm = "bb",     fwindow = 10))

check_bb("CCAI Varimax",
  Varimax(fa.ccai, Tmat = T3_ccai, algorithm = "legacy", fwindow = 1),
  Varimax(fa.ccai, Tmat = T3_ccai, algorithm = "bb",     fwindow = 10))

check_bb("CCAI geominQ",
  geominQ(fa.ccai, Tmat = T3_ccai, algorithm = "legacy", fwindow = 1),
  geominQ(fa.ccai, Tmat = T3_ccai, algorithm = "bb",     fwindow = 10))

check_bb("CCAI tandemI",
  tandemI(fa.ccai, Tmat = T3_ccai, algorithm = "legacy", fwindow = 1),
  tandemI(fa.ccai, Tmat = T3_ccai, algorithm = "bb",     fwindow = 10))

check_bb("CCAI bentlerQ",
  bentlerQ(fa.ccai, Tmat = T3_ccai, algorithm = "legacy", fwindow = 1),
  bentlerQ(fa.ccai, Tmat = T3_ccai, algorithm = "bb",     fwindow = 10))

check_bb("CCAI quartimin",
  quartimin(fa.ccai, Tmat = T3_ccai, algorithm = "legacy", fwindow = 1),
  quartimin(fa.ccai, Tmat = T3_ccai, algorithm = "bb",     fwindow = 10))

# ==============================================================================
# --- NetherlandsTV 2-factor ---
# ==============================================================================

fa.nl <- factanal(factors = 2, covmat = NetherlandsTV,
                  rotation = "none")

check_bb("NetherlandsTV oblimin",
  oblimin(fa.nl, Tmat = T2_nl, algorithm = "legacy", fwindow = 1),
  oblimin(fa.nl, Tmat = T2_nl, algorithm = "bb",     fwindow = 10))

check_bb("NetherlandsTV Varimax",
  Varimax(fa.nl, Tmat = T2_nl, algorithm = "legacy", fwindow = 1),
  Varimax(fa.nl, Tmat = T2_nl, algorithm = "bb",     fwindow = 10))

check_bb("NetherlandsTV geominT",
  geominT(fa.nl, Tmat = T2_nl, algorithm = "legacy", fwindow = 1),
  geominT(fa.nl, Tmat = T2_nl, algorithm = "bb",     fwindow = 10))

check_bb("NetherlandsTV targetT",
  targetT(fa.nl, Tmat = T2_nl,
          Target = matrix(c(1,1,0,0,1,0,0,0,1,1,0,1,0,1), 7, 2),
          algorithm = "legacy", fwindow = 1),
  targetT(fa.nl, Tmat = T2_nl,
          Target = matrix(c(1,1,0,0,1,0,0,0,1,1,0,1,0,1), 7, 2),
          algorithm = "bb",     fwindow = 10))

# ==============================================================================
# --- GriffithMulaik 6-factor ---
# ==============================================================================

fa.gm <- factanal(factors = 6, covmat = GriffithMulaik,
                  n.obs = 523, rotation = "none")

check_bb("GriffithMulaik oblimin",
  oblimin(fa.gm, Tmat = T6_gm, algorithm = "legacy", fwindow = 1),
  oblimin(fa.gm, Tmat = T6_gm, algorithm = "bb",     fwindow = 10))

check_bb("GriffithMulaik geominQ",
  geominQ(fa.gm, Tmat = T6_gm, algorithm = "legacy", fwindow = 1),
  geominQ(fa.gm, Tmat = T6_gm, algorithm = "bb",     fwindow = 10))

check_bb("GriffithMulaik Varimax",
  Varimax(fa.gm, Tmat = T6_gm, algorithm = "legacy", fwindow = 1),
  Varimax(fa.gm, Tmat = T6_gm, algorithm = "bb",     fwindow = 10))

# ==============================================================================
# --- Harman8 2-factor ---
# ==============================================================================

check_bb("Harman8 quartimax",
  quartimax(Harman8, Tmat = T2_harman, algorithm = "legacy", fwindow = 1),
  quartimax(Harman8, Tmat = T2_harman, algorithm = "bb",     fwindow = 10))

check_bb("Harman8 bentlerT",
  bentlerT(Harman8, Tmat = T2_harman, algorithm = "legacy", fwindow = 1),
  bentlerT(Harman8, Tmat = T2_harman, algorithm = "bb",     fwindow = 10))

check_bb("Harman8 geominT",
  geominT(Harman8, Tmat = T2_harman, algorithm = "legacy", fwindow = 1),
  geominT(Harman8, Tmat = T2_harman, algorithm = "bb",     fwindow = 10))

# ==============================================================================
# --- box26 3-factor ---
# ==============================================================================

check_bb("box26 oblimin",
  oblimin(box26, Tmat = T3_box, algorithm = "legacy", fwindow = 1),
  oblimin(box26, Tmat = T3_box, algorithm = "bb",     fwindow = 10))

check_bb("box26 equamax",
  equamax(box26, Tmat = T3_box, algorithm = "legacy", fwindow = 1),
  equamax(box26, Tmat = T3_box, algorithm = "bb",     fwindow = 10))

check_bb("box26 parsimax",
  parsimax(box26, Tmat = T3_box, algorithm = "legacy", fwindow = 1),
  parsimax(box26, Tmat = T3_box, algorithm = "bb",     fwindow = 10))

check_bb("box26 geominQ",
  geominQ(box26, Tmat = T3_box, algorithm = "legacy", fwindow = 1),
  geominQ(box26, Tmat = T3_box, algorithm = "bb",     fwindow = 10))

# ==============================================================================
# --- Cayley vs legacy (orthogonal only) ---
# ==============================================================================

check_bb("Harman8 Varimax cayley",
  Varimax(Harman8, Tmat = T2_harman, algorithm = "legacy", fwindow = 1),
  Varimax(Harman8, Tmat = T2_harman, algorithm = "cayley", fwindow = 10))

check_bb("Harman8 quartimax cayley",
  quartimax(Harman8, Tmat = T2_harman, algorithm = "legacy", fwindow = 1),
  quartimax(Harman8, Tmat = T2_harman, algorithm = "cayley", fwindow = 10))

check_bb("box26 Varimax cayley",
  Varimax(box26, Tmat = T3_ccai, algorithm = "legacy", fwindow = 1),
  Varimax(box26, Tmat = T3_ccai, algorithm = "cayley", fwindow = 10))


cat("legacyVsBB extra tests completed.\n")
if (!all.ok) stop("some tests FAILED")