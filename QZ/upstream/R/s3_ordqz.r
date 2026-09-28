### S3 functions for generalized eigen values.

ordqz <- function(A, B = NULL, cluster = NULL,
    keyword = c("lhp", "rhp", "udi", "udo", "ref", "cef",
                "lhp.fo", "rhp.fo", "udi.fo", "udo.fo"),
    ...){
  ### All keywords for sorting/ording/clustering the generalized eigen values.
  # lhp: Left-half (real(E) < 0)
  # rhp: Right-half (real(E) >= 0)
  # udi: Interior of unit disk (abs(E) < 1)
  # udo: Exterior of unit disk (abs(E) >= 1)
  # ref: Real eigenvalues first (top-left conner)
  # cef: Complex eigenvalues first (top-left conner)
  # lhp.fo: Left-half (real(E) < 0) and finite only
  # rhp.fo: Right-half (real(E) >= 0) and finite only
  # udi.fo: Interior of unit disk (abs(E) < 1) and finite only
  # udo.fo: Exterior of unit disk (abs(E) >= 1) and finite only
  all.keywords <- c("lhp", "rhp", "udi", "udo", "ref", "cef",
                    "lhp.fo", "rhp.fo", "udi.fo", "udo.fo")

  if(!is.null(cluster)){
    select <- cluster
    only.values <- FALSE
    ret <- qz(A = A, B = B, select = select, only.values = only.values, ...)
  } else{
    keyword <- match.arg(tolower(keyword[1]), all.keywords)
    if(!(keyword[1] %in% all.keywords)){
      stop("keyword is unknown.")
    }
  
    # Do usual qz without select.
    select <- NULL
    only.values <- FALSE
    ret <- qz(A = A, B = B, select = select, only.values = only.values, ...)

    # Reordering.
    class.ret <- class(ret)
    tmp.id <- NULL
    select <- rep(FALSE, nrow(A))

    if(keyword == "lhp"){
      # lhp: Left-half (real(E) < 0)
      if(class.ret %in% c("zgges", "dgges")){
        tmp.id <- Re(ret$ALPHA / ret$BETA)
      } else if(class.ret %in% c("zgees", "dgees")){
        tmp.id <- Re(ret$W)
      }
      select[tmp.id < 0] <- TRUE
    } else if(keyword == "rhp"){
      # rhp: Right-half (real(E) >= 0)
      if(class.ret %in% c("zgges", "dgges")){
        tmp.id <- Re(ret$ALPHA / ret$BETA)
      } else if(class.ret %in% c("zgees", "dgees")){
        tmp.id <- Re(ret$W)
      }
      select[tmp.id >= 0] <- TRUE
    } else if(keyword == "udi"){
      # udi: Interior of unit disk (abs(E) < 1)
      if(class.ret %in% c("zgges", "dgges")){
        tmp.id <- Mod(ret$ALPHA / ret$BETA)
      } else if(class.ret %in% c("zgees", "dgees")){
        tmp.id <- Mod(ret$W)
      }
      select[tmp.id < 1] <- TRUE
    } else if(keyword == "udo"){
      # udo: Exterior of unit disk (abs(E) >= 1)
      if(class.ret %in% c("zgges", "dgges")){
        tmp.id <- Mod(ret$ALPHA / ret$BETA)
      } else if(class.ret %in% c("zgees", "dgees")){
        tmp.id <- Mod(ret$W)
      }
      select[tmp.id >= 1] <- TRUE
    } else if(keyword == "lhp.fo"){
      # lhp.fo: Left-half (real(E) < 0) and finite only
      if(class.ret %in% c("zgges", "dgges")){
        tmp.id <- Re(ret$ALPHA / ret$BETA)
      } else if(class.ret %in% c("zgees", "dgees")){
        tmp.id <- Re(ret$W)
      }
      select[tmp.id < 0 & is.finite(tmp.id)] <- TRUE
    } else if(keyword == "rhp.fo"){
      # rhp.fo: Right-half (real(E) >= 0) and finite only
      if(class.ret %in% c("zgges", "dgges")){
        tmp.id <- Re(ret$ALPHA / ret$BETA)
      } else if(class.ret %in% c("zgees", "dgees")){
        tmp.id <- Re(ret$W)
      }
      select[tmp.id >= 0 & is.finite(tmp.id)] <- TRUE
    } else if(keyword == "udi.fo"){
      # udi.fo: Interior of unit disk (abs(E) < 1) and finite only
      if(class.ret %in% c("zgges", "dgges")){
        tmp.id <- Mod(ret$ALPHA / ret$BETA)
      } else if(class.ret %in% c("zgees", "dgees")){
        tmp.id <- Mod(ret$W)
      }
      select[tmp.id < 1 & is.finite(tmp.id)] <- TRUE
    } else if(keyword == "udo.fo"){
      # udo.fo: Exterior of unit disk (abs(E) >= 1) and finite only
      if(class.ret %in% c("zgges", "dgges")){
        tmp.id <- Mod(ret$ALPHA / ret$BETA)
      } else if(class.ret %in% c("zgees", "dgees")){
        tmp.id <- Mod(ret$W)
      }
      select[tmp.id >= 1 & is.finite(tmp.id)] <- TRUE
    } else if(keyword == "ref"){
      # ref: Real eigenvalues first (top-left conner)
      if(class.ret %in% c("zgges", "dgges")){
        tmp.id <- Im(ret$ALPHA / ret$BETA) == 0
      } else if(class.ret %in% c("zgees", "dgees")){
        tmp.id <- Im(ret$W) == 0
      }
      select[tmp.id] <- TRUE
    } else if(keyword == "cef"){
      # cef: Complex eigenvalues first (top-left conner)
      if(class.ret %in% c("zgges", "dgges")){
        tmp.id <- Im(ret$ALPHA / ret$BETA) != 0
      } else if(class.ret %in% c("zgees", "dgees")){
        tmp.id <- Im(ret$W) != 0
      }
      select[tmp.id] <- TRUE
    }

    if(!is.null(B)){
      if(is.complex(A) && is.complex(B)){
        ret <- qz.ztgsen(ret$S, ret$T, ret$Q, ret$Z, select, ijob = 4L,
                         want.Q = TRUE, want.Z = TRUE, ...)
      } else if(is.double(A) && is.double(B)){
        ret <- qz.dtgsen(ret$S, ret$T, ret$Q, ret$Z, select, ijob = 4L,
                         want.Q = TRUE, want.Z = TRUE, ...)
      }
    } else{
      if(is.complex(A)){
        ret <- qz.ztrsen(ret$T, ret$Q, select, job = "B", want.Q = TRUE, ...)
      } else if(is.double(A)){
        ret <- qz.dtrsen(ret$T, ret$Q, select, job = "B", want.Q = TRUE, ...)
      }
    }
  }

  ret
} # End of ordqz().

