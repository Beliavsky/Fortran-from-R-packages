ratetable_add_vals <- function(ratetable, add.years=c(), add.ages=c()){
  # The function takes a ratetable and adds years and ages to the object.

  # New attribute object:
  newat <- attributes(ratetable)
  # Find dates (3 or 4 values):
  wh_date <- newat$type>=3
  # Find ages:
  wh_age <- newat$type==2

  # Find unique years/ages not present in the ratetable:
  add.years <- add.years[!(add.years %in% newat$dimnames$year)]
  add.ages <- add.ages[!(add.ages %in% newat$dimnames$age)]
  # All years:
  nova_leta <- sort(c(add.years, newat$dimnames$year))
  nove_starosti <- as.character(sort(as.integer(c(add.ages, newat$dimnames$age))))
  # Find where the new years appear:
  m_nova_leta <- match(nova_leta, add.years, nomatch = 0)
  m_nove_starosti <- match(nove_starosti, add.ages, nomatch = 0)

  # Adjust newat object:
  newat$dimnames$year <- nova_leta
  newat$dim[wh_date] <- length(nova_leta)
  if(length(add.years)>0) newat$cutpoints[[which(wh_date)]] <- sort(c(newat$cutpoints[[which(newat$dimid=='year')]], as.Date(paste0(add.years, '-01-01'))))
  # newat$cutpoints[[which(wh_date)]] <- as.date(newat$cutpoints[[which(wh_date)]])

  newat$dimnames$age <- nove_starosti
  newat$dim[wh_age] <- length(nove_starosti)
  newat$cutpoints[[which(wh_age)]] <- sort(c(newat$cutpoints[[which(newat$dimid=='age')]], as.integer(add.ages)*365.241))

  # prepare vecs:
  vecs <- vector("list", length=length(newat$type))
  for(i in 1:length(newat$type)){
    vecs[[i]] <- 1:newat$dim[i]
  }

  # Define new ratetable:
  out <- array(NA, dim=newat$dim)

  # Add values for existing columns:
  out[m_nove_starosti==0,m_nova_leta==0,] <- ratetable

  # YEAR: For remaining columns either take earlier or later values:
  wh_diff0 <- which(m_nova_leta!=0)
  wh_eq0 <- which(m_nova_leta==0)
  for(ay in wh_diff0){
    # If first year:
    if(ay==1){
      possible_val <- wh_eq0[ay<wh_eq0]
      choose_year <- min(possible_val)
    } else{
      # If later year take last year before:
      possible_val <- wh_eq0[ay>wh_eq0]
      if(length(possible_val)>0){
        choose_year <- max(possible_val)
        # Except if there is no last year, then take next year:
      } else{
        possible_val <- wh_eq0[ay<wh_eq0]
        choose_year <- min(possible_val)
      }
    }
    # Add values to out:
    vecs_tmp <- vecs
    vecs_tmp[[which(wh_date)]] <- ay
    out[vecs_tmp[[1]],vecs_tmp[[2]],vecs_tmp[[3]]] <- out[vecs_tmp[[1]],choose_year,vecs_tmp[[3]]]
  }
  # AGE: For remaining columns either take earlier or later values:
  wh_diff0 <- which(m_nove_starosti!=0)
  wh_eq0 <- which(m_nove_starosti==0)
  for(ay in wh_diff0){
    # If later age take last age before:
    possible_val <- wh_eq0[ay>wh_eq0]
    if(length(possible_val)>0){
      choose_age <- max(possible_val)
      # Except if there is no last age, then take next age:
    } else{
      possible_val <- wh_eq0[ay<wh_eq0]
      choose_age <- min(possible_val)
    }

    # Add values to out:
    vecs_tmp <- vecs
    vecs_tmp[[which(wh_age)]] <- ay
    out[vecs_tmp[[1]],vecs_tmp[[2]],vecs_tmp[[3]]] <- out[choose_age, vecs_tmp[[2]],vecs_tmp[[3]]]
  }

  # Return:
  attributes(out) <- newat
  out
}

# s2 <- slopop
# add.years <- c('1920', "1950", "2024")
# add.ages <- c(110)
#
# newrt <- ratetable_add_vals(s2, add.years, add.ages)
# # newrt[,"1920",] %>% head()
# # newrt[,"1930",] %>% head()
# identical(newrt[,"1920",], newrt[,"1930",])
# identical(newrt[,"1950",], newrt[,"1948",])
# identical(newrt[,"2024",], newrt[,"2021",])
