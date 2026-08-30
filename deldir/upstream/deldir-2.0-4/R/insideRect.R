insideRect <- function(x,y,rect,rw) {
# Check that rect is a rectangle.
if(!(is.numeric(rect) &  length(rect) == 4)) {
    whinge <- paste0("Argument \"rect\", if not NULL, must",
                     " be a numeric vector of length 4.\n")
    stop(whinge)
}

# Check that rect is a subset of rw and satifies rect[1] <= rect[2]
# and rect[3] <= rect[4].
ok1  <- rw[1] <= rect[1] & rect[1] <= rect[2] & rect[2] <= rw[2]
ok2  <- rw[3] <= rect[3] & rect[3] <= rect[4] & rect[4] <= rw[4]
if(!(ok1 & ok2)) {
    whinge <- paste0("Argument \"rect\" must determine a subset of",
                     " the rectangular\n  window with respect to which",
                     "  the points were triangulated/tessellated.\n")
    stop(whinge)
}

# All is as it should be; determine which points are inside "rect".
okx  <- rect[1] <= x & x <= rect[2]
oky  <- rect[3] <= y & y <= rect[4]
okx & oky
}

