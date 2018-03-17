grd2ind <- function(grd_value, grd){
    ## Projects a grid value of a point to its index
    apply(as.matrix(grd_value), 1, function(x) which.min(abs(grd-x)))
}

## grd2ind(grd_value=c(0,1), grd)
