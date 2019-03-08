crsignal <- function(n, c, l, method = c('a', 'b', 'c')) {
  if (n<10 | n>100) stop("computation only available for n between 10 and 100")
  if (!(c %in% c(0:n-1))) 
    stop(paste("c should be an integer between 0 and",n-1))
  if (!(l %in% c(1:n))) 
    stop(paste("l should be an integer between 1 and",n))
  if (!(method %in% c("a","b","c"))) 
    stop("method should be a, b or c")
  bounds <- data.frame(
    n=10:100,
    ca=c(2,2,3,3,4,4,4,5,5,6,6,6,7,7,8,8,8,9,9,10,10,11,
         11,11,12,12,13,13,14,14,14,15,15,16,16,17,17,17,
         18,18,19,19,20,20,21,21,21,22,22,23,23,24,24,25,
         25,25,26,26,27,27,28,28,29,29,29,30,30,31,31,32,
         32,33,33,34,34,34,35,35,36,36,37,37,38,38,39,39,
         39,40,40,41,41),
    la=c(6,6,7,7,7,7,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,8,8,8,
         8,8,8,8,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9,9,9,9,9,
         9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,
         9,9,9,9,9,9,9,9,9,10,10,10,10,10,10,10,10,10,10),
    cb=c(2,3,3,3,3,4,5,5,5,5,6,7,6,6,6,6,9,9,9,10,11,11,
         11,11,11,12,13,14,13,15,15,15,14,14,17,17,17,17,
         19,19,19,19,19,21,21,21,21,23,23,23,23,23,25,25,
         26,26,27,27,27,28,29,29,29,30,30,31,31,31,32,33,
         33,33,34,33,35,35,35,35,37,37,38,37,39,39,39,39,
         39,41,41,42,41),
    lb=c(6,7,6,6,6,7,8,7,7,7,7,8,7,7,7,7,9,8,8,8,10,9,8,8,
         8,8,9,10,8,11,9,9,8,8,10,9,9,9,12,10,9,9,9,11,10,
         9,9,12,10,10,9,9,11,10,11,10,12,10,10,11,14,11,10,
         11,10,12,11,10,11,13,11,10,11,10,11,11,10,10,12,
         11,12,10,13,11,11,10,10,12,11,12,10),
    cbord=c(3,4,NA,NA,NA,6,6,NA,6,6,NA,NA,7,7,7,NA,10,10,11
            ,NA,12,14,NA,12,13,NA,15,NA,NA,NA,NA,17,NA,NA,
            NA,NA,19,20,20,21,NA,21,21,23,23,NA,23,25,24,26,
            NA,24,27,27,27,27,29,NA,29,29,30,31,30,31,NA,32,
            34,33,33,37,35,NA,36,36,NA,38,36,38,38,39,NA,39,
            41,40,42,NA,41,42,44,43,42),
    lbord=c(5,6,NA,NA,NA,6,7,NA,6,5,NA,NA,6,6,6,NA,7,7,7,NA,
            9,8,NA,7,7,NA,8,NA,NA,NA,NA,8,NA,NA,NA,NA,8,7,11,
            9,NA,8,7,9,8,NA,8,11,9,8,NA,8,9,9,10,9,10,NA,8,
            8,13,9,9,10,NA,9,8,9,8,11,9,NA,10,7,NA,8,9,8,10,
            9,NA,9,12,10,8,NA,8,9,9,10,9))
  res <- FALSE
  if (method=="a") {
    can <- bounds$ca[bounds$n==n]
    lan <- bounds$la[bounds$n==n]
    if (!(c>=can&l<=lan)) res <- TRUE
  } else if (method=="b") {
    cbn <- bounds$cb[bounds$n==n]
    lbn <- bounds$lb[bounds$n==n]
    if (!(c>=cbn&l<=lbn)) res <- TRUE
  } else if (method=="c") {
    cbn <- bounds$cb[bounds$n==n]
    lbn <- bounds$lb[bounds$n==n]
    cbordn <- bounds$cbord[bounds$n==n]
    lbordn <- bounds$lbord[bounds$n==n]
    if (is.na(cbordn)==1 & !(c>=cbn&l<=lbn)) 
      res <- TRUE else if (is.na(cbordn)==0) {
        res <- TRUE
        if (c>=cbn+1&l<=lbn-1) res <- FALSE
        if (c==cbn&l<=lbordn) res <- FALSE
        if (c>=cbordn&l<=lbn) res <- FALSE
      }
  }
  return(res)
} # end function crsignal

crsignal(18,5,7,"a")
crsignal(18,5,8,"a")
crsignal(18,4,7,"a")
crsignal(18,5,7,"b")
crsignal(18,5,8,"b")
crsignal(18,4,7,"b")
crsignal(18,5,7,"c")
crsignal(18,5,8,"c")
crsignal(18,4,7,"c")
crsignal(18,5,6,"c")
crsignal(18,6,7,"c")
