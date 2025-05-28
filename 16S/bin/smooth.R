#!/usr/bin/Rscript
args <- commandArgs(T)
if (length(args) != 2) {
  cat("Rscript smooth.R <input> <output>\n")
  q()
}
data <-
  read.table(
    args[1],
    sep = "\t",
    header = T,
    check.names = F,
    na.strings = ""
  )
sp <- data[, 1]

medi <- function(x) {
  notna <- T #!is.na(x)
  a <- sp[notna]
  b <- x[notna]
  m <- tapply(b, a, median)
  return(m)
}

mea2 <- function(x) {
  notna <- !is.na(x)
  a <- sp[notna]
  b <- x[notna]
  m <- tapply(b, a, mean, na.rm = T)
  return(m)
}

smo0 <- function(x) {
  notna <- T #!is.na(x)
  a <- sp[notna]
  b <- x[notna]
  low <- (a < 10000) & !is.na(b)
  a1 <- a[low]
  b1 <- b[low]
  sm1 <- b1
  sm1[!is.na(sm1)] <- loess(b1 ~ a1, span = .2)$fitted

  a2 <- a[!low]
  b2 <- b[!low]
  if (length(a2) > 0) {
    sm2 <- b2
    sm2[!is.na(sm2)] <- loess(b2 ~ a2, span = .1)$fitted
    sm <- c(sm1, sm2)
  } else {
    sm <- sm1
  }
  out <- tapply(sm, a, mean, na.rm = T)
  return(out)
}

smo1 <- function(x) {
  notna <- T #!is.na(x)
  a <- sp[notna]
  b <- x[notna]
  aa <- c(rev(a), a)
  bb <- c(rev(b), b)
  sm0 <- loess(bb ~ aa, span = .1)$fitted
  len <- length(a)
  sm <- sm0[(len + 1):(2 * len)]
  out <- tapply(sm, a, mean, na.rm = T)
  return(out)
}

smo <- smo1
if (grepl('shannon|simpson', args[1]))
  smo <- smo0

sm_result <- apply(data[, -1, drop = F], 2, smo)
me_result <- apply(data[, -1, drop = F], 2, medi)
new_result <- sm_result
new_result[is.na(sm_result)] <- me_result[is.na(sm_result)]

new <- cbind(unique(sp), new_result)
new[1,] <- 0
write.table(
  new,
  file = args[2],
  row = F,
  sep = "\t",
  quote = F
)
