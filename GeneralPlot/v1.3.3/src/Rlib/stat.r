# cor.test for matrix
cor.mtest <- function (mat, method = "pearson", conf.level = 0.95) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat <- lowCI.mat <- uppCI.mat <-matrix(NA, n, n)
    diag(p.mat) <- 0
    diag(lowCI.mat) <- diag(uppCI.mat) <- 1
    for (i in 1:(n-1)) {
        for (j in (i+1):n) {
            tmp <- cor.test(mat[,i], mat[,j], method = method, conf.level = conf.level)
            p.mat[i,j] <- p.mat[j,i] <-tmp$p.value
            lowCI.mat[i,j] <- lowCI.mat[j,i] <-tmp$conf.int[1]
            uppCI.mat[i,j] <- uppCI.mat[j,i] <-tmp$conf.int[2]
        }
    }
    return(list(p.mat = p.mat, lowCI.mat = lowCI.mat, uppCI.mat = uppCI.mat))
}


