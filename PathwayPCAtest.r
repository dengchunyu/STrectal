#' @description Pathway PCA test
#'
#' @param setenv Specific environment for pathway analysis
#' @param scCounts normalized count matrices. rownames is gene, colnames is cellnames
#' @param min.pathway.size integer Minimum number of observed genes that should be contained in a valid gene set (default=10)
#' @param max.pathway.size integer Maximum number of observed genes in a valid gene set (default=1e3)
#' @param n.randomizations numeric Number of random gene sets (of the same size) to be evaluated in parallel with each gene set (default=5). (This can be kept at 5 or 10, but should be increased to 50-100 if the significance of pathway overdispersion will be determined relative to random gene set models.)
#' @param score.alpha numeric Significance level of the confidence interval for determining upper/lower bounds (default=0.05)
#' @param cells character vector Specific cells to investigate (default=NULL)
#' @param adjusted.pvalues boolean Whether to use adjusted p-values (default=TRUE)
#' @param z.score numeric Z-score to be used as a cutoff for statistically significant patterns (default=qnorm(0.05/2, lower.tail = FALSE))
#'scCounts<-readRDS(system.file("extdata", "sample_BM1_50.rds", package="pagoda2"))
#'scCounts<-as.matrix(scCounts)
#'scCounts<-scCounts[!duplicated(rownames(scCounts)),]
#' @return pathway output

PathwayPCAtest=function(Pathway_list, scCounts,
                                n.randomizations=5,
                                n.cores=1, score.alpha=0.05,
                                adjusted.pvalues=TRUE,
                                randomization=F,
                                z.score = qnorm(0.05/2, lower.tail = FALSE)
                                   ){

   if (any(duplicated(rownames(scCounts)))) {
    stop("Duplicate gene names are not allowed - please reduce")
      }
  if (any(duplicated(colnames(scCounts)))) {
    stop("Duplicate cell names are not allowed - please reduce")
     }
  if (any(is.na(rownames(scCounts)))) {
    stop("NA gene names are not allowed - please fix")
     }
  if (any(is.na(colnames(scCounts)))) {
    stop("NA cell names are not allowed - please fix")
     }

## transpose
  scCounts <- t(scCounts)
  nPcs <- 1
  proper.gene.names <- colnames(scCounts)
    # determine valid pathways
  cm <- Matrix::colMeans(scCounts)
######calculate the pca for each go terms.
  pwpca <- lapply(Pathway_list, function(Pa_id) {
      lab <- proper.gene.names %in% Pa_id
      if (sum(lab)<=1) {
        return(NULL)
       }else{
        result = tryCatch({
        pcs <- irlba(scCounts[,lab], nv=nPcs, nu=0, center=cm[lab])

        pcs$d <- pcs$d/sqrt(nrow(scCounts))
        pcs$rotation <- pcs$v
        pcs$v <- NULL
      # get standard deviations for the random samples
      ngenes <- sum(lab)

      #randomization
      if(randomization){
        z <- do.call(
        rbind,
        lapply(seq_len(n.randomizations), function(i) {
        si <- sample(ncol(scCounts), ngenes)
        result = tryCatch({
          pcs <- irlba(scCounts[,si], nv=nPcs, nu=0, center=cm[si])$d
          },
          error = function(e) {
            return(pcs <- NULL)
             })
      })
      )
      z <- z/sqrt(nrow(scCounts))
      # local normalization of each component relative to sampled PC1 sd
      avar <- pmax(0, (pcs$d^2-mean(z[, 1]^2))/sd(z[, 1]^2))
      }else{
        z<-NA
        avar=1
      }

      if (avar>0.5) {
        # flip orientations to roughly correspond with the means
        pcs$scores <- as.matrix(t(scCounts[,lab] %*% pcs$rotation) - as.numeric((cm[lab] %*% pcs$rotation)))
        cs <- unlist(lapply(seq_len(nrow(pcs$scores)), function(i) sign(cor(pcs$scores[i,], colMeans(t(scCounts[, lab, drop = FALSE])*abs(pcs$rotation[, i]))))))
        pcs$scores <- pcs$scores*cs
        pcs$rotation <- pcs$rotation*cs
        rownames(pcs$rotation) <- colnames(scCounts)[lab]

      } # don't care not significant
      #print(paste0(id," done"))
      #return(list(xp=pcs,z=z,n=ngenes))
      return(list(xp=pcs,z=z,n=ngenes))
      },
      error = function(e) {
        return(NULL)
      })
      }
    })

  vdf <- data.frame(do.call(rbind, lapply(seq_along(pwpca), function(i) {

    result = tryCatch({
      vars <- as.numeric((pwpca[[i]]$xp$d))
      cbind(i = i, var = vars, n = pwpca[[i]]$n, npc = seq(1:ncol(pwpca[[i]]$xp$rotation)))
    },
    error = function(e) {
      return(pcs <- NULL)
    })
    })))
  n.cells <- nrow(scCounts)
  vdf$exp <- RMTstat::qWishartMax(0.5, n.cells, vdf$n, var = 1, lower.tail = FALSE)
  vdf$oe <- (vdf$var)/(vdf$exp)
  df <- data.frame(name = names(pwpca)[vdf$i],score = vdf$oe, stringsAsFactors = FALSE)
  return(df)
}


# BH P-value adjustment with a log option
bh.adjust <- function(x, log = FALSE) {
  nai <- which(!is.na(x))
  ox <- x
  x<-x[nai]
  id <- order(x, decreasing = FALSE)
  if(log) {
    q <- x[id] + log(length(x)/seq_along(x))
  } else {
    q <- x[id]*length(x)/seq_along(x)
  }
  a <- rev(cummin(rev(q)))[order(id)]
  ox[nai]<-a
  ox
}

