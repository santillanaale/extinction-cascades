library(igraph, quietly = TRUE)
library(bipartite, quietly = TRUE)
library(SYNCSA, quietly = TRUE)

calcMetric <- function(dat.web, ...) {
    mets <-  my.networklevel(dat.web, ...)
    mets.group <- grouplevel(dat.web, index=c("mean number of links",
                                         "weighted cluster coefficient"))
    
    mets <- c(mets, mets.group)
    ## the functional redundancy function takes a matrix of sites and
    ## species, and a trait matrix whwere the rownames of the traits
    ## match the column names of the site by species matric. In our
    ## case, the plants and pollinators are the "species"
    ## respectively, and their traits are their interaction partners.

    ## create names for later use in site x species matrix
    rownames(dat.web) <- 1:nrow(dat.web)
    colnames(dat.web) <- 1:ncol(dat.web)

    ## site by species matrix where there is only one "site" since the
    ## networks are site specific, and the columns are the
    ## species.

    ## abundance weighted site by species matrices
    plants <- matrix(rowSums(dat.web),  nrow=1)
    pols <- matrix(colSums(dat.web),  nrow=1)

    colnames(plants) <- rownames(dat.web)
    colnames(pols) <- colnames(dat.web)
    rownames(plants) <- "Plants"
    rownames(pols) <- "Pols"

    ## pull out Functional redundancy score based on: de Bello, F.;
    ## Leps, J.; Lavorel, S. & Moretti, M. (2007). Importance of
    ## species abundance for assessment of trait composition: an
    ## example based on pollinator communities. Community Ecology, 8,
    ## 163:170. and functional complementarity score based on Rao,
    ## C.R. (1982). Diversity and dissimilarity coefficients: a
    ## unified approach. Theoretical Population Biology, 21, 24:43.

    redund.plant <- unlist(rao.diversity(plants,
                                         traits=
                                             dat.web)[c("FunRao",
                                                        "FunRedundancy")])
    redund.pol <- unlist(rao.diversity(pols,
                                       traits=
                                           t(dat.web))[c("FunRao",
                                                         "FunRedundancy")])

    return(c(mets,
             redund.plant,
             redund.pol))

}


## function to simulate 1 null, and calculate statistics on it
calcNullStat <- function(dat.web,
                         null.fun,...) {
    sim.web <- null.fun(dat.web)
    out.met <- calcMetric(sim.web,...)
    return(out.met)
}

##  function that computes summary statistics on simulated null matrices
##  (nulls simulated from web N times)
calcNetworkMetrics <- function(dat.web, N,
                               index= c("niche overlap",
                                        "functional complementarity",
                                        "weighted NODF",
                                        "number of species",
                                        "H2")) {
    ## calculate pvalues
    pvals <- function(stats, nnull){
        rowSums(stats >= stats[, rep(1, ncol(stats))])/(nnull + 1)
    }
    ## calculate zvalues two different ways
    zvals <-function(stats){
        z.sd <- (stats[,1] -
                 apply(stats, 1, mean, na.rm = TRUE))/
            apply(stats, 1, sd, na.rm = TRUE)
        z.sd[is.infinite(z.sd)] <- NA
        return(z.sd)
    }
    ## check that matrix is proper format (no empty row/col and no NAs)
    ## drop empty rows and columns
    dat.web <- as.matrix(bipartite::empty(dat.web))
    ## check to make sure emptied matrix is large enough
    ## to calculate statistics on
    if(all(dim(dat.web) >= 2)) {
        ## calculate null metrics
        null.stat <- replicate(N,
                               calcNullStat(dat.web,
                                            null.fun= vaznull.fast,
                                            index=index,
                                            dist="chao"),
                               simplify=TRUE)
        ## calculate metrics from data
        true.stat <- calcMetric(dat.web,
                                index=index)
        out.mets <- cbind(true.stat, null.stat)
        ## compute z scores
        zvalues <- zvals(out.mets)
        names(zvalues) <- paste("z", names(true.stat), sep="")
        ## compute p-values
        pvalues <- pvals(out.mets, N)
        names(pvalues) <- paste("p", names(true.stat), sep="")
        out <- c(true.stat, zvalues, pvalues)
        return(out)


    }
    return(rep(NA, (length(index) + 6)*3))
}

prepDat <- function(cor.stats, spec.dat,
                    cols.to.keep= c("Site"),
                    net.type){
  dats <- do.call(rbind, cor.stats)
  out <- data.frame(dats)
  out$Site <-  names(cor.stats)
  # out$Year <-  as.factor(sapply(strsplit(names(cor.stats), "\\."),
  #                               function(x) x[2]))
  
  # if(net.type == "YearSR"){
  #     out$SampleRound <-  as.factor(sapply(strsplit(names(cor.stats), "\\."),
  #                                          function(x) x[3]))
  #     temporal.dat <- unique(spec.dat[, c("Site", "Year", "Date", "Doy")])
  #     out <- merge(out, temporal.dat)
  # }
  # site.dats <- unique(spec[, c(cols.to.keep)])
  # site.dats$Site  <- as.character(site.dats$Site)
  # site.dats <- apply(site.dats[, cols.to.keep[-1]], 2, function(x)
  #     tapply(as.numeric(x), site.dats$Site, mean, na.rm=TRUE))
  # site.dats <- as.data.frame(site.dats)
  # site.dats$Site <- rownames(site.dats)
  # out <- merge(out, site.dats, by="Site", all.x=TRUE)
  # 
  rownames(out) <- NULL
  return(out)
}


my.networklevel <- function(web, index = "ALLBUTDD", level = "both", weighted = TRUE,
    ISAmethod = "Bluethgen", SAmethod = "Bluethgen", extinctmethod = "r",
    nrep = 100, CCfun = median, dist = "horn", normalise = TRUE,
    empty.web = TRUE, logbase = "e", intereven = "prod", H2_integer = TRUE,
    fcweighted = TRUE, fcdist = "euclidean", legacy = FALSE)
{
    if (empty.web) {
        web <- empty(web)
    }
    web.e <- empty(web)
    if (NROW(web) < 2 | NCOL(web) < 2)
        warning("Web is really too small to calculate any reasonable index. You will get the values nonetheless, but I wouldn't put any faith in them!")
    allindex <- c("number of species", "connectance", "web asymmetry",
        "links per species", "number of compartments", "compartment diversity",
        "cluster coefficient", "degree distribution", "mean number of shared partners",
        "togetherness", "C score", "V ratio", "discrepancy",
        "nestedness", "NODF", "weighted nestedness", "ISA", "SA",
        "extinction slope", "robustness", "niche overlap", "weighted cluster coefficient",
        "weighted NODF", "partner diversity", "generality", "vulnerability",
        "linkage density", "weighted connectance", "Fisher alpha",
        "interaction evenness", "Alatalo interaction evenness",
        "Shannon diversity", "functional complementarity", "H2")
    index <- unique(index)
    wrong.name <- which(is.na(pmatch(index, c(allindex, "ALL",
        "ALLBUTDD", "info", "quantitative", "binary", "topology",
        "networklevel"))))
    if (length(wrong.name) > 0)
        stop("You selected an index that is not available: ",
            paste(index[wrong.name], collapse = ", "))
    if (length(index) == 1 & !all(index %in% allindex)) {
        index <- switch(index, ALL = allindex, ALLBUTDD = allindex[-which(allindex ==
            "degree distribution")], info = c("number of species",
            "connectance", "web asymmetry", "links per species",
            "number of compartments"), quantitative = c("weighted cluster coefficient",
            "weighted nestedness", "weighted NODF", "functional complementarity",
            "partner diversity", "effective partners", "H2",
            "diversity", "linkage density", "weighted connectance",
            "niche overlap"), binary = c("connectance", "links per species",
            "nestedness", "mean number of partners", "cluster coefficient",
            "C-score", "Fisher alpha"), topology = c("connectance",
            "cluster coefficient", "degree distribution", "togetherness",
            "nestedness", "NODF"), networklevel = c("connectance",
            "web asymmetry", "links per species", "number of compartments",
            "compartment diversity", "cluster coefficient", "nestedness",
            "NODF", "weighted NODF", "ISA", "SA", "linkage density",
            "Fisher alpha", "diversity", "interaction evenness",
            "Alatalo interaction evenness", "H2"), stop("Your index is not recognised! Typo? Check help for options!",
            call. = FALSE))
    }
    if (legacy == FALSE) {
        out <- list()
        if ("connectance" %in% index) {
            suppressWarnings(out$connectance <- sum(web > 0)/prod(dim(web)))
        }
        if ("web asymmetry" %in% index)
            out$"web asymmetry" <- (NCOL(web) - NROW(web))/sum(dim(web))
        if ("links per species" %in% index) {
            L <- sum(web > 0)/sum(dim(web))
            out$"links per species" = L
        }
        if (any(c("number of compartments", "compartment diversity") %in%
            index)) {
            CD <- function(co) {
                if (co$n.compart > 1) {
                  no <- NA
                  for (i in 1:co$n.compart) {
                    comp <- which(abs(co$cweb) == i, arr.ind = TRUE)
                    no[i] <- length(unique(comp[, 1])) + length(unique(comp[,
                      2]))
                  }
                  no <- no/sum(dim(web))
                  CD <- exp(-sum(no * log(no)))
                }
                else {
                  CD <- NA
                }
                CD
            }
            comps <- try(compart(web.e), silent = TRUE)
            if (inherits(comps, "try-error")) {
                ncompart <- compdiv <- NA
            }
            else {
                ncompart <- comps$n.compart
                compdiv <- CD(comps)
            }
            if ("number of compartments" %in% index)
                out$"number of compartments" <- as.integer(ncompart)
            if ("compartment diversity" %in% index)
                out$"compartment diversity" <- compdiv
        }
        if ("cluster coefficient" %in% index) {
            cluster.coef <- function(web, full = FALSE, FUN = mean) {
                web <- as.matrix(web)
                Ci.high <- colSums((web > 0))/nrow(web)
                Ci.low <- rowSums((web > 0))/ncol(web)
                CC <- FUN(Ci.high)
                if (full)
                  out <- list(`cluster coefficient` = CC, `CC values higher` = Ci.high,
                    `CC values lower` = Ci.low)
                else out <- c(`cluster coefficient` = CC)
                out
            }
            out$"cluster coefficient" = as.numeric(cluster.coef(web,
                FUN = CCfun, full = FALSE))
        }
        if ("nestedness" %in% index) {
            nest <- try(nestedtemp(web)$statistic, silent = TRUE)
            out$nestedness <- ifelse(inherits(nest, "try-error"),
                NA, nest)
        }
        if ("NODF" %in% index) {
            NODF <- try(unname(nestednodf(web, order = TRUE,
                weighted = FALSE)$statistic[3]), silent = TRUE)
            out$NODF <- if (inherits(NODF, "try-error"))
                NA
            else NODF
        }
        if ("weighted nestedness" %in% index) {
            wine.res <- try(wine(web.e, nreps = nrep)$wine, silent = TRUE)
            out$"weighted nestedness" <- if (!inherits(wine.res,
                "try-error")) {
                wine.res
            }
            else {
                NA
            }
        }
        if ("weighted NODF" %in% index) {
            wNODF <- try(unname(nestednodf(web, order = TRUE,
                weighted = TRUE)$statistic[3]), silent = TRUE)
            out$"weighted NODF" <- if (inherits(wNODF, "try-error"))
                NA
            else wNODF
        }
        if (any(c("ISA", "interaction strength asymmetry", "dependence asymmetry") %in%
            index)) {
            depL <- web.e/matrix(rowSums(web.e), nrow = NROW(web.e),
                ncol = NCOL(web.e), byrow = FALSE)
            depH <- web.e/matrix(colSums(web.e), nrow = NROW(web.e),
                ncol = NCOL(web.e), byrow = TRUE)
            if (ISAmethod == "Bascompte" & "ISA" %in% index) {
                out$"dependence asymmetry" = mean(abs(depL -
                  depH)/pmax(depL, depH), na.rm = TRUE)
            }
            if (ISAmethod == "Bluethgen" & "ISA" %in% index) {
                web2 <- web
                web2[, which(colSums(web) == 1)] <- 0
                web2[which(rowSums(web) == 1), ] <- 0
                rowsummat <- matrix(rowSums(web2), nrow = NROW(web2),
                  ncol = NCOL(web2), byrow = FALSE)
                colsummat <- matrix(colSums(web2), nrow = NROW(web2),
                  ncol = NCOL(web2), byrow = TRUE)
                depL <- web2/rowsummat
                depH <- web2/colsummat
                depL[depL <= 0] <- NA
                depH[depH <= 0] <- NA
                depLprime <- (depL - 1/rowsummat)/(1 - 1/rowsummat)
                depHprime <- (depH - 1/colsummat)/(1 - 1/colsummat)
                out$"interaction strength asymmetry" = mean(as.matrix(depHprime -
                  depLprime), na.rm = TRUE)
            }
        }
        if ("SA" %in% index) {
            di <- dfun(web)$dprime
            dj <- dfun(t(web))$dprime
            if (SAmethod == "log") {
                lgmeani <- mean(log(di[di > 0]))
                lgmeanj <- mean(log(dj[dj > 0]))
                SA <- (lgmeanj - lgmeani)/sum(lgmeani, lgmeanj)
            }
            if (SAmethod == "Bluethgen") {
                wmeani <- sum(di * rowSums(web.e))/sum(web.e)
                wmeanj <- sum(dj * colSums(web.e))/sum(web.e)
                SA <- (wmeanj - wmeani)/sum(wmeani, wmeanj)
            }
            out$"specialisation asymmetry" <- SA
        }
        if (any(c("linkage density", "weighted connectance") %in%
            index)) {
            preytot.mat <- matrix(rep(colSums(web), NROW(web)),
                NROW(web), byrow = TRUE)
            preyprop.mat <- web/preytot.mat
            predtot.mat <- matrix(rep(rowSums(web), NCOL(web)),
                NROW(web), byrow = FALSE)
            predprop.mat <- web/predtot.mat
            if (logbase == 2 | logbase == "2") {
                H_Nk <- apply(preyprop.mat, 2, function(x) -sum(x *
                  log2(x), na.rm = TRUE))
                H_Pk <- apply(predprop.mat, 1, function(x) -sum(x *
                  log2(x), na.rm = TRUE))
                n_Nk <- ifelse(colSums(web) != 0, 2^H_Nk, 0)
                n_Pk <- ifelse(rowSums(web) != 0, 2^H_Pk, 0)
            }
            if (logbase == "e") {
                H_Nk <- apply(preyprop.mat, 2, function(x) -sum(x *
                  log(x), na.rm = TRUE))
                H_Pk <- apply(predprop.mat, 1, function(x) -sum(x *
                  log(x), na.rm = TRUE))
                n_Nk <- ifelse(colSums(web) != 0, exp(H_Nk),
                  0)
                n_Pk <- ifelse(rowSums(web) != 0, exp(H_Pk),
                  0)
            }
            V <- sum(rowSums(web)/sum(web) * n_Pk)
            G <- sum(colSums(web)/sum(web) * n_Nk)
            LD_q <- 0.5 * (V + G)
            if ("linkage density" %in% index)
                out$"linkage density" <- LD_q
            if ("weighted connectance" %in% index)
                out$"weighted connectance" <- LD_q/sum(dim(web))
        }
        if ("Fisher alpha" %in% index) {
            fish <- try(fisherfit(web)$estimate, silent = TRUE)
            if (inherits(fish, "try-error")) {
                out$"Fisher alpha" <- NA
            }
            else {
                out$"Fisher alpha" <- fish
            }
        }
        if (any(c("interaction evenness", "Alatalo interaction evenness",
            "Shannon diversity") %in% index)) {
            p_i.mat <- web/sum(web)
            SH <- -sum(p_i.mat * log(p_i.mat), na.rm = TRUE)
            if ("Shannon diversity" %in% index)
                out$"Shannon diversity" <- SH
            IE <- ifelse(intereven == "prod", SH/log(prod(dim(web))),
                SH/log(sum(web > 0)))
            if ("interaction evenness" %in% index)
                out$"interaction evenness" <- IE
            if ("Alatalo interaction evenness" %in% index) {
                evenness <- function(web) {
                  pk <- web/sum(web)
                  (Alatalo <- (1/sum(pk^2) - 1)/(exp(-sum(pk *
                    log(pk), na.rm = TRUE)) - 1))
                }
                E <- evenness(web)
                out$"Alatalo interaction evenness" <- E
            }
        }
        if ("H2" %in% index) {
            is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x -
                round(x)) < tol
            if (any(is.wholenumber(web) == FALSE))
                H2_integer <- FALSE
            H2 <- as.numeric(my.H2fun(web, H2_integer = H2_integer)[1])
            out$H2 = ifelse(H2 < 0, 0, H2)
        }
        netw.index <- match(c("connectance", "web asymmetry",
            "links per species", "number of compartments", "compartment diversity",
            "nestedness", "NODF", "weighted nestedness", "weighted NODF",
            "ISA", "SA", "interaction evenness", "Alatalo interaction evenness",
            "Fisher alpha", "H2", "Shannon diversity", "linkage density",
            "weighted connectance"), index)
        exclude.index <- netw.index[!is.na(netw.index)]
        gindex <- if (length(exclude.index) == 0)
            index
        else index[-exclude.index]
        if (length(gindex) > 0)
            outg <- grouplevel(web, index = gindex, level = level,
                weighted = weighted, extinctmethod = extinctmethod,
                nrep = nrep, CCfun = CCfun, dist = dist, normalise = normalise,
                empty.web = empty.web, logbase = logbase, fcweighted = fcweighted,
                fcdist = fcdist)
        if (exists("outg")) {
            if (is.list(outg)) {
                SEQ <- seq(1, 2 * length(outg[[1]]), by = 2)
                sorted.outg <- c(outg[[1]], outg[[2]])
                outg <- sorted.outg[order(c(SEQ, SEQ + 1))]
            }
            out <- c(unlist(out), outg)
        }
        else {
            out <- unlist(out)
        }
    }
    if (legacy == TRUE) {
        out <- .networklevel(web, index = index, ISAmethod = ISAmethod,
            SAmethod = SAmethod, extinctmethod = extinctmethod,
            nrep = nrep, plot.it.extinction = FALSE, plot.it.dd = FALSE,
            CCfun = CCfun, dist = dist, normalise = normalise,
            empty.web = empty.web, logbase = logbase, fcweighted = fcweighted,
            fcdist = fcdist)
    }
    return(out)
}

my.H2fun <- function (web, H2_integer = TRUE)
{
    if (H2_integer & any((web%%1) != 0))
        stop("web does not contain integers! maybe you should set H2_integer to FALSE")
    tot <- sum(web)
    rs <- rowSums(web)
    cs <- colSums(web)
    H2uncorr = -sum(web/tot * log(web/tot), na.rm = TRUE)
    exexpec <- outer(rs, cs/tot)
    if (!H2_integer) {
        newweb <- exexpec
    }
    else {
        expec <- matrix(0, nrow(web), ncol(web))
        difexp <- exexpec - expec
        newweb <- floor(exexpec)
        webfull <- matrix("no", nrow(web), ncol(web))
        while (sum(newweb) < tot) {
            webfull[which(rowSums(newweb) == rs), ] <- "yo"
            webfull[, which(colSums(newweb) == cs)] <- "yo"
            OK <- webfull == "no"
            smallestpossible <- newweb == min(newweb[OK])
            greatestdif <- max(difexp[smallestpossible & OK])
            bestone <- which(OK & smallestpossible & difexp ==
                greatestdif)
            if (length(bestone) > 1)
                bestone <- sample(bestone, 1)
            newweb[bestone] <- newweb[bestone] + 1
            difexp <- exexpec - newweb
        }
        H2_max <- -sum(newweb/tot * log(newweb/tot), na.rm =
                                                             TRUE)
        if (max(exexpec) > 0.3679 * tot) {
            for (tries in 1:500) {
                newmx <- newweb
                difexp <- exexpec - newmx
                greatestdif <- difexp == min(difexp)
                if (length(which(greatestdif)) > 1) {
                  largestvalue = newmx == max(newmx[greatestdif])
                  first <- greatestdif & largestvalue
                }
                else {
                  first = greatestdif
                }
                newmx[first][1] <- newmx[first][1] - 1
                throw = which(rowSums(first) > 0)[1]
                thcol = which(colSums(first) > 0)[1]
                mr = max(difexp[throw, ])
                mc = max(difexp[, thcol])
                if (mr >= mc) {
                  scnd = which(difexp[throw, ] == mr)[1]
                  newmx[throw, scnd] = newmx[throw, scnd] + 1
                  thrd = which(difexp[, scnd] == min(difexp[,
                    scnd]))[1]
                  newmx[thrd, scnd] = newmx[thrd, scnd] - 1
                  newmx[thrd, thcol] = newmx[thrd, thcol] + 1
                }
                else {
                  scnd = which(difexp[, thcol] == mc)[1]
                  newmx[scnd, thcol] = newmx[scnd, thcol] + 1
                  thrd = which(difexp[scnd, ] == min(difexp[scnd,
                    ]))[1]
                  newmx[scnd, thrd] = newmx[scnd, thrd] - 1
                  newmx[throw, thrd] = newmx[throw, thrd] + 1
                }
            }
            newweb <- newmx
        }
    }
    H2_max.improved <- -sum(newweb/tot * log(newweb/tot), na.rm =
                                                              TRUE)
    if(!exists("H2_max")) H2_max <- H2_max.improved
    H2_max <- ifelse(H2_max >= H2_max.improved, H2_max,
                         H2_max.improved)
    newweb <- matrix(0, length(rs), length(cs))
    rsrest = rs
    csrest = cs
    while (round(sum(rsrest), 10) != 0) {
        newweb[which(rsrest == max(rsrest))[1], which(csrest ==
            max(csrest))[1]] = min(c(max(rsrest), max(csrest)))
        rsrest = rs - rowSums(newweb)
        csrest = cs - colSums(newweb)
    }
    Pnew <- newweb/sum(newweb)
    H2_min <- -sum(Pnew * log(Pnew), na.rm = TRUE)
    if (H2uncorr < H2_min)
        H2_min <- H2uncorr
    if (H2_max < H2uncorr)
        H2_max <- H2uncorr
    H_2prime <- (H2_max - H2uncorr)/(H2_max - H2_min)
    c(H2 = H_2prime, H2min = round(H2_min, 3), H2max = round(H2_max,
        3), H2uncorr = round(H2uncorr, 3))
}
