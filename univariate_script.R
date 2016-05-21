univariateF <- function(datMN,
                        samDF,
                        varDF,
                        facC,
                        tesC = c("ttest", "wilcoxon", "anova", "kruskal", "pearson", "spearman")[1],
                        adjC = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")[7],
                        thrN = 0.05) {


    ## Option
    ##---------------

    strAsFacL <- options()$stringsAsFactors
    options(stingsAsFactors = FALSE)
    options(warn = -1)

    if(mode(samDF[, facC]) == "character") {
        facFcVn <- factor(samDF[, facC])
        facLevVc <- levels(facFcVn)
    } else
        facFcVn <- samDF[, facC]

    cat("\nPerforming '", tesC, "'\n", sep="")

    varPfxC <- paste0(make.names(facC), "_", tesC, "_")

    if(tesC %in% c("ttest", "wilcoxon", "pearson", "spearman")) {

        switch(tesC,
               ttest = {
                   staF <- function(y) diff(tapply(y, facFcVn, function(x) mean(x, na.rm = TRUE)))
                   tesF <- function(y) t.test(y ~ facFcVn)[["p.value"]]
               },
               wilcoxon = {
                   staF <- function(y) diff(tapply(y, facFcVn, function(x) median(x, na.rm = TRUE)))
                   tesF <- function(y) wilcox.test(y ~ facFcVn)[["p.value"]]
               },
               pearson = {
                   staF <- function(y) cor(facFcVn, y, method = "pearson", use = "pairwise.complete.obs")
                   tesF <- function(y) cor.test(facFcVn, y, method = "pearson", use = "pairwise.complete.obs")[["p.value"]]
               },
               spearman = {
                   staF <- function(y) cor(facFcVn, y, method = "spearman", use = "pairwise.complete.obs")
                   tesF <- function(y) cor.test(facFcVn, y, method = "spearman", use = "pairwise.complete.obs")[["p.value"]]
               })

        staVn <- apply(datMN, 2, staF)

        fdrVn <- p.adjust(apply(datMN,
                                2,
                                tesF),
                          method = adjC)

        sigVn <- as.numeric(fdrVn < thrN)

        if(tesC %in% c("ttest", "wilcoxon"))
            varPfxC <- paste0(varPfxC, paste(rev(facLevVc), collapse = "-"), "_")

        varDF[, paste0(varPfxC, ifelse(tesC %in% c("ttest", "wilcoxon"), "dif", "cor"))] <- staVn

        varDF[, paste0(varPfxC, adjC)] <- fdrVn

        varDF[, paste0(varPfxC, "sig")] <- sigVn

    } else if(tesC == "anova") {

        ## getting the names of the pairwise comparisons 'class1Vclass2'
        prwVc <- rownames(TukeyHSD(aov(datMN[, 1] ~ facFcVn))[["facFcVn"]])

        aovMN <- t(apply(datMN, 2, function(varVn) {

            aovMod <- aov(varVn ~ facFcVn)
            pvaN <- summary(aovMod)[[1]][1, "Pr(>F)"]
            hsdMN <- TukeyHSD(aovMod)[["facFcVn"]]
            c(pvaN, c(hsdMN[, c("diff", "p adj")]), as.numeric(hsdMN[, "p adj"] < thrN))

        }))

        aovMN[, 1] <- p.adjust(aovMN[, 1], method = adjC)
        sigVn <-  as.numeric(aovMN[, 1] < thrN)
        aovMN <- cbind(aovMN[, 1], sigVn, aovMN[, 2:ncol(aovMN)])
        ## aovMN[which(aovMN[, 2] < 1), (3 + length(prwVc)):ncol(aovMN)] <- NA
        colnames(aovMN) <- paste0(varPfxC,
                                  c(adjC,
                                    "sig",
                                    paste0(prwVc, "_dif"),
                                    paste0(prwVc, "_pva"),
                                    paste0(prwVc, "_sig")))
        aovMN[which(aovMN[, paste0(varPfxC, "sig")] < 1), paste0(varPfxC, c(paste0(prwVc, "_pva"), paste0(prwVc, "_sig")))] <- NA

        varDF <- cbind.data.frame(varDF, as.data.frame(aovMN))

    } else if(tesC == "kruskal") {

        ## getting the names of the pairwise comparisons 'class1Vclass2'
        nemMN <- posthoc.kruskal.nemenyi.test(datMN[, 1], facFcVn, "Tukey")[["p.value"]]
        nemVl <- c(lower.tri(nemMN, diag = TRUE))
        nemClaMC <- cbind(rownames(nemMN)[c(row(nemMN))][nemVl],
                          colnames(nemMN)[c(col(nemMN))][nemVl])
        nemNamVc <- paste0(nemClaMC[, 1], "-", nemClaMC[, 2])
        nemNamVc <- paste0(varPfxC, nemNamVc)

        nemMN <- t(apply(datMN, 2, function(varVn) {

            pvaN <- kruskal.test(varVn ~ facFcVn)[["p.value"]]
            varNemMN <- posthoc.kruskal.nemenyi.test(varVn, facFcVn, "Tukey")[["p.value"]]
            c(pvaN, c(varNemMN))

        }))
        pvaVn <- nemMN[, 1]
        fdrVn <- p.adjust(pvaVn, method = adjC)
        sigVn <- as.numeric(fdrVn < thrN)
        nemMN <- nemMN[, c(FALSE, nemVl)]
        colnames(nemMN) <- paste0(nemNamVc, "_pva")
        nemSigMN <- nemMN < thrN
        mode(nemSigMN) <- "numeric"
        colnames(nemSigMN) <- paste0(nemNamVc, "_sig")
        nemMN[sigVn < 1, ] <- NA
        nemSigMN[sigVn < 1, ] <- NA

        difMN <- sapply(1:nrow(nemClaMC), function(prwI) {
            prwVc <- nemClaMC[prwI, ]
            prwVi <- which(facFcVn %in% prwVc)
            prwFacFc <- factor(as.character(facFcVn)[prwVi], levels = prwVc)
            apply(datMN[prwVi, ], 2, function(varVn) -diff(as.numeric(tapply(varVn, prwFacFc, function(x) median(x, na.rm = TRUE)))))
        })
        colnames(difMN) <- gsub("_sig", "_dif", colnames(nemSigMN))

        nemMN <- cbind(fdrVn, sigVn, difMN, nemMN, nemSigMN)
        colnames(nemMN)[1:2] <- paste0(varPfxC, c(adjC, "sig"))

        varDF <- cbind.data.frame(varDF, as.data.frame(nemMN))

    }

    names(sigVn) <- rownames(varDF)
    sigSumN <- sum(sigVn, na.rm = TRUE)
    if(sigSumN) {
        cat("\nThe following ", sigSumN, " variable", ifelse(sigSumN > 1, "s", ""), " (", round(sigSumN / length(sigVn) * 100), "%) ", ifelse(sigSumN > 1, "were", "was"), " found significant at the ", thrN, " level:\n", sep = "")
        cat(paste(rownames(varDF)[sigVn > 0], collapse = "\n"), "\n", sep = "")
    } else
        cat("\nNo significant variable found at the selected ", thrN, " level\n", sep = "")

    options(stingsAsFactors = strAsFacL)

    return(varDF)

}
