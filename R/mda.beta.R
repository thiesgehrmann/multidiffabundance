###############################################################################
# Run BETA
# Code from Vegan is modified to allow scope to be passed to the internal anova,
#   such that permutations are only performed for the selected variable

`mda.adonis2` <- function(formula, data, permutations = 999, method = "bray",
                          sqrt.dist = FALSE, add = FALSE, by = "terms",
                          parallel = getOption("mc.cores"), na.action = na.fail,
                          strata = NULL, scope = NULL, ...) {
    ## handle missing data
    if (missing(data))
        data <- model.frame(delete.response(terms(formula)),
                            na.action = na.action)
    ## we accept only by = "terms", "margin" or NULL
    if (!is.null(by))
        by <- match.arg(by, c("terms", "margin", "onedf"))
    ## evaluate lhs
    YVAR <- formula[[2]]
    lhs <- eval(YVAR, environment(formula), globalenv())
    environment(formula) <- environment()
    ## Take care that input lhs are dissimilarities
    if ((is.matrix(lhs) || is.data.frame(lhs)) &&
        isSymmetric(unname(as.matrix(lhs))))
        lhs <- as.dist(lhs)
    if (!inherits(lhs, "dist"))
        lhs <- vegdist(as.matrix(lhs), method=method, ...)
    ## adjust distances if requested
    if (sqrt.dist)
        lhs <- sqrt(lhs)
    if (is.logical(add) && isTRUE(add))
        add <- "lingoes"
    if (is.character(add)) {
        add <- match.arg(add, c("lingoes", "cailliez"))
        if (add == "lingoes") {
            ac <- addLingoes(as.matrix(lhs))
            lhs <- sqrt(lhs^2 + 2 * ac)
        }
        else if (add == "cailliez") {
            ac <- addCailliez(as.matrix(lhs))
            lhs <- lhs + ac
        }
    }
    ## adonis0 & anova.cca should see only dissimilarities (lhs)
    if (!missing(data)) # expand and check terms
        formula <- terms(formula, data=data)
    if (is.null(attr(data, "terms"))) # not yet a model.frame?
        data <- model.frame(delete.response(terms(formula)), data,
                            na.action = na.action)
    formula <- update(formula, lhs ~ .)
    sol <- vegan:::adonis0(formula, data = data, method = method)
    ## handle permutations
    perm <- vegan:::getPermuteMatrix(permutations, NROW(data), strata = strata)
    out <- mda.anova.cca(sol, permutations = perm, by = by,
                 parallel = parallel, scope=scope)
    ## attributes will be lost when adding a new column
    att <- attributes(out)
    ## add traditional adonis output on R2
    out <- rbind(out, "Total" = c(nobs(sol)-1, sol$tot.chi, NA, NA))
    out <- cbind(out[,1:2], "R2" = out[,2]/sol$tot.chi, out[,3:4])
    ## Fix output header to show the adonis2() call instead of adonis0()
    att$heading[2] <- deparse(match.call(), width.cutoff = 500L)
    att$names <- names(out)
    att$row.names <- rownames(out)
    attributes(out) <- att
    out
}

`mda.anova.cca` <- function(object, ..., permutations = how(nperm=999), by = NULL,
                            model = c("reduced", "direct", "full"),
                            parallel = getOption("mc.cores"), strata = NULL,
                            cutoff = 1, scope = NULL) {
    EPS <- sqrt(.Machine$double.eps) # for permutation P-values
    model <- match.arg(model)
    ## permutation matrix
    N <- nrow(object$CA$u)
    permutations <- vegan:::getPermuteMatrix(permutations, N, strata = strata)
    seed <- attr(permutations, "seed")
    control <- attr(permutations, "control")
    ## see if this was a list of ordination objects
    dotargs <- list(...)
    ## we do not want to give dotargs to anova.ccalist, but we
    ## evaluate 'parallel' and 'model' here
    if (length(dotargs)) {
        isCCA <- sapply(dotargs, function(z) inherits(z, "cca"))
        if (any(isCCA)) {
            dotargs <- dotargs[isCCA]
            object <- c(list(object), dotargs)
            sol <-
                anova.ccalist(object,
                              permutations = permutations,
                              model = model,
                              parallel = parallel)
            attr(sol, "Random.seed") <- seed
            attr(sol, "control") <- control
            return(sol)
        }
    }
    ## We only have a single model: check if it is empty
    if (is.null(object$CA) || is.null(object$CCA) ||
        object$CCA$rank == 0 || object$CA$rank == 0)
        return(anova.ccanull(object))
    ## by cases
    if (!is.null(by)) {
        by <- match.arg(by, c("terms", "margin", "axis", "onedf"))
        if (is.null(object$terms))
            stop("model must be fitted with formula interface")
        sol <- switch(by,
                      "terms" = vegan:::anova.ccabyterm(object,
                      permutations = permutations,
                      model = model, parallel = parallel),
                      "margin" = mda.anova.ccabymargin(object,
                      permutations = permutations,
                      model = model, parallel = parallel,
                      scope = scope),
                      "axis" = vegan:::anova.ccabyaxis(object,
                      permutations = permutations,
                      model = model, parallel = parallel,
                      cutoff = cutoff),
                      "onedf" = vegan:::anova.ccaby1df(object,
                       permutations = permutations,
                       model = model, parallel = parallel)
                      )
        attr(sol, "Random.seed") <- seed
        attr(sol, "control") <- control
        return(sol)
    }
    ## basic overall test: pass other arguments except 'strata'
    ## because 'permutations' already is a permutationMatrix
    tst <- permutest.cca(object, permutations = permutations,
                         model = model, parallel = parallel, ...)
    Fval <- c(tst$F.0, NA)
    Pval <- (sum(tst$F.perm >= tst$F.0 - EPS) + 1)/(tst$nperm + 1)
    Pval <- c(Pval, NA)
    table <- data.frame(tst$df, tst$chi, Fval, Pval)
    if (inherits(object, c("capscale", "dbrda")) && object$adjust == 1)
        varname <- "SumOfSqs"
    else if (inherits(object, "rda"))
        varname <- "Variance"
    else
        varname <- "ChiSquare"
    colnames(table) <- c("Df", varname, "F", "Pr(>F)")
    head <- paste0("Permutation test for ", tst$method, " under ",
                  tst$model, " model\n", vegan:::howHead(control))
    mod <- paste("Model:", c(object$call))
    structure(table, heading = c(head, mod), Random.seed = seed,
              control = control, F.perm = tst$F.perm,
              class = c("anova.cca", "anova", "data.frame"))
}

`mda.anova.ccabymargin` <- function(object, permutations, scope, ...){
    EPS <- sqrt(.Machine$double.eps)
    nperm <- nrow(permutations)
    ## We need term labels but without Condition() terms
    if (!is.null(scope) && is.character(scope))
        trms <- scope
    else
        trms <- drop.scope(object)
    ## Condition() not considered marginal
    alltrms <- intersect(attr(terms(object$terminfo), "term.labels"),
                         attr(terms(object), "term.labels"))
    trmlab <- intersect(alltrms, trms)
    if (length(trmlab) == 0)
        stop("the scope was empty: no available marginal terms")
    ## baseline: all terms
    big <- permutest(object, permutations, ...)
    dfbig <- big$df[2]
    chibig <- big$chi[2]
    scale <- big$den/dfbig
    ## Collect all marginal models. This differs from old version
    ## (vegan 2.0) where other but 'nm' were partialled out within
    ## Condition(). Now we only fit the model without 'nm' and compare
    ## the difference against the complete model.
    Y <- ordiYbar(object, "init")
    X <- model.matrix(object)
    ## we must have Constraints to get here, but we may also have
    ## Conditions
    if (!is.null(object$pCCA)) {
        Z <- X$Conditions
        X <- X$Constraints
    } else {
        Z <- NULL
    }
    ass <- object$terminfo$assign
    if (is.null(ass))
        stop("old style result object: update() your model")
    ## analyse only terms of 'ass' thar are in scope
    scopeterms <- which(alltrms %in% trmlab)
    mods <- lapply(scopeterms, function(i, ...)
           permutest(vegan:::ordConstrained(Y, X[, ass != i, drop=FALSE], Z, "pass"),
                     permutations, ...), ...)
    ## Chande in df
    Df <- sapply(mods, function(x) x$df[2]) - dfbig
    ## F of change
    Chisq <- sapply(mods, function(x) x$chi[2]) - chibig
    Fstat <- (Chisq/Df)/(chibig/dfbig)
    ## Simulated F-values
    Fval <- sapply(mods, function(x) x$num)
    ## Had we an empty model we need to clone the denominator
    if (length(Fval) == 1)
        Fval <- matrix(Fval, nrow = nperm)
    Fval <- sweep(-Fval, 1, big$num, "+")
    Fval <- sweep(Fval, 2, Df, "/")
    Fval <- sweep(Fval, 1, scale, "/")
    ## Simulated P-values
    Pval <- (colSums(sweep(Fval, 2, Fstat - EPS, ">=")) + 1)/(nperm + 1)
    ## Collect results to anova data.frame
    out <- data.frame(c(Df, dfbig), c(Chisq, chibig),
                      c(Fstat, NA), c(Pval, NA))
    if (inherits(object, c("capscale", "dbrda")) && object$adjust == 1)
        varname <- "SumOfSqs"
    else if (inherits(object, "rda"))
        varname <- "Variance"
    else
        varname <- "ChiSquare"
    dimnames(out) <- list(c(trmlab, "Residual"),
                          c("Df", varname, "F", "Pr(>F)"))
    head <- paste0("Permutation test for ", object$method, " under ",
                   big$model, " model\n",
                   "Marginal effects of terms\n",
                   vegan:::howHead(attr(permutations, "control")))
    mod <- paste("Model:", c(object$call))
    attr(out, "heading") <- c(head, mod)
    attr(out, "F.perm") <- Fval
    class(out) <- c("anova.cca", "anova", "data.frame")
    out
}
                   
###############################################################################
                   
mda.beta <- function(mda.D, beta.permutations=999, beta.parallel=8, beta.hack=TRUE, ...){
    
    suppressPackageStartupMessages(require(vegan))
    suppressPackageStartupMessages(require(tidyverse))
    D <- mda.D

    do <- function(f_idx){
        f <- D$formula$norand[[f_idx]]
        mainvar <- D$formula$main_var[f_idx]
        
        if ( length(D$formula$rand_slope[[f_idx]]) > 0 ){
            message(paste0(c("[MDA] mda.beta. Formula on ", f_idx, " contains random slope effects. Adonis2 can not handle random slopes. Run will continue without random slopes")))
        }
        
        block <- NULL
        if ( length(D$formula$rand_intercept[[f_idx]]) > 0 ){
            rblock <- D$formula$rand_intercept[[f_idx]][1]
            block <- trimws(gsub(")", "", unlist(strsplit(rblock, split="\\|"))[[2]]))
            message(paste0(c("[MDA] mda.beta: Formula on ", f_idx, " contains a random intercept effect. Adonis2 cannot handle random effects. However, we will translate this into a block design. This means that permutations will be performed within these blocks. If you do not desire this, replace the random intercept with a fixed effect. Will block with (1|", block, ").")))
            if (length(D$formula$rand_intercept[[f_idx]]) > 1){
                message(paste0(c("[MDA] mda.beta: Formula on ", f_idx, " contains more than one random intercept effect. Adonis2 can only handle one random effect. Will continue with (1|", block, ").")))
            }
        }
        f.cache <- if (is.null(block)){f}else{paste0(c(mda.deparse(f), " + (1|", block, ")"), collapse="")}
        f.cache <- as.formula(f.cache)
        adj_var <- unlist(lapply(strsplit(D$formula$adj_vars[[f_idx]],split='[+]'), str_trim))
        variables <- c(mainvar, adj_var, block)

        meta_data.nona <- na.omit(D$meta_data[,variables,drop=FALSE])
        count_data.nona <- D$count_data[rownames(meta_data.nona),]

        Dcopy <- D
        Dcopy$count_data <- count_data.nona
        Dcopy$meta_data <- meta_data.nona
        
        out <- mda.cache_load_or_run_save(Dcopy, "beta", f.cache, {
            #meta_data.nona <- na.omit(D$meta_data[,variables,drop=FALSE])
            #count_data.nona <- D$count_data[rownames(meta_data.nona),]
            dist <- mda.cache_load_or_run_save(Dcopy, "beta_dist", y~x, {vegdist(count_data.nona)})
            strata <- if (is.null(block)){NULL}else{meta_data.nona[,block]}
            f <- update(f, dist ~ .)
            assign("dist",dist,envir=environment(f)) # Dumbest shit I ever saw
            
            mda.adonis2(f, meta_data.nona, na.action = 'na.exclude', strata=strata,
                        by = "margin", permutations=beta.permutations, parallel=beta.parallel,
                        scope=if(beta.hack){mainvar}else{NULL})
        }, order_invariant=!beta.hack)

        res.full <- out
        res.full <- res.full[,c("R2","F", "Pr(>F)")]


        names(res.full)[names(res.full)=="R2"] <- "effectsize"
        names(res.full)[names(res.full)=="Std. Error"] <- "se"
        names(res.full)[names(res.full)=="F"] <- "stat"
        names(res.full)[names(res.full)=="Pr(>F)"] <- "pvalue"
        
        res.full$variable <- rownames(res.full)
        res.full$se <- rep(NA, dim(res.full)[1])
        res.full$taxa <- rep("mda.beta", dim(res.full)[1])
        res.full$qvalue.withinformula <- res.full$pvalue
        
        res.full$formula <- rep(mda.deparse(f), dim(res.full)[1])
        res.full$method <- rep("beta", dim(res.full)[1])
        res.full$n <- rep(mda.meta.n(D, mainvar), dim(res.full)[1])
        res.full$freq <- rep(mda.meta.freq(D, mainvar), dim(res.full)[1])

        res <- res.full[startsWith(res.full$variable, mainvar),]
        res$variable <- rep(mainvar, dim(res)[1])

        return(list(res=res, res.full=res.full))
    }

    R <- lapply(1:length(D$formula$main_var), do)

    res <- bind_rows(lapply(R, function(x){x$res}))
    res$qvalue <- p.adjust(res$pvalue, "fdr")
    res.full <- bind_rows(lapply(R, function(x){x$res.full}))

    ###############################################################################
    # Format Output

    column.order <- c("taxa","variable","effectsize","se","stat","pvalue","qvalue.withinformula","qvalue","formula","method","n","freq")
    res <- res[,column.order]
    
    return(list(res=res, res.full=res.full, summary=mda.summary(res)))
}