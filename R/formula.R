formula.rhs <- function(fn){
    fn.char <- as.character(fn)
    rhs <- paste0(c('~', fn.char[if (length(fn.char) == 3){3} else {2}]), collapse=' ')
    as.formula(rhs)
}

formula.ismixed <- function(fn){
    any(grepl('[|]', as.character(fn)))
}

formula.ismixed.intercept <- function(fn){
    length(formula.parts.random.intercept(fn)) > 0
}

formula.ismixed.slope <- function(fn){
    length(formula.parts.random.slope(fn)) > 0
}

formula.norand <- function(fn){
    if (!formula.ismixed(fn)){
        return(fn)
    }
    
    fe.parts <- formula.parts.fixed(fn)
    fn.fixed <- paste0(c('~', paste0(unlist(fe.parts), collapse=' + ')), collapse=' ')
    return(as.formula(fn.fixed))
}

formula.parts <- function(fn){
    attributes(terms(formula.rhs(fn), keep.order=TRUE))$term.labels
}

formula.parts.fixed <- function(fn){
    parts <- formula.parts(fn)
    
    intercept <- attributes(terms(formula.rhs(fn)))$intercept
    fe.parts <- parts[!grepl('[|]', parts)]

    if (intercept == 0){
        fe.parts <- c(fe.parts, ' - 1')
    }
    fe.parts
}

formula.parts.random <- function(fn){
    parts <- formula.parts(fn)
    parts[grepl('[|]', parts)]
}

formula.parts.random.slope <- function(fn){
    parts <- formula.parts.random(fn)
    if (length(parts) == 0){
        return(parts)
    }
    parts[sapply(parts, function(p){
        ps <-  strsplit(p, split='[|]')[[1]]
        effect <- trimws(ps[1])
        effect != '1'
        })]
}

formula.parts.random.intercept <- function(fn){
    parts <- formula.parts.random(fn)
    if (length(parts) == 0){
        return(parts)
    }
    parts[sapply(parts, function(p){
        ps <-  strsplit(p, split='[|]')[[1]]
        effect <- trimws(ps[1])
        effect == '1'
        })]
}

formula.model.matrix <- function(fn, data){
    data.new <- model.matrix(terms(fn, keep.order=TRUE), data)
    return(data.new)
}

formula.reformulate <- function(fn, data){
    f <- if (formula.ismixed(fn)) {formula.reformulate.mixed} else {formula.reformulate.fixed}
    new <- f(fn, data)
    
    parts <- formula.parts.fixed(fn)
    
    nfreq <- mda.nfreq(new$data)
    nfreq$variable <- sapply(rownames(nfreq), function(v){new$map[v]})
    nfreq['(Intercept)',] <- c(NA, NA, '(Intercept)')
    nfreq$variable.mda <- rownames(nfreq)
    
    fn.orig <- Reduce(function(fn, i){
                    as.formula(stringr::str_replace_all(mda.deparse(fn), names(new$map)[i], new$map[i]))
                }, 1:length(new$map), init=new$fn)
    
    list(
        fn.raw=fn,
        fn.orig=fn.orig,
        fn=new$fn,
        data=as.data.frame(new$data),
        map=new$map,
        nfreq=nfreq,
        
        
        parts.fixed=formula.parts.fixed(new$fn),
        parts.random=formula.parts.random(new$fn),
        parts.random.slope=formula.parts.random.slope(new$fn),
        parts.random.intercept=formula.parts.random.intercept(new$fn),
        
        checksum.order_invariant = digest(c(sort(formula.parts(fn.orig)), new$data), algo="md5"),
        checksum.order_variant   = digest(c(formula.parts(fn.orig), new$data), algo="md5")
        )
}

formula.reformulate.fixed <- function(fn, data){
    intercept <- attributes(terms(formula.rhs(fn)))$intercept
    
    data.new <- formula.model.matrix(fn, data)
    
    if (intercept == 1){
        data.new <- subset(data.new, select=-`(Intercept)`)
    }
    
    data.names.map <- colnames(data.new)
    names(data.names.map) <- sapply(1:dim(data.new)[2], function(x){sprintf('mda_p%05d', x)})
    colnames(data.new) <- names(data.names.map)
    
    fn.new  <- paste0(c('~', paste0(names(data.names.map), collapse=' + ')), collapse=' ')
    
    if (intercept == 0){
        fn.new <- paste0(c(fn.new, ' - 1'), collapse='')
    }
    
    return(list(
        fn=as.formula(fn.new),
        data=data.new,
        map=data.names.map
        ))
}

formula.reformulate.mixed <- function(fn, data){

    parts <- formula.parts(fn)
    
    intercept <- attributes(terms(formula.rhs(fn)))$intercept
    fe.parts <- parts[!grepl('[|]', parts)]
    re.parts <- parts[grepl('[|]', parts)]

    fn.fixed <- as.formula(paste0(c('~', paste0(unlist(fe.parts), collapse=' + ')), collapse=' '))
    fe.data <- subset(formula.model.matrix(fn.fixed, data), select=-`(Intercept)`)
    
    fe.data.missing_rows <- setdiff(rownames(data), rownames(fe.data))
    fe.data.missing_data <- as.data.frame(matrix(nrow=length(fe.data.missing_rows), ncol=length(colnames(fe.data))), rownames=fe.data.missing_rows)
    colnames(fe.data.missing_data) <- colnames(fe.data)

    fe.data <- rbind(fe.data, fe.data.missing_data)

    # Get data for just the random effect parts
    re <- lapply(re.parts, function(part){

        parts <-  strsplit(part, split='[|]')[[1]]
        effect <- trimws(parts[1])
        group <-  trimws(parts[2])

        effect.form <- as.formula(paste0(c('~', paste0(effect, collapse=' + ')), collapse=" "))
        effect.data <- subset(model.matrix(effect.form, data=data), select=-`(Intercept)`)
        effect.data <- cbind(data[,group,drop=FALSE], effect.data)

        return(list(
            data=effect.data,
            group=group))
    })

    # Merge that data to a single matrix
    re.data <- lapply(re, function(x){x$data})
    re.data <- re.data[unlist(lapply(re.data, function(x){!is.null(colnames(x))}))]
    re.data <- Reduce(
        function(x, y){
            newcols <- setdiff(colnames(y), colnames(x))
            bound <- cbind(x, y[,newcols])
            colnames(bound) <- c(colnames(x), newcols)
            bound
        }, re.data)
    
    re.data.missing_rows <- setdiff(rownames(data), rownames(re.data))
    re.data.missing_data <- as.data.frame(matrix(nrow=length(re.data.missing_rows), ncol=length(colnames(re.data))), rownames=re.data.missing_rows)
    colnames(re.data.missing_data) <- colnames(re.data)
    re.data <- rbind(re.data, re.data.missing_data)
    
    # Merge the fixed and random effects into one table
    data.new <- cbind(fe.data, re.data[,setdiff(colnames(re.data), colnames(fe.data)), drop=FALSE])
    data.names.map <- colnames(data.new)
    names(data.names.map) <- sapply(1:dim(data.new)[2], function(x){sprintf('mda_p%05d', x)})

    data.names.pam <- names(data.names.map)
    names(data.names.pam) <- data.names.map

    colnames(data.new) <- names(data.names.map)

    fe.form.parts <- sapply(colnames(fe.data), function(x){data.names.pam[x]})

    re.form.parts <- sapply(re, function(x){
        effect <- paste0(sapply(tail(colnames(x$data), n=-1), function(x){data.names.pam[x]}), collapse=' + ')
        group <- data.names.pam[x$group]
        if (nchar(effect) == 0){
            effect <- '1'
        }
        mixed <- paste0(c('(', effect, '|', group, ')'), collapse='')
    })

    fn.new <- paste0(c('~', paste0(c(fe.form.parts, re.form.parts), collapse=' + ')), collapse=' ')
    
    if (intercept == 0){
        fn.new <- paste0(c(fn.new, ' -1'), collapse='')
    }

    list(
        fn=as.formula(fn.new),
        data=data.new,
        map=data.names.map
        )
    
}

mda.nfreq <- function(data, maxunique=100){
    n <- function(var){
        sum(unlist(lapply(data[,var],function(x){!(is.na(x))})))
    }

    freq <- function(var){
        if (is.numeric(data[,var])){
            ""
        }
        else if (length(unique(data[,var])) >= maxunique){
            ""
        } else {
            paste0(mapply(function(x,y){paste0(c(x,y), collapse=":")}, names(table(data[,var])), as.character(table(data[,var]))), collapse=", ")
        }
    }
    
    as.data.frame(list(
        n=sapply(colnames(data), n),
        freq=sapply(colnames(data), freq)
    ))
}
