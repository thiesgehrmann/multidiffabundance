 ###############################################################################
# custom
mda.custom <- function(mda.D, custom.file, custom.name="custom.method", ...){
    D <- mda.D

    # modified from https://stackoverflow.com/a/58880544
    tmp <- new.env(parent=parent.frame())
    source(custom.file, local = tmp)
    print(names(parent.frame()))
    for(x in names(tmp)) {
        if( x %in% names(parent.frame()) ){
            mda.message(paste0(c("Overwriting existing environment variable `", x, "`"), collapse=""), type="warning")
        }
        assign(x, tmp[[x]], envir = parent.frame())
    }

    method <- parent.frame()[[custom.name]]

    r <- method(D, ...)

    return(r)
}