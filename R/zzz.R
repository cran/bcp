.onLoad <- function(libname, pkgname) {
    library.dynam("bcp", pkgname, libname);
}

.noGenerics <- TRUE

.onUnload <- function(libpath) {
    library.dynam.unload("bcp", libpath);
}