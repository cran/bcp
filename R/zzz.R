.First.lib <- function(lib, pkg) {
	library.dynam("bcp", pkg, lib)
}