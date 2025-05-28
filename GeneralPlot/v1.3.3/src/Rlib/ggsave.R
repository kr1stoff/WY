plot_dev <- function (device, filename = NULL, dpi = 300) 
{
    force(filename)
    force(dpi)
    if (is.function(device)) 
        return(device)
    eps <- function(filename, ...) {
        grDevices::postscript(file = filename, ..., onefile = FALSE, 
            horizontal = FALSE, paper = "special")
    }
    devices <- list(eps = eps, ps = eps, 
	    tex  = function(filename, ...) grDevices::pictex(file = filename, ...),
		pdf  = function(filename, ...) Cairo::CairoPDF(file = filename, ...),
        svg  = function(filename, ...) svglite::svglite(file = filename, ...), 
		emf  = function(...) grDevices::win.metafile(...), 
        wmf  = function(...) grDevices::win.metafile(...), 
		png  = function(...) grDevices::png(..., res = dpi, units = "in"), 
		jpg  = function(...) grDevices::jpeg(..., res = dpi, units = "in"), 
		jpeg = function(...) grDevices::jpeg(..., res = dpi, units = "in"), 
		bmp  = function(...) grDevices::bmp(..., res = dpi, units = "in"), 
		tiff = function(...) grDevices::tiff(..., res = dpi, units = "in")
	)
    if (is.null(device)) {
        device <- to_lower_ascii(tools::file_ext(filename))
    }
    if (!is.character(device) || length(device) != 1) {
        abort("`device` must be NULL, a string or a function.")
    }
    dev <- devices[[device]]
    if (is.null(dev)) {
        abort(glue("Unknown graphics device '{device}'"))
    }
    dev
}
environment(plot_dev) <- asNamespace("ggplot2")
assignInNamespace("plot_dev", plot_dev, ns = "ggplot2")

