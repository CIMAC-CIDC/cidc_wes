library(pdftools)

convert_pdf2png <- function(arg_in, arg_out, arg_pg){
   pdf_convert(arg_in, format='png', pages = arg_pg:arg_pg, filenames = c(arg_out), dpi = 72, antialias = TRUE, opw = "", upw = "", verbose = FALSE)
}

args <- commandArgs( trailingOnly = TRUE )
arg_in= args[1]
arg_out = args[2]
arg_pg = args[3]

convert_pdf2png(arg_in, arg_out, arg_pg)
