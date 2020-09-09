library(maftools)

generateMafPlots <- function(mafs, summary_png_f, vaf_png_f) {
   png(summary_png_f);
   plotmafSummary(maf = mafs, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE);
   dev.off();

   png(vaf_png_f);
   plotVaf(maf=mafs);
   dev.off();

}

args <- commandArgs( trailingOnly = TRUE )
arg_in= args[1]
arg_summary_png_f = args[2]
arg_vaf_png_f = args[3]

mafs = read.maf(maf=arg_in)
generateMafPlots(mafs, arg_summary_png_f, arg_vaf_png_f)
