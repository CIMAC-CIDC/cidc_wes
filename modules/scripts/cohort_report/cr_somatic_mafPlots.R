library(maftools)

generateMafPlots <- function(mafs, summary_png_f, onco_png_f, titv_png_f, vaf_png_f, tcga_png_f, interact_png_f) {
   png(summary_png_f);
   plotmafSummary(maf = mafs, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE);
   dev.off();

   png(onco_png_f);
   oncoplot(maf = mafs, top=25);
   dev.off();

   png(titv_png_f);
   tmp.titv = titv(maf=mafs, plot=F, useSyn=T);
   plotTiTv(res= tmp.titv);
   dev.off();

   png(vaf_png_f);
   plotVaf(maf=mafs);
   dev.off();

   png(tcga_png_f);
   tmp.mutLoad = tcgaCompare(maf=mafs, cohortName="Cohort")
   dev.off();

   png(interact_png_f);
   somaticInteractions(maf = mafs, top = 25, pvalue = c(0.05, 0.1))
   dev.off();
}

args <- commandArgs( trailingOnly = TRUE )
arg_in= args[1]
arg_summary_png_f = args[2]
arg_onco_png_f = args[3]
arg_titv_png_f = args[4]
arg_vaf_png_f = args[5]
arg_tcga_png_f = args[6]
arg_interact_png_f = args[7]

mafs = read.maf(maf=arg_in)
generateMafPlots(mafs, arg_summary_png_f, arg_onco_png_f, arg_titv_png_f, arg_vaf_png_f, arg_tcga_png_f, arg_interact_png_f)
