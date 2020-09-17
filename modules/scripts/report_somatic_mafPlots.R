library(maftools)

generateMafPlots <- function(mafs, summary_png_f, vaf_png_f, lolli1_png_f, lolli2_png_f, lolli3_png_f, lolli4_png_f) {
   png(summary_png_f);
   plotmafSummary(maf = mafs, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE);
   dev.off();

   png(vaf_png_f);
   plotVaf(maf=mafs);
   dev.off();
   
   #lollipop plots for top 4 genes
   geneSummary <- getGeneSummary(mafs);
   topgenes <- head(geneSummary[,1], n=4);
   lolli_f_names= c(lolli1_png_f, lolli2_png_f, lolli3_png_f, lolli4_png_f);
   for (tp in topgenes) {
       for (i in 1:length(tp)) {
          png(lolli_f_names[i]);
	  lollipopPlot(maf=mafs, gene=tp[i], showMutationRate=T);
	  dev.off();
       }
   }

}

args <- commandArgs( trailingOnly = TRUE )
arg_in= args[1]
arg_summary_png_f = args[2]
arg_vaf_png_f = args[3]

arg_lolli1_png_f = args[4]
arg_lolli2_png_f = args[5]
arg_lolli3_png_f = args[6]
arg_lolli4_png_f = args[7]


mafs = read.maf(maf=arg_in)
generateMafPlots(mafs, arg_summary_png_f, arg_vaf_png_f, arg_lolli1_png_f, arg_lolli2_png_f, arg_lolli3_png_f, arg_lolli4_png_f)
