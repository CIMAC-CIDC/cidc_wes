library(maftools)

generateMafPlots <- function(mafs, summary_png_f, vaf_png_f, lolli1_png_f, lolli2_png_f, lolli3_png_f, lolli4_png_f, lolli5_png_f) {
   png(summary_png_f);
   plotmafSummary(maf = mafs, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE);
   dev.off();

   png(vaf_png_f);
   plotVaf(maf=mafs);
   dev.off();
   
   #lollipop plots for top 5 genes
   geneSummary <- getGeneSummary(mafs);
   lolli_f_names= c(lolli1_png_f, lolli2_png_f, lolli3_png_f, lolli4_png_f, lolli5_png_f);
   #topgenes <- head(geneSummary[,1], n=5);
   topgenes <- head(geneSummary[,1], n=length(lolli_f_names));
   for (tp in topgenes) {
       for (i in 1:length(tp)) {
          g <- tp[i];
          png(lolli_f_names[i]);
	  out <- tryCatch(
             { #TRY to generate the lollipop plot
	       lollipopPlot(maf=mafs, gene=g, showMutationRate=T);
	     }, error=function(cond) {
	        #Error state: print out helpful message and generate blank png
	        message(paste("WARNING: Lollipop plot call failed for gene", g));
	        message(cond);
		#Generate blank png
		par(mar=c(0,0,0,0), xpd=NA, mgp=c(0,0,0), oma=c(0,0,0,0), ann=F)
		plot.new()
		plot.window(0:1, 0:1)
	     });
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
arg_lolli5_png_f = args[8]


mafs = read.maf(maf=arg_in)
generateMafPlots(mafs, arg_summary_png_f, arg_vaf_png_f, arg_lolli1_png_f, arg_lolli2_png_f, arg_lolli3_png_f, arg_lolli4_png_f, arg_lolli5_png_f)
