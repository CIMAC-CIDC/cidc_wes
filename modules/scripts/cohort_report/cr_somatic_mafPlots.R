library(maftools)

generateMafPlots <- function(mafs, cancerGeneList_f, summary_png_f, onco_png_f, titv_png_f, vaf_png_f, tcga_png_f, interact_png_f, lolli1_png_f, lolli2_png_f, lolli3_png_f, lolli4_png_f, lolli5_png_f, lolli6_png_f) {
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

   #lollipop plots for top 6 genes
   #read in geneList
   cncrGnList <- read.table(file=cancerGeneList_f, sep = '\t', header = TRUE);
   geneSummary <- getGeneSummary(mafs);
   topgenes <- unlist(geneSummary[,1]);
   topCnrGenes <- head(topgenes[topgenes %in% cncrGnList$Hugo.Symbol], n=6);
   #print(topCnrGenes);

   lolli_f_names= c(lolli1_png_f, lolli2_png_f, lolli3_png_f, lolli4_png_f, lolli5_png_f, lolli6_png_f);
   i = 1;
   for (g in topCnrGenes) {
      png(lolli_f_names[i]);
      lollipopPlot(maf=mafs, gene=g, showMutationRate=T);
      dev.off();
      i = i + 1;
   }
}

args <- commandArgs( trailingOnly = TRUE )
arg_in= args[1]
args_geneList = args[2]

arg_summary_png_f = args[3]
arg_onco_png_f = args[4]
arg_titv_png_f = args[5]
arg_vaf_png_f = args[6]
arg_tcga_png_f = args[7]
arg_interact_png_f = args[8]

arg_lolli1_png_f = args[9]
arg_lolli2_png_f = args[10]
arg_lolli3_png_f = args[11]
arg_lolli4_png_f = args[12]
arg_lolli5_png_f = args[13]
arg_lolli6_png_f = args[14]

mafs = read.maf(maf=arg_in)

#this is starting to get ugly!
generateMafPlots(mafs, args_geneList, arg_summary_png_f, arg_onco_png_f, arg_titv_png_f, arg_vaf_png_f, arg_tcga_png_f, arg_interact_png_f, arg_lolli1_png_f, arg_lolli2_png_f, arg_lolli3_png_f, arg_lolli4_png_f, arg_lolli5_png_f, arg_lolli6_png_f)