library(TcellExTRECT)

#bam_file = "data/CBUPUKJZF.01.sorted.dedup.bam"
#bam_file = "analysis/align/CBUPRG5WT.01/CBUPRG5WT.01.sorted.dedup.bam"

tcell_est <- function(bam_file, capture_bed, out_dir, out_name) {
   #Load twist custom exons
   #twist_capture_bed = './ref_files/hg38/target_beds/twist.broad.mdacc.mocha.liftover.hg38.sorted.merged.bed';

   TCRA_exons_hg38_twist <- createExonDFBed(capture_bed, 'hg38')

   #load segment file
   data("tcra_seg_hg38")
   data(TCRA_exons_hg38)


   #get coverage data- generates a .txt file
   cov.file <- getCovFromBam(bamPath = bam_file,
                             outPath = paste0(out_dir, "/"),
                             vdj.seg = tcra_seg_hg38)
			  
   cov_df <- loadCov(cov.file) #load the .txt file that was just generated

   #Run
   TCRA.out <- runTcellExTRECT(cov_df, TCRA_exons_hg38_twist, tcra_seg_hg38, 'hg38')
   #Alt w/ purity and copy number:
   #TCRA.out <- adjustTcellExTRECT(TCRA.out, purity = 0.5, TCRA.cn = 3)


   #write to file
   out_f <- paste0(out_dir, "/", out_name, "_tcellextrect.txt")
   write.csv(TCRA.out, out_f, row.names=F, quote = F)

   #PLOT:
   plot_f <- paste0(out_dir, "/", out_name, "_tcellplot.pdf")
   pdf(plot_f)
   plotTcellExTRECT(cov_df, TCRA_exons_hg38_twist, tcra_seg_hg38,'hg38', sample_name = out_name)
   dev.off()
}

args <- commandArgs( trailingOnly = TRUE )
arg_bam= args[1]
arg_bed= args[2] #path to twist capture regionds
arg_out_dir = args[3]
arg_out_name = args[4]

tcell_est(arg_bam, arg_bed, arg_out_dir, arg_out_name);

