#Script to generate depth coverage by certain thresholds
# input- mosdepth.region.dist.txt files, output path, sample name
# output- csv file with % reads with 20x coverage, 50x, 100x, 150x etc

depth_summary<-function(arg_in, arg_sample_name, arg_out){
   df<-read.table(arg_in);
   df_ag<-aggregate(df$V3, by=list(Category=df$V2), FUN=sum)
   df_ag$x<-df_ag$x/(sum(df_ag$x))
   depth20=(1-cumsum(df_ag[,2]))[21]
   depth50=(1-cumsum(df_ag[,2]))[51]
   depth100=(1-cumsum(df_ag[,2]))[101]
   depth150=(1-cumsum(df_ag[,2]))[151]

   results<-matrix(c(depth20,depth50,depth100,depth150), nrow=1)
   colnames(results) <-c('x20','x50','x100','x150')
   write.csv(results, arg_out, row.names=F, quote = F)
   
   tmp<-as.data.frame(data.matrix(results))
}

args <- commandArgs( trailingOnly = TRUE )
arg_in= args[1]
arg_sample_name = args[2]
arg_out = args[3]

#print(arg_in)
#print(arg_out)
#print(arg_sample_name)

depth_summary(arg_in, arg_sample_name, arg_out)
