suppressMessages(library("optparse"))
suppressMessages(library("tools"))
suppressMessages(library("Mfuzz"))
suppressMessages(library("GenomicFeatures"))
suppressMessages(library("DESeq2"))
suppressMessages(library("edgeR"))
suppressMessages(library("ggplot2"))

#nohup Rscript time_series.R 4T 3 /disk1/xilu/data/project-hypoxia/hypoxia_plus_minus/bedtools_deseq2/4T.sig.q.0.01.2 cluster_with_4T_3 /disk1/xilu/data/project-hypoxia/hypoxia_plus_minus/time_series_deseq2_q_0.01 > 4T.time.series.3.o &
set.seed("1234567")

args<-commandArgs(TRUE)
input1<-args[1]
input3<-args[2]
output<-args[3]

#load("/disk1/xilu/collaborate/niklas/5.1.differential_peak_vs_d0/atacDDS_deseq.RData")
#normalized_counts_n <- counts(atacDDS_deseq,normalized=FALSE)
normalized_counts<-read.table("/disk1/xilu/collaborate/niklas/12.30m_bam_peakmapped_with_10m_bam/differential_peaks/58.sample.peak.remove.chry.random.chry.xls",header=T)
peak_name<-paste(normalized_counts[,1],normalized_counts[,2],normalized_counts[,3],sep="_")
normalized_counts_n<-normalized_counts[,c(4:dim(normalized_counts)[[2]])]
rownames(normalized_counts_n)<-peak_name

gene_length<-read.table("/disk1/xilu/collaborate/niklas/12.30m_bam_peakmapped_with_10m_bam/differential_peaks/length.txt",sep="\t",header=FALSE)
rownames(gene_length) <- gene_length[,1]
gene_length2<-data.frame(length = gene_length[,2])
rownames(gene_length2)<-rownames(gene_length) 

comparison<-read.table("/disk1/xilu/collaborate/niklas/4.differential_peaks/final.comparison",header=T,sep="\t")

a<-NULL
for ( i in comparison[,1]){
    b<-which(colnames(normalized_counts_n) == i)
    a<-c(a,b)
}
normalized_counts<-normalized_counts_n[,a]

#data_norm<-merge(normalized_counts,gene_length2,by="row.names",all.x=T)

#data_norm_tmp = data_norm[,-1]
#rownames(data_norm_tmp)<-data_norm[,1]
#data_norm<-data_norm_tmp
#rpkm<-rpkm(data_norm[,1:dim(data_norm)[2]-1],gene.length=data_norm$length,normalized.lib.sizes=TRUE, log=FALSE)
rpkm<-rpkm(normalized_counts,gene.length=gene_length2$length,normalized.lib.sizes=TRUE, log=FALSE)

#for each time value containing replicates, we calculate the fpkm count means ###
exprs<-as.matrix(rpkm)
time<-factor(comparison[,2])
#time<-factor(c("experiment_4T","6h_RAW","6h_365","6h_365","12h_4T","12h_4T","experiment_4T","experiment_RAW","experiment_RAW","experiment_365","experiment_365","6h_4T","6h_4T","6h_RAW","12h_RAW","24h_365","36h_4T","36h_4T","36h_RAW","36h_RAW","12h_RAW","12h_365","12h_365","24h_4T","24h_4T","24h_RAW","24h_RAW","24h_365","36h_365","72h_4T","72h_RAW","72h_RAW","72h_365","72h_365","36h_365","48h_4T","48h_4T","48h_RAW","48h_RAW","48h_365","48h_365","72h_4T"))
count=1
for ( i in unique(time) ){
  if ( dim(as.data.frame(exprs[,which(time==i)]))[2] == 1 ){
    mean_rpkm=data.frame(exprs[,which(time==i)])
  } else {
    mean_rpkm=data.frame(rowMeans(exprs[,which(time==i)]))
  }
  colnames(mean_rpkm)=i
  if (count == 1){
    mean_rpkm_ok=mean_rpkm
  } else {
    mean_rpkm_ok=merge(mean_rpkm_ok,mean_rpkm,by="row.names")
    rownames(mean_rpkm_ok)=mean_rpkm_ok[,1]
    mean_rpkm_ok=mean_rpkm_ok[,-1]
  }
  count=count+1
}
exprs_with_time=as.matrix(mean_rpkm_ok, header=TRUE, sep="\t",row.names=1,as.is=TRUE)

sig<-read.table(input1,header=FALSE)
exprs_with_time_4T<-exprs_with_time[,c("experiment_d0","experiment_d19","experiment_d35","experiment_d65")]
sig_4T<- exprs_with_time_4T[rownames(exprs_with_time_4T) %in% as.matrix(sig),]

nb_clusters <- input3
nb_clusters <- as.numeric(nb_clusters)
filename<-output
exprSet=ExpressionSet(assayData=sig_4T)
exprSet.r=filter.NA(exprSet, thres=0.25)
exprSet.f=fill.NA(exprSet.r,mode= "mean")
tmp=filter.std(exprSet.f,min.std=0, visu=FALSE)
exprSet.s=standardise(tmp)
# clustering
m1=mestimate(exprSet.s)
print (nb_clusters)
print (m1)
cl=mfuzz(exprSet.s,c=nb_clusters,m=m1,ite=1000)
print ("ok_______")    
output = getwd()
membership_cutoff = c(0.4,0.5,0.6,0.7)
#filename = output
data<-sig_4T

for (membership in membership_cutoff){
  membership=as.numeric(membership)
  # create one output folder per membership
  dir=paste(output,paste(filename,membership, sep="_"),sep="/")
  dir.create(dir, showWarnings = FALSE)
  #--------------------------------------------------------------------------------------------------------------------------
  # membership cut-off part and plot clusters
  pdf(paste(dir,paste(paste("clusters_Mfuzz_membership_equals_",membership,sep=""),".pdf",sep=""), sep="/"))
  #mfuzz.plot2(exprSet.s,cl=cl,time.labels=colnames(data),min.mem=membership, colo="fancy", x11=FALSE)
  mfuzz.plot(exprSet.s,cl=cl,mfrow=c(1,1),time.labels=c("experiment_d0","experiment_d19","experiment_d35","experiment_d65"),min.mem=membership,new.window=FALSE)

  dev.off()
  setEPS()
  postscript(paste(dir,paste(paste("clusters_Mfuzz_membership_equals_",membership,sep=""),".eps",sep=""), sep="/"))
  #--------------------------------------------------------------------------------------------------------------------------

  #--------------------------------------------------------------------------------------------------------------------------
  # generates one genes list per cluster
  acore.list=acore(exprSet.s,cl=cl,min.acore=membership)
  print("----")
  print(paste("Membership",membership,sep=" : "))
  for (cluster in 1:nb_clusters){
    print(paste(paste("Number of genes in cluster", cluster, sep=" "),dim(acore.list[[cluster]])[1], sep=" : "))
    cluster_table=merge(data,acore.list[[cluster]][2], by="row.names", all.y=TRUE)
    write.table(cluster_table,paste(dir,paste(paste("list_of_genes_in_cluster",cluster,sep="_"),".txt",sep=""),sep="/"), sep="\t",row.names=F, dec=".")
  }
  print("----")
  #--------------------------------------------------------------------------------------------------------------------------
}

fcm_centroids <- cl[1]
fcm_centroids_df <- data.frame(fcm_centroids)
colnames(fcm_centroids_df) <- as.matrix(c("experiment_d0","experiment_d19","experiment_d35","experiment_d65"))
fcm_centroids_df$cluster <- row.names(fcm_centroids_df)
centroids_long <- tidyr::gather(fcm_centroids_df,"sample",'value',1:4)
centroids_long$sample<-factor(centroids_long$sample,levels=c("experiment_d0","experiment_d19","experiment_d35","experiment_d65"))
for (membership in membership_cutoff){
  membership=as.numeric(membership)
dir=paste(output,paste(filename,membership, sep="_"),sep="/")
pdf(paste(dir,paste(paste("clusters_",nb_clusters,membership,sep=""),".pdf",sep=""), sep="/"))
ggplot(centroids_long, aes(x=sample,y=value, group=cluster, colour=as.factor(cluster))) +
  geom_line() +
  xlab("Time") +
  ylab("Expression") +
  labs(title= "Cluster Expression by Time",color = "Cluster")
dev.off()
}
