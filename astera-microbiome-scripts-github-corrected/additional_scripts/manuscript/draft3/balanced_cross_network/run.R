library(igraph);library(tibble)
o<-'manuscript/draft3/balanced_cross_network';p<-'manuscript/draft3/carbon_paragraph_stats'
for(e in parse('run_nt_vs_till_balance_sensitivity.R'))if(is.call(e)&&as.character(e[[1]])=='<-'&&as.character(e[[2]]) %in% c('clr_transform','build_corr_network','network_stats'))eval(e)
markers<-c('16S','ITS','AMF');mats<-list();metas<-list()
for(k in markers){mats[[k]]<-as.matrix(read.csv(file.path(p,paste0(k,'_counts.csv')),row.names=1,check.names=FALSE));metas[[k]]<-read.csv(file.path(p,paste0(k,'_metadata.csv')))}
ids<-Reduce(intersect,lapply(mats,rownames));m<-metas[['16S']];m<-m[match(ids,m$sample_id),];nt<-m$sample_id[m$management_group=='NT'];ti<-m$sample_id[m$management_group=='till'];n<-min(length(nt),length(ti));set.seed(20261003);res<-list();draws<-list()
for(i in 1:20)for(group in c('NT','till')){
 keep<-sample(if(group=='NT')nt else ti,n);draws[[length(draws)+1]]<-data.frame(iter=i,group,sample_id=keep)
 parts<-lapply(markers,function(k){x<-mats[[k]][keep,,drop=FALSE];x<-x[,colMeans(x>0)>=.1,drop=FALSE];if(ncol(x)>120)x<-x[,order(colSums(x),decreasing=TRUE)[1:120],drop=FALSE];colnames(x)<-paste0(k,'|',colnames(x));x})
 x<-do.call(cbind,parts);g<-build_corr_network(x);s<-as.data.frame(network_stats(g,group));s$iter<-i;s$n_each<-n;res[[length(res)+1]]<-s
}
r<-do.call(rbind,res);write.csv(r,file.path(o,'cross_marker_balanced_iterations.csv'),row.names=FALSE);write.csv(do.call(rbind,draws),file.path(o,'resampled_ids.csv'),row.names=FALSE)
s<-aggregate(r[,c('nodes','edges','mean_degree','density','modularity','transitivity')],list(Group=r$label),median,na.rm=TRUE);write.csv(s,file.path(o,'cross_marker_balanced_medians.csv'),row.names=FALSE);print(s)
writeLines(c(paste('20 balanced resamplings;',n,'samples per group from',length(nt),'NT and',length(ti),'till common samples.'),'Seed 20261003. Sampling without replacement; draw identities saved. All three markers aligned by sample ID. Original cross-network selection: >=10% prevalence and top120 counts per marker, followed by original balanced-network global top250 filter, CLR with pseudocount1, Spearman abs(rho)>=0.6 and BH q<0.05 over all candidate edges. Isolates excluded. Cross-marker network includes within- and between-marker edges. Louvain uses absolute correlation weights; global transitivity is unweighted. Original six single-marker median rows retained; their historical draws differ from these newly computed cross-marker draws. These are descriptive sample-balanced networks and do not control profile dependence.',paste('igraph',packageVersion('igraph'),'Hmisc',packageVersion('Hmisc'))),file.path(o,'README.txt'))
