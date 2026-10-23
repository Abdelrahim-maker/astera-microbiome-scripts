library(vegan)
o <- 'manuscript/draft3/carbon_paragraph_stats'
markers <- c('16S','ITS','AMF'); mats<-list(); metas<-list();out<-list();checks<-list()
for (marker in markers) {
 x<-as.matrix(read.csv(file.path(o,paste0(marker,'_counts.csv')),row.names=1,check.names=FALSE));m<-read.csv(file.path(o,paste0(marker,'_metadata.csv')));stopifnot(identical(rownames(x),m$sample_id));m$management<-factor(m$management_group,levels=c('NT','till'));mats[[marker]]<-x;metas[[marker]]<-m
}
common<-Reduce(intersect,lapply(mats,rownames)); parts<-list()
for (marker in markers) {
 x<-mats[[marker]][common,,drop=FALSE];x<-x[,colMeans(x>0)>=.05,drop=FALSE]
 if(ncol(x)>200)x<-x[,order(colSums(x),decreasing=TRUE)[1:200],drop=FALSE]
 z<-log(x+1);parts[[marker]]<-z-rowMeans(z)
}
joint<-do.call(cbind,parts);joint<-scale(joint);stopifnot(all(is.finite(joint)))
m<-metas[['16S']];m<-m[match(common,m$sample_id),];metas[['Integrated']]<-m
for (marker in c(markers,'Integrated')) {
 m<-metas[[marker]];d<-if(marker=='Integrated')dist(joint) else vegdist(decostand(mats[[marker]],'total'),'bray')
 for (variant in c('recorded','exclude_conflicting')) {
  keep<-if(variant=='recorded')rep(TRUE,nrow(m)) else tolower(as.character(m$carbon_conflict))!="true"
  md<-m[keep,];dd<-as.dist(as.matrix(d)[keep,keep]);set.seed(20261003)
  fit<-adonis2(dd ~ management + carbon + management:carbon,data=md,permutations=9999,by='terms')
  t<-as.data.frame(fit);t$term<-rownames(t);t$marker<-marker;t$variant<-variant;t$n<-nrow(md);out[[length(out)+1]]<-t
 }
}
r<-do.call(rbind,out);names(r)[names(r)=='Pr(>F)']<-'p';r$q_BH<-NA_real_
for(v in unique(r$variant)){ix<-which(r$variant==v & !is.na(r$p));r$q_BH[ix]<-p.adjust(r$p[ix],'BH')}
write.csv(r,file.path(o,'sequential_permanova_rerun.csv'),row.names=FALSE)
legacy<-read.csv('analysis_outputs/carbon_analysis/carbon_permanova_by_domain.csv',check.names=FALSE)
for(marker in markers){a<-subset(legacy,domain==marker & test_type=='sequential');b<-r[r$marker==marker & r$variant=='recorded',];stopifnot(max(abs(a$R2-b$R2))<1e-8);checks[[marker]]<-data.frame(marker,max_R2_difference=max(abs(a$R2-b$R2)),passed=TRUE)}
a<-read.csv('analysis_outputs/carbon_analysis/carbon_joint_permanova.csv');a<-subset(a,test_type=='sequential');b<-r[r$marker=='Integrated' & r$variant=='recorded',];checks[['Integrated']]<-data.frame(marker='Integrated',max_R2_difference=max(abs(a$R2-b$R2)),passed=max(abs(a$R2-b$R2))<1e-8)
write.csv(do.call(rbind,checks),file.path(o,'legacy_R2_validation.csv'),row.names=FALSE)
writeLines(c(paste('vegan',packageVersion('vegan')),'9999 unrestricted permutations; seed 20261003 reset for each model. Integrated carbon is from the 16S records.','Integrated input rows explicitly matched by sample ID. Sensitivity removes conflicting samples after fixing feature selection and integrated scaling on the common cohort.'),file.path(o,'run_notes.txt'))
print(r[,c('marker','variant','term','R2','p','q_BH')]);print(do.call(rbind,checks))
