library(vegan)
o <- 'manuscript/draft3/carbon_adjusted'
t <- read.csv(file.path(o,'community_carbon.csv'))
z <- list()
for (marker in c('16S','ITS','AMF')) for (variant in c('recorded','exclude_conflicting')) {
 m <- read.csv(file.path(o,paste0(marker,'_',variant,'_paired_metadata.csv')))
 y <- as.matrix(read.csv(file.path(o,paste0(marker,'_',variant,'_coordinates.csv'))))
 fit <- rda(y ~ carbon + Condition(Site_Name + Depth_cm), data=m)
 f <- (fit$CCA$tot.chi/fit$CCA$rank)/(fit$CA$tot.chi/fit$CA$rank)
 # Residual df uses n minus design rank, not residual matrix rank.
 dfres <- nrow(y)-qr(model.matrix(~Site_Name+Depth_cm+carbon,m))$rank
 f <- fit$CCA$tot.chi/(fit$CA$tot.chi/dfres)
 partial <- fit$CCA$tot.chi/(fit$CCA$tot.chi+fit$CA$tot.chi)
 a <- t[t$marker==marker & t$variant==variant & t$model=='within_profile_carbon_adjusted_depth', ]
 stopifnot(nrow(a)==1, abs(f-a$pseudo_F)<1e-8,abs(partial-a$partial_R2)<1e-8)
 z[[length(z)+1]]<-data.frame(marker,variant,F_vegan=f,partial_R2_vegan=partial,passed=TRUE)
}
write.csv(do.call(rbind,z),file.path(o,'vegan_validation.csv'),row.names=FALSE)
cat('PASS: all six partial models independently reproduced in vegan',as.character(packageVersion('vegan')),'\n')
