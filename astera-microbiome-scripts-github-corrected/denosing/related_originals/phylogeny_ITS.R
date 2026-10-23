####phylogeny
##load libraries
# Set the R library location
.libPaths(new = "/fs/project/PAS1212/abdelrahim")
install.packages("phangorn")
library(phangorn)
install.packages("DECIPHER")
library(DECIPHER)
library(dada2)
indir <- "/fs/ess/PAS1212/Afaf_rich/ITS/results/results/dada"
seqtab <- readRDS(file.path(indir, "seqtab.rds"))
# Extract the ASV sequences from the seqtab object
# This propagates to the tip labels of the tree.
# At this stage ASV labels are full ASV sequence
seqs <- getSequences(seqtab)
names(seqs) <- seqs 
alignment <- AlignSeqs(DNAStringSet(seqs),
                       anchor = NA,
                       iterations = 5,
                       refinements = 5)

print('Computing pairwise distances from ASVs...')
Sys.time()
phang.align <- phangorn::phyDat(as(alignment, 'matrix'), type = 'DNA')
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order is not sequence order
fit = pml(treeNJ, data = phang.align)
print('...Done.')
Sys.time()

print('Fit GTR model...')
Sys.time()
fitGTR <- update(fit, k = 4, inv = 0.2)
print('...Done.')
Sys.time()

print('Computing likelihood of tree...')
Sys.time()
fitGTR <- optim.pml(fitGTR,
                    model = 'GTR',
                    optInv = TRUE,
                    optGamma = TRUE,
                    rearrangement = 'NNI',
                    control = pml.control(trace = 0))
print('...Done')
Sys.time()


# Root the tree
rooted_tree <- root_by_longest_edge(fitGTR)
write.tree(fitGTR$tree, file = "GTR.phy")

