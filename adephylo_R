# Load necessary libraries
library(ape)
library(adephylo)
library(phylobase)

# Read the tree from your file
tree <- read.tree("topology.tre")
rooted_tree <- root(tree, outgroup = tree$tip.label[1], resolve.root = TRUE)
data <- data.frame(
  TipLabel = rooted_tree$tip.label,
  PCA1 = runif(length(rooted_tree$tip.label), -1.5, 1.5),  # Assuming range -1.5 to 1.5
  PCA2 = runif(length(rooted_tree$tip.label), -1.5, 1.5)
)
rownames(data) <- data$TipLabel
phylo4d_obj <- phylo4d(rooted_tree, data[, -1])  # Exclude the TipLabel column
par(mar = rep(0.1, 4))
pdf("phylogenetic_tree_custom_traits_all_plots.pdf", width = 10, height = 10)
table.phylo4d(phylo4d_obj, var.lab = c("PCA1", "PCA2"), show.node = FALSE, cex.lab = 0.8, ratio.tree = 0.5)
ppca_result <- ppca(phylo4d_obj, scale = FALSE, scannf = FALSE, nfposi = 1, nfnega = 1, method = "Abouheif")
eig_values <- ppca_result$eig
significant_eigs <- which(eig_values > mean(eig_values))
tempcol <- rep("grey", length(eig_values))
tempcol[significant_eigs] <- "black"
barplot(eig_values, main = "pPCA eigenvalues", cex.main = 1.8, col = tempcol)
par(mar = rep(0.1, 4))
plot(ppca_result, ratio.tree = 0.7)

# Dotchart of contributions to PCA1 and PCA2
dotchart(data$PCA1, lab = data$TipLabel, main = "Contributions of PCA1")
abline(v = 0, lty = 2)
dotchart(data$PCA2, lab = data$TipLabel, main = "Contributions of PCA2")
abline(v = 0, lty = 2)

# Prepare data for phylogenetic autocorrelation tests
obj.ppca <- phylo4d_obj
tdata(obj.ppca, type = "tip") <- ppca_result$li

# Figure 1: Plot with labels
myLab <- paste(" ", rownames(ppca_result$li), sep = "")
par(mar = c(0.1, 2.4, 2.1, 1))
table.phylo4d(obj.ppca, ratio = 0.7, var.lab = c("1st global PC", "1st local\n   PC"), tip.label = myLab, box = FALSE, cex.lab = 1.4, cex.sym = 1.2, show.node.label = TRUE)
add.scatter.eig(ppca_result$eig, 1, 1, 1, csub = 1.2, posi = "topleft", ratio = 0.23)

# Figure 2: Plot arrows
s.arrow(ppca_result$c1, xlim = c(-1, 1), clab = 1.3, cgrid = 1.3)

# Final plot - phylogenetic tree with trait labels
plot(phylo4d_obj, type = "phylogram", edge.color = "blue", tip.color = "red", cex.leg = 0.6, font = 2, use.edge.length = TRUE)

# Close the PDF device
dev.off()
