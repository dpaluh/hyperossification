#geomorph shape analyses for Paluh et al.: Evolution of hyperossification expands skull diversity in frogs

#load packages
library(geomorph) #version 3.0.3
library(ape)
library(geiger)
library(phytools)
library(beeswarm)
library(RevGadgets)

#read in classifier file
classifier  <- read.csv("classifier.csv", header=T, row.names=1)

#read in NTS files
filelist <- list.files(pattern = ".nts")
mydata <- readmulti.nts(filelist)
readmulti.nts(filelist)
dim(mydata)
data(mydata)

#Procrustes analysis
Y.gpa<-gpagen(mydata)
PCA <- plotTangentSpace(Y.gpa$coords, axis1 = 1, axis2 = 2, label=TRUE)

#Save PC scores
PCscores <- PCA$pc.scores
write.csv(PCscores, file = "PCscores.csv")
pc1 <- read.csv("PCscores.csv", header=TRUE, row.names=1)
pc1 <- as.matrix(pc1)[,1]
pc1
pc2 <- read.csv("PCscores.csv", header=TRUE, row.names=1)
pc2 <- as.matrix(pc2)[,2]
pc2
pc3 <- read.csv("PCscores.csv", header=TRUE, row.names=1)
pc3 <- as.matrix(pc3)[,3]
pc3

#read in Jetz and Pyron 2018 tree
tree <- read.tree("amph_shl_new_Posterior_7238.1.tre")

#prune tips
species.to.keep <-read.csv("classifier.csv", header = TRUE)
pruned.tree<-drop.tip(tree,tree$tip.label[-na.omit(match(species.to.keep[,1],tree$tip.label))])
write.tree(pruned.tree, file = "pruned_anura_final.phy")
mytree <- read.tree("pruned_anura_final.phy")
name.check(mytree, classifier)

#phylogenetic signal for skull shape
PS.shape<-physignal(Y.gpa$coords, mytree, iter=10000)
summary(PS.shape)

#phylogenetic signal for skull size
PS.size<-physignal(Y.gpa$Csize, mytree, iter=10000)
summary(PS.size)

#geomorph dataframe and sorting tips & traits
gdf <- geomorph.data.frame(Y.gpa, hyperossified = classifier$Hyperossified, diet = classifier$feeding_biology, diet2 = classifier$feeding_biology2, microhabitat = classifier$microhabitat, phragmosis = classifier$phragmosis, phy = mytree)
compare <- treedata(mytree, classifier, sort=TRUE)
classifier_sorted <- as.data.frame(compare$data)
tree_sorted <- compare$phy
write.tree(tree_sorted, "tree_sorted.tre")
tree <- read.tree("tree_sorted.tre")

#morphological disparity
md <-morphol.disparity(coords ~ hyperossified, data = gdf, iter = 10000)
summary(md)

#compare evolutionary rates; "method = simulation" and "method = permutation"
coords<-Y.gpa$coords
ER<-compare.evol.rates(coords, phy=tree, gp=classifier_sorted$Hyperossified, iter = 10000, method = c("simulation"), print.progress = TRUE)
summary(ER)
ER_perm<-compare.evol.rates(coords, phy=tree, gp=classifier_sorted$Hyperossified, iter = 10000, method = c("permutation"), print.progress = TRUE)
summary(ER_perm)

#phylogenetic MANOVA: hyperossification
procD.pgls(coords ~ hyperossified, phy = phy, data = gdf, iter=10000)

#Phylogenetic regression
all <- procD.allometry(coords ~ Csize,logsz = TRUE, data=gdf, iter=10000, RRPP=TRUE, print.progress = FALSE)
PhyloReg<-procD.pgls(coords ~ Csize, phy = phy, data = gdf, iter=10000)
summary(PhyloReg)

#size * hyperossification factor interaction and homogeneity of slopes test 
HOS <- advanced.procD.lm(f1= coords ~ log(Csize) + hyperossified,
                         f2= ~ log(Csize) * hyperossified, groups = ~ hyperossified, phy = tree,
                         slope = ~ log(Csize), angle.type = "deg", iter = 10000, data = gdf)
summary(HOS, formula = FALSE) 


#phylogenetic MANOVAs for microhabitat, feeding biology, & phragmosis

#microhabitat only
micro <- procD.pgls(coords ~ microhabitat, phy = phy, data = gdf, iter=10000)
summary(micro, formula = TRUE)

#microhabitat only & post hoc comparisons
micro_PW.means.test <- advanced.procD.lm(f1= coords ~ 1, f2= ~ microhabitat, 
                                               groups = ~microhabitat, phy=mytree, data=gdf, iter=10000)
summary(micro_PW.means.test, formula = TRUE)

#microhabitat*hyperossification & post hoc comparisons
micro_hyperos_PW.means.test <- advanced.procD.lm(f1= coords ~ microhabitat*hyperossified, f2= ~ microhabitat + hyperossified, 
                                         groups = ~microhabitat*hyperossified, phy=mytree, data=gdf, iter=10000)
summary(micro_hyperos_PW.means.test, formula = TRUE)


#diet only
diet <- procD.pgls(coords ~ diet, phy = phy, data = gdf, iter=10000)
summary(diet, formula = TRUE)

#diet*hyperossification & post hoc comparisons
diet_hyperos_PW.means.test <- advanced.procD.lm(f1= coords ~ diet*hyperossified, f2= ~ diet + hyperossified, 
                                        groups = ~diet*hyperossified, phy=mytree, data=gdf, iter=10000)
summary(diet_hyperos_PW.means.test, formula = TRUE)

#diet only with unknown category
diet_2 <- procD.pgls(coords ~ diet2, phy = phy, data = gdf, iter=10000)
summary(diet_2, formula = TRUE)

#diet only & post hoc comparisons with unknown category
diet_2_PW.means.test <- advanced.procD.lm(f1= coords ~ 1, f2= ~ diet2, 
                                         groups = ~diet2, phy=mytree, data=gdf, iter=10000)
summary(diet_2_PW.means.test, formula = TRUE)

#diet*hyperossification & post hoc comparisons with unknown category
diet2_hyperos_PW.means.test <- advanced.procD.lm(f1= coords ~ diet2*hyperossified, f2= ~ diet2 + hyperossified, 
                                        groups = ~diet2*hyperossified, phy=mytree, data=gdf, iter=10000)
summary(diet2_hyperos_PW.means.test, formula = TRUE)


#phragmosis only
phragmosis <- procD.pgls(coords ~ phragmosis, phy = phy, data = gdf, iter=10000)

#phragmosis*hyperossification
phrag_hyperos<-procD.pgls(coords ~ phragmosis * hyperossified, phy = phy, data = gdf, iter=10000)
summary(phrag_hyperos, formula = TRUE)



##PLOT FIGURE 2
compare <- treedata(mytree, classifier, sort=TRUE)
classifier_sorted <- as.data.frame(compare$data)
tree_sorted <- compare$phy

par(mfrow=c(2,2))

diet.sort <- classifier_sorted$feeding_biology
diet.colors=c('skyblue','sienna1')
microhabitatsort <- classifier_sorted$microhabitat
microhabitat.colors=c('#6A5ACD','#20B2AA','#8B4513','#a6dba0')

plotGMPhyloMorphoSpace(tree_sorted,Y.gpa$coords, tip.labels = FALSE, node.labels = FALSE, ancStates = FALSE, xaxis=1, yaxis=2,
                       plot.param=list(t.bg=microhabitat.colors[as.numeric(microhabitatsort)],t.pch=c(21,22)[as.numeric(classifier_sorted$Hyperossified)],n.cex=0,l.col="gray", lwd=1, txt.cex=0.5))
legend("topleft", legend= levels(microhabitatsort), pch=c(16), col = microhabitat.colors2, 
       bty = "n", cex = 1, pt.cex=1.5)
plotGMPhyloMorphoSpace(tree_sorted,Y.gpa$coords, tip.labels = FALSE, node.labels = FALSE, ancStates = FALSE, xaxis=3, yaxis=2,
                       plot.param=list(t.bg=microhabitat.colors[as.numeric(microhabitatsort)],t.pch=c(21,22)[as.numeric(classifier_sorted$Hyperossified)],n.cex=0.0,l.col="gray", lwd=1, txt.cex=0.5))

plotGMPhyloMorphoSpace(tree_sorted,Y.gpa$coords, tip.labels = FALSE, node.labels = FALSE, ancStates = FALSE, xaxis=1, yaxis=2,
                       plot.param=list(t.bg=diet.colors[as.numeric(diet.sort)],t.pch=c(21,22)[as.numeric(classifier_sorted$Hyperossified)],n.cex=0.0,l.col="gray", lwd=1,txt.cex=0.2))
legend("topleft", legend= levels(diet.sort), pch=c(16), col = diet.colors, 
       bty = "n", cex = 1, pt.cex=1.5)
plotGMPhyloMorphoSpace(tree_sorted,Y.gpa$coords, tip.labels = FALSE, node.labels = FALSE, ancStates = FALSE, xaxis=3, yaxis=2,
                       plot.param=list(t.bg=diet.colors[as.numeric(diet.sort)],t.pch=c(21,22)[as.numeric(classifier_sorted$Hyperossified)],n.cex=0.0,l.col="gray", lwd=1, txt.cex=0.2))

dev.off()
#diet with unknown category
par(mfrow=c(1,2))
diet.sort2 <- classifier_sorted$feeding_biology2
diet.colors2 =c('grey33','white', 'sienna1')

plotGMPhyloMorphoSpace(tree_sorted,Y.gpa$coords, tip.labels = FALSE, node.labels = FALSE, ancStates = FALSE, xaxis=1, yaxis=2,
                       plot.param=list(t.bg=diet.colors2[as.numeric(diet.sort2)],t.pch=c(21,22)[as.numeric(classifier_sorted$Hyperossified)],n.cex=0.0,l.col="gray", lwd=1,txt.cex=0.2))
legend("topleft", legend= levels(diet.sort2), pch=c(16), col = diet.colors2, 
       bty = "n", cex = 1, pt.cex=1.5)
plotGMPhyloMorphoSpace(tree_sorted,Y.gpa$coords, tip.labels = FALSE, node.labels = FALSE, ancStates = FALSE, xaxis=3, yaxis=2,
                       plot.param=list(t.bg=diet.colors2[as.numeric(diet.sort2)],t.pch=c(21,22)[as.numeric(classifier_sorted$Hyperossified)],n.cex=0.0,l.col="gray", lwd=1))


##PLOT FIGURE 3
regscore <- all$Reg.proj
logCS<- log(Y.gpa$Csize)
regression.colors <- c("white","black")
dev.off()
par(mfrow=c(1,1))
plot(logCS,regscore,
     xlab="log(Centroid Size)", ylab ="Shape (Regression Score)",
     col = c("black"),
     bg = regression.colors[as.numeric(classifier$Hyperossified)],
     pch=21,cex=1, axes=TRUE)
legend("topleft", legend= levels(classifier$Hyperossified), pch=c(19,15), col = c("white","black"), 
       bty = "n", cex = 1)
text(logCS,regscore, labels = rownames(classifier), cex=0.3, pos = 4)

#treadline for hyperossified taxa
combined <- cbind(logCS,classifier)
CSyes <- subset(combined,combined$Hyperossified %in% "yes")
x1 <- CSyes$logCS
combined <- cbind(regscore,classifier)
RegscoreYes <- subset(combined,combined$Hyperossified %in% "yes")
y1 <- RegscoreYes$regscore
(slopeOLS <- cov(x1,y1)/var(x1))								
(interceptOLS <- mean(y1)-slopeOLS*mean(x1))			
(rOLS <- slopeOLS*sd(x1)/sd(y1))						
(r_sq <- rOLS^2)
O <- cbind(slopeOLS,interceptOLS,r_sq)
OLS <- as.vector(O)
table <- cbind(OLS)
rownames(table) <- c('Slope', 'Intercept', 'r-squared')
table
abline(a=interceptOLS, b=slopeOLS, lty=1, col="black", lwd=2)
#treadline for non-hyperossified taxa
combined <- cbind(logCS,classifier)
CSNo <- subset(combined,combined$Hyperossified %in% "no")
x1 <- CSNo$logCS
combined <- cbind(regscore,classifier)
RegscoreNo <- subset(combined,combined$Hyperossified %in% "no")
y1 <- RegscoreNo$regscore
(slopeOLS <- cov(x1,y1)/var(x1))								
(interceptOLS <- mean(y1)-slopeOLS*mean(x1))			
(rOLS <- slopeOLS*sd(x1)/sd(y1))						
(r_sq <- rOLS^2)
O <- cbind(slopeOLS,interceptOLS,r_sq)
OLS <- as.vector(O)
table <- cbind(OLS)
rownames(table) <- c('Slope', 'Intercept', 'r-squared')
table
abline(a=interceptOLS, b=slopeOLS, lty=1, col="gray80", lwd=2)

#beeswarm plots
dev.off()
par(mfrow=c(1,3))
library(beeswarm)
microhabitat.colors=c('#6A5ACD','#20B2AA','#8B4513','#a6dba0')
swarm.colors <- c("white","black")
feeding.colors=c('lightblue','orange')
phrag.colors=c('lightblue','red')
outline.colors=c('gray','black')

beeswarm(pc1, method = "hex", vertical = TRUE, cex = 1.5, pch =21, pwcol = outline.colors[as.numeric(classifier$Hyperossified)], pwbg  = microhabitat.colors[as.numeric(classifier$microhabitat)], ylim=c(-0.3,0.3))
bxplot(pc1, probs = c(0.05, 0.5, 0.95), add = TRUE)
beeswarm(-pc2, method = "hex", vertical = TRUE, cex = 1.5, pch = 21, pwcol = outline.colors[as.numeric(classifier$Hyperossified)], pwbg  = feeding.colors[as.numeric(classifier$feeding_biology)],  ylim=c(-0.3,0.3))
bxplot(-pc2, probs = c(0.05, 0.5, 0.95), add = TRUE)
beeswarm(-pc3, method = "hex", vertical = TRUE, cex = 1.5, pch = 21, pwcol = outline.colors[as.numeric(classifier$Hyperossified)], pwbg  = phrag.colors[as.numeric(classifier$phragmosis)], ylim=c(-0.3,0.3))
bxplot(-pc3, probs = c(0.05, 0.5, 0.95), add = TRUE)



## PLOT FIGURE 1 - ASR & trait mapping
dev.off()
par(mfrow=c(1,1))

#plotting contMap of skull shape (PC2)
cols = c("white","black")
obj<-contMap(mytree, pc2, type = "fan", fsize = 0.01, res=1)
obj<-setMap(obj,colors=c("white","black"))
par(oma=c(1.5,1.5,1.5,1.5))
plot(obj, type="fan", fsize = 0.01, res=1000)
#plot(obj, type="fan", fsize = 0.5, res=1000) # larger tip labels


#get node states from RevGadgets
library(RevGadgets)
freeK_RJ_tree_file = "ase_freeK_RJ.tree"
#plot RevBayes ancestral states
freeK_RJ <- plot_ancestral_states(freeK_RJ_tree_file, summary_statistic="MAP",
                                  tip_label_size=0.5,
                                  xlim_visible=NULL,
                                  node_label_size=0,
                                  show_posterior_legend=TRUE,
                                  node_size_range=c(1,3),
                                  alpha=0.75)
plot(freeK_RJ)

#extract node states and posterior probabilities 
node_states<-freeK_RJ$data$anc_state_1[159:314]
pp<-freeK_RJ$data$anc_state_1_pp[159:314]

#set node colors
mycolnodes<-character(length(node_states))
mycolnodes[node_states=="1"]<-"white"
mycolnodes[node_states=="2"]<-"black"

#plot nodes w/ posterior probability
nodelabels(pch=21, col="black", adj=0, bg=mycolnodes, cex = 1.5*as.numeric(pp))

#set tip colors
hyperossified <- classifier_sorted$Hyperossified
mycol<-character(length(classifier_sorted))
mycol[hyperossified=="yes"]<-"black"
mycol[hyperossified=="no"]<-"white"
diet <- classifier_sorted$feeding_biology
microhabitat <- classifier_sorted$microhabitat
mycol2<-character(length(classifier_sorted))
mycol2[diet=="invertebrate_predator"]<-"skyblue"
mycol2[diet=="vertebrate_predator"]<-"sienna1"
mycol3<-character(length(classifier_sorted))
mycol3[microhabitat=="aquatic"]<-"#6A5ACD"
mycol3[microhabitat=="arboreal"]<-"#20B2AA"
mycol3[microhabitat=="fossorial"]<-"#8B4513"
mycol3[microhabitat=="terrestrial"]<-"#a6dba0"

#plot tip labels
tiplabels(pch=21, col="black", adj=1, bg=mycol, cex=1.5)
par(oma=c(0.75,0.75,0.75,0.75))
tiplabels(pch=21, col="black", adj=1, bg=mycol2, cex=1.5)
par(oma=c(0,0,0,0))
tiplabels(pch=21, col="black", adj=1, bg=mycol3, cex=1.5)


###PLOT APPENDIX FIGURES
##PLOT FIGURE S2: contmaps
par(mfrow=c(1,3))
obj<-contMap(mytree, pc1, fsize = 0.5, res=1, outline=FALSE)
obj<-contMap(mytree, pc2, fsize = 0.5, res=1, outline=FALSE)
obj<-contMap(mytree, pc3, fsize = 0.5, res=1, outline=FALSE)

#PLOT FIGURE S5: labeled PCA plots with microhabitat
par(mfrow=c(1,2))
plot(-pc1,-pc2, bg=microhabitat.colors[as.numeric(classifier$microhabitat)], pch=c(21,22)[as.numeric(classifier$Hyperossified)])
text(-pc1,-pc2, labels = rownames(classifier), cex=0.3, pos = 4)
plot(-pc3,-pc2, bg=microhabitat.colors[as.numeric(classifier$microhabitat)], pch=c(21,22)[as.numeric(classifier$Hyperossified)])
text(-pc3,-pc2, labels = rownames(classifier), cex=0.3, pos = 4)

#PLOT FIGURE S7: labeled PCA plots with diet (vertebrate predator, invertebrate predator, unknown diet)
par(mfrow=c(1,2))
diet.sort2 <- classifier$feeding_biology2
diet.colors2 =c('grey33','white', 'sienna1')
plot(-pc1,-pc2, bg=diet.colors2[as.numeric(classifier$feeding_biology2)], pch=c(21,22)[as.numeric(classifier$Hyperossified)])
text(-pc1,-pc2, labels = rownames(classifier), cex=0.3, pos = 4)
plot(-pc3,-pc2, bg=diet.colors2[as.numeric(classifier$feeding_biology2)], pch=c(21,22)[as.numeric(classifier$Hyperossified)])
text(-pc3,-pc2, labels = rownames(classifier), cex=0.3, pos = 4)
