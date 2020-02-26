##hp
##05.31.2016

setwd('B:/')
source('fcns/config.r')
source('fcns/my_heat_pwys.r')
setwd('kahn/emrah/emrah_stoolProteomics/host')
source('../../my_barplot.r')
library('plyr')

##Parse data
dat <- read.csv('150427KEASAM00828_MUSTMT10PLX-FTK2_UNiprotMouse_Annotation_128Cnormaliza....csv')
dat <- dat[dat$Protein.FDR.Confidence %in% c('High', 'Medium', 'Low') & dat$Contaminant == FALSE, ]
rownames(dat) <- dat$Accession
rownames(dat) <- gsub('-.*', '', rownames(dat))

pheno <- data.frame(Group = c(rep('Chow', 3), rep('HFD', 3), rep('Vanc', 2), rep('Metr', 2)),
                    Replicate = c(1:3, 1:3, 1:2, 1:2),
                    Label = c('126', '127C', '127N', '128N', '129C', '128C', '129N', '130C', '130N', '131'))
rownames(pheno) <- paste(pheno$Group, pheno$Replicate, sep = '.')

annot <- dat[, c('Accession', 'Gene.ID', 'Description')]
annot$Gene.ID[annot$Gene.ID == ''] <- annot$Accession[annot$Gene.ID == '']
dat.tissue <- read.csv("B:/Annotations/encode_mouse_19_tissues_expr/mouse_tissue_expr.csv", row.names = 1)
annot <- cbind(annot, dat.tissue[match(annot$Gene.ID, dat.tissue$Symbol), paste0("Tissue.", 1:3)])
mf <- read.csv("B:/Annotations/msigdb/msigdb_mf_table.csv", row.names = 1)
annot$Molecular.Function <- mf[toupper(annot$Gene.ID), ]
signalIP <- read.delim("B:/Annotations/protein_secretion/hpa_protein-class_SignalPPredictedSecreted.tab")
signalIP <- signalIP[grep('Predicted secreted proteins', signalIP$Protein.class), ]
annot$SignalPPredictedSecreted <- toupper(annot$Gene.ID) %in% signalIP$Gene
annot$SignalPPredictedSecreted[grep(';', annot$Gene.ID)] <- NA

mat <- data.matrix(dat[, grep('Abundance..F1', colnames(dat))])
colnames(mat) <- gsub('.*F1\\.\\.|\\.\\.(Sam|Con).*', '', colnames(mat))
mat <- mat[, match(pheno$Label, colnames(mat))]
colnames(mat) <- rownames(pheno)
# filter NAs
hist(rowMeans(is.na(mat)))
mat <- mat[rowMeans(is.na(mat)) < 0.8, ]
sum(is.na(mat))
# normalize to column sum and log transform
mat <- log2(apply(mat, 2, function(v) v*max(colSums(mat))/sum(v)))
boxplot(mat, las=2)

##Qc
multi.pca(mat, pheno, color = 'Group', name = 'stoolproteomics_pca')
# see trend, flat
voom(2^mat, plot = TRUE)

##Stats
contr.v = c(HFDvsChow = 'HFD - Chow', VancvsHFD = 'Vanc - HFD',
            MetrvsHFD = 'Metr - HFD', MetrvsVanc='Metr - Vanc')
mtt <- limma.contrasts(mat, pheno$Group, contr.v)
mtt.df <- data.frame(signif(mtt, 3), annot[rownames(mtt), ])
sig.hist(mtt, name = 'stoolproteomics_sig_hist')
write.csv(mtt.df, 'stoolproteomics_protein_stats.csv', na = '')

##Plots
# volcanoe plots
mtt.df.vol <- mtt.df
mtt.df.vol$Gene.ID <- gsub(';.*', '', mtt.df.vol$Gene.ID)
pdf('stoolproteomics_volcanoes.pdf', 4, 4)
multi.volcano(mtt.df.vol, sym.col = 'Gene.ID', n.top = 5, name = NA)
dev.off()

# heat maps
topp1 <- rownames(mtt)[1:50]
my.heat(mat[topp1, ], symbols = annot[topp1, 'Gene.ID'], name = 'stoolproteomics_all_group',
        Colv = FALSE, margins = c(7, 10), labCol = pheno$Group)
topp2 <- rownames(mtt)[order(combine.pvalues(mtt[, -grep('Chow', colnames(mtt))]))][1:50]
my.heat(mat[topp2, pheno$Group != 'Chow'], symbols = annot[topp2, 'Gene.ID'], 
        name = 'stoolproteomics_no_chow', Colv = FALSE, margins = c(7, 10), 
        labCol = pheno$Group[pheno$Group != 'Chow'])
topp3 <- rownames(mtt)[order(combine.pvalues(mtt[, -grep('Vanc|Metr', colnames(mtt))]))][1:50]
my.heat(mat[topp3, !pheno$Group %in% c('Vanc', 'Metr')], symbols = annot[topp3, 'Gene.ID'], 
        name = 'stoolproteomics_no_antibiotics', Colv = FALSE, margins = c(7, 10), 
        labCol = pheno$Group[!pheno$Group %in% c('Vanc', 'Metr')])

top4 <- rownames(mtt)[order(mtt$HFDvsChow.p)][1:50]
top5 <- rownames(mtt)[order(mtt$VancvsHFD.p)][1:10]
top6 <- rownames(mtt)[order(mtt$MetrvsHFD.p)][1:10]
topp <- unique(c(top4, top5, top6))
my.heat2(mat[topp, ], symbols = mtt.df.vol[topp, 'Gene.ID'], name = 'stoolproteomics_all_groups_2', Colv = FALSE, labCol = pheno$Group, sc = 'z')

# bar plots
pdf('stoolproteomics_topprotein_barplots.pdf')
for(p in topp1) {
    dat2p <- data.frame(Protein = 2^(mat[p, ])/mean(2^mat[p, 4:6]), Group = pheno$Group)
    dat2p <- ddply(dat2p, .(Group), summarise, N = length(Protein), Mean = mean(Protein),
                   SE = sd(Protein)/sqrt(N))
    
    g <- ggplot(dat2p, aes(Group, Mean, fill = Group)) 
    g <- g + geom_bar(stat = 'identity', width = 0.7) + guides(fill=FALSE) 
    g <- g + theme_bw() + geom_errorbar(aes(ymax = Mean + SE, ymin = Mean - SE, width = 0.35))
    g <- g + labs(y = 'Relative abundance', title = paste(annot[p, 'Gene.ID'], p, sep = '/'))
    plot(g)
}
dev.off()

grp <- factor(pheno$Group)
levels(grp) <- list(C = "Chow", H = "HFD", M = "Metr", V = "Vanc")
pdf("stoolproteomics_topprotein_barplots_2.pdf", 2.5, 3)
my.barplot(2^mat/rowMeans(2^mat[, 4:6]), grp, 
           main.v = annot[rownames(mat), "Gene.ID"], fill = c("#CCCCCC","#FF6666", "#6699FF","#66CC66"), 
           xlab = "", ylab="Relative Abundance")
dev.off()

##Pwys
# pdb.files <- c(c2.kegg = 'c2.cp.kegg.v4.0', c2.reactome = 'c2.cp.reactome.v4.0', c5.cc = 'c5.cc.v5.1')
# for (pdb.ind in seq_along(pdb.files)){
#     G <- gmtToG(paste0('B:/Annotations/gene_sets/', pdb.files[pdb.ind], '.symbols.gmt'))
#     G <- processG(G)
#     G <- my.chip2chip(G = G, annot = annot, split='; ', symbol.col = 'Gene.ID')
#     sp <- sigpwy.contrasts(G=G, tab=mat, grp=pheno$Group, contrasts.v=contr.v, minNPS=2,
#                            annot=annot[rownames(mat), ], prefix = names(pdb.files)[pdb.ind],
#                            add.comb=FALSE)
#     write.csv(sp, paste0(names(pdb.files)[pdb.ind], '_pwys.csv'), row.names=FALSE)
# }

# fry
pdb.files <- c(cp = 'c2.cp.v6.0', go = 'c5.all.v6.0')
for (i in seq_along(pdb.files)){
    G <- gmtToG(paste0('B:/Annotations/gene_sets/', pdb.files[i], '.symbols.gmt'))
    G <- processG(G)
    G <- my.chip2chip(G = G, annot = annot, split='; ', symbol.col = 'Gene.ID')
    pwys <- fry.contrasts(dat=mat, G=G, stats=mtt.df, name=names(pdb.files)[i], grp=pheno$Group, contrasts.v=contr.v, min.ngenes=2)
    
    my.heat.pwys(dat=pwys[1:min(50, nrow(pwys)), ], contrast.v=contr.v, p.cutoff=0.05, name=paste0(names(pdb.files)[i], '_pwys'), Colv=FALSE)
}

##Check
all(rownames(mat)%in% rownames(annot))
all(colnames(mat) == rownames(pheno))
mtt.df[1, 'Accession'] == 'Q7TSG5-2'
mtt.df[1, 'HFDvsChow.FC'] < -2.5
mtt.df[which(mtt.df$Gene.ID == 'Scgb1b2'), 'VancvsHFD.logFC'] > 0.5
tapply(mat['Q6DYE8', ], pheno$Group, mean)[c('Chow', 'HFD')] #17.9, 16.9
