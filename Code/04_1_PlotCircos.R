library(circlize)
library(data.table) 

res <- read.csv('../Results/meta females.csv')

bool_significant <- res$difference_type != 'None'
res <- res[bool_significant,]

df1 = setDT(res[,c(11,4)])[, .N, by = c(names(res[,c(11,4)]))]
# Make the circular plot
pdf(file='../Results/Plots/Fig2_CirclosPlot.pdf', width=10, height=10)

#chordDiagram(df1, transparency = 0.25, annotationTrack = c("grid"), preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(df1))))), big.gap = 20, small.gap = 3)
chordDiagram(df1, transparency = 0.25, annotationTrack = c("grid"), preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(df1))))))

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(cex=0.75,CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))}, 
              bg.border = NA) # here set bg.border to NA is important
circos.clear()
dev.off()