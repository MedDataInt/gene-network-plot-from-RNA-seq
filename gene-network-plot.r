# library
library(igraph)
library(network)
library(sna)
library(ndtv)
library(visNetwork)
library(extrafont)


install.packages('sna')
install.packages('ndtv')
install.packages('visNetwork')
install.packages('extrafont')

par(mfrow = c(1,1))
nodes <- read.csv("node2.csv", header=T, as.is=T)
links <- read.csv("edge2.csv", header=T, as.is=T)

## here are the two dataframe information 
# cols of edge2: "from", "to", "weight", and "type"
# cols of node2: "id", "gene","gene.type",  "type.label", "fold.size"

# "from" and "to" are the "id"s
# "gene.type": 1,2,3
# "type.label": transcription factor, up-regulated, down-regulated
# "fold.size": fold change KO vs.WT
# "weight": -log10(p-value)
# "type": 'hyperlink'

head(nodes)
head(links)
nrow(nodes); length(unique(nodes$id))
nrow(links); nrow(unique(links[,c("from", "to")]))

links <- aggregate(links[,3], links[,-3], sum)
#links <- links[order(links$from, links$to),]

colnames(links)[4] <- "weight"
rownames(links) <- NULL
library('igraph')
net <- graph_from_data_frame(d=links, vertices=nodes, directed=T)
net

E(net) # The edges of the "net" object
V(net) # The vertices of the "net" object
E(net)$type # Edge attribute "type"
V(net)$gene # Vertex attribute "gene"

plot(net) # not a pretty picture!

net <- simplify(net, remove.multiple = F, remove.loops = T)
plot(net, edge.arrow.size=.4,vertex.label=NA)

plot(net, edge.arrow.size=.4, edge.curved=.1)
plot(net, edge.arrow.size=.2, edge.color="orange",
     vertex.color="orange", vertex.frame.color="#ffffff",
     vertex.label=V(net)$gene, vertex.label.color="black", edge.curved=.1)


# Generate colors based on media type:
colrs <- c("gray", "tomato", "gold")
V(net)$color <- colrs[V(net)$gene.type]
# Compute node degrees (#links) and use that to set node size:
deg <- degree(net, mode="all")
V(net)$size <- deg*3
# We could also use the audience size value:
V(net)$size <- V(net)$fold.size*4
# The labels are currently node IDs.
# Setting them to NA will render no labels:
V(net)$label <- NA
# Setting the frame color
V(net)$frame.color <- "white"
# Set edge width based on weight:
E(net)$width <- E(net)$weight/6
#change arrow size and edge color:
E(net)$arrow.size <- .5
E(net)$edge.color <- "gray80"
E(net)$width <- 1+E(net)$weight/8
plot(net)
plot(net, edge.curved=.3)

# legend(x=-1.5, y=-1.1, c("PU.1","upregulated gene", "downregulated gene"), pch=21,
       #col="#777777", pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1)

## alternatively: label gene
plot(net, vertex.label=V(net)$gene,
     vertex.label.font=2, vertex.label.color="gray40",
     vertex.label.cex=1.2, edge.color='gray', edge.curved=.3)

# change edge colors
edge.start <- ends(net, es=E(net), names=F)[,1]
edge.col <- V(net)$color[edge.start]
plot(net, edge.color=edge.col, edge.curved=.1)


# Circle layout
l <- layout_in_circle(net)
plot(net, layout=l)

# 3D sphere layout
l <- layout_on_sphere(net)
plot(net, layout=l)

# force-directed layouts
l <- layout_with_fr(net)
plot(net, layout=l)

#  Kamada Kawai
l <- layout_with_kk(net)
plot(net, layout=l)

# LGL algorithm is meant for large, connected graph
plot(net, layout=layout_with_lgl)


# MDS (multidimensional scaling) algorithm
plot(net, layout=layout_with_mds)

# random 
l <- layout_randomly(net.bg)

## others 
layouts <- grep("^layout_", ls("package:igraph"), value=TRUE)[-1]
# Remove layouts that do not apply to our graph.
layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama|tree", layouts)]
par(mfrow=c(3,3), mar=c(1,1,1,1))
for (layout in layouts) {
  print(layout)
  l <- do.call(layout, list(net))
  plot(net, edge.arrow.mode=0, layout=l, main=layout) }

# end

plot(net, layout=l)
