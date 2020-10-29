
library(tidyverse)
library(viridis)
library(patchwork)
library(ggraph)
library(igraph)

# create a data frame giving the hierarchical structure of your individuals. 
my_hierarchy <- read.csv('cluster_478_high_hierarchy.csv')
my_connect <- read.csv('cluster_478_high_connections.csv')

# create a vertices data.frame. One line per object of our hierarchy, giving features of nodes.
my_vertices <- data.frame(name = unique(c(as.character(my_hierarchy$from), as.character(my_hierarchy$to))) )

# Preparation to draw labels properly:
my_vertices$id=NA
my_nleaves=length(my_vertices$name)
my_myleaves=seq(1:my_nleaves)
my_vertices$id[ my_myleaves ] = seq(1:my_nleaves)
my_vertices$angle= 90 - 360 * my_vertices$id / my_nleaves
my_vertices$hjust<-ifelse( my_vertices$angle < -90, 1, 0)
my_vertices$angle<-ifelse(my_vertices$angle < -90, my_vertices$angle+180, my_vertices$angle)

# Create a graph object with the igraph library
my_mygraph <- graph_from_data_frame( my_hierarchy, vertices=my_vertices )

# The connection object must refer to the ids of the leaves:
my_from <- match( my_connect$from, my_vertices$name)
my_to <- match( my_connect$to, my_vertices$name)


# plot
#ggraph(my_mygraph, layout = 'dendrogram', circular = TRUE) + 
#  geom_node_text(aes(x = x*1.01, y=y*1.01, filter = leaf, label=my_vertices$name, angle = my_vertices$angle), size=1.5, alpha=1) +
#  geom_conn_bundle(data = get_con(from = my_from, to = my_to), alpha=0.2, colour="skyblue", tension = 0) + 
#  geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05)) +
#  theme_void()
ggraph(my_mygraph, layout = 'dendrogram', circular = TRUE) + 
  geom_conn_bundle(data = get_con(from = my_from, to = my_to), alpha = 0.2, colour="skyblue", tension = 0) + 
  geom_node_text(aes(x = x*1.01, y=y*1.01, filter = leaf, label=id, angle = angle, hjust=hjust), size=1.5, alpha=1) +
  coord_fixed() +
  theme_void() +
  theme(
    legend.position="none",
    plot.margin=unit(c(0,0,0,0),"cm"),
  ) +
  expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2))
