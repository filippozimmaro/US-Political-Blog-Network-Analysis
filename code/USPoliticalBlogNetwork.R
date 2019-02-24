library(igraph)
setwd("/Users/shichaoma/Google Drive/Job Market/Non-academic/InsightInterview/")

# US political polblogs network 
# Adamic and Glance 2005: The political blogosphere and the 2004 US Election
rm(list=ls())
load("polblogs.Rdata")
V(polblogs)$color = c('blue','red')[(V(polblogs)$LeftRight + 1)]
# V(polblogs)$size <- log(degree(polblogs))

plot(polblogs, 
     vertex.label = NA,
     vertex.size = 4,
     vertex.frame.color = 'gray',
     edge.arrow.size = 0.1,
     edge.width = .1,
     margin = 0)

matrix.polblogs = get.adjacency(polblogs,sparse=F)
diag(matrix.polblogs) = 0
LeftRight = V(polblogs)$LeftRight

# delete those with no incoming or outgoing relations
matrix.polblogs.small = matrix.polblogs
LeftRight.small = LeftRight
threshold = 1
while(any(apply(matrix.polblogs.small,1,"sum") < threshold) || any(apply(matrix.polblogs.small,2,"sum") < threshold)){
  inf.Index = (apply(matrix.polblogs.small,1,"sum") < threshold) | (apply(matrix.polblogs.small,2,"sum") < threshold)
  matrix.polblogs.small = matrix.polblogs.small[!inf.Index,!inf.Index]
  LeftRight.small = LeftRight.small[!inf.Index]
}

cat("n =",nrow(matrix.polblogs.small))
polblogs.small = graph.adjacency(matrix.polblogs.small)
V(polblogs.small)$LeftRight = LeftRight.small
V(polblogs.small)$color = c('blue','red')[(V(polblogs.small)$LeftRight + 1)]

plot(polblogs.small, 
     vertex.label = NA,
     vertex.size = 4,
     vertex.frame.color = 'gray',
     edge.arrow.size = 0.1,
     edge.width = .1,
     margin = 0)


save(matrix.polblogs.small, polblogs.small, LeftRight.small, 
     file = "polblogsSmall.Rdata")

# summary statistics
outgoing = apply(matrix.polblogs.small,1,"sum")
incoming = apply(matrix.polblogs.small,2,"sum")
mean(outgoing)
mean(incoming)
min(outgoing)
min(incoming)
max(outgoing)
max(incoming)
median(outgoing)
median(incoming)
sum(outgoing)
sum(incoming)
sum(outgoing>=200)
sum(incoming>=200)
sum(outgoing>=100)
sum(incoming>=100)
sum(outgoing>=50)
sum(incoming>=50)