```library("geepack")

Metadata = read.table("for_stats.csv", header = TRUE, sep= ",") 

fit3 <- geeglm(persistant ~ Sex+AgeGroup+Size+PEN+TET+CTRIX+CIPRO+AZITH,
               data = Metadata,
               id = Cluster,
               family = binomial,
               corstr = "independence")
summary(fit3)```
