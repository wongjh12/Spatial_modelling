# plots for publishing
# 1) 4 groups of cell types and each subsystems by regions build heatmap based on posthoc kruskal wallis test
# 2) boxplot flux for each region 

# PRE-PROCESSING in EXCEL (muresults)
**** add subsystem col for missing ones in excel
Exchange
Transport, intercellular
Biomass
**** Change NA to 0 and remove all EX_ rxns and comm reactions

##### PLOTS IN R
# First, manipulate data in R
# upload muresultsS_ as data
values <-data[c(3:length(data))]
allvalues <- data.frame(col1 = unlist(values, use.names = FALSE))

#length(data)-2
allsubsys <- data.frame(lapply(data[2], function(x) paste(x,'_R1')))
allsubsys <- rbind(allsubsys,data.frame(lapply(data[2], function(x) paste(x,'_R2'))))
allsubsys <- rbind(allsubsys,data.frame(lapply(data[2], function(x) paste(x,'_R3'))))
allsubsys <- rbind(allsubsys,data.frame(lapply(data[2], function(x) paste(x,'_R4'))))
allsubsys <- rbind(allsubsys,data.frame(lapply(data[2], function(x) paste(x,'_R5'))))
allsubsys <- rbind(allsubsys,data.frame(lapply(data[2], function(x) paste(x,'_R6'))))
#allsubsys <- rbind(allsubsys,data.frame(lapply(data[2], function(x) paste(x,'_R7'))))
#allsubsys <- rbind(allsubsys,data.frame(lapply(data[2], function(x) paste(x,'_R8'))))

# edit no. of R_ based on how many regions in slide
region <- data.frame(rep(c('R1','R2','R3','R4','R5','R6'),each=nrow(data)))
subsys<- data.frame(rep(c(data[2]),each=length(data)-2))
subsys <- data.frame(col1 = unlist(subsys, use.names = FALSE))

library(stringr)
rxns <- str_remove_all(data[[1]], 'model[0-9]_')

#### change the names of celltype for each region
rxns <-str_replace_all(data[[1]], c('model1' = 'CD19','model2'='CD3','model3'='Tumor'))
reactions<- data.frame(rep(c(rxns),2)) #same cell types for 2 regions
rxns <-str_replace_all(data[[1]], c('model1' = 'CD3','model2'='NE','model3'='Tumor'))
reactions<- rbind(reactions,data.frame(rep(c(rxns),2)))
rxns <-str_replace_all(data[[1]], c('model1' = 'CD3','model2'='Tumor'))
temp <- data.frame(c(rxns))
names(temp) <- names(reactions) 
reactions<- rbind(reactions,temp)
rxns <-str_replace_all(data[[1]], c('model1' = 'Tumor'))
temp <- data.frame(c(rxns))
names(temp) <- names(reactions) 
reactions<- rbind(reactions,temp)
rxns <-str_replace_all(data[[1]], c('model1' = 'CD3','model2'='Tumor'))
temp <- data.frame(c(rxns))
names(temp) <- names(reactions) 
reactions<- rbind(reactions,temp)
####

# combbine into dataframe
allvalues <- cbind(allsubsytems = c(allsubsys), allvalues);
allvalues <- cbind(regions = c(region), allvalues);
allvalues <- cbind(subsytems = subsys, allvalues);
allvalues <- cbind(reactions = reactions, allvalues);

# add headers
colnames(allvalues)[1] <- "Reactions"
colnames(allvalues)[2] <- "Subsystems"
colnames(allvalues)[3] <- "Regions"
colnames(allvalues)[4] <- "Both"

allvalues['col1']<- abs(allvalues['col1'])
X<-split(allvalues, allvalues$Subsystems)

# compare each subsystem for all regions 
kresults = list()
for (i in 1:length(X)){
ktest = kruskal.test(col1 ~ Regions,
  data = X[[i]]
)
kresults[[i]] <- ktest['p.value']
}

kresults <- unlist(kresults, recursive = FALSE)
temp <- which(kresults<0.05)

names(X) = gsub('/',', ',names(X))
result <- data.frame(name = names(X)[temp],pvalue = unlist(kresults)[temp],index=temp)
result <- result[order(result$pvalue,decreasing=FALSE),]

############ plot heatmap
library(ggstatsplot)
library(ggplot2)
library(RColorBrewer)
source("heatmap.2.mod.R")

for (name in result$name[1:10]) {
tiff(filename=paste(name,"heatmap.tiff"), units="in", width=10, height=10, res=600)
temp <- X[[which(names(X)==name)]]
temp <- temp[, !colnames(temp) %in% c("Subsystems", "Both")]
temp[temp==0]<-NA  ### add in bc i changed NA to 0 in excel muresults

X2<-split(temp, temp$Reactions)
test<-do.call("rbind",X2)
test<-reshape(test, idvar = "Reactions", timevar = "Regions", direction = "wide")
rownames(test)=test[,1]
test<- test[,-1]

# edit no. of R_ based on how many regions in slide
colnames(test)=c('R1','R2','R3','R4','R5','R6')
test<-test[rowSums(is.na(test)) != ncol(test), ]
test[is.na(test)] <- 0

temp<-sapply(str_split(rownames(test),"_"), `[`, 1)
temp<-table(temp)
rowcol = rep(c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072"), temp)
colcol = brewer.pal(7, 'Set1')
par(cex.main=2, mar=c(0,0,10,0.5))

heatmap.2.mod(data.matrix(test), scale="row",
              ,trace="none", Rowv=FALSE, dendrogram="column",col=brewer.pal(11,"PRGn")[0:11], margins=c(5,10)
              , RowSideColors=rowcol,lhei=c(1.5,10),lwid = c(1.5,8),labRow = FALSE
              ,cexCol=2,key.title="",key.ylab="",keysize=10,cex.lab=2,key.par=list(mgp=c(1.5, 0.5, 0)))
legend("right", title = "Cell Type",legend=c("CD19","CD3","NE","Tumor"), 
       fill=c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072"), cex=1.2, box.lty=0)
dev.off()
}

############ Plot boxplot
# replace nan with 0
Xnew <- lapply(X, function(x) replace(x,is.na(x),0))

# manually change the i (dont know why cannot loop :<)
i=
index = result$index[i]
name = result$name[i]
tiff(filename=paste(name,"boxplot.tiff"), units="in", width=18, height=28, res=600)
ggbetweenstats(
    data = Xnew[[index]],
    x = Regions,
    y = col1,
    type = "nonparametric", # ANOVA or Kruskal-Wallis
    plot.type = "box",
    pairwise.comparisons = TRUE,
    pairwise.display = "significant",
    centrality.plotting = FALSE,
    bf.message = FALSE,
    palette='Set1',
    ylab='Flux',
    xlab='',
    ggsignif.args = list(textsize = 10, tip_length = 0.01)
)+theme(axis.title = element_text(size = 40),axis.text = element_text(size = 35),plot.title = element_text(size = 40),plot.subtitle = element_text(size = 20),axis.title.y = element_text(vjust=3),plot.margin = unit(c(0,0,0,1.5), "cm"),axis.text.x = element_text(face="bold"))+ ggplot2::labs(subtitle = NULL)
dev.off() 
