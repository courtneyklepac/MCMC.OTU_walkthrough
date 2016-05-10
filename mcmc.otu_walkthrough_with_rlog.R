#####################################
# installing the package: execute this line just once when first using this script:
install.packages("MCMC.OTU")
install.packages("MCMC.qpcr") # for trellisByGene function

#####################################
# Looking for differentially represented OTUs 

library(MCMC.OTU)
library(MCMC.qpcr)

# loading Green et al data (comes with the package). 
# To load your own data, use functions read.table or read.csv
# rows - samples, columns - OTUs
# MAKE SURE THERE IS COLUMN NAMED sample LISTING ALL SAMPLE NAMES!
data(green.data)

# removing low-count samples and OTUs (?purgeOutliers to see how to adjust cutoffs)
# also adjust otu.columns value for your dataset
goods=purgeOutliers(green.data,count.columns=4:156,otu.cut=0.001)

# >>>>>>> this is how the format of the data table should look 
# >>>>>>> note the 'sample' column giving a unique name for each sample - this one is required
# >>>>>>> (row names are not required)
head(goods)

# what is the proportion of samples with data for these OTUs?
withData=apply(goods[,4:length(goods[1,])],2,function(x){sum(x>0)/length(x)})
hist(withData)

# what percentage of total counts each OTU represents?
props=apply(goods[,4:length(goods[1,])],2,function(x){sum(x)/sum(goods[,4:length(goods[1,])])})
barplot(sort(props,decreasing=T),xaxt="n",log="y")

# stacking the data; adjust otu.columns and condition.columns values for your data
gss=otuStack(goods,count.columns=c(4:length(goods[1,])),condition.columns=c(1:3))

# fitting the model. Replace the formula specified in 'fixed' with yours, add random effects if present. 
# See ?mcmc.otu for these and other options. 
mm=mcmc.otu(
	fixed="bank+species+bank:species",
	data=gss,
	nitt=55000,thin=50,burnin=5000 # a long MCMC chain to improve modeling of rare OTUs
	)

# selecting the OTUs that were modeled reliably
# (OTUs that are too rare for confident parameter estimates are discarded) 
acpass=otuByAutocorr(mm,gss)

# calculating differences and p-values between all pairs of factor combinations
smm0=OTUsummary(mm,gss,otus=acpass,summ.plot=FALSE) 

# adjusting p-values for multiple comparisons:
smmA=padjustOTU(smm0)

# significant OTUs at FDR<0.05:
sigs=signifOTU(smmA)
sigs

# plotting the significant ones
smm1=OTUsummary(mm,gss,otus=sigs)

# now plotting them by species
smm1=OTUsummary(mm,gss,otus=sigs,xgroup="species")

# trellis by OTU (some massaging to the summary object first to make trellisByGene eat it):
smg=smm1
names(smg$summary)[1]="gene"
trellisByGene(smg,xFactor="bank",groupFactor="species",nrow=1)

# table of log10-fold changes and p-values: this one goes into supplementary info in the paper
smmA$otuWise[sigs]

#############################################
# Principal coordinate analysis - technocally not part of MCMC.OTU package except for the logLin function
 
library(vegan)
library(MCMC.OTU)
data(green.data)

# purging under-sequenced samples and OTUs represented in less than 10% of all samples 
goods2=purgeOutliers(green.data,count.columns=4:156,zero.cut=0.1,otu.cut=0)

# creating a log-transfromed normalized dataset for PCoA:
goods.log=logLin(data=goods2,count.columns=4:length(names(goods2)))

# computing Manhattan distances (sum of all log-fold-changes) and performing PCoA:
goods.dist=vegdist(goods.log,method="manhattan")
goods.pcoa=pcoa(goods.dist)

# plotting by bank:
scores=goods.pcoa$vectors
conditions=goods2[,1:3]
margin=0.01
quartz()
plot(scores[,1], scores[,2],type="n",
	xlim=c(min(scores[,1])-margin,max(scores[,1])+margin),
	ylim=c(min(scores[,2])-margin,max(scores[,2])+margin),
	mgp=c(2.3,1,0),
	xlab=paste("Axis1 (",round(goods.pcoa$values$Relative_eig[1]*100,1),"%)",sep=""),
	ylab=paste("Axis2 (",round(goods.pcoa$values$Relative_eig[2]*100,1),"%)",sep=""),
	main="PCoA colored by bank")
points(scores[conditions$bank=="east",1],scores[conditions$bank=="east",2])
points(scores[conditions$bank=="west",1],scores[conditions$bank=="west",2],pch=19)
legend("bottomright", c("East","West"), pch=c(1, 19), cex=0.8)

# plotting by species:
margin=0.01
conditions=goods2[,1:3]
plot(scores[,1], scores[,2],type="n",
	xlim=c(min(scores[,1])-margin,max(scores[,1])+margin),
	ylim=c(min(scores[,2])-margin,max(scores[,2])+margin),
	mgp=c(2.3,1,0),
	xlab=paste("Axis1 (",round(goods.pcoa$values$Relative_eig[1]*100,1),"%)",sep=""),
	ylab=paste("Axis2 (",round(goods.pcoa$values$Relative_eig[2]*100,1),"%)",sep=""),
	main="PCoA colored by species")
points(scores[conditions$species=="franksi",1],scores[conditions$species=="franksi",2])
points(scores[conditions$species=="faveolata",1],scores[conditions$species=="faveolata",2],pch=19)
legend("bottomright", c("franksi","faveolata"), pch=c(1, 19), cex=0.8)


#############################################
# exploring correlations between OTUs 

library(MCMC.OTU)
data(green.data)
goods=purgeOutliers(green.data,count.columns=4:156,otu.cut=0.001,zero.cut=0.2)

# creating a log-transformed normalized dataset ignoring zero counts:
nl=startedLog(data=goods,count.columns=4:length(names(goods)),logstart=0)
names(nl)=c("I","II","III","IV","V","VI","VII")

# displaying a matrix of scatterplots and p-values of OTU correlations
# (onlu p-values better than 0.1 are displayed)
pairs(nl,lower.panel=panel.smooth,upper.panel=panel.cor.pval)

#----------------------------------------


