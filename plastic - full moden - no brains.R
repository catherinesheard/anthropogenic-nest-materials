x<-read.csv("NestMaterialsData.csv")
x$PhyloName<-gsub(" ", "_",x$PhyloName)
x<-x[-which(x$duplicate==1),]
library(MCMCglmm)

trees<-read.tree("100birdtrees.tre") #100 trees from the Jetz et al. 2012 Global Bird Tree
t100<-trees[1:100] #this step isn't necessary if your code runs perfectly, but it's nice to save a copy and then trim is rather than have to re-load the tree if something goes wrong

tree<-t100[[1]] #select one tree for trimming purposes

#prepare the data
if(sum(is.na(x$PhyloName)>0)){x<-x[-which(is.na(x$PhyloName)),]}
if(sum(is.na(x$Plastic)>0)){x<-x[-which(is.na(x$Plastic)),]}
if(sum(is.na(x$Number.Mat)>0)){x<-x[-which(is.na(x$Number.Mat)),]}
if(sum(is.na(x$BodyMass)>0)){x<-x[-which(is.na(x$BodyMass)),]}
if(sum(is.na(x$str.none)>0)){x<-x[-which(is.na(x$str.none)),]}
if(sum(is.na(x$loc.artificial)>0)){x<-x[-which(is.na(x$loc.artificial)),]}
if(sum(is.na(x$RangeSize)>0)){x<-x[-which(is.na(x$RangeSize)),]}
if(sum(is.na(x$IUCN.simplified)>0)){x<-x[-which(is.na(x$IUCN.simplified)),]}
if(sum(is.na(x$HFI)>0)){x<-x[-which(is.na(x$HFI)),]}
if(sum(is.na(x$Biome.simplified)>0)){x<-x[-which(is.na(x$Biome.simplified)),]}
if(sum(is.na(x$WoS.Papers)>0)){x<-x[-which(is.na(x$WoS.Papers)),]}

t100<-lapply(t100,drop.tip,tip=setdiff(tree$tip.label,x$PhyloName)) #trim out everything from the tree that's not in the dataset

#t100 and x should now match


######
#prepare the data
######

x$zMat<-scale(x$Number.Mat)
x$zMass<-scale(x$BodyMass)
x$zRange<-scale(x$RangeSize)
x$zHFI<-scale(x$HFI)
x$zWoS<-scale(sqrt(x$WoS.Papers))

######
#set up a dummy run
######

i=1 #this is arbitrary

tree<-t100[[i]]  

animalA<-inverseA(tree)$Ainv 

#set up priors for this model
gelmanprior<-list(B=list(mu=rep(0,26), #mu has to be the number of parameters to be estimated (k-1 for each category of size k, 1 for each continuous variable, 1 for the intercept) 
                         V=gelman.prior(~zWoS+zMat+zMass+zRange+zHFI+as.factor(IUCN.simplified)+as.factor(Biome.simplified)
                                        +as.factor(str.scrapenone)+as.factor(str.cup)+as.factor(str.platform)+as.factor(str.dome)+as.factor(str.ex)
                                        +as.factor(loc.artificial)+as.factor(loc.ground.water)+as.factor(loc.cavity)+as.factor(loc.rock)+as.factor(loc.veg), #you'll have to adjust this equation to match your linear model, and mu above may have to change
                                        data = x,  scale=1+pi^2/3)), 
                  R=list(V=1,fix=1),G=list(G1=list(V=1E-10,nu=-1)))

mod<-MCMCglmm(as.factor(Plastic)~zWoS+zMat+zMass+zRange+zHFI+as.factor(IUCN.simplified)+as.factor(Biome.simplified)
              +as.factor(str.scrapenone)+as.factor(str.cup)+as.factor(str.platform)+as.factor(str.dome)+as.factor(str.ex)
              +as.factor(loc.artificial)+as.factor(loc.ground.water)+as.factor(loc.cavity)+as.factor(loc.rock)+as.factor(loc.veg),
              random=~PhyloName, 
              ginverse=list(PhyloName=animalA), 
              prior = gelmanprior, 
              verbose=TRUE, 
              family="categorical", 
              data = x,
              nitt=5500*2,
              thin=10,
              burnin=500*2,
              pl=TRUE, 
              pr=TRUE, 
              slice=TRUE) 



Final.disp<-mod #set up a structure that we'll populate with the full model
Final.disp$VCV[((i-1)*10+1):(i*10), ]<-mod$VCV[1:10,] 
Final.disp$Sol[((i-1)*10+1):(i*10), ]<-mod$Sol[1:10,] 
Final.disp$Liab[((i-1)*10+1):(i*10), ]<-mod$Liab[1:10,] 

nsamp.l<-nrow(mod$VCV)
start1.l=list(R=mod$VCV[nsamp.l,"units"], G=list(G1=mod$VCV[nsamp.l,"PhyloName"]))

save(Final.disp,file="plastic-full.Rdata")

######
#run the full model
######

for(i in 1:100){ #loop through 100 trees
  tree<-t100[[i]] #select the ith tree 
  
  animalA<-inverseA(tree)$Ainv 
  
  mod<-MCMCglmm(as.factor(Plastic)~zWoS+zMat+zMass+zRange+zHFI+as.factor(IUCN.simplified)+as.factor(Biome.simplified)
                +as.factor(str.scrapenone)+as.factor(str.cup)+as.factor(str.platform)+as.factor(str.dome)+as.factor(str.ex)
                +as.factor(loc.artificial)+as.factor(loc.ground.water)+as.factor(loc.cavity)+as.factor(loc.rock)+as.factor(loc.veg),
                random=~PhyloName, 
                ginverse=list(PhyloName=animalA), 
                prior = gelmanprior, 
                verbose=FALSE, 
                family="categorical", 
                start= start1.l,
                data = x,
                nitt=30000, #this combo of nitt/thin/burnin will result in 10 samples / tree
                thin=2000,
                burnin=10000,
                pl=TRUE,
                pr=TRUE,
                slice=TRUE)
  
  print(i) #print which tree you're on (for your enjoyment as the loop runs)
  
  Final.disp$VCV[((i-1)*10+1):(i*10), ]<-mod$VCV[1:10,]  #put this tree's run into the overall structure
  Final.disp$Sol[((i-1)*10+1):(i*10), ]<-mod$Sol[1:10,] 
  Final.disp$Liab[((i-1)*10+1):(i*10), ]<-mod$Liab[1:10,] 
  
  nsamp.l<-nrow(mod$VCV)
  start1.l=list(R=mod$VCV[nsamp.l,"units"], G=list(G1=mod$VCV[nsamp.l,"PhyloName"]))
  
  #save(Final.disp,file="plastic-full.Rdata") 
  
}

save(Final.disp,file="plastic-full.Rdata")

