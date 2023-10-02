### Authors: Angelo Soto-Centeno & Pedro Ivo Monico
# date: 14 of September, 2023
##########################################################################################

#Run packages 
library(dismo)
library(tidyverse)
library(spThin)
library(maptools)
library(sf)
library(rgdal)
library(raster)
library(rJava)
library(ENMeval)
library(tidyverse)
library(sf)
library(ecospat)
library(maxnet)
library(viridis)
library(ENMTools)

#Set working directory for the relevant geographical partition dataset 
setwd("/Users/Pedro/Desktop/Eptesicus_fuscus/data/PAPER.MODELS/caribbean")

# Upload bio climatic variables at the appropriated resolution
raster_files <- list.files("/Users/Pedro/Projects/PhD/WorldClimateData/wc2-5", full.names = T, pattern = ".bil") 

predictors<-stack(raster_files) 

sp1 <- read.csv(file = "./thin_data.csv", header = T)  

#check dataset information
head(sp1)
summary(sp1)

# get map data
data("wrld_simpl")
plot(wrld_simpl, xlim = c(-116, -65), ylim = c(17,33), axes = TRUE, col = NA)

# add species observation localities
points(x = sp1$lon, y = sp1$lat, col = "red", pch = 21, cex = 0.75)

#Check if the predictors are working (in this case predicting just BIO1)
plot(predictors$bio1) 


##########################################################################################

#clipping the predictors 
geo.extent <- extent(-108, -63, 03, 32) 

#Check if the predictors are working
#(in this case predicting just BIO1)
predictors <- crop(predictors, geo.extent)
plot(predictors$bio1) 


##########################################################################################

# create the modeling xy data.frame & verify it
xy <- sp1[c("lon", "lat")] 
summary(xy)


##########################################################################################
# NOTE: the maxent.jar program MUST be within the java directory of the package dismo
##to verify/get the directory of dismo use: 
system.file("java", package = "dismo")

sp.sf <- st_as_sf(xy, coords = c("lon","lat"))
summary(sp.sf) 
class(sp.sf)   
crs(predictors) <- raster::crs(sp.sf) # match the predictor & species CRS.
sp.buf <- sf::st_buffer(sp.sf, dist = 2.25) %>% sf::st_union() %>% sf::st_sf()
plot(predictors[[1]], main = names(predictors)[1])
points(xy)
plot(sp.buf, border = "blue", lwd = 1.75, add = TRUE)
predictors1 <- crop(predictors, sp.buf)
predictors1 <- raster::mask(predictors1, sp.buf)
plot(predictors1)
geo.ext.sqbuff <- extent(min(xy$lon)-25, max(xy$lon)+25, min(xy$lat)-25, max(xy$lat)+25) 

# crop original predictors
predictors2 <- crop(predictors, geo.ext.sqbuff)
plot(predictors2)


##########################################################################################

#create the default model folder & produce the model

dir.create("./maxentdef")

setwd("/Users/Pedro/Desktop/Eptesicus_fuscus/data/PAPER.MODELS/caribbean/maxentdef")

mxnt.dflt <- maxent(predictors1, xy, path = "/Users/Pedro/Desktop/Eptesicus_fuscus/data/PAPER.MODELS/caribbean/maxentdef") #this step will take some time; 
dflt.dist <- predict(mxnt.dflt, predictors2, progress = 'text') #this step will take some time; 
plot(dflt.dist) 

ecospat.boyce(dflt.dist, xy, window.w = "default", res = 100, PEplot = T)
writeRaster(dflt.dist,"ena_default.tif",format="GTiff",NAflag=-99999,overwrite=T)  #exporting the data from default model

getFCs <- function(html) {
  htmlRead <- readLines(html)
  featureTypes <- htmlRead[grep("Feature types", htmlRead)]
  substr(featureTypes, start=21, stop=nchar(featureTypes)-4)
}

def.results <- getFCs(paste("/Users/Pedro/Desktop/Eptesicus_fuscus/data/PAPER.MODELS/caribbean/maxentdef", "/maxent.html", sep = "")) # "MxntDflt/" needs to point to the directory where your default model lives

def.results <- strsplit(def.results, " ")[[1]]
def.results <- lapply(def.results, function(x) gsub("hinge", "H", x))
def.results <- lapply(def.results, function(x) gsub("linear", "L", x))
def.results <- lapply(def.results, function(x) gsub("product", "P", x))
def.results <- lapply(def.results, function(x) gsub("threshold", "T", x))
def.results <- lapply(def.results, function(x) gsub("quadratic", "Q", x))
def.results <- lapply(def.results, function(x) paste(x, collapse = ""))
def.results <- paste(unlist(def.results),collapse = "")
# print the Feature Classes used in the default model
def.results 

##########################################################################################
#create the custom model folder & produce the model

dir.create("/Users/Pedro/Desktop/Eptesicus_fuscus/data/PAPER.MODELS/caribbean/enm_val") 

setwd("/Users/Pedro/Desktop/Eptesicus_fuscus/data/PAPER.MODELS/caribbean/enm_val")


#WHICH ALGORITM FITS BETTER ?
eval.results <- ENMevaluate(occs = xy, envs = predictors,  
                            algorithm = 'maxent.jar', partitions = 'block', 
                            tune.args = list(fc = c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"), rm = 1:3))
#this step will take some time;
### if dismo::maxent is preferred, change algorithm = "maxent.jar" or "maxnet" 


head(eval.results@results)
write.csv(eval.results@results, file = "eval.results.csv", row.names = F)

AICmods <- which(eval.results@results$AICc == min(na.omit(eval.results@results$AICc)))
eval.results@results[AICmods, ] #The result of evaluation is here
def.results #this is default results (just for comparison)

eval.results@results[["fc"]]
eval.results@results[["AICc"]] 
evalplot.stats(e = eval.results, stats = "or.mtp", color = "fc", x.var = "rm") #regularization multiplier  ##GRAPHIC HERE
evalplot.stats(e = eval.results, stats = "auc.val", color = "fc", x.var = "rm") #regularization multiplier  #GRAPHIC HERE


#If "ERROR", try quit R and opening again & DONT run the evaluation

#Switch these arguments with the ones generated by evaluation

setwd("/Users/Pedro/Desktop/Eptesicus_fuscus/data/PAPER.MODELS/caribbean/")
mxnt.best <- maxent(predictors1, xy, args = c("linear=true", "quadratic=true", "product=true", 
                                              "hinge=true", "threshold=false", "betamultiplier=1"), path = "Mxnt_Cust/") 

plot(mxnt.best) # use this to build variable contribution plot; 

response(mxnt.best) ### this is the line to builde the response curves, from bio 1 to 19; 

mxnt.best.dist.log <- predict(mxnt.best, predictors2, args = c("outputformat=cloglog"), progress = "text")   #BEST MODEL (CUSTOM)
  
ecospat.boyce(mxnt.best.dist.log, xy, window.w = "default", res = 100, PEplot = T)  

plot(dflt.dist, main = "Default Model", xlab = "longitude", ylab = "latitude")

plot(mxnt.best.dist.log, main = "Best Model", xlab = "longitude", ylab = "latitude")

writeRaster(mxnt.best.dist.log,"enm_maxbest.tif",format="GTiff",NAflag=-99999,overwrite=T)




####################################### ENM TOOLS

setwd("/Users/Pedro/Desktop/Eptesicus_fuscus/data")

install.packages("devtools")
library(devtools)
install_github("danlwarren/ENMTools")
library(ENMTools)


#rm(predictors)

#creating the geographical projection' 
geo.extent<- extent(-120, -62, 4, 35)

env<- crop(predictors, geo.extent)

env<- check.env(env) #removing NA's; 

plot(env$bio1) # to see how the projection looks like 

summary(is.na(env$bio1)) #here you can confirm if there is any NA; 


############# AMERICA CENTRAL ###############

#creating the object list 
fuscus.central<- read.csv(file = "/Users/Pedro/Desktop/Eptesicus_fuscus/data/ORIGINAL.MODELS/miradorensis/thin_miradorensis/sp1_thin1_mira.csv", header = T)
colnames(fuscus.central)

central <- enmtools.species()
central #object for ENMTools

#presence points
central <- enmtools.species(species.name = "central", 
                            presence.points = fuscus.central[,2:3])

head(central[["presence.points"]])

#range
central$range <- background.raster.buffer(central$presence.points, 150000, mask = env) #this 10000 number i have to play with it; buffers; 
central

#background
central$background.points <- background.points.buffer(points = central$presence.points,
                                                      radius = 120000, n = 1000, mask = env[[1]])    # usually should be closer and highly present? 


central<-check.species(central) #if you find nothing, keep going

#visual
interactive.plot.enmtools.species(central)




############# EAST NORTH AMERICA ###############

#creating the object list 
fuscus.east<- read.csv(file = "/Users/Pedro/Desktop/Eptesicus_fuscus/data/ORIGINAL.MODELS/osceola/thin_osceola/sp1_thin1.csv", header = T)
colnames(fuscus.east)

east <- enmtools.species() # creating the object for ENM_Tools
east 

#records
east <- enmtools.species(species.name = "east", 
                         presence.points = fuscus.east[,2:3])


head(east[["presence.points"]])

#range
east$range <- background.raster.buffer(east$presence.points, 270000, mask = env) 
east

#background 
east$background.points <- background.points.buffer(points = east$presence.points,
                                                   radius = 220000, n = 1200, mask = env[[1]]) #23k reaches The Bahamas; 
# >1200 reaches The Bahamas

#last-check
east<-check.species(east) #if you find nothing, keep going

#visual
interactive.plot.enmtools.species(east) #data organized in a visual way 




############# CARIBBEAN ###############

#creating the object list 
fuscus.caribbean<- read.csv(file = "/Users/Pedro/Desktop/Eptesicus_fuscus/data/ORIGINAL.MODELS/caribbean/caribbean.csv", header = T)
colnames(fuscus.caribbean)
caribbean<- enmtools.species() # creating the object for ENM_Tools
caribbean

#fulfilling the list with the species information:
#records
caribbean <- enmtools.species(species.name = "caribbean", 
                              presence.points = fuscus.caribbean[,1:2])
head(caribbean[["presence.points"]])

#range
caribbean$range <- background.raster.buffer(caribbean$presence.points, 190000, mask = env) 
caribbean

#background
caribbean$background.points <- background.points.buffer(points = caribbean$presence.points,
                                                        radius = 50000, n = 3500, mask = env[[1]]) 

caribbean<-check.species(caribbean) #if you find nothing, keep going; error message could indicate that something is missing or that the data is not correctly wrangled' 

#visual
interactive.plot.enmtools.species(caribbean) #data organized in a visual way 




##################### FOR THE ENMTOOL MODELS #####################

## Those are the raster for the produced ENM models; Upload the produced .tif files for each .mx object. 

east.mx<-raster("/Users/Pedro/Desktop/Eptesicus_fuscus/data/ORIGINAL.MODELS/osceola/eo_maxbestBROAD.tif")
#plot(east.mx)

caribbean.mx<-raster("/Users/Pedro/Desktop/Eptesicus_fuscus/data/ORIGINAL.MODELS/caribbean/ec_maxbest.tif")
#plot(caribbean.mx)

central.mx<-raster("/Users/Pedro/Desktop/Eptesicus_fuscus/data/ORIGINAL.MODELS/miradorensis/em_maxbestNEW.tif")

#plot if you wanna check 
plot(central.mx)
plot(caribbean.mx)
plot(east.mx)

######################################################################################################## not necessary to run because the models were produced already in the previous steps

#To build ENM models using Maxent algorithm; target are the three populations of the study 
# 2 min for each 
#central.maxent<- enmtools.maxent(central, env, test.prop = 0.2)
#east.maxent<- enmtools.maxent(east, env, test.prop = 0.2)
#caribbean.maxent<- enmtools.maxent(caribbean, env, test.prop = 0.2)

#response curves
#caribbean.maxent$response.plots #to see individual variables, type "$bio1" after the end 
#east.maxent$response.plots
#central.maxent$response.plots

#visualization
#visualize.enm(caribbean.maxent, env, layers = c("bio1", "bio2"), plot.test.data = TRUE)

#bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19"

#Metrics: breadth, correlation and overlap
#raster.breadth(caribbean.maxent)
#raster.breadth(east.maxent)
#raster.breadth(central.maxent)

#visuals
#plot(caribbean.maxent$suitability,main="ENM CARIBBEAN EPTESICUS")
#plot(central.maxent$suitability,main="ENM CENTRAL EPTESICUS")
#plot(east.maxent$suitability,main="ENM EAST EPTESICUS")

#plot(caribbean.mx$suitability,main="ENM CARIBBEAN EPTESICUS")



#overlap of the produced ENM maxent models

#CENTRAL X CARIBBEAN
#raster.overlap(central.maxent,caribbean.maxent)
#env.overlap(central.maxent, caribbean.maxent, env, tolerance = .001) #2 minutes to run 

#CENTRAL X EAST 
#raster.overlap(central.maxent,east.maxent)
#env.overlap(central.maxent, east.maxent, env, tolerance = .001) #2 minutes to run

#EAST X CARIBBEAN 
#raster.overlap(east.maxent,caribbean.maxent)
#env.overlap(east.maxent, caribbean.maxent, env, tolerance = .001) #2 minutes to run


# Niche identity / Equivalency (this part requires some time of running) (In this script,is 5 replications, but for the analysis it was used *100 replications) 

# (caribbean x central' caribbean x east'  east x central)

comp1 <- identity.test(species.1 = caribbean, species.2 = central, env = env, type = "mx", nreps = 5) # 5 to 10 min to run 

comp2 <- identity.test(species.1 = caribbean, species.2 = east, env = env, type = "mx", nreps = 5) # 5 to 10 min to run 

comp1

#comp1 RESULTS
#Identity test caribbean vs. central

#Identity test p-values:
# D          I   rank.cor      env.D      env.I    env.cor 
# 0.00990099 0.00990099 0.00990099 0.01980198 0.00990099 0.21782178 


# Replicates:

#|          |         D|         I|   rank.cor|     env.D|     env.I|   env.cor|
#|:---------|---------:|---------:|----------:|---------:|---------:|---------:|
#|empirical | 0.2160605| 0.4873301| -0.0762752| 0.5597426| 0.7685762| 0.1802140|
#|rep 1     | 0.6749714| 0.9120120|  0.6253389| 0.7789051| 0.9246243| 0.4166723|
# |rep 2     | 0.5911022| 0.8522744|  0.4489585| 0.6592615| 0.8818413| 0.2632105|
#|rep 3     | 0.5353831| 0.8361496|  0.3871761| 0.7513423| 0.9117484| 0.3206016|
#|rep 4     | 0.5826590| 0.8383930|  0.3189648| 0.7232946| 0.8946683| 0.3107636|
# |rep 5     | 0.7401286| 0.9260603|  0.6894539| 0.7127105| 0.9037471| 0.4160483|



comp2
#Identity test caribbean vs. east

#Identity test p-values:
# D          I   rank.cor      env.D      env.I    env.cor 
#0.00990099 0.00990099 0.00990099 0.01980198 0.01980198 0.12871287 


#Replicates:


# |          |         D|         I|   rank.cor|     env.D|     env.I|    env.cor|
#|:---------|---------:|---------:|----------:|---------:|---------:|----------:|
#|empirical | 0.0535375| 0.1848503| -0.4596765| 0.3475679| 0.6231987| -0.0560213|
#|rep 1     | 0.4547916| 0.7248457|  0.3202735| 0.6589281| 0.8850187|  0.0817281|
#|rep 2     | 0.3813446| 0.6331652| -0.0512680| 0.5522964| 0.8229799| -0.1929562|
#|rep 3     | 0.4637200| 0.7502731|  0.3239378| 0.7644635| 0.9326222|  0.2848177|
#|rep 4     | 0.4447625| 0.6933133|  0.3228111| 0.6916062| 0.9058630|  0.2369352|
#|rep 5     | 0.5163069| 0.7718138|  0.3655167| 0.7472372| 0.9344983|  0.2688439|




comp3 <- identity.test(species.1 = central, species.2 = east, env = env, type = "mx", nreps = 5) # 5 to 10 min to run 

comp3
plot(comp3)

#Identity test central vs. east

#Identity test p-values:
#D          I   rank.cor      env.D      env.I    env.cor 
#0.00990099 0.00990099 0.00990099 0.00990099 0.00990099 0.00990099 


#Replicates:


# |          |         D|         I|   rank.cor|     env.D|     env.I|    env.cor|
# |:---------|---------:|---------:|----------:|---------:|---------:|----------:|
# |empirical | 0.1870447| 0.4136682| -0.0635895| 0.0946166| 0.2512275| -0.2522867|
# |rep 1     | 0.6230172| 0.8213258|  0.4269868| 0.3287460| 0.5873707|  0.2241143|
# |rep 2     | 0.5933268| 0.8225532|  0.4648996| 0.4110708| 0.6756012|  0.3827927|
# |rep 3     | 0.6459920| 0.8316392|  0.4337313| 0.6647849| 0.8929727|  0.2925966|
# |rep 4     | 0.5955370| 0.8316453|  0.5417124| 0.3381972| 0.6139958|  0.0176340|
# |rep 5     | 0.6965089| 0.9046119|  0.6171723| 0.6636025| 0.8969380|  0.2805391|




#background / similarity test
#1
background.mx.sym.1 <- background.test(species.1 = central, species.2 = caribbean, env = env, type = "mx", nreps = 5, test.type = "symmetric" )
background.mx.asym.1 <- background.test(species.1 = central, species.2 = caribbean, env = env, type = "mx", nreps = 5, test.type = "asymmetric" )

background.mx.sym.1

#D          I   rank.cor      env.D      env.I    env.cor 
#0.01980198 0.01980198 0.19801980 0.06930693 0.09900990 0.11881188 


#Replicates:



# |          |         D|         I|   rank.cor|     env.D|     env.I|    env.cor|
# |:---------|---------:|---------:|----------:|---------:|---------:|----------:|
# |empirical | 0.1739977| 0.4275451| -0.2565609| 0.5227043| 0.7430521| -0.0134625|
#|rep 1     | 0.4345747| 0.7023205| -0.2186021| 0.4105546| 0.7054067|  0.0325962|
#|rep 2     | 0.3526325| 0.6608741| -0.6326760| 0.3582365| 0.6374075| -0.5222464|
#|rep 3     | 0.2307149| 0.4954960| -0.7033062| 0.1619794| 0.4103553| -0.2929877|
# |rep 4     | 0.2980523| 0.5992369| -0.0526177| 0.3995110| 0.6805335| -0.1905117|
# |rep 5     | 0.3299054| 0.6289511| -0.6281788| 0.3366177| 0.6271979| -0.2547544|


background.mx.asym.1    

# background test p-values:

# D         I   rank.cor      env.D      env.I    env.cor 
#0.02970297 0.02970297 0.19801980 0.06930693 0.08910891 0.14851485 


#|          |         D|         I|   rank.cor|     env.D|     env.I|    env.cor|
#|:---------|---------:|---------:|----------:|---------:|---------:|----------:|
#|empirical | 0.1739977| 0.4275451| -0.2565609| 0.5220755| 0.7424080| -0.0158906|
#|rep 1     | 0.2144904| 0.4639217| -0.2605428| 0.3429457| 0.6119267| -0.2333038|
#|rep 2     | 0.2430835| 0.5177826| -0.3007870| 0.5486727| 0.7789067| -0.1003547|
#|rep 3     | 0.2269828| 0.5105974| -0.2170224| 0.3144568| 0.6002117| -0.1833148|
#|rep 4     | 0.2305245| 0.5151442| -0.1427881| 0.3103530| 0.5978643| -0.2020002|
#|rep 5     | 0.2852135| 0.5562241| -0.1194714| 0.4004422| 0.6390998| -0.2157515|

#2
background.mx.sym.2 <- background.test(species.1 = central, species.2 = east, env = env, type = "mx", nreps = 5, test.type = "symmetric" )
background.mx.asym.2 <- background.test(species.1 = central, species.2 = east, env = env, type = "mx", nreps = 5, test.type = "asymmetric" )


background.mx.sym.2
#D          I   rank.cor      env.D      env.I    env.cor 
#0.01980198 0.01980198 0.00990099 0.19801980 0.17821782 0.34653465 


#Replicates:



#|          |         D|         I|   rank.cor|     env.D|     env.I|    env.cor|
#|:---------|---------:|---------:|----------:|---------:|---------:|----------:|
#|empirical | 0.1379039| 0.3183597| -0.1556778| 0.1184480| 0.2447277| -0.0877904|
#|rep 1     | 0.3043413| 0.5161515| -0.4885719| 0.1506112| 0.3180192| -0.0892467|
#|rep 2     | 0.2465369| 0.4972257| -0.6623449| 0.1571045| 0.3540361|  0.1164738|
# |rep 3     | 0.2205061| 0.4491638| -0.6845214| 0.2059674| 0.3374234| -0.0028609|
#|rep 4     | 0.2666387| 0.4923478| -0.6525493| 0.1742112| 0.2977730| -0.0462902|
#|rep 5     | 0.3182757| 0.5321066| -0.4468737| 0.1855077| 0.3471382| -0.0342600|

background.mx.asym.2

#Asymmetric background test
#central vs. east background


#background test p-values:

#  D          I   rank.cor      env.D      env.I    env.cor 
#0.07920792 0.07920792 0.00990099 0.18811881 0.12871287 0.28712871 


#Replicates:



# |          |         D|         I|   rank.cor|     env.D|     env.I|    env.cor|
# |:---------|---------:|---------:|----------:|---------:|---------:|----------:|
# |empirical | 0.1379039| 0.3183597| -0.1556778| 0.1168912| 0.2431486| -0.0848198|
# |rep 1     | 0.2190552| 0.4353208| -0.6003600| 0.1714065| 0.3543908|  0.0142442|
# |rep 2     | 0.2190505| 0.4296139| -0.5474365| 0.1556447| 0.3431080| -0.1397644|
# |rep 3     | 0.2137359| 0.4299325| -0.3899859| 0.1499636| 0.3161488| -0.0453452|
# |rep 4     | 0.1746739| 0.3808549| -0.6661985| 0.2086208| 0.4200127| -0.0251583|
# |rep 5     | 0.2200806| 0.4281426| -0.4890823| 0.2214843| 0.4325599|  0.2616175|




#3
background.mx.sym.3 <- background.test(species.1 = east, species.2 = caribbean, env = env, type = "mx", nreps = 5, test.type = "symmetric" )
background.mx.asym.3 <- background.test(species.1 = east, species.2 = caribbean, env = env, type = "mx", nreps = 5, test.type = "asymmetric" )

background.mx.sym.3
#Symmetric background test
#east background vs. caribbean background


#background test p-values:

#  D          I   rank.cor      env.D      env.I    env.cor 
#0.11881188 0.06930693 0.09900990 0.15841584 0.13861386 0.04950495 


#Replicates:

x

#|          |         D|         I|   rank.cor|     env.D|     env.I|    env.cor|
#|:---------|---------:|---------:|----------:|---------:|---------:|----------:|
#|empirical | 0.0347382| 0.1259763| -0.5328273| 0.1354889| 0.3197832|  0.0807271|
# |rep 1     | 0.0240692| 0.1418159| -0.8977824| 0.2886840| 0.5550705| -0.5176677|
#|rep 2     | 0.0424812| 0.1613581| -0.6354659| 0.2272195| 0.4636820| -0.0355909|
#|rep 3     | 0.0272811| 0.1161109| -0.8367145| 0.3393143| 0.5926517| -0.2529813|
#|rep 4     | 0.0689137| 0.2030075| -0.5703760| 0.1809296| 0.3704872|  0.0103422|
#|rep 5     | 0.0443866| 0.1542828| -0.8443531| 0.1453527| 0.3474098| -0.2101570|


background.mx.asym.3  

#Asymmetric background test
#east vs. caribbean background


#background test p-values:

# D          I   rank.cor      env.D      env.I    env.cor 
#0.36633663 0.48514851 0.24752475 0.11881188 0.14851485 0.08910891 



#Replicates:



# |          |         D|         I|   rank.cor|     env.D|     env.I|    env.cor|
# |:---------|---------:|---------:|----------:|---------:|---------:|----------:|
# |empirical | 0.1053021| 0.2700551| -0.3202557| 0.2187810| 0.4465181|  0.0413645|
# |rep 1     | 0.0568942| 0.1991765| -0.4899973| 0.2109871| 0.4449248|  0.0210511|
# |rep 2     | 0.1139559| 0.2904439| -0.4639158| 0.1597514| 0.3589417| -0.2221129|
# |rep 3     | 0.0817081| 0.2731677| -0.2226270| 0.2304331| 0.4786799|  0.1271186|
# |rep 4     | 0.0550264| 0.1902155| -0.3809163| 0.1279015| 0.3135561| -0.2900639|
# |rep 5     | 0.0905294| 0.2571146| -0.4849795| 0.1792491| 0.4137912| -0.1642217|