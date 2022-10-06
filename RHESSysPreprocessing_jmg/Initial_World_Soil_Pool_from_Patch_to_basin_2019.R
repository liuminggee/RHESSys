#install.packages('tidyverse')
#install.packages('tokenizers')
library(dplyr)
library('tidyverse')
# read the txt file

# after spin up the patches you use the second part to initalize the whole basin
# then run 500 years of basin scale to get the final state
# then you can run the first part code to initialize veg to zero prepare the lai spin up file
# but first you need to run the awk file to add extra row spinup_defualt_id 1 to veg_parm_Id for the world file
# the change all the veg value to zero


dir <-'C:/Work_in_WSU2018/R_Project/RHESSysPreprocessing/'

#dir<-'C:/Users/PETBUser/Documents'
# c:\Users\PETBUser\Dropbox\ecohydro\for_cluster\final_state_backup_WS\
setwd(dir)
dir()
world0<-read.table("./Data/trail_creek_twolayer.world",skip=0,fill=TRUE)
world0<-read.table("./trail_creek_subbasin_k.world",skip=0,fill=TRUE)
head(world0)

names(world0)<-c("value","name")

world0<-world0 %>% mutate (index=1:nrow(world0))
#write.table(world0,"cf_maps.txt",sep="\t\t",row.names=FALSE,quote = FALSE,col.names=FALSE)

world0_bk<-world0
# make all the vegtation initial value is zero


str(world0)
dim(world0)



# this is for initial the soil pools based on different vegetation and soil types



####soil9<-world0 %>% filter(name=="soil_parm_ID") %>% filter(value==12)

veg1<-world0 %>% filter(name=="veg_parm_ID") %>% filter(value==1) # pine
veg2<-world0 %>% filter(name=="veg_parm_ID") %>% filter(value==2) # decidous
veg3<-world0 %>% filter(name=="veg_parm_ID") %>% filter(value==3)  # grass
veg4<-world0 %>% filter(name=="veg_parm_ID") %>% filter(value==4) # noveg
veg5<-world0 %>% filter(name=="veg_parm_ID") %>% filter(value==5)  # shrub
veg49<-world0 %>% filter(name=="veg_parm_ID") %>% filter(value==49) # underdeci
veg50<-world0 %>% filter(name=="veg_parm_ID") %>% filter(value==50) # underpine
veg51<-world0 %>% filter(name=="veg_parm_ID") %>% filter(value==51) # underurban



# get the index of these vegitation from the world file
# 9 means the loam soil and 12 means the sandy-loam soil
# here we consider

index91<-veg1$index
index92<-veg2$index
index93<-veg3$index
index94<-veg4$index
index95<-veg5$index
# may not need below
index949<-veg49$index
index950<-veg50$index
index951<-veg51$index




dir()

# start reading the values needed to be assigned
patch91<-read.table('./Data/patch91.world2.state') # 1- evergreen
patch92<-read.table('./Data/patch91.world2.state')  # 2- decidous
patch93<-read.table('./Data/patch93.world2.state') # grass
patch94<-read.table('./Data/patch91.world2.state')  # no-veg
patch95<-read.table('./Data/patch125.world2.state.Y2979M12D30H24.state')  # shrub
#patch91<-patch91 %>% mutate(index=1:nrow(patch91))
#names(patch91)<-c("value","name","index")

patch91=patch91[c(55:65),] # here is to get the soil data form litter_cs_litr1c to soil_cs_soil4c
patch92=patch92[c(55:65),]
patch93=patch93[c(55:65),]
patch94=patch94[c(55:65),]
patch95=patch95[c(55:65),]


patch91[c(7,8),1]=0  # this is to set the soil_ns.sminn and soil_ns.nitrate to zero
patch92[c(7,8),1]=0 # deci
patch93[c(7,8),1]=0  # grass
patch94[,1]=0   # urban so all the states of soil are zero
patch95[c(7,8),1]=0 # this is the shrub



## now initialize the world file for soil

###patch91[,1]=0


#index91<- na.omit(index91)

for( i in 1:length(index91)){

  temp<-c((index91[i]-14): (index91[i]-4))
  world0[temp,]<-patch91[,]


}

for( i in 1:length(index92)){

  temp<-c((index92[i]-14): (index92[i]-4))
  world0[temp,]<-patch92[,]


}

for( i in 1:length(index93)){

  temp<-c((index93[i]-14): (index93[i]-4))
  world0[temp,]<-patch93[,]


}

for( i in 1:length(index94)){

  temp<-c((index94[i]-14): (index94[i]-4))
  world0[temp,]<-patch94[,]


}

for( i in 1:length(index95)){

  temp<-c((index95[i]-14): (index95[i]-4))
  world0[temp,]<-patch95[,]


}


#### now change the vegetation canopy_strata_ID same to


patchID<-world0 %>% filter(name=="patch_ID")

index0 <- patchID$index
index1 <- index0+34
index2 <- index0+89

test1 <- world0[index1,]
test2 <- world0[index2,]

# the max patch id is 16256 so start the understory id from 30000 to make less confusion

# assign the upperstory ID

world0[index1,]$value <- patchID$value
# assign the understory
world0[index2,]$value <- patchID$value + max(patchID$value) +10000 ## this number should be the largest patch ID in your


################################
# change the understory cover fraction
# all the upper story cover fraction set to 0.6 instead of using 1
# all the understory cover fraction change based on different vegetation type


# then set the cover fraction for overstory

overid<-world0 %>% filter(name=="num_stratum") %>% filter(value==2) # all the overstory is 0.6
index.overid <-overid$index+3
world0[index.overid, 1]=0.6

# set up the understory cover fraction
# the pine 50 is 0.55
pine.id <- world0 %>% filter(name=="veg_parm_ID") %>% filter(value ==50) #underpine is 50
index.pin <-pine.id$index+1 # skip the spinup_default _id
world0[index.pin, 1] =0.55

# check the distribution of value

#pine.id <- pine.id %>% mutate(index2 = index+1)
#result1 <- world0[index.pin, ]
#summary(result1)

# set up the other not 50 understory is 0.1

grass.id <- world0 %>% filter(name=="veg_parm_ID") %>% filter( value==51) # grass and shrub is 0.1

summary(grass.id)
index.grass <- grass.id$index+1

head(world0[index.grass, ])
world0[index.grass,1] =0.1

# check the distribution of value
#grass.id <- grass.id %>% mutate(index2 = index+1)

#id <-grass.id$index2

#result <- world0[id, ]
#summary(result)

# set up the dedicous understory is 0.6

decid.id <- world0 %>% filter(name=="veg_parm_ID") %>% filter( value==49) # decidous is 0.6
index.decid <- decid.id$index+1

head(world0[index.decid, ])
world0[index.decid,1] =0.6

# check the distribution of value
#decid.id <- decid.id %>% mutate(index2 = index+1)

#id <-decid.id$index2

#result <- world0[id, ]
#summary(result)








### output the world files after change the soil pool and understory



# controlt the result not showing e
options("scipen"=999, "digits"=8) # seems like important



world0_backup2<-world0
#world0 <- world0_backup2




world0<- world0 %>% select(value,name)
write.table(world0,"TC_world_R_20190212.txt",sep="\t\t",row.names=FALSE,quote = FALSE,col.names=FALSE)




















###below is select both the soil type and vegetation type

t=1
rm(index91)
index91<-vector()
for (i in 1:nrow(soil9)){
  for (j in 1:nrow(veg1)){
    tem<-(soil9$index[i]-veg1$index[j])
    if (tem==(-31)){
      index91[t]<-soil9$index[i]
      t=t+1
    }
  }


}


