set.seed(10)

#load libraries and read important functions
library('Rcpp')
library(raster)
library(sf)
library(ggplot2)
library(dplyr)
library(viridis)
library(progress)
library(rnaturalearth)
library(rnaturalearthdata)


sourceCpp('aux1.cpp')
source('gibbs sampler.R')
source('gibbs functions.R')
source('helper functions.R')



#################
### Load data ###
#################


#load data
dat<- read.csv("Snail Kite Gridded Data.csv", header = T, sep = ",")
obs<- get.summary.stats_obs(dat)  #frequency of visitation in each location (column) for each time segment (row)

#geographical coordinates of locations
utm.crs<- CRS('+init=epsg:32617')
extent<- extent(min(dat$utmlong), max(dat$utmlong), min(dat$utmlat), max(dat$utmlat))
res<- 5000
buffer<- 10000
grid.coord<- grid.summary.table(dat=dat, crs=utm.crs, extent=extent, res=res, buffer=buffer)
dat<- left_join(dat, grid.coord, by="grid.cell") #add gridded locs to DF

#Define initial activity centers (top 50 by # of obs)
tmp<- colSums(obs[,-1]) %>% data.frame(grid.cell = colnames(obs[,-1]), nobs = .) %>%
      arrange(desc(nobs)) %>% slice(n=1:50) %>% select(grid.cell)
ind<- sample(as.numeric(tmp[,1]), size = 20, replace = F)
# nobs<- colSums(obs[,ind+1])
ac.coord.init<- grid.coord[ind,]

#top 50
ac.coord.init2<- grid.coord[as.numeric(tmp[,1]),]

#potential locations for activity centers (AC)
possib.ac=grid.coord #these don't have to be identical (i.e., we can define AC's on a coarser grid)


#########################
### Run Gibbs sampler ###
#########################

#basic setup
ngibbs=1000
nburn=ngibbs/2
n.ac=20
gamma1=0.1

#run gibbs sampler
options(warn=2)

pb <- progress_bar$new(
  format = " iteration (:current/:total) [:bar] :percent [Elapsed: :elapsed, Remaining: :eta]",
  total = ngibbs, clear = FALSE, width= 100)

res=gibbs.activity.center(dat=obs[-1,],grid.coord=grid.coord[,-3],n.ac=n.ac,
                          ac.coord.init=ac.coord.init[,-3],gamma1=gamma1,
                          possib.ac=possib.ac[,-3])


#plot output and look at frequency of AC visitation
plot(res$logl,type='l')
plot(res$phi,type='l')

#######################
### Load Saved Data ###
#######################

# ac<- data.frame(ac=res$z[ngibbs,])
# ac.coords<- matrix(NA, length(unique(ac$ac)), 2)
# colnames(ac.coords)<- c("x","y")
# tmp<- res$coord[ngibbs,]
# 
# for (i in 1:length(unique(ac$ac))) {
#   ac.coords[i,]<- round(c(tmp[i], tmp[i+n.ac]), 0)
# } 
# 
# ac.coords<- data.frame(ac.coords, ac=1:length(unique(ac$ac)))

ac.coords<- read.csv("Activity Center Coordinates.csv", header = T, sep = ',')
ac<- read.csv("ac.csv", header = T, sep = ',')
table(ac$ac)
obs<- cbind(ac, obs)


############################
### Add ACs to Dataframe ###
############################


obs1<- obs %>% filter(id == 1) %>% mutate(time.seg = 1:nrow(.)) %>% dplyr::select(ac, time.seg)
obs12<- obs %>% filter(id == 12) %>% mutate(time.seg = 1:nrow(.)) %>% dplyr::select(ac, time.seg)
obs19<- obs %>% filter(id == 19) %>% mutate(time.seg = 1:nrow(.)) %>% dplyr::select(ac, time.seg)
obs27<- obs %>% filter(id == 27) %>% mutate(time.seg = 1:nrow(.)) %>% dplyr::select(ac, time.seg)

dat1<- dat %>% filter(id == 1) %>% left_join(obs1, by = "time.seg")
dat12<- dat %>% filter(id == 12) %>% left_join(obs12, by = "time.seg")
dat19<- dat %>% filter(id == 19) %>% left_join(obs19, by = "time.seg")
dat27<- dat %>% filter(id == 27) %>% left_join(obs27, by = "time.seg")

dat<- rbind(dat1, dat12, dat19, dat27)

#Calculate number of obs per AC
dat %>% filter(id==1) %>% dplyr::select(ac) %>% table()
dat %>% filter(id==12) %>% dplyr::select(ac) %>% table()
dat %>% filter(id==19) %>% dplyr::select(ac) %>% table()
dat %>% filter(id==27) %>% dplyr::select(ac) %>% table()



## Map

#Load world map data
usa <- ne_states(country = "United States of America", returnclass = "sf")
fl<- usa %>% filter(name == "Florida")
fl<- sf::st_transform(fl, crs = "+init=epsg:32617") #change projection to UTM 17N

nests<- dat %>% group_by(id) %>% select(c(id, utmlong, utmlat)) %>% slice(n=1)


# ACs and initial values
ggplot() +
  geom_sf(data = fl) +
  coord_sf(xlim = c(min(dat$utmlong-20000), max(dat$utmlong+20000)),
           ylim = c(min(dat$utmlat-20000), max(dat$utmlat+20000)), expand = FALSE) +
  geom_point(data = ac.coord.init2, aes(x, y, color = "A (highest: n=50)"), size = 5, shape = 1) +
  geom_point(data = ac.coord.init, aes(x, y, color = "B (initial: n=20)"), size = 2) +
  geom_point(data = ac.coords, aes(x, y, color = "C (model: n=20)"), size = 3, alpha = 0.5) +
  geom_point(data = nests, aes(utmlong, utmlat, color = "Nests"), shape = 17, size = 2) +
  labs(x="Longitude", y="Latitude") +
  scale_color_manual("", values = c("grey40",viridis(n=5)[c(3,5)],"red")) +
  guides(color = guide_legend(override.aes = list(shape = c(1, 16, 16, 17)))) +
  theme_bw()

# ACs and snail kite locs
ggplot() +
  geom_sf(data = fl) +
  coord_sf(xlim = c(min(dat$utmlong-20000), max(dat$utmlong+20000)),
           ylim = c(min(dat$utmlat-20000), max(dat$utmlat+20000)), expand = FALSE) +
  geom_point(data = dat, aes(utmlong, utmlat, color=ac), size=1, alpha = 0.6) +
  geom_point(data = ac.coords, aes(x, y, color = ac), size = 4, pch = 1, stroke = 1) +
  scale_color_viridis_c("Activity Center") +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  facet_wrap(~id)


###################
### Save Output ###
###################

setwd("~/Documents/Snail Kite Project/Data/R Scripts/cluster_tsegments_loc")

# write.csv(ac.coords, "Activity Center Coordinates.csv", row.names = F)
# write.csv(dat, "Snail Kite Gridded Data_AC.csv", row.names = F)
