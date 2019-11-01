rm(list=ls(all=TRUE))
set.seed(10)

#read important functions
setwd('U:\\GIT_models\\activcenter_subset_locations')
library('Rcpp')
sourceCpp('aux1.cpp')
source('gibbs sampler.R')
source('gibbs functions.R')

#get data
dat=read.csv('fake data.csv',as.is=T)
grid.coord=read.csv('fake data grid.csv',as.is=T)

#basic setup
ngibbs=1000
nburn=ngibbs/2
n.ac=20
gamma1=0.1

#initial coordinates (define this based on data instead of uninformative as below)
ind=sample(nrow(grid.coord),size=n.ac)
ac.coord.init=grid.coord[ind,]
possib.ac=grid.coord #these don't have to be identical (i.e., we can define AC's on a coarser grid)

#run gibbs sampler
options(warn=2)
dat=data.matrix(dat)
res=gibbs.activity.center(dat=dat,grid.coord=grid.coord,n.ac=n.ac,
                          ac.coord.init=ac.coord.init,gamma1=gamma1,
                          possib.ac=possib.ac)

