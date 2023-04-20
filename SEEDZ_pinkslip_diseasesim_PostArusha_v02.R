# program to simulate stochastic SIR model on simulated movement network 
# Paul Johnson 
# 24 January 2019

rm(list = ls())
graphics.off()

# simulate fast and slow disease (with no spatial variation in risk)
# target vaccination at high mp.geomean wards (top 5%)
# target movement ban at high betweeness wards (top 5%)

root.path <- 
  if(Sys.info()["nodename"] == "boydorr-ws4-ibls-gla-ac-uk") 
    "~/projects/SEEDZ/" else
      "~/OneDrive - University of Glasgow/Projects/SEEDZ/"
setwd(paste0(root.path, "movement/analysis/modelling/disease.sim/"))

# packages
library(igraph)
library(SimInf)
# vignette("SimInf")
library(parallel)
library(scales)
library(RColorBrewer)

# functions

# load sim.hurdle.move function
source('market.ward.sim.v3.R')


# adapted from http://www.r-bloggers.com/great-circle-distance-calculations-in-r/
# Calculates the geodesic distance between two points specified by radian latitude/longitude using the
# Spherical Law of Cosines (slc)
gcd.slc <- function(long, lat) {
  R <- 6371 # Earth mean radius [km]
  latr <- lat*pi/180
  longr <- long*pi/180
  x <- prod(sin(latr)) + prod(cos(latr)) * cos(diff(longr))
  x <- min(x, 1)
  x <- max(x, -1)
  d <- acos(x) * R
  return(d) # Distance in km
}

# function to find nearest location 
which.nearest <-
  function(latlon, latlon.mat, show.dist = FALSE) {
    d <- 
      apply(latlon.mat, 1, 
            function(x) 
              gcd.slc(long = c(x[2], latlon[2]), lat = c(x[1], latlon[1])))
    i <- which.min(d)
    if(show.dist) print(paste(round(d[i], 2), "km"))
    attributes(i)$min.dist.km <- d[i]
    i
  }

# function to split nodes into subnodes to create within-node heterogeneity
split.nodes.sir<- 
  function(n.subnodes, u0, adj = NULL) {
    compartments <- c("S", "I", "R")
    sub.comp <- 
      paste0(rep(compartments, each = n.subnodes), rep(1:n.subnodes, length(compartments)))
    n.sub.comp <- length(sub.comp)
    if(is.null(adj)) {
      d <- sqrt(n.subnodes)
      adj <- 
        structure(rep(0, n.subnodes^2), dim = c(n.subnodes, n.subnodes), 
                  dimnames = list(1:n.subnodes, 1:n.subnodes))
      for(i in 1:nrow(adj)) {
        for(j in 1:ncol(adj)) {
          if(((i - j) == 1 && j%%d != 0) | ((i - j) == -1 && i%%d != 0) | abs(i - j) == d) adj[i, j] <- 1
        }; rm(j)
      }; rm(i)
    }
    stopifnot(all(dim(adj) == n.subnodes))
    transitions.tab <-
      sapply(1:n.subnodes, function(sn) {
        nbrs <- rownames(adj)[adj[sn, ] == 1]
        env.contam <- 
          if(length(nbrs) > 0) paste0(" + (", paste0("I", nbrs, collapse = "+"), ")*coupling") else
            NULL
        N <- paste0("(S", sn, "+I", sn, "+R", sn, ")")
        si.trans <-
          paste0("S", sn, 
                 " -> ", N, " > 0 ? beta*S", sn, "*(I", sn, 
                 env.contam,
                 ")/", N,
                 " : 0 -> I", sn, " + Ic")
        ir.trans <-
          paste0("I", sn, " -> gamma*I", sn, " -> R", sn)
        c(si.trans, ir.trans)
      })
    transitions <- c(t(transitions.tab))
    
    u0.list <- 
      lapply(compartments, 
             function(cmp) {
               out <- t(sapply(1:nrow(u0), function(i) rmultinom(1, u0[i, cmp], rep(1, n.subnodes))))
               if(n.subnodes == 1) out <- t(out)
               colnames(out) <- paste0(cmp, 1:n.subnodes)
               out
             })
    u0.out <- data.frame(do.call("cbind", u0.list))
    stopifnot(all(names(u0.out) %in% sub.comp) & all(sub.comp %in% names(u0.out)))
    u0.out <- u0.out[, sub.comp]
    u0.out$Ic <- u0$Ic
    
    E <- matrix(rep(0, n.sub.comp * 2),
                nrow = n.sub.comp,
                dimnames = list(sub.comp))
    E[substr(rownames(E), 1, 1) == "S", 1] <- 1
    E[, 2] <- 1
    E <- 
      rbind(
        cbind(E, rbind(diag(n.subnodes), diag(n.subnodes), diag(n.subnodes))), 
        Ic = 0)
    
    stopifnot(all(rownames(E) %in% names(u0.out)) & all(names(u0.out) %in% rownames(E)))
    
    list(compartments = rownames(E), transitions = transitions, E = E, u0 = u0.out, adj = adj)
  }



# the date the epidemic starts
start.date <- as.Date("2015-01-01")

# runs through which months?
n.months <- 12
n.days <- n.months/12 * 365
all.months <- formatC(1:n.months, width = 2, flag = "0")

# which day of the year (1:365) is the first of each month?
all.months.day1 <-
  sapply(all.months, function(mm) {
    print(mm)
    yyyy <- 2014 + ceiling(as.numeric(mm)/12)
    mm <- ((as.numeric(mm) - 1) %% 12) + 1
    as.Date(paste(yyyy, mm, "01", sep = "-")) - start.date + 1
  })

# select species to model
sp <- c("animals", "cattle", "caprine")[2]


# import wards data - required to identify the nearest primary market
wards <- read.csv(paste0(root.path, "mapping/WardSpatialData.csv"), stringsAsFactors = FALSE, quote = "")

# drop tanga
wards <- droplevels(wards[wards$region %in% c("Arusha", "Manyara", "Kilimanjaro"), ])
all.ward.names <- wards$fullname
dim(wards)
dim(unique(wards))
head(wards)
head(wards)
table(wards$region)
length(table(wards$district))

# wards with <1000 animals are boosted to 1000
wards[wards[, paste0("pop.", sp)] < 1000, paste0("pop.", sp)] <- 1000
sort(wards[, paste0("pop.", sp)])
dim(wards)
ward.names <- wards$fullname
rownames(wards) <- ward.names

# which wards contain secondary markets?
is.mkt <- ward.names %in% c("arusha/monduli/meserani", "arusha/arusha/bwawani", "kilimanjaro/hai/machamekusini")
mkts <- ward.names[is.mkt]
not.mkts <- ward.names[!is.mkt]
not.mkts.tab <- matrix(sample(not.mkts), ncol = 5)

# import list of contiguous wards, to allow spread via local diffusion as well as 
# via the movement network
contig.wards.all <- 
  read.csv(paste0(root.path, "mapping/contiguous.wards.csv"), 
           stringsAsFactors = FALSE, quote = "")
rownames(contig.wards.all) <- paste(contig.wards.all$from, contig.wards.all$to, sep = ".")


# choose ward(s) to seed epidemic in (not used in this version of the program)
seed.wards.list <- # chosen to have high outdegree and RVF risk
  sample(c(
    "arusha/karatu/karatu",       
    "manyara/babati/bashneti",
    "manyara/mbulu/imboru"))

# exlude seed ward from intervention
exclude.seed <- FALSE

# how far apart are the seed wards?
gcd.slc(wards[seed.wards.list[c(1, 2)], "long"], 
        wards[seed.wards.list[c(1, 2)], "lat"])
gcd.slc(wards[seed.wards.list[c(1, 3)], "long"], 
        wards[seed.wards.list[c(1, 3)], "lat"])
gcd.slc(wards[seed.wards.list[c(2, 3)], "long"], 
        wards[seed.wards.list[c(2, 3)], "lat"])


# import network measures for wards
netmeasures <-
  read.csv(paste0("networks/simmoves.monthly.2015.RVF.gemma/", sp,"/AllWards.csv"), 
           stringsAsFactors = FALSE, quote = "", row.names = 1)
head(netmeasures)

# create ranks for selected measures

# sets of closely correlated measures:
#   year.out.deg, mp.outdegree, mp.degree, cy.outdegree, cy.degree
#   year.in.degree, mp.indegree, cy.indegree
#   betweenness, cy.betweenness


meas.to.rank <-
  c("mp.eigval", "mp.indegree",	"mp.outdegree", "mp.betweenness", "betweenness", "year.in.deg", "year.out.deg")
for(x in meas.to.rank) netmeasures[, paste0(x, ".rank")] <- rank(-netmeasures[, x])
rm(x)

# make geometric mean of mp.indegree and mp.outdegree
netmeasures$mp.geomean <- 
  exp(apply(log(netmeasures[, c("mp.indegree", "mp.outdegree")]), 1, mean))

# create sum of ranks for selected measures 
netmeasures$mp.comb.rank <- rank(rowSums(netmeasures[, paste0(meas.to.rank[1:4], ".rank")]))

#tw <- c("arusha/monduli/meserani", "kilimanjaro/hai/machamekusini", "arusha/arusha/bwawani", "kilimanjaro/mwanga/mgagao")
#netmeasures[tw, ]
#plot(netmeasures[, sort(names(netmeasures)[grep("rank", names(netmeasures))])])


# restrict to included wards
netmeasures <- netmeasures[ward.names, ]
head(netmeasures)
dim(netmeasures)
setdiff(rownames(netmeasures), wards$fullname)
setdiff(wards$fullname, rownames(netmeasures))
#plot(netmeasures[, c("betweenness.rank", "year.in.deg.rank", "mp.comb.rank")])

# add net measures to the wards data frame
wards <- cbind(wards, netmeasures[rownames(wards), ])

# import movement parameters

# load predicted probability of movement
mu.z.list <-
  lapply(all.months, function(month) {
    m <- formatC(1 + (as.numeric(month) - 1) %% 12, width = 2, flag = "0")
    mu.z <- 
      as.matrix(
        read.csv(paste0("networks/simmoves.monthly.2015/", sp,"/", sp, ".month", m, ".mvt.matrix.prob.2018-06-07-00.csv"), 
                 row.names = 1))
    # make rownames and colnames consistent                  
    colnames(mu.z) <- rownames(mu.z)
    # keep movements to and from included wards
    mu.z[ward.names, ward.names]
  })


# load predicted number moved given movement
mu.c.list <-
  lapply(all.months, function(month) {
    m <- formatC(1 + (as.numeric(month) - 1) %% 12, width = 2, flag = "0")
    mu.c <- 
      as.matrix(
        read.csv(paste0("networks/simmoves.monthly.2015/", sp,"/", sp, ".month", m, ".mvt.matrix.rate.2018-06-07-00.csv"), 
                 row.names = 1))
    # make rownames and colnames consistent                  
    colnames(mu.c) <- rownames(mu.c)
    # keep movements to and from included wards
    mu.c[ward.names, ward.names]
  })

# load dispersion (theta) parameter for the ZTNB distribution
theta.list <- 
  lapply(all.months, function(month) {
    m <- formatC(1 + (as.numeric(month) - 1) %% 12, width = 2, flag = "0")
    read.csv(paste0("networks/simmoves.monthly.2015/", sp,"/", sp, ".month", m, ".mvt.matrix.rate.theta.2018-06-07-00.csv"))$theta
  })

names(mu.z.list) <- names(mu.c.list) <- names(theta.list) <- all.months

# calculate expected livestock flow for the whole time perion (all.months)
mu.yr <- Reduce("+", lapply(all.months, function(m) mu.z.list[[m]] * mu.c.list[[m]]))

# define a market as a ward with expected outflow >0 for the year
mkts.all <- unique(c(mkts, rownames(mu.yr)[rowSums(mu.yr) != 0]))

# estimate betweenness
between <- wards$betweenness
names(between) <- rownames(wards)

# estimate outdegree
outdeg <- wards$year.out.deg
names(outdeg) <- rownames(wards)


# estimate indegree
indeg <- wards$year.in.deg
names(indeg) <- rownames(wards)


# select wards for targetting interventions

# those with highest rift risk
rvf.des <- ward.names[order(-wards[ward.names, "rvf.risk.dummy"])]
mkts %in% rvf.des[1:80]
# check they are ordered by RVF risk
wards[rvf.des[1:5], "rvf.risk.dummy"]

# those with most animals (change to density, once the model includes density?)
n.animals.des <- ward.names[order(-wards[ward.names, paste0("pop.", sp)])]
mkts %in% n.animals.des[1:80]
# check they are ordered by no of animals
wards[n.animals.des[1:5], paste0("pop.", sp)]

# highest betweenness 
between.des <- names(rev(sort(between[ward.names])))
mkts %in% between.des[1:80]
# check they are ordered by betweenness
wards[between.des[1:5], "betweenness"]

# highest betweenness combined with highest RVF risk
between.rvf.des <-
  ward.names[order(rank(-wards$rvf.risk.dummy) + wards$betweenness.rank)]
#  c(between.des[between.des %in% wards$fullname[wards$rvf.risk.dummy.bin == 1]],
#    between.des[between.des %in% wards$fullname[wards$rvf.risk.dummy.bin == 0]])
# check they are ordered by betweenness and rvf risk
wards[between.rvf.des[1:5], c("betweenness", "rvf.risk.dummy")]

# highest outdegree 
outdegree.des<- names(rev(sort(outdeg[ward.names])))
mkts %in% outdegree.des[1:80]
# check they are ordered by outdegree
wards[outdegree.des[1:5], "year.out.deg"]

# highest outdegree combined with highest RVF risk
outdegree.rvf.des <-
  ward.names[order(rank(-wards$rvf.risk.dummy) + wards$year.out.deg.rank)]
# check they are ordered by outdegree and rvf risk
wards[outdegree.rvf.des[1:5], c("year.out.deg", "rvf.risk.dummy")]

# highest indegree
indegree.des <- names(rev(sort(indeg[ward.names])))
mkts %in% indegree.des[1:80]
# check they are ordered by indegree
wards[indegree.des[1:5], "year.in.deg"]

# highest indegree combined with highest RVF risk
indegree.rvf.des <-
  ward.names[order(rank(-wards$rvf.risk.dummy) + wards$year.in.deg.rank)]
# check they are ordered by indegree and rvf risk
wards[indegree.rvf.des[1:5], c("year.in.deg", "rvf.risk.dummy")]

# highest combined multiplex rank
mp.comb.des <- rownames(wards)[order(wards$mp.comb.rank)]
mkts %in% mp.comb.des[1:80]
# check they are ordered by multiplex rank
wards[mp.comb.des[1:5], "mp.comb.rank"]

# highest combined multiplex rank combined with highest RVF risk
mp.comb.rvf.des <-
  ward.names[order(rank(-wards$rvf.risk.dummy) + wards$mp.comb.rank)]
# check they are ordered by combined multiplex rank and rvf risk
wards[mp.comb.rvf.des[1:5], c("mp.comb.rank", "rvf.risk.dummy")]

# highest MP geometric mean degree
mp.geomean.des <- rownames(wards)[order(-wards$mp.geomean)]
mkts %in% mp.geomean.des[1:80]
mkts %in% mp.geomean.des[1:20]
table(mkts.all %in% mp.geomean.des[1:20])
# check they are ordered by geometric mean degree
wards[mp.geomean.des[1:5], "mp.geomean"]

# highest MP betweenness
mp.between.des <- rownames(wards)[order(-wards$mp.betweenness)]
mkts %in% mp.between.des[1:80]
mkts %in% mp.between.des[1:20]
table(mkts.all %in% mp.between.des[1:20])

# check they are ordered by MP betweenness
wards[mp.between.des[1:5], "mp.betweenness"]


# highest MP eigenvalue centrality
mp.eigval.des <- rownames(wards)[order(-wards$mp.eigval)]
mkts %in% mp.eigval.des[1:80]
mkts %in% mp.eigval.des[1:20]
table(mkts.all %in% mp.eigval.des[1:20])
# check they are ordered by eigenvalue centrality
wards[mp.eigval.des[1:5], "mp.eigval"]



# how similar are the top wards selected by
#   mp eigenvalue centrality?
#   mp betweenness?
#   mp geomean degree?

n.compare <- 20
unique(c(mp.between.des[1:n.compare], mp.geomean.des[1:n.compare], mp.eigval.des[1:n.compare]))

# between and geomean degree differ on 9 wards in top 20
setdiff(mp.between.des[1:n.compare], between.des[1:n.compare])
setdiff(between.des[1:n.compare], mp.between.des[1:n.compare])

# between and geomean degree differ on 7 wards in top 20
setdiff(mp.between.des[1:n.compare], mp.geomean.des[1:n.compare])
setdiff(mp.geomean.des[1:n.compare], mp.between.des[1:n.compare])

# between and eigval differ on 15 wards in top 20
setdiff(mp.between.des[1:n.compare], mp.eigval.des[1:n.compare])
setdiff(mp.eigval.des[1:n.compare], mp.between.des[1:n.compare])

# geomean and eigval differ on 10 wards in top 20
setdiff(mp.geomean.des[1:n.compare], mp.eigval.des[1:n.compare])
setdiff(mp.eigval.des[1:n.compare], mp.geomean.des[1:n.compare])


###################
###################
###################
# GLOBAL SETTINGS #
###################
###################
###################


# how many subnodes to split the wards into, to simulate heterogeneity within wards 
n.subnodes <- 1 # 1 means don't split

# strength of coupling between subnodes (used if n.subnodes > 1)
coupling <- 0.1

# SIR parameters
inf.pars <- 
  list(
    beta = c(1.5/7, 3/7),      # S -> I rate
    gamma = 1/7,       # I -> R rate
    i0 = 10,           # number of infectious animals to seed into ward(s) at the start of the epidemic
    var.beta = FALSE)  # should the transmission rate vary between wards       

# R0
(R0 <- inf.pars[["beta"]] / inf.pars[["gamma"]])
inf.pars <- 
  c(inf.pars, 
    cov = 0.7) #max(0.5, min(0.1 + 1 - 1/R0, 0.9)))

# reporting threshold
rep.thresh <- 0.01

# allow beta to vary with RVF risk
#if(inf.pars[["var.beta"]] > 0) {
#  wards$beta <-
#    inf.pars["beta"] * wards[ward.names, "rvf.risk.dummy"] / mean(wards[ward.names, "rvf.risk.dummy"])
#} else wards$beta <- inf.pars["beta"]

# R0
#wards$R0 <- wards$beta / inf.pars["gamma"]
#hist(wards$R0)
#mean(wards$R0)
#wards$cov <- 1 - 1/wards$R0
#wards$cov[wards$cov < 0] <- 0
#wards$cov[wards$cov > 1] <- 1

# how many vaccines are available? Inf means unlimited
n.vax.dose <- c(Inf, round(sum(wards[, paste0("pop.", sp)]) / 10))[1]

# how many replicates of each scenario to run?
nrep <- 3 * nrow(not.mkts.tab)

# initial node states
# these inital states are common to all simulations (outside the loop)
n <- length(ward.names)
u0.outer <- 
  data.frame(
    S = rep(0, n),
    I = rep(0, n),
    R = rep(0, n))
rownames(u0.outer) <- ward.names

mkts.i <- which(ward.names %in% mkts)
ward.names[mkts.i]

# add susceptible animals
u0.outer$S[!is.mkt] <- wards[ward.names[!is.mkt], paste0("pop.", sp)]

# add cumulative count of infecteds and vax column
u0.outer$Ic <- u0.outer$I

# inspect u0.outer
head(u0.outer)
dim(u0.outer)


# movement ban starts where and when? (use day=Inf for no ban, start=1 for no movement ever)
# note that the movement ban happens before contiguous movements, so will not include them
# if contiguous movements are also being banned, set contig.move.n <- 0

# target how many wards?
n.ward <- c(few = 20, many = 80, all = n)

mban.months <- all.months[-1]
mban.list <-
  list(
    none = list(ward = "none", mban.months = "none"),
    all = list(ward = ward.names[ward.names %in% mkts.all], mban.months = mban.months),
    rand = list(ward = sample(mkts.all), mban.months = mban.months),
    between = list(ward = between.des[between.des %in% mkts.all], mban.months = mban.months),
    indegree = list(ward = indegree.des[indegree.des %in% mkts.all], mban.months = mban.months),
    outdegree =  list(ward = outdegree.des[outdegree.des %in% mkts.all], mban.months = mban.months),
    mp.comb = list(ward = mp.comb.des[mp.comb.des %in% mkts.all], mban.months = mban.months),
    mp.eigval = list(ward = mp.eigval.des[mp.eigval.des %in% mkts.all], mban.months = mban.months),
    mp.between = list(ward = mp.between.des[mp.between.des %in% mkts.all], mban.months = mban.months),
    mp.geomean = list(ward = mp.geomean.des[mp.geomean.des %in% mkts.all], mban.months = mban.months),
    rvf = list(ward = rvf.des[rvf.des %in% mkts.all], mban.months = mban.months),
    between.rvf = list(ward = between.rvf.des[between.rvf.des %in% mkts.all], mban.months = mban.months),
    indegree.rvf = list(ward = indegree.rvf.des[indegree.rvf.des %in% mkts.all], mban.months = mban.months),
    outdegree.rvf = list(ward = outdegree.rvf.des[outdegree.rvf.des %in% mkts.all], mban.months = mban.months),
    mp.comb.rvf = list(ward = mp.comb.rvf.des[mp.comb.rvf.des %in% mkts.all], mban.months = mban.months),
    n.animals = list(ward = n.animals.des[n.animals.des %in% mkts.all], mban.months = mban.months))

# where, when and what proportion to vaccinate? use day=Inf or cov=0 for never
vax.list <-
  list(
    none = 
      data.frame(
        ward = "none",
        day = Inf,
        cov = 0,
        stringsAsFactors = FALSE),
    all.lo = 
      data.frame(
        ward = ward.names, 
        day = 0, 
        cov = 0.1, 
        stringsAsFactors = FALSE,
        row.names = NULL),
    all.hi = 
      data.frame(
        ward = ward.names, 
        day = 0, 
        cov = inf.pars["cov"], 
        stringsAsFactors = FALSE,
        row.names = NULL),
    rand = 
      data.frame(
        ward = sample(ward.names), 
        day = 0, 
        cov = inf.pars["cov"], 
        stringsAsFactors = FALSE,
        row.names = NULL),
    between = 
      data.frame(
        ward = between.des, 
        day = 0, 
        cov = inf.pars["cov"], 
        stringsAsFactors = FALSE,
        row.names = NULL),
    mp.between = 
      data.frame(
        ward = mp.between.des, 
        day = 0, 
        cov = inf.pars["cov"], 
        stringsAsFactors = FALSE,
        row.names = NULL),
    indegree = 
      data.frame(
        ward = indegree.des, 
        day = 0, 
        cov = inf.pars["cov"], 
        stringsAsFactors = FALSE,
        row.names = NULL),
    outdegree = 
      data.frame(
        ward = outdegree.des, 
        day = 0, 
        cov = inf.pars["cov"], 
        stringsAsFactors = FALSE,
        row.names = NULL),
    mp.comb = 
      data.frame(
        ward = mp.comb.des, 
        day = 0, 
        cov = inf.pars["cov"], 
        stringsAsFactors = FALSE,
        row.names = NULL),
    mp.geomean = 
      data.frame(
        ward = mp.geomean.des, 
        day = 0, 
        cov = inf.pars["cov"], 
        stringsAsFactors = FALSE,
        row.names = NULL),
    mp.eigval= 
      data.frame(
        ward = mp.eigval.des, 
        day = 0, 
        cov = inf.pars["cov"], 
        stringsAsFactors = FALSE,
        row.names = NULL),
    rvf = 
      data.frame(
        ward = rvf.des, 
        day = 0, 
        cov = inf.pars["cov"], 
        stringsAsFactors = FALSE,
        row.names = NULL),
    between.rvf = 
      data.frame(
        ward = between.rvf.des, 
        day = 0, 
        cov = inf.pars["cov"], 
        stringsAsFactors = FALSE,
        row.names = NULL),
    indegree.rvf = 
      data.frame(
        ward = indegree.rvf.des, 
        day = 0, 
        cov = inf.pars["cov"], 
        stringsAsFactors = FALSE,
        row.names = NULL),
    outdegree.rvf = 
      data.frame(
        ward = outdegree.rvf.des, 
        day = 0, 
        cov = inf.pars["cov"], 
        stringsAsFactors = FALSE,
        row.names = NULL),
    mp.comb.rvf = 
      data.frame(
        ward = mp.comb.rvf.des, 
        day = 0, 
        cov = inf.pars["cov"], 
        stringsAsFactors = FALSE,
        row.names = NULL),
    n.animals = 
      data.frame(
        ward = n.animals.des, 
        day = 0, 
        cov = inf.pars["cov"], 
        stringsAsFactors = FALSE,
        row.names = NULL),
    mkts.monthly = 
      list(
        ward = ward.names[is.mkt], 
        day = all.months.day1 + 3, 
        cov = inf.pars["cov"], 
        stringsAsFactors = FALSE))


# allow vaccine coverage to vary among wards according to RVF risk
if(FALSE) {
  for(vn in names(vax.list)[!names(vax.list) %in% c("none", "all.lo")]) {
    vax.list[[vn]]$cov <- wards[vax.list[[vn]]$ward, "cov"] 
  }; rm(vn)
}


vax.sets <-
  list(
    c("none"),
    c("none", "all.hi"),
    c("none", "all.hi", "mp.between", "mp.geomean", "mp.eigval", "rand"),
    c("none", "all.hi", "mp.between", "mp.geomean", "mp.eigval", "rand", "n.animals"),
    c("none", "all.hi", "all.lo", "rand", "between", "indegree", "mp.comb"),
    c("none", "all.hi", "all.lo", "rand"),
    c("none", "all.hi", "between", "indegree", "mp.comb"),
    c("none", "all.hi", "rvf"),
    c("none", "all.hi", "between.rvf", "indegree.rvf", "mp.comb.rvf", "rvf"))

vax.plans <- unique(unlist(vax.sets))

mban.sets <- 
  list(
    c("none"),
    c("none", "all"),
    c("none", "all", "mp.between", "mp.geomean", "mp.eigval", "rand"),
    c("none", "all", "mp.between", "mp.geomean", "mp.eigval", "rand", "n.animals"),
    c("none", "all", "rand", "between", "indegree", "mp.comb"),
    c("none", "all", "rand", "rvf"),
    c("none", "all", "rand", "between.rvf", "indegree.rvf", "mp.comb.rvf", "rvf"),
    "all")

mban.plans <- unique(unlist(mban.sets))

# intervention table
intervention.vax <- 
  expand.grid(vax = vax.sets[[4]],
              mban = "none", 
              n.vax = n.ward["few"],
              n.mban = 1,
              beta = inf.pars[["beta"]],
              gamma = inf.pars[["gamma"]],
              rep = 1:nrep, 
              stringsAsFactors = FALSE)

intervention.mban <- 
  expand.grid(vax = "none", 
              mban = mban.sets[[4]], 
              n.vax = 1,
              n.mban = n.ward["few"],
              beta = inf.pars[["beta"]],
              gamma = inf.pars[["gamma"]],
              rep = 1:nrep,
              stringsAsFactors = FALSE)
intervention <- rbind(intervention.vax, intervention.mban)
if(all(intervention$vax == "none")) 
  intervention$plan <- intervention$mban else
    if(all(intervention$mban == "none")) intervention$plan <- intervention$vax else
      intervention$plan <- paste0("VX: ", intervention$vax, "; MB: ", intervention$mban)
dim(intervention)

intervention$n.vax[intervention$vax %in% "none"] <- 0
intervention$n.vax[intervention$vax %in% c("all.lo", "all.hi")] <- n
intervention$n.mban[intervention$mban %in% "none"] <- 0
intervention$n.mban[intervention$mban %in% "all"] <- length(mkts.all)
intervention <- intervention[!duplicated(intervention), ]
intervention$k <- 1:nrow(intervention)
intervention$contig.move.n <- 1


# choose ward to seed with infected animals
#intervention$seed.ward <- seed.wards.list[(intervention$rep %% length(seed.wards.list)) + 1]
#intervention$seed.ward <- "rand05"
intervention$seed.ward <- "balanced"

# check all chosen interventions have been coded
stopifnot(
  c(all(intervention$vax %in% names(vax.list)), 
    all(intervention$mban %in% names(mban.list))))

# set colours for plotting if only looking at vax
plans <- 
  if(all(intervention$mban == "none")) vax.plans else
    if(all(intervention$vax == "none")) mban.plans else
      unique(intervention$plan)
cols <- 
  c("black",
    brewer.pal(min(12, max(3, length(plans) - 1)), "Set1"))[1:length(plans)]

if(!is.null(cols)) {
  if(any(is.na(cols))) cols[is.na(cols)] <- brewer.pal(8, "Dark2")[1:sum(is.na(cols))]
}

names(cols) <- plans
cols[cols == "#FFFF33"] <- "#1E90FF"

head(intervention)


######################################
# loop over different policy options
start.time <- Sys.time()
input.list <- mclapply(intervention$k, function(k) {
  print(paste("Simulating events for scenario", k, "of", nrow(intervention)))
  print(intervention[k, ])
  
  ward.names.ns <-
    if(exclude.seed) ward.names[!ward.names %in% intervention$seed.ward[k]] else
      ward.names
  
  # select movement ban
  mban.all <- mban.list[[intervention$mban[k]]]
  if(intervention$mban[k] == "rand") {
    mban.all$ward <- sample(mban.all$ward)
  }
  
  # load parameters for simulating monthly movements
  par.list <-
    lapply(all.months, function(m) {
      # predicted probability of movement
      mu.z <- mu.z.list[[m]]
      
      # predicted number moved given movement
      mu.c <- mu.c.list[[m]]
      
      # dispersion (theta) parameter for the ZTNB distribution
      theta <- theta.list[[m]]
      
      # problem:
      # because the pink slips record only onward journeys from markets (n = 111), 
      # wards with no market (n = 287) appear to produce no animals. 
      # a simple solution is to assume that they send all their animals to
      # the nearest market. this could be improved by taking account of Gemma's
      # survey data, and also the trade-off between distance and size.
      # go with the simple solution for now.
      # identify wards that sent out no animals
      zero.out <- data.frame(ward = rownames(mu.z)[rowSums(mu.z) == 0], stringsAsFactors = FALSE)
      rowSums(mu.z[zero.out$ward, ]) # zero probability of outflow
      rowSums(mu.c[zero.out$ward, ]) # but nonzero rate given any outflow
      # for those wards that send out no animals 
      # assume that this is due to missing data and assume that all
      # animals from this ward went to the nearest active primary wards (active = outflow > 0) 
      zero.out$nearest.active.primary <- 
        sapply(zero.out$ward, function(w) {
          nap.index <- 
            which.nearest(
              unlist(wards[w, c("lat", "long")]),
              wards[!wards$fullname %in% zero.out$ward & wards$market == "Primary", 
                    c("lat", "long")])
          names(nap.index)
        })
      # how many primary markets have outflow now, but no inflow?
      rownames(mu.z)[rowSums(mu.z) > 0 & !rownames(mu.z) %in% zero.out$nearest.active.primary]
      wards[rownames(mu.z)[rowSums(mu.z) > 0 & !rownames(mu.z) %in% zero.out$nearest.active.primary], "prodsys"]
      # for now, just assume that all of these wards are supplied internally.
      # now distribute the outflow from each "nearest active primary" over 
      # the upstream wards with no recorded outflow
      for(nap in unique(zero.out$nearest.active.primary)) {
        mu.z[zero.out$ward[zero.out$nearest.active.primary == nap], nap] <- 1/sum(zero.out$nearest.active.primary == nap)
      }; rm(nap)
      sort(rowSums(mu.z))
      # diagonals should all be zero
      all(diag(mu.z) == 0)
      
      # apply movement ban via the movement probability matrix mu.z
      if(length(mban.all$ward) > 0 & intervention$mban[k] != "none" & m %in% mban.all$mban.months) {
        keep.wards <- mban.all$ward[mban.all$ward %in% ward.names.ns]
        keep.wards <- keep.wards[1:(min(length(ward.names.ns), max(1, intervention$n.mban[k])))]
        mu.z[keep.wards, ] <- 0
        mu.z[, keep.wards] <- 0
      }
      
      # multiply matrices to give expected flow
      mu <- mu.z * mu.c
      # check all wards now have outputs
      print(table(rowSums(mu) == 0))
      # output
      list(mu.z = mu.z, mu.c = mu.c, mu = mu, theta = theta)
    })
  names(par.list) <- all.months
  
  # remove secondary market wards from the contiguous wards network
  # because these wards are assumed to have no standing population of animals
  contig.wards <- contig.wards.all[contig.wards.all$from %in% not.mkts & contig.wards.all$to %in% not.mkts, ]
  
  # simulate balanced moves
  sim.moves.list <- 
    lapply(par.list, 
           function(mpar) 
             sim.hurdle.move(mu.z = mpar$mu.z, mu.c = mpar$mu.c, alpha = mpar$theta, is.mkt = is.mkt))
  
  # create events
  # start with event type "extTrans": moves between nodes
  
  # turn the move matrix list into a table of pairwise moves
  contig.move.n <- intervention$contig.move.n[k]
  evts.list <-
    lapply(names(sim.moves.list), 
           function(mth) {
             print(mth)
             #             date1 <- as.Date(paste("2015", mth, "01", sep = "-"))
             date1 <- all.months.day1[mth]
             move <- sim.moves.list[[mth]]
             
             # make a data frame of movements, one row per movement
             move.pair <- 
               expand.grid(fr.ward = rownames(move), 
                           to.ward = colnames(move), 
                           time = NA, n.moved = NA,
                           stringsAsFactors = FALSE)
             rownames(move.pair) <- paste(move.pair$fr.ward, move.pair$to.ward, sep = ".")
             
             
             # remove rows where no animals were moved
             # plus contiguous ward pairs, if any neighbour-neighbour movement 
             dim(move.pair)
             keep.pairs <- 
               apply(which(move > 0, arr.ind = TRUE), 1, 
                     function(i) paste(rownames(move)[i[1]], colnames(move)[i[2]], sep = "."))
             if(contig.move.n > 0) keep.pairs <- unique(c(keep.pairs, rownames(contig.wards)))
             move.pair <- move.pair[keep.pairs, ]
             dim(move.pair)
             
             # separate into the four different types of movement
             move.order <- c("wm", "mm", "mw", "ww")
             move.pair$move.type <- 
               factor(tolower(paste0(c("w", "m")[(move.pair$fr.ward %in% mkts) + 1], 
                                     c("w", "m")[(move.pair$to.ward %in% mkts) + 1])), 
                      move.order)
             
             for(m in move.order) {
               move.pair$time[move.pair$move.type == m] <- match(m, move.order) + date1
               move.pair$n.moved[move.pair$move.type == m] <- 
                 diag(
                   as.matrix(
                     move[move.pair$fr.ward[move.pair$move.type == m], move.pair$to.ward[move.pair$move.type == m]]))
               move.pair
             }; rm(m)
             move.pair <- move.pair[order(move.pair$time), ]
             
             # check that the same total number of animals is being moved
             sum(move.pair$n.moved) == sum(move)
             
             # simulate movements between contiguous wards
             # by moving an extra contig.move.n animals between all contiguous wards each month
             if(contig.move.n > 0) {
               move.pair[rownames(contig.wards), "n.moved"] <- 
                 move.pair[rownames(contig.wards), "n.moved"] + contig.move.n
             }
             
             #print(paste(nrow(contig.wards) * contig.move.n, sp, "moved between neighbours in month", mth))
             
             # output data frame of movement events
             if(nrow(move.pair) > 0.5) {
               move.df <- 
                 data.frame(
                   event = "extTrans", # options: exit (death), enter (birth), intTrans (e.g. S -> I), extTrans (movement)
                   # Events scheduled at same time processed inorder: exit, enter, internal trans, external trans
                   time = move.pair$time,
                   node = match(move.pair$fr.ward, ward.names),
                   dest = match(move.pair$to.ward, ward.names),
                   n = move.pair$n.moved,
                   proportion = 0,
                   select = 2, # chooses the column (compartment) of E the event operates on
                   shift = 0,
                   stringsAsFactors = FALSE)
             } else {
               move.df <- 
                 data.frame(event = character(0), time = integer(0), node = integer(0),
                            dest = integer(0), n = integer(0), proportion = integer(0),
                            select = integer(0), shift = integer(0)) 
             }
             
             
             # what proportion of movements are of 1 animal 
             print(paste0(round(100 * mean(move.df$n == 1)), "% of market movements have batch size = 1 in month ", mth))
             # measure changes in ward populations for each node
             active.nodes <- sort(unique(c(move.df$node, move.df$dest)))
             active.nodes <- active.nodes[!active.nodes %in% which(is.mkt)]
             inout.df.list <- 
               lapply(active.nodes, function(j) {
                 # n out
                 n.out <- sum(move.df$n[move.df$node == j & move.df$event == "extTrans"])
                 # n in
                 n.in <- sum(move.df$n[move.df$dest == j & move.df$event == "extTrans"])
                 imbalance <- n.in - n.out
                 # create births or deaths to correct imbalance in in/out-flow
                 if(n.out != n.in) {
                   demog.df <- 
                     data.frame(
                       event = ifelse(imbalance > 0, "exit", "enter"), # options: exit (death), enter (birth), intTrans (e.g. S -> I), extTrans (movement)
                       time = abs(max(sign(imbalance) * move.pair$time)) + sign(imbalance),
                       node = j,
                       dest = 0,
                       n = abs(imbalance),
                       proportion = 0,
                       select = ifelse(imbalance > 0, 2, 1), # chooses the column (compartment) of E the event operates on
                       shift = 0,
                       stringsAsFactors = FALSE)
                   return(demog.df)
                 } else return(NULL)
               })
             output.df <- rbind(move.df, do.call("rbind", inout.df.list))
             # check balance between inflow, outflow, births and deaths 
             imbalances <- 
               sapply(1:n, function(i) {
                 sum(output.df$n[output.df$event == "extTrans" & output.df$dest == i]) -
                   sum(output.df$n[output.df$event == "extTrans" & output.df$node == i]) +
                   sum(output.df$n[output.df$event == "enter" & output.df$node == i]) - 
                   sum(output.df$n[output.df$event == "exit" & output.df$node == i])
               })
             stopifnot(all(imbalances == 0))
             if(nrow(output.df) > 0.5) output.df$month <- mth
             output.df
           })
  sapply(evts.list, dim)
  
  # rbind the movements into a data frame
  evts <- do.call("rbind", evts.list)
  table(evts$time)
  
  # select vaccination schedule
  vax.df <- vax.list[[intervention$vax[k]]]
  
  # add vaccination event
  if(!all(vax.df$day == 0) & any(vax.df$day <= max(evts$time)) & any(vax.df$cov > 0)) {
    if(intervention$vax[k] == "rand") vax.df <- vax.df[sample(nrow(vax.df)), ]
    vax.df <- vax.df[vax.df$ward %in% ward.names.ns, ]
    vax.df <- vax.df[1:min(intervention$n.vax[k], nrow(vax.df)), ]
    vaccination <- 
      data.frame(event = "intTrans", 
                 time = vax.df$day, 
                 node = match(vax.df$ward, ward.names), 
                 dest = 0, 
                 n = 0, 
                 proportion = vax.df$cov,
                 select = 1, 
                 shift = 1,
                 month = NA,  # don't need month for vax events
                 stringsAsFactors = FALSE)
    evts <- rbind(evts, vaccination)
  }
  evts <- evts[order(evts$time), ]
  evts
  
  # inspect events table
  head(evts)
  dim(evts)
  table(evts$event) # 0 = exit; 1 = enter; 2 = internal transfer; 3 = external transfer
  
  # makes sure integer columns really are integer
  evts$time <- as.integer(evts$time)
  evts$node <- as.integer(evts$node)
  evts$dest <- as.integer(evts$dest)
  evts$n <- as.integer(evts$n)
  evts$select <- as.integer(evts$select)
  evts$shift <- as.integer(evts$shift)
  
  # output events table
  evts 
  
}, mc.cores = detectCores(), mc.preschedule = TRUE)

# check all input events objects are data frames 
stopifnot(sapply(input.list, class) == "data.frame")

# run model
sim.res.list <- lapply(intervention$k, function(k, ...) {
  print(paste("Simulating disease transmission for scenario", k, "of", nrow(intervention)))
  print(intervention[k, ])
  
  # which ward(s) to put infected animals into  
  seed.ward <- 
    if(substr(intervention$seed.ward[k], 1, 4) == "rand") {
      sample(not.mkts, gsub("rand", "", intervention$seed.ward[k]))
    } else
      if(intervention$seed.ward[k] == "balanced") {
        not.mkts.tab[intervention$rep[k] %% nrow(not.mkts.tab) + 1, ]
      } else
        intervention$seed.ward[k]
  
  # add infected animals to u0
  u0 <- u0.outer
  u0[seed.ward, "I"] <- u0[seed.ward, "Ic"] <- inf.pars["i0"] 
  u0$S <- u0$S - u0$I
  
  # get vaccination plan for scenario k 
  vax.df <- vax.list[[intervention$vax[k]]]
  vax.per.ward <- 0
  if(all(vax.df$day == 0) & !(intervention$vax[k] %in% c("none"))) {
    ward.names.ns <- ward.names[!ward.names %in% intervention$seed.ward[k]]
    vax.df <- vax.df[vax.df$ward %in% ward.names.ns, ]
    vax.df <- vax.df[1:min(intervention$n.vax[k], nrow(vax.df)), ]
    if(intervention$vax[k] == "rand") vax.df <- vax.df[sample(nrow(vax.df)), ]
    if(!(intervention$vax[k] %in% c("none", "all.lo", "all.hi"))) {
      vax.df <- vax.df[cumsum(u0[vax.df$ward, "S"] * vax.df$cov) < n.vax.dose, ]
    }
    vax.per.ward <- round(u0[vax.df$ward, "S"] * vax.df$cov)
    u0[vax.df$ward, "R"] <- u0[vax.df$ward, "R"] + vax.per.ward
    u0[vax.df$ward, "S"] <- u0[vax.df$ward, "S"] - vax.per.ward
  }
  n.ward.vax <- nrow(vax.df[vax.df$ward != "none", ])
  vax.used <- sum(vax.per.ward)
  if(is.na(vax.used)) vax.used <- 0
  print(intervention$vax[k])
  print(paste(n.ward.vax, "wards vaccinated.", vax.used, "doses used."))
  
  # set up SIR model
  mod.input <- split.nodes.sir(n.subnodes = n.subnodes, u0 = u0)
  
  model <- mparse(transitions = mod.input$transitions,
                  compartments = rownames(mod.input$E),
                  gdata = c(beta = intervention$beta[k], gamma = intervention$gamma[k], coupling = coupling),
                  u0 = mod.input$u0,
                  tspan = 1:n.days,
                  events = input.list[[k]],
                  #ldata = matrix(rep(intervention$beta[k], nrow(wards)), nrow = 1, dimnames = list("beta", NULL)),
                  E = mod.input$E,
                  N = cbind(mod.input$E[, 1]) * n.subnodes * 2)  # determine how vaccination (col 1) happens: it shifts susceptibles down to R
  
  # run the simulation
  result <- run(model, threads = 1)
  
  #plot(result, node = which(u0$Ic > 0)[1], compartments = "Ic", main = n.subnodes)
  #plot(result, node = sample(nrow(u0), 1), compartments = "Ic", main = n.subnodes)
  #plot(result, node = which(u0$Ic > 0), compartments = "Ic", main = n.subnodes, range = FALSE)
  #u0[which(u0$Ic > 0)[1], ]
  
  #plot(run(SIR(u0 = u0[which(u0$Ic > 0), ], tspan = 1:n.days, beta = inf.pars$beta[2], gamma = intervention$gamma[k])), range = FALSE)
  
  # plot result (only use for debugging, hence FALSE)
  result.prop <- result
  if(FALSE) {
    if(paste(compartments, collapse = "") == "SEIIcRV") {
      for(i in 6*(1:n)) {
        result.prop@U[(i-5):i, ][-4, ] <- prop.table(result.prop@U[(i-5):i, ][-4, ], 2)
        result.prop@U[(i-5):i, ][4, ] <- 0
      }; rm(i)
    }
    if(paste(compartments, collapse = "") == "SIIcR") {
      for(i in 4*(1:n)) {
        result.prop@U[i-(3:1), ] <- prop.table(result.prop@U[i-(3:1), ], 2)
        result.prop@U[i, ] <- 0
      }; rm(i)
    }
    plot(result.prop, range = FALSE, node = which(!is.mkt), lwd = 1)
  }
  
  event.tab <- table(model@events@event)
  print(paste("There have been", event.tab[1], "deaths,", event.tab[2], "births, and", 
              event.tab[3], "movements"))
  
  # extract the total cumulative incidence across all wards and the total number of animals
  use.wards <- match(not.mkts, rownames(u0))
  sir.comp <- mod.input$compartments[mod.input$compartments != "Ic"]
  i.comp <- sir.comp[substring(sir.comp, 1, 1) == "I"]
  cumInc.full <- sapply(use.wards, function(w) trajectory(result, node = w)$Ic)
  cumInc.full.prop <- 
    sapply(use.wards, function(w) {
      trj <- trajectory(result, node = w)
      out <- trj$Ic / rowSums(trj[, sir.comp])
      out[is.na(out) | is.infinite(out)] <- 0
      out
    })
  Inc.full <- sapply(use.wards, function(w) rowSums(cbind(trajectory(result, node = w)[, i.comp])))
  Inc.full.prop <- 
    sapply(use.wards, function(w) {
      trj <- trajectory(result, node = w)
      out <- rowSums(cbind(trj[, i.comp])) / rowSums(trj[, sir.comp])
      out[is.na(out) | is.infinite(out)] <- 0
      out
    })
  Inc <- rowSums(Inc.full)
  cumInc <- rowSums(cumInc.full)
  N <- 
    rowSums(sapply(use.wards, function(w)
      rowSums(trajectory(result, node = w)[, sir.comp])))
  out <-
    cbind(
      intervention[k, ],
      day = 1:length(cumInc),
      cumInc = cumInc, 
      Inc = Inc, 
      N = N, 
      cumInc.ward = rowSums(cumInc.full.prop > rep.thresh),
      Inc.ward = rowSums(Inc.full.prop > rep.thresh),
      N.ward = n,
      cumInc.district = sapply(1:n.days, function(d) length(unique(wards$district[cumInc.full[d, ] > 0]))),
      N.district = length(unique(wards$district)),
      row.names = NULL)
  attr(out, "cumInc.full.prop") <- cumInc.full.prop
  attr(out, "Inc.full.prop") <- Inc.full.prop
  attr(out, "n.ward.vax") <- n.ward.vax
  attr(out, "vax.used") <- vax.used
  out
}, mc.cores = detectCores(), mc.preschedule = TRUE)


# identify and remove failed runs
intervention$failed <- sapply(sim.res.list, class) != "data.frame"
intervention[intervention$failed, ]
if(any(intervention$failed)) warning(paste(sum(intervention$failed), "of", nrow(intervention), "disease simulation failed"))

# add vaccination stats to the intervention table
intervention$n.ward.vax <- sapply(sim.res.list, attr, "n.ward.vax")
intervention$vax.used <- sapply(sim.res.list, attr, "vax.used")

summary(intervention$vax.used[!intervention$vax %in% c("none", "all.hi", "all.lo")])

sim.res.list <- sim.res.list[which(!intervention$failed)]
intervention <- intervention[!intervention$failed, ]

for(use.beta in inf.pars$beta) {
  
  #####################################################
  # rbind the simulation results list into a data frame
  sim.res <- do.call("rbind.data.frame", sim.res.list[intervention$k[intervention$beta == use.beta]])
  unique(sim.res$k)
  unique(intervention$k)
  
  # estimate reduction due to interventions relative to do-nothing
  # stratified by rep (because each rep used the same set of starting wards)
  res.tab.red <- 
    sapply(1:nrep, function(r) {
      sr <- sim.res[sim.res$rep == r & sim.res$day == n.days, ]
      out <- log(sr$cumInc / sr$cumInc[sr$vax == "none" & sr$mban == "none"])
      names(out) <- paste0("VX: ", sr$vax, "; MB: ", sr$mban, "; R0: ", round(sr$beta/sr$gamma, 2))
      out
    })
  
  # get mean cumulative incidence as a fraction on the total N cattle
  res.tab.cuminc <- 
    sapply(1:nrep, function(r) {
      sr <- sim.res[sim.res$rep == r & sim.res$day == n.days, ]
      out <- sr$cumInc / sr$N
      names(out) <- paste0("VX: ", sr$vax, "; MB: ", sr$mban, "; R0: ", round(sr$beta/sr$gamma, 2))
      out
    })
  
  # make table of results for output
  res.tab <- 
    data.frame(
      cumInc = apply(res.tab.cuminc, 1, mean),
      cumInc.rel.log = apply(res.tab.red, 1, mean))
  res.tab$cumInc.red <- 1 - exp(res.tab$cumInc.rel.log)
  res.tab$cumInc.red.se <- exp(res.tab$cumInc.rel.log) * apply(res.tab.red, 1, function(x) sd(x)/sqrt(length(x) - 1))
  res.tab$cumInc.rel.log <- NULL
  print(res.tab)
  
  # average across replicate simulations
  sim.res.mean <- 
    do.call("rbind",
            lapply(plans[plans %in% sim.res$plan], function(p) {
              dat <- sim.res[sim.res$plan == p, ]
              dat.mean <- dat[dat$rep == min(dat$rep), ]
              for(nm in names(dat.mean)[!names(dat.mean) %in% c("day", names(intervention))]) {
                dat.mean[, nm] <- apply(sapply(unique(dat$rep), function(r) dat[dat$rep == r, nm]), 1, mean)
              }
              dat.mean$k <- dat$k[dat$rep == min(dat$rep)]
              stopifnot(dat.mean$day == 1:n.days)
              dat.mean
            }))
  
  sim.res <- sim.res.mean # including this line hides the individual results and only plots the mean 
  
  if(is.null(cols)) {
    cols <- 
      c("black", 
        brewer.pal(min(12, max(3, length(plans) - 1)), "Set1"))[1:length(plans)]
    if(any(is.na(cols))) cols[is.na(cols)] <- brewer.pal(8, "Dark2")[1:sum(is.na(cols))]
    names(cols) <- plans
  }
  
  ltys <- rep(1, length(plans))
  ltys[-grep("MB none", plans)] <- 2
  names(ltys) <- plans
  lwd <- c(0.4, 2)
  alpha.trans <- c(1, 1)
  
  pdf(paste0("inc.tab.beta", round(use.beta, 3) ,".pdf"), height = 11.69, width = 8.27)
  old.par <- par(mar = c(3, 4, 1.5, 1) + 0.1, par(bg = grey(1)))
  plot.layout <- layout(rbind(1, 2))
  plot(1, type = "n", 
       ylim = c(0, min(100, 100 * max(sim.res$cumInc.ward/sim.res$N.ward))), 
       #ylim = c(0, 21), #c(0, min(100, 100 * max(sim.res$cumInc.ward/sim.res$N.ward))), 
       xlim = c(1, n.days),
       ylab = paste0("Proportion (%) of wards with > ", 100 * rep.thresh, "% animals infected"), 
       xlab = "")
  mtext("Day", 1, 2)
  legend("topleft", col = cols[plans], lwd = lwd[2], legend = plans, 
         lty = 1, bty = "n")
  sim.res.plot.list <- list(sim.res, sim.res.mean)
  lapply(1:length(sim.res.plot.list), function(i) {
    x <- sim.res.plot.list[[i]]
    invisible(lapply(unique(x$k), function(k) {
      print(intervention[intervention$k == k, ])
      cumInc.tab <- x[x$k == k,]
      print(cumInc.tab$cumInc.ward[n.days])
      lines(1:nrow(cumInc.tab), 100 * cumInc.tab$cumInc.ward/cumInc.tab$N.ward, 
            col = alpha(cols[intervention$plan[intervention$k == k]], alpha.trans[i]), 
            lwd = lwd[i], 
            lty = 1)
      lines(1:nrow(cumInc.tab), 100 * cumInc.tab$Inc.ward/cumInc.tab$N.ward, 
            col = alpha(cols[intervention$plan[intervention$k == k]], alpha.trans[i]), 
            lwd = lwd[i], 
            lty = 2)
      
    }))
  })
  
  
  plot(1, type = "n", 
       #ylim = c(0, 11), # c(0, min(100, 100 * max(sim.res$cumInc/sim.res$N))), 
       ylim = c(0, min(100, 100 * max(sim.res$cumInc/sim.res$N))), 
       xlim = c(1, n.days),
       ylab = "Proportion (%) of animals infected", xlab = "")
  mtext("Day", 1, 2)
  #title("Population cumulative incidence")
  legend("topleft", col = cols[plans], lwd = lwd[2], legend = plans, 
         lty = 1, bty = "n")
  
  lapply(1:length(sim.res.plot.list), function(i) {
    x <- sim.res.plot.list[[i]]
    invisible(lapply(unique(x$k), function(k) {
      print(intervention[intervention$k == k, ])
      cumInc.tab <- x[x$k == k,]
      # how many districts and cases after 6 months
      time.points <- c(halfway = ceiling(n.days/2), end = n.days)
      lapply(names(time.points), 
             function(tp) {
               print(paste("After", time.points[tp], "days", cumInc.tab$cumInc.district[time.points[tp]], 
                           "districts have at least one case, and there were a total of", 
                           cumInc.tab$cumInc[time.points[tp]], "cases"))
             })
      print(cumInc.tab$cumInc.ward[365])
      lines(1:nrow(cumInc.tab), 100 * cumInc.tab$cumInc/cumInc.tab$N, 
            col = alpha(cols[intervention$plan[intervention$k == k]], c(alpha.trans, 1)[i]), 
            lwd = lwd[i], 
            lty = 1)
      lines(1:nrow(cumInc.tab), 100 * cumInc.tab$Inc/cumInc.tab$N, 
            col = alpha(cols[intervention$plan[intervention$k == k]], c(alpha.trans, 1)[i]), 
            lwd = lwd[i], 
            lty = 2)
    }))
  })
  
  dev.off()
  
  # assess efficacy of interventions
  # export example simulations
  
  # which simulations to export
  # pick a run that gives close to the average ward prevalence on the last day
  sim.res.end <- sim.res[sim.res$day == n.days & sim.res$mban == "none" & sim.res$vax == "none", ]
  
  export.k <- 
    sim.res.end$k[which.min(abs(sim.res.end$cumInc.ward - round(mean(sim.res.end$cumInc.ward))))[1]]
  
  # incidence (wards stay "infected") or infectious status?
  output.incidence.tab <- FALSE
  
  # cut points for converting day to month
  cut.month <- c(all.months.day1 - 0.5, n.days + 0.5)
  
  lapply(export.k, function(k) {
    # export movements
    evts <- input.list[[k]]
    lapply(sort(unique(na.omit(evts$month))), function(mth) {
      evts <- evts[evts$event %in% "extTrans" & evts$month %in% mth, ]
      mvmat <- 
        matrix(rep(0, length(all.ward.names)^2), 
               ncol = length(all.ward.names),
               dimnames = list(all.ward.names, all.ward.names))
      for(w in ward.names) {
        node <- match(w, ward.names)
        dest <- evts$dest[evts$node == node]
        n.moved <- evts$n[evts$node == node]
        mvmat[w, ward.names[dest]] <- n.moved
      }
      outfile <-
        paste("networks/simmoves.monthly.2015.RVF", 
              sp, paste0(sp, ".month", mth, ".mvt.matrix.simmoves.", k, ".csv"),
              #sp, paste0(sp, ".month", mth, ".mvt.matrix.simmoves.spatial.csv"),
              sep = "/")
      print(table(rowSums(mvmat) == 0))
      write.csv(mvmat, outfile)
    })
    
    # export binary infected/not = 1/0 matrix of wards X months
    cumInc.full.prop <-   attr(sim.res.list[[k]], "Inc.full.prop")
    colnames(cumInc.full.prop) <- not.mkts
    # add back on the wards that were dropped due to too few animals
    for(w in all.ward.names[!all.ward.names %in% colnames(cumInc.full.prop)]) {
      cumInc.full.prop <- cbind(cumInc.full.prop, 0)
      colnames(cumInc.full.prop)[colnames(cumInc.full.prop) == ""] <- w
    }; rm(w)
    cumInc.full.prop <- cumInc.full.prop[, all.ward.names]
    start.day <- apply(cumInc.full.prop > rep.thresh, 2, function(x) min(c(Inf, which(x))))
    
    inf.tab <- 
      t(apply(cumInc.full.prop > rep.thresh, 2, 
              function(x) all.months %in% 
                unique(cut(which(x),  breaks = cut.month, labels = all.months))) + 0)
    colnames(inf.tab) <- all.months
    if(output.incidence.tab) {
      inf.tab <- t(apply(inf.tab, 1, cumsum) > 0) + 0
    }
    colnames(inf.tab) <- paste0("month", colnames(inf.tab))
    outfile <-
      paste("networks/simmoves.monthly.2015.RVF", 
            sp, paste0(sp, ".infected.ward.tab.", k, ".csv"),
            sep = "/")
    write.csv(inf.tab, outfile)
  })
  
  
  
  if(FALSE) {
    invisible(lapply(unique(sim.res$k), function(k) {
      print(intervention[k, ])
      cumInc.tab <- sim.res[sim.res$k == k,]
      # how many districts and cases after 6 months
      time.points <- c(six = 183, twelve = 365)
      lapply(names(time.points), 
             function(tp) {
               print(paste("After", tp, "months", cumInc.tab$cumInc.district[time.points[tp]], 
                           "districts have at least one case, and there were a total of", 
                           cumInc.tab$cumInc[time.points[tp]], "cases"))
             })
      print(cumInc.tab$cumInc.ward[365])
      lines(1:nrow(cumInc.tab), 100 * cumInc.tab$cumInc/cumInc.tab$N, 
            col = alpha(cols[intervention$plan[k]], alpha.trans), lwd = lwd, 
            lty = ltys[intervention$plan[k]])
      
    }))
  }
  round(res.tab, 2)
  write.csv(cbind(intervention[intervention$rep == 1 & intervention$beta == use.beta, ], round(res.tab, 5)), paste0("phil.trans.results.beta", round(use.beta, 3), ".csv"),
            row.names = FALSE)
}

finish.time <- Sys.time()
print(finish.time - start.time)

#rm(input.list)
#save.image("SEEDZ_pinkslip_diseasesim_PostArusha_v02DRAFT.RData")
