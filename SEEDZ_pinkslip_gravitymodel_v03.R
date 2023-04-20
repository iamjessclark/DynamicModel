# gravity-hurdle model of preliminary pink slip data
# started 25th July 2017

# clear objects from memory, reset graphics devices
rm(list = ls())
graphics.off()
Sys.setenv(TZ="Europe/London")

library(pscl)
library(splines)
library(piecewiseSEM)
library(parallel)
library(GLMMmisc)
library(scales)
library(ggmap)
library(MASS)
library(glmmTMB)

# functions

# show 20 random rows of a table
sample.rows <- function(x, n = 20) x[sort(sample(nrow(x), n)), ]


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



# wald z-test function 
# (adapted from http://stijnruiter.nl/blog/?p=309 - broken link)
# useful for GLMMs when you want to test for inclusion in the model of a 
# categorical variable that has more than 2 levels. in a linear regression model 
# you could use a likelihood ratio test, but this is not recommended for GLMMs.
# However this Wald test method is also dubious at low n because it doesn't take account of the 
# degrees of freedom.
# (see Bolker et al. Trends in Ecology & Evolution, Volume 24, Issue 3, 127-135, 29 January 2009).
# 2017-07-02: adapted to deal with glmmTMB fits

Wald<-
  function(object, R, q = NULL, gmmTMB.model = "cond")
  {
    require(lme4)
    if (!is.matrix(R)) stop("Restrictions must be a matrix")
    if(is.null(q)) q <- rep(0, nrow(R))
    b <- fixef(object)
    if(class(b) == "fixef.glmmTMB") b <- b[[gmmTMB.model]]
    vc <- vcov(object)
    if("vcov.glmmTMB" %in% class(vc)) vc <- vc[[gmmTMB.model]]
    w <- t(R %*% b - q) %*% solve(R %*% vc %*% t(R)) %*% (R %*% b - q)
    pw <- pchisq(w[1], length(q), lower.tail = FALSE)
    cat("*************\n* Wald Test *\n*************\n")
    cat("lme4 fixed effects:\n")
    print(fixef(object))
    cat("\nRestrictions:\n")
    print(R)
    cat("\nq = ",q)
    cat("\nChi-square:", round(w[1],3), " df = ", length(q))
    cat("\nProb x>chisq:", round(pw, 5), "\n")
    return(pw)
  }


# funtion to adjust binomial GLMM for jensen's inequality 
jensen.adjust <-
  function(p, V, method = "mcculloch", inverse = FALSE) {
    stopifnot(!(method == "mcculloch" & inverse))
    Beta <- qlogis(p)
    if(method == "mcculloch") {
      return(plogis(Beta - 0.5 * V * tanh(Beta * (1 + 2 * exp(-0.5 * V))/6)))
    }
    if(method == "zeger") {
      if(inverse) 
        return(plogis(Beta * sqrt(256 * V / (75 * pi^2) + 1))) else
          return(plogis(Beta/sqrt(1 + ((16 * sqrt(3))/(15 * pi))^2 * V)))
    }
  }


# function to capitalise initial letters of words.
# this function is adapted from the examples in the R Help for toupper()
capwords<-
  function(s,strict=FALSE,first.only=FALSE)
  {
    cap<-
      function(s)
        paste(
          toupper(substring(s,1,1)),
          {s <- substring(s,2); if(strict) tolower(s) else s},
          sep = "", collapse = " " )
    if(first.only)
      cap(s)
    else
      sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
  }


# read in ward data
wards <- 
  read.csv("~/OneDrive - University of Glasgow/Projects/SEEDZ/mapping/WardSpatialData.csv", 
           stringsAsFactors = FALSE, quote = "")
dim(wards)
dim(unique(wards))

head(wards)
sum(is.na(wards))

for(n in c("name", "district", "region")) {
  wards[, n] <- gsub(" ", "", gsub("[^[:alpha:] ]", "", tolower(wards[, n])))
  wards[!grepl('^[A-Za-z]+$', wards[, n]), n] <- NA
} ; rm(n)

sum(is.na(wards))

wards$name <- paste(wards$region, wards$district, wards$name, sep = "/")

head(wards)

table(wards$region)

# select study regions
include.regions <- c("arusha", "manyara", "kilimanjaro")
wards <- droplevels( wards[wards$region %in% include.regions, ])
length(wards$name)
length(unique(wards$name))

# load movement data
file.info("movt.analysis.data.csv")
dat.all <- read.csv("movt.analysis.data.csv", stringsAsFactors = FALSE)
dim(dat.all)
head(dat.all)
apply(is.na(dat.all), 2, sum)
dat.all[is.na(dat.all$date_yyyy), ]

# year
dat.all$year <- factor(dat.all$date_yyyy)

# year and month
dat.all$year.month <- 
  factor(paste(dat.all$date_yyyy, formatC(dat.all$date_mm, flag = "0", width = 2), sep = "-"))
grep("NA", dat.all$year.month)
dat.all$month <- factor(dat.all$date_mm)

# half-year and year
dat.all$half <- factor(dat.all$date_mm < 6.5, c(TRUE, FALSE), c("1st", "2nd"))
table(dat.all$half, exclude = NULL)
dat.all$year.half <- factor(paste(dat.all$year, dat.all$half, sep = "-"))
table(dat.all$year.half, exclude = NULL)

# quarter-year and year
dat.all$quarter <- 
  cut(dat.all$date_mm, c(0, 3.5, 6.5, 9.5, 12.5), 
      labels = c("1st", "2nd", "3rd", "4th"))
table(dat.all$quarter)
plot(table(dat.all$date_mm))
table(dat.all$quarter, exclude = NULL)
dat.all$year.quarter <- factor(paste(dat.all$year, dat.all$quarter, sep = "-"))
table(dat.all$year.quarter, exclude = NULL)

# clear up NAs
dat.all$year[grep("NA", dat.all$year)] <- NA
dat.all$year <- droplevels(dat.all$year)
levels(dat.all$year)
grep("NA", dat.all$year.month)
dat.all$year.month[grep("NA", dat.all$year.month)] <- NA
dat.all$year.month <- droplevels(dat.all$year.month)
grep("NA", dat.all$year.month)
levels(dat.all$year.month)
dat.all$year.half[grep("NA", dat.all$year.half)] <- NA
dat.all$year.half <- droplevels(dat.all$year.half)
levels(dat.all$year.half)
dat.all$year <- droplevels(dat.all$year)
grep("NA", dat.all$year.half)
table(dat.all$year.quarter, exclude = NULL)
dat.all$year.quarter[grep("NA", dat.all$year.quarter)] <- NA
dat.all$year.quarter <- droplevels(dat.all$year.quarter)
levels(dat.all$year.quarter)
grep("NA", dat.all$year.quarter)
table(dat.all$year.quarter, exclude = NULL)

head(dat.all)


# choose time period to aggregate over
dat.all$period <- dat.all$year.month


plot(table(dat.all$year.month), las = 3, cex.axis = 0.7)

dim(dat.all)
dat <- 
  droplevels(
    dat.all[dat.all$fr.match.type == "exact" & 
              dat.all$to.match.type == "exact" &
              !is.na(dat.all$period), ])
dat[sample(nrow(dat), 10), ]
rm(dat.all)
dim(dat)
apply(is.na(dat), 2, sum)
apply(is.na(dat), 2, mean)

levels(dat$year.month)

# type of transport
table(dat$to.match.type)
dat$transport <- 
  factor(dat$transport, exclude = NULL, labels = c("Foot", "Motor", "Not rec."))
dat$species <- factor(dat$species, c("Cattle", "Sheep/goats", "Mixed"))
dat[sample(nrow(dat), 10), ]

table(is.na(dat$date))
table(is.na(dat$period))

# impute zero distances as 0.5 km
sort(unique(dat$distance))[1:10]
dat[dat$distance == 0, ]
sort(dat$distance)[1:200]
dat$distance[dat$distance == 0] <- 0.5


# create data set with all possible movements from the origin wards (complete list)
# to the destinations
fr.wards <- unique(dat$fr.ward)
to.wards <- unique(wards$name)

setdiff(fr.wards, to.wards) # should return character(0)

# how large (max) will the data set be?
length(fr.wards) * length(to.wards) * nlevels(dat$period)

# make template data set of all possible ward movements in each period
all(fr.wards %in% to.wards)

datw.full <- 
  expand.grid(fr.ward = to.wards, to.ward = to.wards, period = levels(dat$period))
datw.full$month <- as.numeric(substr(datw.full$period, 6, 7)) - 1 # -1 so that Jan is ref
datw.full$year <- factor(substr(datw.full$period, 1, 4))
str(datw.full)
dim(datw.full)

# remove rows where origin = destination; we're not interested in within ward movement
datw.full <- 
  droplevels(datw.full[as.character(datw.full$fr.ward) != as.character(datw.full$to.ward), ])
dim(datw.full)

# give rows unique names including origin, destination and period
datw.full$fr.to.date <- 
  factor(paste(datw.full$fr.ward, datw.full$to.ward, datw.full$period, sep = "."))
rownames(datw.full) <- datw.full$fr.to.date

# make a factor for market-period too
datw.full$fr.date <- 
  factor(paste(datw.full$fr.ward, datw.full$period, sep = "."))

# create the same variable in the source data set to allow matching
dat$fr.to.date <- 
  paste(dat$fr.ward, dat$to.ward, dat$period, sep = ".")


# add origin and destination labels (these are the names that are presented in tables and figures, 
# e.g. "kilimanjaroairport" rather than "kilimanjarointernationalairport")
#datw.full[1, ]
#match.label <- 
#  unique(cbind(c(dat$fr.ward, dat$to.ward), c(dat$fr.match.label, dat$to.match.label)))
#dim(match.label)
#head(match.label)
#rownames(match.label) <- match.label[, 1]

#datw.full$fr.match.label <- dat$fr.match.label[match(as.character(datw.full$fr.to.date), dat$fr.to.date)]
#table(is.na(datw.full$fr.match.label))
#dat$fr.match[dat$fr.match.label == "kilimanjaroairport"]
#dat$to.match.label

# number of animals moved of each species
dat$n.permits <- 1
spp <- c("animals", "cattle", "caprine")
animals.sum.tab <- 
  do.call("rbind", 
          by(dat[dat$fr.to.date %in% datw.full$fr.to.date, c(spp, "n.permits")], 
             dat$fr.to.date[dat$fr.to.date %in% datw.full$fr.to.date], 
             colSums))
head(animals.sum.tab)
datw.full[, c(spp, "n.permits")] <- 0
datw.full[rownames(animals.sum.tab), c(spp, "n.permits")] <- animals.sum.tab[, c(spp, "n.permits")]
rm(animals.sum.tab)

# note that now not all fr.wards produce movemnents in datw.full, because datw.full
# excludes origins that only send animals out of the 3 study regions
# e.g.
dat[dat$fr.ward == "manyara/hanang/simbay", ]
sum(datw.full$animals[datw.full$fr.ward == "manyara/hanang/simbay"])

# examine animals variable in datw.full
summary(datw.full$animals)
sum(is.na(datw.full$animals))
sum(dat$animals)
sum(datw.full$animals) # fewer, as movements into and out of the 3 regions excluded
plot(tapply(datw.full$animals, datw.full$period, sum))
mean(datw.full$animals)
mean(datw.full$animals[datw.full$animals > 0.5])
mean(datw.full$cattle[datw.full$cattle > 0.5])
sum(datw.full$cattle)/sum(datw.full$animals)
table(datw.full$n.permits)

sum(datw.full$n.permits == 1) / sum(datw.full$n.permits > 0)
# ...so if we break down movements into months, 58% of movements are single permits 

# add human and livestock population size, population density, area, status

# define continuous variables
cont.vars <- 
  c("area.km2", "pop2012", "popdens2012",
    paste0("pop.", spp), paste0("popdens.", spp))
wards[1:3, cont.vars]

# add ward status from the census
datw.full$fr.status <- 
  factor(wards$status[match(datw.full$fr.ward, wards$name)], 
         c("Rural Ward", "Mixed Ward", "Urban Ward"))
datw.full$to.status <- 
  factor(wards$status[match(datw.full$to.ward, wards$name)], 
         c("Rural Ward", "Mixed Ward", "Urban Ward"))

table(datw.full$fr.status, datw.full$to.status, exclude = NULL)

# add Will de Glanville's production system classification
datw.full$fr.prodsys <- 
  factor(wards$prodsys[match(datw.full$fr.ward, wards$name)], 
         c("Pastoral", "Agropastoral", "Smallholder", "Urban"))
datw.full$to.prodsys <- 
  factor(wards$prodsys[match(datw.full$to.ward, wards$name)], 
         c("Pastoral", "Agropastoral", "Smallholder", "Urban"))

table(datw.full$fr.prodsys, datw.full$to.prodsys, exclude = NULL)
table(datw.full$fr.status, datw.full$fr.prodsys, exclude = NULL)
table(datw.full$to.status, datw.full$to.prodsys, exclude = NULL)
table(wards$prodsys, wards$status)

# identify wards containing primary and secondary markets

# binary indicators for presence of 1ary and 2ary markets in origin and destination wards
datw.full$fr.primary <- as.integer(wards$market[match(datw.full$fr.ward, wards$name)] == "Primary")
datw.full$fr.secondary <- as.integer(wards$market[match(datw.full$fr.ward, wards$name)] == "Secondary")
table(datw.full$fr.primary, datw.full$fr.secondary)
datw.full$to.primary <- as.integer(wards$market[match(datw.full$to.ward, wards$name)] == "Primary")
datw.full$to.secondary <- as.integer(wards$market[match(datw.full$to.ward, wards$name)] == "Secondary")
table(datw.full$to.primary, datw.full$to.secondary)

# factors for presence of 1ary and 2ary markets in origin and destination wards
datw.full$to.market.f <- wards$market[match(datw.full$to.ward, wards$name)]
datw.full$fr.market.f <- wards$market[match(datw.full$fr.ward, wards$name)]
table(datw.full$fr.market.f)
table(datw.full$to.market.f)
table(datw.full$fr.market.f, datw.full$to.market.f)

# factor for type of journey: primary to primary, primary to secondary, etc
datw.full$market.f <- factor(paste(datw.full$fr.market.f, datw.full$to.market.f, sep = "-"))
levels(datw.full$market.f)
table(datw.full$market.f)

# factors for presence of 1ary and 2ary markets in origin and destination wards
tapply(datw.full$animals, datw.full$market.f, mean)
tapply(datw.full$animals > 0, datw.full$market.f, sum)
round(100 * tapply(datw.full$animals > 0, datw.full$market.f, mean), 2)


# add variable with distance to primary and secondary market
datw.full$fr.primary5 <- as.integer(wards$mindist.primary[match(datw.full$fr.ward, wards$name)] < 5)
mean(datw.full$fr.primary5)
datw.full$fr.primary5.f <- factor(datw.full$fr.primary5, 0:1, c("<5km", ">5km"))
datw.full$to.primary5 <- as.integer(wards$mindist.primary[match(datw.full$to.ward, wards$name)] < 5)
mean(datw.full$to.primary5)
datw.full$to.primary5.f <- factor(datw.full$to.primary5, 0:1, c("<5km", ">5km"))
datw.full$fr.secondary10 <- as.integer(wards$mindist.secondary[match(datw.full$fr.ward, wards$name)] < 10)
mean(datw.full$fr.secondary10)
datw.full$fr.secondary10.f <- factor(datw.full$fr.secondary10, 0:1, c("<10km", ">10km"))
datw.full$to.secondary10 <- as.integer(wards$mindist.secondary[match(datw.full$to.ward, wards$name)] < 10)
mean(datw.full$to.secondary10)
datw.full$to.secondary10.f <- factor(datw.full$to.secondary10, 0:1, c("<10km", ">10km"))
sum(datw.full$to.secondary10)

table(datw.full$market.f, datw.full$fr.primary5.f)
table(datw.full$market.f, datw.full$to.primary5.f)
table(datw.full$market.f, datw.full$fr.secondary10.f)
table(datw.full$market.f, datw.full$to.secondary10.f)

# how do total cattle numbers compare with LINKS total volume?
# quick check of secondary markets (meserani, themi, weruweru) from LINKS (http://www.lmistz.net).
# also plot distribution of data over the period to identify missing chunks of time

# make a secondary markets plot comparing LINKS and pink slip data over years
secmkt <-
  expand.grid(
    year = levels(dat$year),
    market = c("meserani", "themi", "weruweru", "dosidosi", "mgagao"))

secmkt$cattle.permits.est <- 
  sapply(1:nrow(secmkt), function(i) {
    sum(dat$cattle[dat$fr.match.label == secmkt$market[i] & dat$year == secmkt$year[i]])
  }) * 20 / length(unique(dat$batch))

secmkt$permits.months <- 
  sapply(1:nrow(secmkt), function(i) {
    length(unique(na.omit(dat$date_mm[dat$fr.match.label == secmkt$market[i] & dat$year == secmkt$year[i]])))
  })



secmkt$cattle.links[secmkt$year == "2009" & secmkt$market == "meserani"] <- 16271 # 10
secmkt$cattle.links[secmkt$year == "2011" & secmkt$market == "meserani"] <- 26905 # 12
secmkt$cattle.links[secmkt$year == "2013" & secmkt$market == "meserani"] <- 27748 # 12
secmkt$cattle.links[secmkt$year == "2015" & secmkt$market == "meserani"] <- 18635 # 11

secmkt$cattle.links[secmkt$year == "2009" & secmkt$market == "themi"] <- NA
secmkt$cattle.links[secmkt$year == "2011" & secmkt$market == "themi"] <- 6949 # 12
secmkt$cattle.links[secmkt$year == "2013" & secmkt$market == "themi"] <- 6870 # 12
secmkt$cattle.links[secmkt$year == "2015" & secmkt$market == "themi"] <- 3449 # 10

secmkt$cattle.links[secmkt$year == "2009" & secmkt$market == "weruweru"] <- 1209 # 3 months
secmkt$cattle.links[secmkt$year == "2011" & secmkt$market == "weruweru"] <- 2206 # 4 months
secmkt$cattle.links[secmkt$year == "2013" & secmkt$market == "weruweru"] <- NA
secmkt$cattle.links[secmkt$year == "2015" & secmkt$market == "weruweru"] <- NA

secmkt$cattle.links[secmkt$year == "2009" & secmkt$market == "dosidosi"] <- 8436 # 9
secmkt$cattle.links[secmkt$year == "2011" & secmkt$market == "dosidosi"] <- 3681 # 4
secmkt$cattle.links[secmkt$year == "2013" & secmkt$market == "dosidosi"] <-  748 # 1 month
secmkt$cattle.links[secmkt$year == "2015" & secmkt$market == "dosidosi"] <- 8510 # 6 months

secmkt$cattle.links[secmkt$year == "2009" & secmkt$market == "mgagao"] <- 5120 # 12
secmkt$cattle.links[secmkt$year == "2011" & secmkt$market == "mgagao"] <- 5080 # 12
secmkt$cattle.links[secmkt$year == "2013" & secmkt$market == "mgagao"] <- 3775 # 11 months
secmkt$cattle.links[secmkt$year == "2015" & secmkt$market == "mgagao"] <- 4734 # 12 months of records


# missingness in links -- how many months have any total cattle numbers recorded?

secmkt$links.months <- NA

secmkt$links.months[secmkt$year == "2009" & secmkt$market == "meserani"] <- 10
secmkt$links.months[secmkt$year == "2011" & secmkt$market == "meserani"] <- 12
secmkt$links.months[secmkt$year == "2013" & secmkt$market == "meserani"] <- 12
secmkt$links.months[secmkt$year == "2015" & secmkt$market == "meserani"] <- 11

secmkt$links.months[secmkt$year == "2009" & secmkt$market == "themi"] <- NA
secmkt$links.months[secmkt$year == "2011" & secmkt$market == "themi"] <- 12
secmkt$links.months[secmkt$year == "2013" & secmkt$market == "themi"] <- 12
secmkt$links.months[secmkt$year == "2015" & secmkt$market == "themi"] <- 10

secmkt$links.months[secmkt$year == "2009" & secmkt$market == "weruweru"] <- 3 
secmkt$links.months[secmkt$year == "2011" & secmkt$market == "weruweru"] <- 4 
secmkt$links.months[secmkt$year == "2013" & secmkt$market == "weruweru"] <- NA
secmkt$links.months[secmkt$year == "2015" & secmkt$market == "weruweru"] <- NA

secmkt$links.months[secmkt$year == "2009" & secmkt$market == "dosidosi"] <- 9
secmkt$links.months[secmkt$year == "2011" & secmkt$market == "dosidosi"] <- 4
secmkt$links.months[secmkt$year == "2013" & secmkt$market == "dosidosi"] <- 1
secmkt$links.months[secmkt$year == "2015" & secmkt$market == "dosidosi"] <- 6

secmkt$links.months[secmkt$year == "2009" & secmkt$market == "mgagao"] <- 12
secmkt$links.months[secmkt$year == "2011" & secmkt$market == "mgagao"] <- 12
secmkt$links.months[secmkt$year == "2013" & secmkt$market == "mgagao"] <- 11
secmkt$links.months[secmkt$year == "2015" & secmkt$market == "mgagao"] <- 12


secmkt$permits.months
plot(jitter(secmkt$links.months, factor = 0.5), secmkt$cattle.links/secmkt$cattle.permits.est)
plot(jitter(secmkt$links.months, factor = 0.5), jitter(secmkt$permits.months, factor = 0.5))
cor.test(secmkt$links.months, secmkt$permits.months, method = "spearman", use = "pairwise.complete.obs")

secmkt <- droplevels(na.omit(secmkt[secmkt$permits.months >= 8 & secmkt$links.months >= 8, ]))


hist(secmkt$cattle.links/secmkt$cattle.permits.est, nclass = 50)

plot(cattle.permits.est ~ cattle.links, data = secmkt, 
     pch = as.numeric(year), col = as.numeric(market) * 2)#, cex = sqrt(links.months/permits.months) + 0.2)
abline(0, 1)
abline(0, 0.5, lty = 2)
with(secmkt, 
     legend("topleft", paste0(market, " ", year, " (",links.months, "/", permits.months, ")"),
            pch = as.numeric(year), col = as.numeric(market) * 2))
text(x = 20000, y = 20000, expression(permits==links))
text(x = 20000, y = 20000/2, expression(permits==links/2))


# primary markets, 2015

# dosidosi, LINKS total volume: 8510
sum(dat$cattle[dat$fr.match.label == "dosidosi" & dat$year == "2015"])/8510
sum(dat$cattle[dat$to.match.label == "dosidosi" & dat$year == "2015"])/8510
hist(dat$date.dec[dat$fr.match.label == "dosidosi" & dat$year == "2015"], xlim = c(2015, 2016))
# NB permit data for dosidosi covers whole of 2015, LINKS only 2015-06-15 to 2015-11-16



# mgagao, LINKS total volume: 4734
sum(dat$cattle[dat$fr.match.label == "mgagao" & dat$year == "2015"])/4734
sum(dat$cattle[dat$to.match.label == "mgagao" & dat$year == "2015"])/4734
hist(dat$date.dec[dat$fr.match.label == "dosidosi" & dat$year == "2015"], xlim = c(2015, 2016))

# mbulu - no LINKS data for 2015

sort(table(dat$fr.match.label))
sort(table(dat$to.match.label))

round(c(0.01, 0.02, 0.5) * c(70, 5, 15))


# back to the permit data... 

# add coordinates
datw.full$fr.lat <- wards$lat[match(datw.full$fr.ward, wards$name)]
datw.full$fr.long <- wards$long[match(datw.full$fr.ward, wards$name)]
datw.full$to.lat <- wards$lat[match(datw.full$to.ward, wards$name)]
datw.full$to.long <- wards$long[match(datw.full$to.ward, wards$name)]

# calculate distance in units of km
datw.full$dist.km <- 
  round(
    apply(datw.full[, c("fr.long", "to.long", "fr.lat", "to.lat")], 1,
          function(x) gcd.slc(long = x[1:2], lat = x[3:4])), 
    2)

summary(datw.full$dist.km)
sort(unique(datw.full$dist.km))[1:10] # the shortest distances

# are all the zero distances within-ward journeys and vice versa?
table(as.character(datw.full$fr.ward) == as.character(datw.full$to.ward), datw.full$dist.km < 0.00001)
# yes

any(is.na(datw.full$dist.km))

# impute zero distances as the minimum between ward distance
datw.full$dist.km[datw.full$dist.km < 0.0000001] <- NA
datw.full$dist.km[is.na(datw.full$dist.km)] <- min(datw.full$dist.km, na.rm = TRUE)
datw.full$log10.dist.km <- log10(datw.full$dist.km) 

datw.full$dist.km.cat <- 
  cut(datw.full$dist.km, 
      c(-Inf, quantile(datw.full$dist.km, seq(0.1, 0.9, by = 0.1)), Inf))
table(datw.full$dist.km.cat, exclude = NULL)

head(datw.full)
hist(datw.full$dist.km)
hist(datw.full$animals, nclass = 1000, ylim = c(0, 100))

datw.full[sample(nrow(datw.full), 10), ]

# a way to deal with missing data
# fit a zero-inflated model where the response is the number of permits produced by
# an origin in a given month. significant zero-inflation suggests that the zeroes are not
# due to low and fluctuating numbers of outward movements, but some other process 
# such as loss of permit books.

# first, assume that origins that never produced a permit are obligate zeroes and remove them
fr.wards.study <- as.character(sort(unique(datw.full$fr.ward[datw.full$animals > 0])))
origins <- 
  data.frame(
    fr.date = as.character(unique(datw.full$fr.date[datw.full$fr.ward %in% fr.wards.study])),
    stringsAsFactors = FALSE)
origins$year <-
  factor(substr(origins$fr.date, nchar(origins$fr.date) - 6, nchar(origins$fr.date) - 3))
origins$month <-
  factor(substr(origins$fr.date, nchar(origins$fr.date) - 6, nchar(origins$fr.date)))
origins$ward <-
  factor(substr(origins$fr.date, 1, nchar(origins$fr.date) - 8))
dim(origins)
origins$n.permits <- 
  with(droplevels(datw.full[datw.full$fr.ward %in% fr.wards.study, ]), 
       tapply(n.permits, fr.date, sum))[origins$fr.date]
table(origins$n.permits)
origins <- origins[order(origins$ward, origins$month), ]

# find wards that have any months with no permits
fr.wards.study.z <- levels(origins$ward)[!tapply(origins$n.permits != 0, origins$ward, all)]

# fit negative binomial GLMs with and without zero-inflation, and test for zero-inflation
zi.res <-   
  t(sapply(fr.wards.study.z, 
         function(w) {
           # w <- fr.wards.study.z[15]
           # w <- "arusha/meru/kikatiti"
           # w <- "manyara/babati/mwada"
           
           print(w)
           y <- as.vector(droplevels(origins[origins$ward == w, ])$n.permits)
           print(table(y))
           zero.prop <- mean(y == 0)
           if(sum(y != 0) <= 4) 
             return(c(zero.prop = zero.prop, zi.prop = NA, p.val = 1, mean.y = mean(y), 
                      n.nonzero = sum(y != 0), exp.true.zeroes = NA))
           fit.zi <- zeroinfl(y ~ 1, dist = ifelse(all(y %in% 0:1), "poisson", "negbin"))
           zi.prop <- as.vector(plogis(coef(fit.zi)["zero_(Intercept)"]))
           fit.nb <-
             if(all(y %in% 0:1))
               glm(y ~ 1, family = "poisson") else
                 fit.nb <- glm.nb(y ~ 1, init.theta = 1) # control = glm.control(maxit = 25, trace = 3), 
           p.val <- 
             pchisq(2 * (logLik(fit.zi) - logLik(fit.nb)), 
                    1, lower.tail = FALSE)
           c(zero.prop = zero.prop, 
             zi.prop = zi.prop, 
             p.val = p.val, 
             mean.y = mean(y), 
             n.nonzero = sum(y != 0),
             exp.true.zeroes = round((zero.prop - zi.prop) * length(y), 1))
         }))
zi.res <- data.frame(zi.res)

# which wards are significantly zero-inflated (using P < 0.05)?
table(zi.res$p.val < 0.05)
zi.res[zi.res$p.val < 0.05, ]
sort(zi.res$zi.prop[zi.res$p.val < 0.05] / zi.res$zero.prop[zi.res$p.val < 0.05])
range(zi.res$exp.true.zeroes[zi.res$p.val < 0.05])
# for the 16 significantly ZI wards, the proportion of zeroes that are "FALSE"
# zeroes (suggesting that these zeroes are due to, e.g., missing forms rather than fewer 
# movements) the vast majority of zeroes are false (87-100%), so that the number of the 48 months with 
# true zeroes is never >2. 
# proceed by treating all zeroes from these 16 wards as false (so delete these rows).

drop.zero.fr.date <- 
  origins$fr.date[origins$ward %in% rownames(zi.res)[zi.res$p.val < 0.05] & origins$n.permits == 0]

datw <- 
  droplevels(
    datw.full[
      !(datw.full$fr.date %in% drop.zero.fr.date) & 
        (datw.full$fr.ward %in% fr.wards.study), ])
dim(datw.full)
dim(datw)
table(datw$fr.ward)

# OPTIONALLY temporarily subsample data to allow multiple models to be fitted using laptop
#set.seed(1234)
#datw <- droplevels(datw[sample(nrow(datw), sort(ceiling(nrow(datw)/10))), ])
#table(datw$cattle > 0, datw$fr.prodsys)

# centre log10.dist.km and other continuous variables on mean. 
# this needs to be done on the reduced data set that will be used for modelling
# then used to centre the full data set on the mean of the reduced data set

for(cv in cont.vars) {
  datw[, paste0("fr.", cv)] <- wards[match(datw$fr.ward, wards$name), cv]
  datw[, paste0("to.", cv)] <- wards[match(datw$to.ward, wards$name), cv]
  datw.full[, paste0("fr.", cv)] <- wards[match(datw.full$fr.ward, wards$name), cv]
  datw.full[, paste0("to.", cv)] <- wards[match(datw.full$to.ward, wards$name), cv]
  datw[, paste0("log10.fr.", cv, ".cent")] <- scale(log10(datw[, paste0("fr.", cv)]), scale = FALSE)
  datw[, paste0("log10.to.", cv, ".cent")] <- scale(log10(datw[, paste0("to.", cv)]), scale = FALSE)
  datw.full[, paste0("log10.fr.", cv, ".cent")] <- 
    log10(datw.full[, paste0("fr.", cv)]) - attr(datw[, paste0("log10.fr.", cv, ".cent")], "scaled:center")
  datw.full[, paste0("log10.to.", cv, ".cent")] <- 
    log10(datw.full[, paste0("to.", cv)]) - attr(datw[, paste0("log10.to.", cv, ".cent")], "scaled:center")
}; rm(cv)

datw$log10.dist.km.cent <- scale(datw$log10.dist.km, scale = FALSE)
datw.full$log10.dist.km.cent <- datw.full$log10.dist.km - attr(datw$log10.dist.km.cent, "scaled:center")

####### MODEL FITTING ########

# select species to model

(sp <- spp[2])
any.sp <- paste0("any.", sp)

# create a binary indicator of "any movement"
datw[, any.sp] <- as.integer(datw[, sp] > 0.5)

# create data set for zero-truncated count model by dropping pairs of wards with no traffic
datw1 <- datw[datw[, any.sp] == 1, ]
nrow(datw1)
hist(datw1[, sp], nclass = 1000, xlim = c(1, max(datw1$animals)))

datw1[datw1$animals == max(datw1$animals), ]

sum(datw1$cattle)


# select the fixed effects for both models

# list the effects
# these will be used to fit drop1 models, so care must be taken
# with interactions
# the main effect can't be dropped without also dropping the interaction 

# create spline fixed effects
cont.fe <-
  c("log10.dist.km.cent",
    "log10.fr.pop2012.cent", 
    "log10.to.pop2012.cent", 
    paste0("log10.fr.pop.", sp,".cent"), 
    paste0("log10.to.pop.", sp,".cent"), 
    "log10.fr.area.km2.cent", 
    "log10.to.area.km2.cent", 
    "month")


# get knots and boundary knots for both parts of hurdle model, z and c

spline.z <- 
  lapply(datw[, cont.fe], ns, df = 3)
spline.c <- 
  lapply(datw1[, cont.fe], ns, df = 3)
spline.z.attr <-
  lapply(spline.z, function(x) attributes(x)[c("knots", "Boundary.knots")])
spline.c.attr <-
  lapply(spline.c, function(x) attributes(x)[c("knots", "Boundary.knots")])


spline.fe.z <- 
  paste0("ns(", cont.fe, 
         ", knots = c(", 
         apply(round(sapply(spline.z.attr[cont.fe], "[[", 1), 5), 2, paste, collapse = ", "),
         "), Boundary.knots = c(",
         apply(round(sapply(spline.z.attr[cont.fe], "[[", 2), 3), 2, paste, collapse = ", "), 
         "))") 
spline.fe.c <- 
  paste0("ns(", cont.fe, 
         ", knots = c(", 
         apply(round(sapply(spline.c.attr[cont.fe], "[[", 1), 5), 2, paste, collapse = ", "),
         "), Boundary.knots = c(",
         apply(round(sapply(spline.c.attr[cont.fe], "[[", 2), 3), 2, paste, collapse = ", "), 
         "))") 

fac.fe <- 
  c("year", 
    "market.f", 
    "fr.prodsys", "to.prodsys")

fixed.eff.z <- c(spline.fe.z, fac.fe)

fixed.eff.c <- c(spline.fe.c, fac.fe)

fixed.eff.z.fac <- 
  fixed.eff.z[fixed.eff.z %in% names(datw[, sapply(datw, is.factor)])]
fixed.eff.c.fac <- 
  fixed.eff.c[fixed.eff.c %in% names(datw1[, sapply(datw1, is.factor)])]

# make up model formulae for the full model and all "drop1" models
form.rhs.z <-
  sapply((length(fixed.eff.z) + 1):1, 
         function(i) {
           if(i > length(fixed.eff.z)) return(paste(fixed.eff.z[-i], collapse = " + "))
           drop.fe <- fixed.eff.z[i]
           keep.fe <- fixed.eff.z[-grep(gsub("\\)", "\\\\)", gsub("\\(", "\\\\(", drop.fe)), fixed.eff.z)]
           paste(keep.fe, collapse = " + ")
         })
length(form.rhs.z)

form.rhs.c <-
  sapply((length(fixed.eff.c) + 1):1, 
         function(i) {
           if(i > length(fixed.eff.c)) return(paste(fixed.eff.c[-i], collapse = " + "))
           drop.fe <- fixed.eff.c[i]
           keep.fe <- fixed.eff.c[-grep(gsub("\\)", "\\\\)", gsub("\\(", "\\\\(", drop.fe)), fixed.eff.c)]
           paste(keep.fe, collapse = " + ")
         })
length(form.rhs.c)


(form.z <- paste(any.sp, "~ ", form.rhs.z))
(form.c <- paste(sp, "~ ", form.rhs.c))
names(form.rhs.z) <- names(form.z) <- c("full", rev(fixed.eff.z))
names(form.rhs.c) <- names(form.c) <- c("full", rev(fixed.eff.c))


# any factor levels that have no variation in the response? 
invar.lev.list.z <-
  lapply(fixed.eff.z.fac, 
         function(fac) {
           out <- by(datw[, sp], datw[, fac], var)
           invar.lev <- names(out)[out == 0]
           which(datw[, fac] %in% invar.lev)
         })
names(invar.lev.list.z) <- fixed.eff.z.fac
invar.lev.list.z

invar.lev.list.c <-
  lapply(fixed.eff.c.fac, 
         function(fac) {
           out <- by(datw1[, sp], datw1[, fac], var)
           invar.lev <- names(out)[out == 0]
           which(datw1[, fac] %in% invar.lev)
         })
names(invar.lev.list.c) <- fixed.eff.c.fac
invar.lev.list.c

nlevels(datw$period)

# fit the "zero" (binary) part of the hurdle model

# fit full ZI model 
start.time <- Sys.time()
fit.z <-
  glmmTMB(formula(paste(form.z["full"], "+ (1 | fr.ward) + (1 | to.ward) + (1 | period)")),
              family = binomial, data = datw)
print(Sys.time() - start.time)
    

# fit all drop1 ZI models 
dim(datw)
start.time <- Sys.time()
fit.z.stats <-
  t(sapply(form.z[fixed.eff.z], function(form) { 
    print(paste("Fitting model", match(form, form.z), "of", length(form.z)))
    fit <-  
      glmmTMB(formula(paste(form, "+ (1 | fr.ward) + (1 | to.ward) + (1 | period)")),
              family = binomial, data = datw)
    n <- nrow(model.matrix(fit.z))
    fe.var.full <- var(model.matrix(fit.z) %*% fixef(fit.z)$cond) * (n - 1) / n
    fe.var.red <- var(model.matrix(fit) %*% fixef(fit)$cond) * (n - 1) / n
    print(Sys.time() - start.time)
    c(p = anova(fit.z, fit)[2, "Pr(>Chisq)"],
      fe.var.red = fe.var.red, re.var.red = sum(sapply(VarCorr(fit)$cond, sum)),
      eff.size = 1 - fe.var.red/fe.var.full, 
      aic.d1 = AIC(fit) - AIC(fit.z), bic.d1 = BIC(fit) - BIC(fit.z))
  }))


# fit full count model 
start.time <- Sys.time()
fit.c <-
  glmmTMB(formula(paste(form.c["full"], "+ (1 | fr.ward) + (1 | to.ward) + (1 | period)")),
          family = "truncated_nbinom2",
          data = datw1)
print(Sys.time() - start.time)


# fit all drop1 count models
start.time <- Sys.time()
fit.c.stats <-
  t(sapply(form.c[fixed.eff.c], function(form) {  
    print(paste("Fitting model", match(form, form.c), "of", length(form.c)))
    fit <-  
      glmmTMB(formula(paste(form, "+ (1 | fr.ward) + (1 | to.ward) + (1 | period)")),
              family = "truncated_nbinom2",
              data = datw1)
    n <- nrow(model.matrix(fit.c))
    fe.var.full <- var(model.matrix(fit.c) %*% fixef(fit.c)$cond) * (n - 1) / n
    fe.var.red <- var(model.matrix(fit) %*% fixef(fit)$cond) * (n - 1) / n
    print(Sys.time() - start.time)
    c(p = anova(fit.c, fit)[2, "Pr(>Chisq)"],
      fe.var.red = fe.var.red, re.var.red = sum(sapply(VarCorr(fit)$cond, sum)),
      eff.size = 1 - fe.var.red/fe.var.full, 
      aic.d1 = AIC(fit) - AIC(fit.c), bic.d1 = BIC(fit) - BIC(fit.c))
  }))



# make table reporting p-values and summary stats for both zero and count models
fit.z.stats[fixed.eff.z, ]

names(fixed.eff.z) <- sapply(strsplit(gsub("ns\\(", "", fixed.eff.z), ","), "[", 1)
names(fixed.eff.c) <- sapply(strsplit(gsub("ns\\(", "", fixed.eff.c), ","), "[", 1)

res.tab <- 
  data.frame(
    predictor.z = names(fixed.eff.z),
    eff.size.z = paste0(as.vector(round(100 * fit.z.stats[fixed.eff.z, "eff.size"])), "%"),
    p.z = as.vector(round(fit.z.stats[fixed.eff.z, "p"], 3)),
    predictor.c = names(fixed.eff.c),
    eff.size.c = paste0(as.vector(round(100 * fit.c.stats[fixed.eff.c, "eff.size"])), "%"),
    p.c = as.vector(round(fit.c.stats[fixed.eff.c, "p"], 3))
  )
all(res.tab$predictor.z == res.tab$predictor.c)
res.tab[, -4]


table(datw$caprine > 0, datw$cattle > 0)

by(datw$dist.km, datw$caprine > 0, mean)
by(datw$dist.km, datw$cattle > 0, mean)


sum(datw$dist.km * datw$cattle) / sum(datw$cattle)
sum(datw$dist.km * datw$caprine) / sum(datw$caprine)



# check correlations between predictor variables
#cormat <- abs(cor(model.matrix(fit.z)[, -1]))
#cormat <- abs(cor(model.matrix(fit.c)[, -1]))
#diag(cormat) <- 0
#summary(c(cormat))
#hist(cormat[lower.tri(cormat)], nclass = 200)
#which(cormat > 0.8, arr.ind = TRUE)
#which(cormat > 0.5, arr.ind = TRUE)
# ...there are a few over 0.5, but nothing concerning, spline sections and
# productions system levels -- we expect strongish negative correlations among 
# factor levels

# wald test for the fixed effects
eff.names.z <- names(fixef(fit.z)$cond)
p.vals.z <-
  sapply(fixed.eff.z, function(fe) {
    print(fe)
    p.names <- grep(fe, eff.names.z, value = TRUE, fixed = TRUE)
    R <- 
      t(sapply(p.names, 
               function(en) 
                 grep(en, eff.names.z, value = TRUE, fixed = TRUE) == eff.names.z)) + 0
    Wald(fit.z, R)
  })
p.vals.z
round(p.vals.z, 3)


# wald test for the fixed effects
eff.names.c <- names(fixef(fit.c)$cond)
p.vals.c <-
  sapply(fixed.eff.c, function(fe) {
    print(fe)
    p.names <- grep(fe, eff.names.c, value = TRUE, fixed = TRUE)
    R <- 
      t(sapply(p.names, 
               function(en) 
                 grep(en, eff.names.c, value = TRUE, fixed = TRUE) == eff.names.c)) + 0
    Wald(fit.c, R)
  })
p.vals.c
round(p.vals.c, 3)

round(cbind(p.vals.z, p.vals.c), 3)

# how do the Wald and LRT p-values compare?
plot(fit.z.stats[, "p"], p.vals.z)
abline(0, 1)
plot(fit.c.stats[, "p"], p.vals.c)
abline(0, 1)

# where sp == "animals"
# model proportion that are sheep and goats
if(sp == "animals" & FALSE) {
  
  library(lme4)
  datw1$more.caprine <- as.integer(datw1$caprine/datw1$animals > 0.5)
  nrow(datw1)
  length(unique(datw1$fr.to.date))
  fit.p <-
    glmmTMB(cbind(caprine, cattle) ~ year + (1 | fr.ward) + (1 | to.ward), 
            family = binomial,
            data = datw1)
  summary(fit.p)
  
  fit.p <-
    glmer(more.caprine ~ year + (1 | fr.ward) + (1 | to.ward), 
          family = binomial,
          data = datw1)
  summary(fit.p)
  
  fit.p <-
    glmer(paste("more.caprine ~", form.rhs.c["full"],  "+ (1 | fr.ward) + (1 | to.ward)"), 
          family = binomial,
          data = datw1, control = glmerControl(optimizer = "bobyqa"))
  summary(fit.p)
  
}


# save point -- saved workspace here    
# save.image(file = paste0("SEEDZ_pinkslip_gravitymodel_v03DRAFT_", sp, "_", Sys.Date(), ".RData"))

# save and load here

# load("SEEDZ_pinkslip_gravitymodel_v03DRAFT_cattle_2018-09-15.RData")

### load("SEEDZ_pinkslip_gravitymodel_v03DRAFT_cattle_2018-06-14.RData")
### load("SEEDZ_pinkslip_gravitymodel_v03DRAFT_caprine_2018-06-16.RData")

#ranef(fit.z)
summary(fit.z)
ranef.z.fr.ward <- ranef(fit.z)$cond$fr.ward[, "(Intercept)"]
names(ranef.z.fr.ward) <- rownames(ranef(fit.z)$cond$fr.ward)

qqnorm(ranef.z.fr.ward)
qqline(ranef.z.fr.ward)
quiet.origins <- names(sort(ranef.z.fr.ward)[1:3])
wards[wards$fullname %in% quiet.origins, 
      c("fullname", "lat", "long", "area.km2", "status", "popdens2012", "popdensglw.cattle", "prodsys", "market")]

busy.origins <- names(sort(ranef.z.fr.ward)[length(ranef.z.fr.ward) - (0:2)])
wards[wards$fullname %in% busy.origins, 
      c("fullname", "lat", "long", "area.km2", "status", "popdens2012", "popdensglw.cattle", "prodsys", "market")]

sapply(quiet.origins, function(qo) mean(datw[datw$fr.ward == qo, any.sp]))
sapply(busy.origins, function(qo) mean(datw[datw$fr.ward == qo, any.sp]))


ranef.z.to.ward <- ranef(fit.z)$cond$to.ward[, "(Intercept)"]
names(ranef.z.to.ward) <- rownames(ranef(fit.z)$cond$to.ward)

qqnorm(ranef.z.to.ward)
qqline(ranef.z.to.ward)
quiet.destinations <- names(sort(ranef.z.to.ward)[1:3])
wards[wards$fullname %in% quiet.destinations, 
      c("fullname", "area.km2", "status", "popdens2012", "popdensglw.cattle", "prodsys", "market")]
sapply(quiet.destinations, function(qo) mean(datw[datw$to.ward == qo, any.sp]))

busy.destinations <- names(sort(ranef.z.to.ward)[length(ranef.z.to.ward) - (0:2)])
wards[wards$fullname %in% busy.destinations, 
      c("fullname", "long", "area.km2", "status", "popdens2012", "popdensglw.cattle", "prodsys", "market")]

sapply(busy.destinations, function(qo) mean(datw[datw$to.ward == qo, any.sp]))


# calculate an R-squared on the liability scale for the zero model
# (latent-scale marginal R2_GLMM)
n <- nrow(model.matrix(fit.z))
re.var.full <- sum(sapply(VarCorr(fit.z)$cond, sum))
fe.var.full <- var(model.matrix(fit.z) %*% fixef(fit.z)$cond) * (n - 1) / n
fe.var.full / (re.var.full + re.var.full)
# cattle: 40%

# calculate an R-squared on the latent scale for the count model
# (latent-scale marginal R2_GLMM)
n <- nrow(model.matrix(fit.c))
re.var.full <- sum(sapply(VarCorr(fit.c)$cond, sum))
fe.var.full <- var(model.matrix(fit.c) %*% fixef(fit.c)$cond) * (n - 1) / n
fe.var.full / (re.var.full + re.var.full)
# cattle: 24%


datw.full[
  intersect(
    intersect(grep("kikatiti", datw.full$fr.ward), grep("kiusa", datw.full$to.ward)),
    grep("2015-01", datw.full$period)), "cattle"]

datw.full[
  intersect(
    intersect(grep("poli", datw.full$fr.ward), grep("kiusa", datw.full$to.ward)),
    grep("2015-01", datw.full$period)), "cattle"]

datw.full["arusha/arusha/bwawani.arusha/arushaurban/elerai.2009-01", "cattle"]
datw.full["arusha/monduli/meserani.arusha/arushaurban/elerai.2009-01", "cattle"]
datw.full["manyara/simanjiro/naberera.arusha/arushaurban/elerai.2009-01", "cattle"]
datw.full["arusha/arusha/bwawani.arusha/meru/kikatiti.2009-01", "cattle"]
datw.full["arusha/monduli/meserani.arusha/meru/kikatiti.2009-01", "cattle"]
datw.full["manyara/simanjiro/naberera.arusha/meru/kikatiti.2009-01", "cattle"]
datw.full["arusha/arusha/bwawani.arusha/arusha/kisongo.2009-01", "cattle"]
datw.full["arusha/monduli/meserani.arusha/arusha/kisongo.2009-01", "cattle"]
datw.full["manyara/simanjiro/naberera.arusha/arusha/kisongo.2009-01", "cattle"]



head(datw1)

# plot predictions from the models
pred.plot <-
  function(fit, data, newdata, cent, minmax = min, pos = -minmax(-2, -4), offset = 0,
           xcat.col = NULL, xcat.lty = NULL, 
           x.at = sort(c(10^(-6:6), 5 * 10^(-6:6))), jit.fac = 50, rug.lwd = 0.5, x.fn = log10,
           ...) {
    nlev <- sapply(newdata, nlevels)
    xcat <- names(nlev)[nlev > 1.5]
    xvar <- sapply(newdata[, names(nlev)[nlev < 0.5]], var)
    xcont <- names(xvar)[xvar > 0]
    
    Var <- sum(unlist(VarCorr(fit)$cond))
    X <- model.matrix(lm(terms(fit), data = data), data = newdata)
    newdata$eta <- X %*% fixef(fit)$cond
    
    if(family(fit)$link == "logit" & family(fit)$family == "binomial") {
      newdata$pred.resp <- 
        plogis(newdata$eta - 0.5 * Var * tanh(newdata$eta * (1 + 2 * exp(-0.5 * Var))/6))
    }
    if(family(fit)$link == "log" & family(fit)$family == "truncated_nbinom2") {
      lambda <- exp(newdata$eta + 0.5 * Var)
      newdata$pred.resp <- lambda / (1 - dnbinom(0, mu = lambda, size = sigma(fit)))
    }
    ylim <- c(0, max(newdata$pred.resp)) * 1.1
    print(ylim)
    xlim <- range(data[, xcont]) + c(-1, 1) * diff(range(data[, xcont])) / 40
    plot(formula(paste("pred.resp ~", xcont)), data = newdata, type = "n", 
         axes = FALSE, ylim = ylim, xlim = xlim, ...)
    rug(jitter(data[fit$frame[[1]] > 0.5, xcont], factor = jit.fac), 
        col = alpha("black", 0.1), side = 3, lwd = rug.lwd)
    rug(jitter(data[fit$frame[[1]] < 0.5, xcont], factor = jit.fac), 
        col = alpha("black", 0.1), side = 1, lwd = rug.lwd)
    axis(2)
    axis(1, labels = x.at, at = x.fn(x.at) - cent)
    box()
    
    lapply(1:nlevels(newdata[, xcat]), function(i) {
      xcat.lev <- levels(newdata[, xcat])[i]
      points(formula(paste("pred.resp ~", xcont)), 
             data = newdata[newdata[, xcat] == xcat.lev, ], 
             type = "l", 
             lty = ifelse(is.null(xcat.lty), i, xcat.lty[i]), 
             col = ifelse(is.null(xcat.col), "black", xcat.col[i]) )
      text(x = minmax(newdata[, xcont]),
           y = newdata$pred.resp[newdata[, xcont] == minmax(newdata[, xcont]) & 
                                   newdata[, xcat] == xcat.lev], 
           labels = xcat.lev, pos = pos, offset = offset)
    })
  }

# distance & market.f
newdata <-
  expand.grid(
    log10.dist.km.cent = seq(min(datw$log10.dist.km.cent), max(datw$log10.dist.km.cent), length.out = 100), 
    log10.fr.pop2012.cent = 0, 
    log10.to.pop2012.cent = 0, 
    log10.fr.pop.animals.cent = 0,
    log10.to.pop.animals.cent = 0,
    log10.fr.area.km2.cent = 0,
    log10.to.area.km2.cent = 0,
    month = 0,
    year = "2015",
    fr.prodsys = "Agropastoral", # 57% of origin wards are agropastoral
    to.prodsys = "Agropastoral", # 40% of destination wards are agropastoral
    market.f = levels(datw$market.f),
    fr.ward = fr.wards[1], 
    to.ward = to.wards[1],
    period = levels(datw$period)[1]) 
newdata[, any.sp] <- 1
names(newdata) <- gsub("animals", sp, names(newdata))

fr.to.mkt <- do.call("rbind", strsplit(as.character(newdata$market.f), "-"))
xcat.col <- c("black", "black", "black", "red", "red", "red", "blue", "blue", "blue")
xcat.lty <- rep(1:3, 3)
levels(newdata$market.f)

pred.plot(fit = fit.z, 
          data = datw, 
          cent = attr(datw$log10.dist.km.cent, "scaled:center"),
          ylab = "Probability of livestock movement", 
          xlab = expression(Distance ~ (km)),
          rug.lwd = 0.05,
          minmax = min,
          newdata = newdata,
          xcat.col = xcat.col,
          xcat.lty = xcat.lty,
          col = match(fr.to.mkt[, 1], unique(fr.to.mkt[, 1])),
          offset = 1000)
title("Gravity-model-predicted probability of cattle\nmovement, by straight-line distance and market type")
legend("topright", legend = c("", levels(newdata$market.f)), lty = c(0, xcat.lty), col = c("white", xcat.col), 
       bty = "n", cex = 0.8)

rm(newdata)



# distance & fr.prodsys
newdata <-
  expand.grid(
    log10.dist.km.cent = seq(min(datw$log10.dist.km.cent), max(datw$log10.dist.km.cent), length.out = 100), 
    log10.fr.pop2012.cent = 0, 
    log10.to.pop2012.cent = 0, 
    log10.fr.pop.animals.cent = 0,
    log10.to.pop.animals.cent = 0,
    log10.fr.area.km2.cent = 0,
    log10.to.area.km2.cent = 0,
    month = 0,
    year = "2015",
    fr.prodsys = levels(datw$to.prodsys), # 57% of origin wards are agropastoral
    to.prodsys = "Agropastoral", # 40% of destination wards are agropastoral
    market.f = levels(datw$market.f)[4],
    fr.ward = fr.wards[1], 
    to.ward = to.wards[1],
    period = levels(datw$period)[1]) 
newdata[, any.sp] <- 1
names(newdata) <- gsub("animals", sp, names(newdata))

pred.plot(fit = fit.z, 
          data = datw, 
          cent = attr(datw$log10.dist.km.cent, "scaled:center"),
          ylab = "Probability of livestock movement", 
          xlab = expression(Distance ~ (km)),
          minmax = min,
          newdata = newdata) 
title("Gravity-model-predicted probability of cattle\nmovement, by straight-line distance and origin production system")
rm(newdata)

# distance & to.prodsys
newdata <-
  expand.grid(
    log10.dist.km.cent = seq(min(datw$log10.dist.km.cent), max(datw$log10.dist.km.cent), length.out = 100), 
    log10.fr.pop2012.cent = 0, 
    log10.to.pop2012.cent = 0, 
    log10.fr.pop.animals.cent = 0,
    log10.to.pop.animals.cent = 0,
    log10.fr.area.km2.cent = 0,
    log10.to.area.km2.cent = 0,
    month = 0,
    year = "2015",
    fr.prodsys = "Agropastoral", # 57% of origin wards are agropastoral
    to.prodsys = levels(datw$to.prodsys), # 40% of destination wards are agropastoral
    market.f = levels(datw$market.f)[4],
    fr.ward = fr.wards[1], 
    to.ward = to.wards[1],
    period = levels(datw$period)[1]) 
newdata[, any.sp] <- 1
names(newdata) <- gsub("animals", sp, names(newdata))

pred.plot(fit = fit.z, 
          data = datw, 
          cent = attr(datw$log10.dist.km.cent, "scaled:center"),
          ylab = "Probability of livestock movement", 
          xlab = expression(Distance ~ (km)),
          minmax = min,
          newdata = newdata) 
title("Gravity-model-predicted probability of cattle\nmovement, by straight-line distance and destination production system")
rm(newdata)


# destination population density & market type
newdata <-
  expand.grid(
    log10.dist.km.cent = 0,
    log10.fr.pop2012.cent = 0, 
    log10.to.pop2012.cent = 
      seq(min(datw$log10.to.pop2012.cent), max(datw$log10.to.pop2012.cent), length.out = 100), 
    log10.fr.pop.animals.cent = 0,
    log10.to.pop.animals.cent = 0,
    log10.fr.area.km2.cent = 0,
    log10.to.area.km2.cent = 0,
    month = 0,
    year = "2015",
    secondary.f = levels(datw$secondary.f),
#    fr.status = "Rural Ward", 
#    to.status = "Urban Ward",
    fr.ward = fr.wards[1], 
    to.ward = to.wards[1])
newdata[, any.sp] <- 1
names(newdata) <- gsub("animals", sp, names(newdata))

pred.plot(fit = fit.z, 
          data = datw, 
          cent = attr(datw$log10.to.pop2012.cent, "scaled:center"),
          ylab = "Probability of livestock movement", 
          xlab = expression(Ward ~ population ~ density ~ (km^-2)),
          minmax = max,
          x.at = c(500, 1000, 2000, 5000, 10000, 20000, 50000, 100000),
          newdata = newdata) 
title("Gravity-model-predicted probability of livestock movement,\nby destination population density and origin market type")
rm(newdata)

# month and year

newdata <-
  expand.grid(
    log10.dist.km.cent = -0.391,   # -0.391 equates to 50 km
    log10.fr.pop2012.cent = 0, 
    log10.to.pop2012.cent = 0,
    log10.fr.pop.animals.cent = 0,
    log10.to.pop.animals.cent = 0,
    log10.fr.area.km2.cent = 0,
    log10.to.area.km2.cent = 0,
    market.f = levels(datw$market.f)[5],
    month = 0:11,
    year = levels(datw$year),
    fr.prodsys = "Agropastoral", # 57% of origin wards are agropastoral
    to.prodsys = "Agropastoral", # 40% of destination wards are agropastoral
    fr.ward = fr.wards[1], 
    to.ward = to.wards[1],
    period = levels(datw$period)[1]) 
newdata[, any.sp] <- 1
names(newdata) <- gsub("animals", sp, names(newdata))

pred.plot(fit = fit.z, 
          data = datw, 
          cent = 0,
          ylab = "Probability of livestock movement", 
          xlab = "Month",
          minmax = max,
          x.at = 0:11,
          x.fn = I,
          jit.fac = 2,
          rug.lwd = 0.05,
          newdata = newdata) 
title("Gravity-model-predicted probability of livestock movement,\nby calendar month and year")

rm(newdata)



# plot predictions from count model

# fr.prodsys and distance
newdata <-
  expand.grid(
    log10.dist.km.cent = seq(min(datw1$log10.dist.km.cent), max(datw1$log10.dist.km.cent), length.out = 100), 
    log10.fr.pop2012.cent = 0, 
    log10.to.pop2012.cent = 0, 
    log10.fr.pop.animals.cent = 0,
    log10.to.pop.animals.cent = 0,
    log10.fr.area.km2.cent = 0,
    log10.to.area.km2.cent = 0,
    month = 0,
    year = "2015",
    fr.prodsys = "Agropastoral", # 57% of origin wards are agropastoral
    to.prodsys = levels(datw$to.prodsys), # 40% of destination wards are agropastoral
    market.f = levels(datw$market.f)[4],
    fr.ward = fr.wards[1], 
    to.ward = to.wards[1],
    period = levels(datw$period)[1]) 
newdata[, sp] <- 1
names(newdata) <- gsub("animals", sp, names(newdata))

pred.plot(fit = fit.c, 
          data = datw1, 
          cent = attr(datw$log10.dist.km.cent, "scaled:center"),
          ylab = "Number of livestock moved", 
          xlab = expression(Distance ~ (km)),
          minmax = max,
          pos = 2,
          newdata = newdata) 
title("Gravity-model-predicted number of livestock moved,\nby straight-line distance and destination production system")
rm(newdata)


# to.prodsys and distance
newdata <-
  expand.grid(
    log10.dist.km.cent = seq(min(datw1$log10.dist.km.cent), max(datw1$log10.dist.km.cent), length.out = 100), 
    log10.fr.pop2012.cent = 0, 
    log10.to.pop2012.cent = 0, 
    log10.fr.pop.animals.cent = 0,
    log10.to.pop.animals.cent = 0,
    log10.fr.area.km2.cent = 0,
    log10.to.area.km2.cent = 0,
    month = 0,
    year = "2015",
    fr.prodsys = levels(datw$to.prodsys), # 57% of origin wards are agropastoral
    to.prodsys = "Agropastoral", # 40% of destination wards are agropastoral
    market.f = levels(datw$market.f)[4],
    fr.ward = fr.wards[1], 
    to.ward = to.wards[1],
    period = levels(datw$period)[1]) 
newdata[, sp] <- 1
names(newdata) <- gsub("animals", sp, names(newdata))

pred.plot(fit = fit.c, 
          data = datw1, 
          cent = attr(datw$log10.dist.km.cent, "scaled:center"),
          ylab = "Number of livestock moved", 
          xlab = expression(Distance ~ (km)),
          minmax = max,
          pos = 2,
          newdata = newdata) 
title("Gravity-model-predicted number of livestock moved,\nby straight-line distance and origin production system")
rm(newdata)


# fr.status & distance
newdata <-
  expand.grid(
    log10.dist.km.cent = seq(min(datw1$log10.dist.km.cent), max(datw1$log10.dist.km.cent), length.out = 100), 
    log10.fr.popdens2012.cent = 0, 
    log10.to.popdens2012.cent = 0, 
    log10.fr.pop.cattle.cent = 0,
    log10.to.pop.cattle.cent = 0,
    log10.fr.area.km2.cent = 0,
    log10.to.area.km2.cent = 0,
    month = 0,
    year = "2015",
    secondary.f = levels(datw1$secondary.f)[1],
    fr.status = levels(datw1$fr.status), 
    to.status = "Urban Ward",
    period = rev(levels(datw1$period))[1], 
    fr.ward = fr.wards[1], 
    to.ward = to.wards[1]) 
newdata[, sp] <- 1
names(newdata) <- gsub("animals", sp, names(newdata))

pred.plot(fit = fit.c, 
          data = datw1, 
          cent = attr(datw$log10.fr.area.km2.cent, "scaled:center"),
          ylab = "Number of livestock moved", 
          xlab = expression(Distance ~ (km)),
          minmax = max,
          pos = 2,
          newdata = newdata) 

title("Gravity-model-predicted number of livestock moved,\nby straight-line distance and type of origin")
rm(newdata)

# to.status & distance
newdata <-
  expand.grid(
    log10.dist.km.cent = seq(min(datw1$log10.dist.km.cent), max(datw1$log10.dist.km.cent), length.out = 100), 
    log10.fr.popdens2012.cent = 0, 
    log10.to.popdens2012.cent = 0, 
    log10.fr.pop.cattle.cent = 0,
    log10.to.pop.cattle.cent = 0,
    log10.fr.area.km2.cent = 0,
    log10.to.area.km2.cent = 0,
    year = "2015",
    month = 0,
    secondary.f = levels(datw$secondary.f)[1],
    fr.status = "Rural Ward",
    to.status = levels(datw1$to.status), 
    period = rev(levels(datw1$period))[1], 
    fr.ward = fr.wards[1], 
    to.ward = to.wards[1]) 
newdata[, sp] <- 1
names(newdata) <- gsub("animals", sp, names(newdata))

pred.plot(fit = fit.c, 
          data = datw1, 
          cent = attr(datw$log10.fr.area.km2.cent, "scaled:center"),
          ylab = "Number of livestock moved", 
          xlab = expression(Distance ~ (km)),
          minmax = max,
          pos = 2,
          newdata = newdata) 

title("Gravity-model-predicted number of livestock moved,\nby straight-line distance and type of destination")
rm(newdata)


# month and year

newdata <-
  expand.grid(
    log10.dist.km.cent = -0.391,   # -0.391 equates to 50 km
    log10.fr.pop2012.cent = 0, 
    log10.to.pop2012.cent = 0,
    log10.fr.pop.animals.cent = 0,
    log10.to.pop.animals.cent = 0,
    log10.fr.area.km2.cent = 0,
    log10.to.area.km2.cent = 0,
    market.f = levels(datw1$market.f)[5],
    month = 0:11,
    year = levels(datw1$year),
    fr.prodsys = "Agropastoral", # 57% of origin wards are agropastoral
    to.prodsys = "Agropastoral", # 40% of destination wards are agropastoral
    fr.ward = fr.wards[1], 
    to.ward = to.wards[1],
    period = levels(datw1$period)[1]) 
newdata[, sp] <- 1
names(newdata) <- gsub("animals", sp, names(newdata))

pred.plot(fit = fit.c, 
          data = datw1, 
          cent = 0,
          ylab = "Number of livestock moved", 
          xlab = "Month",
          minmax = max,
          x.at = 0:11,
          x.fn = I,
          pos = 1,
          jit.fac = 1,
          newdata = newdata) 
title("Gravity-model-predicted number of livestock moved,\nby calendar month and year")
rm(newdata)


#boxplot(popdens2012 ~ factor(status, levels(datw$fr.status)), data = wards, log = "y")
#boxplot(pop2012 ~ factor(status, levels(datw$fr.status)), data = wards, log = "y")
#boxplot(area.km2 ~ factor(status, levels(datw$fr.status)), data = wards, log = "y")
#boxplot(fr.popdens2012 ~ fr.status, data = datw1, log = "y")
#boxplot(to.popdens2012 ~ to.status, data = datw1, log = "y")
#boxplot(fr.pop2012 ~ fr.status, data = datw1, log = "y")
#boxplot(to.pop2012 ~ to.status, data = datw1, log = "y")
cor(datw1$fr.popdens2012, as.numeric(datw1$fr.status), method = "spearman")
cor(datw1$to.popdens2012, as.numeric(datw1$to.status), method = "spearman")
cor(wards$popdens2012, as.numeric(factor(wards$status, levels(datw$fr.status))), method = "spearman")
cor(wards$pop2012, as.numeric(factor(wards$status, levels(datw$fr.status))), method = "spearman")
cor(wards$area.km2, as.numeric(factor(wards$status, levels(datw$fr.status))), method = "spearman")
cor(datw1$fr.pop2012, as.numeric(datw1$fr.status), method = "spearman")
cor(datw1$to.pop2012, as.numeric(datw1$to.status), method = "spearman")
#plot(animals ~ to.popdens2012, data = datw1, log = "xy")
#plot(animals ~ factor(round(log10.to.popdens2012.cent)), data = datw1, log = "y")
table(factor(round(datw1$log10.to.popdens2012.cent)))


###############
# simulations #
###############

# random effects and dispersion parameter of negative binomial
Sigma2.z <- VarCorr(fit.z)$cond
Sigma2.c <- VarCorr(fit.c)$cond
theta <- sigma(fit.c)

# date-hour for output files
time.stamp <- gsub(" ", "-", substr(Sys.time(), 1, 13))

#set.seed(round(as.numeric(Sys.time()) * 1000 - 1528300000000))

# output movement network parameters from hurdle model for 12 months in 2015
sim.monthly2015 <-
  lapply(sort(unique(datw$month)), function(m) {
    # m <- 0
    
    # make subset of full data set
    datw.p <- (datw.full[datw.full$month == m & datw.full$year == "2015", ])
    dim(datw.p)
    
    # condition on the origin and destination wards and the period
    ranef.cond <- c("fr.ward", "to.ward", "period")
    
    # define response variables (any animal and number of animals given >= 1 animal)
    any.sp.sim <- paste0(any.sp, ".sim")
    sp.sim <- paste0(sp, ".sim")
    
    sum(datw.p[, sp])
    datw.p[, any.sp] <- as.integer(datw.p[, sp] > 0.5)
    table(datw.p[, any.sp])
    
    # binomial part
    # fixed effects part
    MM.z <- model.matrix(formula(form.z), data = datw.p)
    rownames(MM.z) <- rownames(datw.p)
    B.z <- fixef(fit.z)$cond         #[-grep("period", names(fixef(fit.z)$cond))]
    eta.f <- c(MM.z %*% B.z)
    names(eta.f) <- rownames(MM.z)
    sum(plogis(eta.f))

    # random effects part
    re.tab <-
      sapply(names(Sigma2.z), 
             function(re) {
               if(!re %in% ranef.cond) {
                 out <- rnorm(nlevels(datw.p[, re]), sd = sqrt(Sigma2.z[[re]]))
                 names(out) <- levels(datw.p[, re])
               } else {
                 out <- ranef(fit.z)$cond[[re]][, "(Intercept)"]
                 names(out) <- levels(fit.z$frame[, re])
               } 
               out[as.character(datw.p[, re])]
             })
    
    rownames(re.tab) <- rownames(datw.p)
    dim(re.tab)
    colSums(is.na(re.tab))
    head(re.tab)
    re.tab[is.na(re.tab[, "fr.ward"]), "fr.ward"] <- -Inf  # because plogis(-Inf) returns 0
    
    eta <- rowSums(cbind(eta.f, re.tab))
    summary(eta)
    
    #eta.adj <- sum(jensen.adjust(plogis(eta), sum(unlist(Sigma2.z[-match(ranef.cond, names(Sigma2.z))])))    )
#    sum(plogis(eta[rownames(datw)[rownames(datw) %in% names(eta)]]))
#    sum(rbinom(length(eta[rownames(datw)[rownames(datw) %in% names(eta)]]), 1, 
#               plogis(eta[rownames(datw)[rownames(datw) %in% names(eta)]])))
#    sum(rbinom(length(eta.f[rownames(datw)[rownames(datw) %in% names(eta.f)]]), 1, 
#               plogis(eta.f[rownames(datw)[rownames(datw) %in% names(eta.f)]])))
#    sum(datw.p[, any.sp])

    datw.p$pred.prob.mvt <- plogis(eta)
    
    rm(eta, eta.f, re.tab)
    
    # negative binomial part
    # simulate animal numbers moved for all rows
    MM.c <- model.matrix(formula(form.c), data = datw.p)#[1:30,]
    eta.f <- c(MM.c %*% fixef(fit.c)$cond)
    
    # random part of model
    re.tab <-
      sapply(names(Sigma2.c), 
             function(re) {
               if(!re %in% ranef.cond) {
                 out <- rnorm(nlevels(datw.p[, re]), sd = sqrt(Sigma2.c[[re]]))
                 names(out) <- levels(datw.p[, re])
               } else {
                 out <- ranef(fit.c)$cond[[re]][, "(Intercept)"]
                 names(out) <- rownames(ranef(fit.c)$cond[[re]])
               }
               out <- out[as.character(datw.p[, re])]
               if(any(is.na(out))) {
                 names(out) <- as.character(datw.p[, re])
                 imputed.re.names <- unique(names(out[is.na(out)]))
                 imputed.re <- rnorm(length(imputed.re.names), sd = sqrt(Sigma2.c[[re]]))
                 names(imputed.re) <- imputed.re.names
                 out[names(out) %in% imputed.re.names] <- 
                   imputed.re[names(out[names(out) %in% imputed.re.names])]
               }
               out
             })
    
    rownames(re.tab) <- NULL
    dim(re.tab)
    colSums(is.na(re.tab))
    head(re.tab)
    re.tab[rowSums(is.na(re.tab)) > 0, ]
    
    eta <- rowSums(cbind(eta.f, re.tab))
    
    datw.p$pred.n.mvt <- exp(eta)

    # scale up the probilities and rates to account for the
    length(unique(dat$batch)) / 20
    # subsample
    (infl.fac <- 20 / length(unique(dat$batch)))
    # how would increasing the number of permits affect
    # (a) the probability of seeing any permits, p?
    # (b) the number of permits found, given any permits are found?
    # at the extremes, all the effect will be in either
    # (a) infl.fac times as many directed edges will have permits (assuming p is very small)
    # (b) the same number of edges will have permits but the number of animals will be times infl.fac.
    
    # if probability of detecting movement with in t time is p, the probability of 
    # detecting it in k*t time is:
    # probability of not detecting mvt in t is 1 - p
    # probability of not detecting mvt in k*t is (1 - p)^k
    # probability of detecting mvt in k*t time is 1 - (1 - p)^k
    
    # start by inflating the probability of finding any permit
    datw.p$pred.prob.mvt.adj <- 1 - (1 - datw.p$pred.prob.mvt)^infl.fac
    
    # how much has the number of permits increased by?
    mean(datw.p$pred.prob.mvt.adj) / mean(datw.p$pred.prob.mvt)
    
    # make up the difference
    infl.fac.c <- infl.fac / (mean(datw.p$pred.prob.mvt.adj) / mean(datw.p$pred.prob.mvt))
    datw.p$pred.n.mvt.adj <- datw.p$pred.n.mvt * infl.fac.c

    # overall inflation:
    mean(datw.p$pred.prob.mvt.adj) / mean(datw.p$pred.prob.mvt) * infl.fac.c
    
    
    # simulate any movement from binomial
    datw.p[, any.sp.sim] <- rbinom(nrow(datw.p), 1, datw.p$pred.prob.mvt.adj)
    table(datw.p[, any.sp.sim], exclude = NULL)
    mean(datw.p[, any.sp.sim])
    mean(datw.p[, any.sp])
    mean(datw.p[rownames(datw)[rownames(datw) %in% rownames(datw.p)], any.sp.sim])
    mean(datw.p[rownames(datw)[rownames(datw) %in% rownames(datw.p)], any.sp])
    
    table(datw.p[, any.sp.sim], datw.p[, any.sp])
    cor(datw.p[, any.sp.sim], datw.p[, any.sp])
    
    
    # simulate from zero-truncated negative binomial
    datw.p$temp.count <-
      rnbinom(n = nrow(datw.p), size = sigma(fit.c), mu = datw.p$pred.n.mvt.adj)
    
    while(any(datw.p$temp.count == 0)) {
      zero.counts <- datw.p$temp.count == 0
      n.zero.counts <- sum(zero.counts)
      print(n.zero.counts)
      datw.p$temp.count[zero.counts] <- 
        rnbinom(n = n.zero.counts, size = sigma(fit.c), mu = datw.p$pred.n.mvt[zero.counts])
    }
    
    datw.p[, sp.sim] <- datw.p[, any.sp.sim] * datw.p$temp.count

    cor(datw.p[, sp], datw.p[, sp.sim], method = "spearman")
    #abline(0, 1)
    mean(datw.p[, sp.sim])
    #boxplot(list(sim = datw$animals.sim, real = datw$animals))
    
    
    if(FALSE) {
      out.degree <- tapply(datw.p[, any.sp], datw.p$fr.ward, sum)
      in.degree <- tapply(datw.p[, any.sp], datw.p$to.ward, sum)
      out.degree.sim <- tapply(datw.p[, any.sp.sim], datw.p$fr.ward, sum)
      in.degree.sim <- tapply(datw.p[, any.sp.sim], datw.p$to.ward, sum)
      
      par(mfrow = c(3, 2))
      xlim.out <- c(0, max(out.degree,out.degree.sim))
      xlim.in <- c(0, max(in.degree, in.degree.sim))
      hist(out.degree, nclass = 50, xlim = xlim.out)
      hist(in.degree, nclass = 50, xlim = xlim.in)
      hist(out.degree.sim, nclass = 50, xlim = xlim.out)
      hist(in.degree.sim, nclass = 50, xlim = xlim.in)
      cor(out.degree.sim, out.degree, method = "spearman")
      plot(out.degree, out.degree.sim)
      cor(in.degree.sim, in.degree, method = "spearman")
      plot(in.degree, in.degree.sim)
      
      summary(in.degree)
      summary(in.degree.sim)
      summary(out.degree)
      summary(out.degree.sim)
      
    }
    

    # create empty matrix of wards
    nlevels(datw.p$fr.ward)
    nlevels(datw.p$to.ward)
    length(unique(c(levels(datw.p$fr.ward), levels(datw.p$to.ward))))
    
    mu1.z <-
      structure(
        rep(NA, nlevels(datw.p$to.ward)^2), 
        dim = c(nlevels(datw.p$to.ward), nlevels(datw.p$to.ward)),
        dimnames = list(outmove = levels(datw.p$to.ward), inmove = levels(datw.p$to.ward)))
    mu1.z[1:10, 1:10]
    dim(mu1.z)
    sim.moves <- mu1.c <- mu1.z
    
    
    # fill probabilities, movement rates and simulated moves
    # (takes < 30 s)
    for(i in 1:nrow(datw.p)) {
      mu1.z[as.character(datw.p$fr.ward[i]), as.character(datw.p$to.ward[i])] <- 
        datw.p$pred.prob.mvt.adj[i]
      mu1.c[as.character(datw.p$fr.ward[i]), as.character(datw.p$to.ward[i])] <-
        datw.p$pred.n.mvt.adj[i]
#      sim.moves[as.character(datw.p$fr.ward[i]), as.character(datw.p$to.ward[i])] <-
#        datw.p[i, sp.sim]
    }; rm(i)
    mu1.z[is.na(mu1.z)] <- 0
    mu1.c[is.na(mu1.c)] <- 0
    #sim.moves[is.na(sim.moves)] <- 0

    # write the parameter matrices to disc
    mtxt <- formatC(m + 1, flag = "0", width = 2)
    write.csv(round(mu1.z, 6), paste0("networks/simmoves.monthly.2015/", sp, "/", sp, ".month", mtxt, ".mvt.matrix.prob.", time.stamp,".csv"))
    write.csv(round(mu1.c, 4), paste0("networks/simmoves.monthly.2015/", sp, "/", sp, ".month", mtxt, ".mvt.matrix.rate.", time.stamp,".csv"))
    write.csv(cbind(theta = round(theta, 4)), 
              paste0("networks/simmoves.monthly.2015/", sp, "/", sp, ".month", mtxt, ".mvt.matrix.rate.theta.", time.stamp,".csv"),
              row.names = FALSE)
    # write the simulated movements to disc
    #write.csv(sim.moves, paste0("networks/simmoves.monthly.2015/", sp, "/", sp, ".month", mtxt, ".mvt.matrix.simmoves.", time.stamp,".csv"))
    
    list(
      sp = sp, year = "2015", calendar.month = mtxt, time.stamp = time.stamp,
      mu.z = mu1.z, mu.c = mu1.c, theta = theta) #, sim.moves = sim.moves)
    
  })



# calculate expected number of animals for each row
sp.pred <- paste0(sp, ".pred")
datw$pred.z <- predict(fit.z)
datw$pred.c <- predict(fit.c, newdata = datw)
datw[, sp.pred] <- datw$pred.z * datw$pred.c

boxplot(formula(paste(sp.pred, "~ factor(fr.ward)")) , 
        data = datw[datw$fr.secondary == 1 & (datw[, any.sp] == 1), ], log = "y",
        ylim = c(min(datw[datw[, any.sp] == 1, sp.pred]), max(datw[, sp])))
boxplot(formula(paste(sp, "~ factor(fr.ward)")) , 
        data = datw[datw$fr.secondary == 1 & (datw[, any.sp] == 1), ], log = "y",
        ylim = c(min(datw[datw[, any.sp] == 1, sp.pred]), max(datw[, sp])))
cor(datw[, sp.pred], datw[, sp], method = "spearman")
sum(datw[, sp.pred])
sum(datw[, sp])


tapply(datw[datw$fr.secondary == 1, sp.pred], factor(datw$fr.ward[datw$fr.secondary == 1]), sum)
tapply(datw[datw$fr.secondary == 1, sp], factor(datw$fr.ward[datw$fr.secondary == 1]), sum)


# draw a couple of maps of the simulated data
map.fr <- 
  get_googlemap(
    center = 
      apply(datw[datw[, any.sp] == 1, c("fr.long", "fr.lat")], 2, 
            function(x) mean(range(x, na.rm = TRUE))), 
    zoom = 7, size = c(350, 410),
    maptype = "terrain", color = "bw", scale = 2)


sim.or.real <- sample(c(any.sp, rep(any.sp.sim, 1)))
lapply(1:length(sim.or.real), function(i) {
  movemapzoom <-
    ggmap(map.fr) +
    geom_segment(data = datw[datw[, sim.or.real[i]] == 1, ], 
                 mapping = aes(x = jitter(fr.long, factor = 5), 
                               y = jitter(fr.lat, factor = 5), 
                               xend = jitter(to.long, factor = 5), 
                               yend = jitter(to.lat, factor = 5),
                               colour = period),
                 arrow = arrow(length = unit(0.01, "npc")),
                 size = 0.1) 
  
  ggsave(paste0("network.map.", sim.or.real[i], i, ".jpg"), movemapzoom)
})

table(datw[, any.sp])
table(datw[, any.sp.sim])

# fit models to simulated data, compare estimates and residuals to real data


fit.z.sim <-
  update(fit.z, paste0(any.sp.sim, "~ ."), data = datw)
summary(fit.z.sim)

plot(fixef(fit.z)$cond, fixef(fit.z.sim)$cond, 
     cex = 2 * (1 - coef(summary(fit.z.sim))$cond[, "Pr(>|z|)"]))
abline(0, 1)

# residuals plot, slow (< 1 min) 
plot(fitted(fit.z), resid(fit.z))
points(fitted(fit.z.sim), resid(fit.z.sim), col = "red", pch = 3, cex = 0.6)



fit.c.sim <-
  update(fit.c, paste0(sp.sim, "~ ."), data = datw[datw$any.animals.sim == 1, ])
summary(fit.c.sim)

plot(fixef(fit.c)$cond, fixef(fit.c.sim)$cond, 
     cex = 2 * (1 - coef(summary(fit.c.sim))$cond[, "Pr(>|z|)"]))
abline(0, 1)

# residuals plot
plot(fitted(fit.c), resid(fit.c))
points(fitted(fit.c.sim), resid(fit.c.sim), col = "red", pch = 3, cex = 0.6)



######################################################################
# plots, maps, other non-essential code below this point             #
######################################################################



if(F) {
  fit.hurdle <- hurdle(form.c, data = datw, dist = "negbin")
  AIC(fit.hurdle)
  summary(fit.hurdle)
  
  
  # compute an Rsq for the count part
  beta.c <- coef(fit.hurdle, model = "count")
  X.c <- model.matrix(fit.hurdle, model = "count")
  S2f.c <- var(X.c %*% beta.c)
  lambda <- mean(X.c %*% beta.c)
  theta <- fit.hurdle$theta
  S2e.c <- trigamma(1/(1/lambda + 1/theta))
  (Rsq.c <- S2f.c / (S2e.c + S2f.c))
  
  # compute an Rsq for the zero part
  beta.z <- coef(fit.hurdle, model = "zero")
  X.z <- model.matrix(fit.hurdle, model = "zero")
  S2f.z <- var(X.z %*% beta.z)
  p <- plogis(mean(X.c %*% beta.c))
  # compare with
  mean(datw$any.animals)
  S2e.z <- 1/(p*(1-p))
  (Rsq.z <- S2f.z / (S2e.z + S2f.z))
  
  
  mean(predict(fit.hurdle, type = "prob"))
  
  # try fitting the binomial part in MCMCglmm
  
  
  prior <- 
    list(
      R = list(V = 1, fix = 1),
      G = list(
        G1 = list(V = 1, nu = 0.002),
        G2 = list(V = 1, nu = 0.002)))
  
  library(MCMCglmm)
  fit.mcmc.z <- 
    MCMCglmm(form.z,
             random = ~ fr.ward + to.ward,
             family = "categorical", prior = prior,
             data = datw,
             burnin = 30000, nitt = 130000, thin = 100)
  plot(fit.mcmc.z)
  est <- apply(fit.mcmc$Sol, 2, median)
  ci95 <- HPDinterval(fit.mcmc$Sol)
  data.frame(ss = i, b.name = names(unlist(B))[-2], b = unlist(B)[-2], est = est, ci95)
  
  prior <- 
    list(
      R = list(V = 1, nu = 0.002))
  fit.mcmc <- 
    MCMCglmm(
      animals ~
        log10.dist.km + factor(period),
      #    random = ~ fr.ward + to.ward,
      family = "poisson", prior = prior,
      data = datw,
      burnin = 30000, nitt = 130000, thin = 100)
  plot(fit.mcmc)
  est <- apply(fit.mcmc$Sol, 2, median)
  ci95 <- HPDinterval(fit.mcmc$Sol)
  data.frame(ss = i, b.name = names(unlist(B))[-2], b = unlist(B)[-2], est = est, ci95)
  
}

# slow, so skip for now, but might be useful for partial R^2:
if(FALSE) {
}



if(FALSE) {
  AIC(fit.glmmadmb.c)
  
  start.time <- Sys.time()
  fit.glmmadmb.c.ixn <-
    glmmadmb(animals ~ 
               ns(log10.dist.km.cent, 3) + 
               ns(log10.fr.popdens2012.cent, 3) * fr.status + 
               ns(log10.to.popdens2012.cent, 3) * to.status + 
               period + 
               (1 | fr.ward) + (1 | to.ward),
             family = "truncnbinom", data = datw1)
  Sys.time() - start.time
  
  AIC(fit.glmmadmb.c.ixn)
  AIC(fit.glmmadmb.c.ixn)
  
  fit.glmmadmb.c <- fit.glmmadmb.c.ixn
  
}

summary(fit.glmmTMB.c)

plot(animals ~ to.status, data = datw1, log = "y")
plot(animals ~ period, data = datw1, log = "y")
table(datw1$to.status)

if(FALSE) {
  # check results against ZI-poisson in MCMCglmm
  prior <- 
    list(
      R = list(V = 1, nu = 0.002),
      G = list(
        G1 = list(V = 1, nu = 0.002),
        G2 = list(V = 1, nu = 0.002)))
  library(MCMCglmm)
  fit.mcmc <- 
    MCMCglmm(formula(form.c["full"]),
             random = ~ fr.ward + to.ward,
             family = "ztpoisson", prior = prior,
             data = datw1,
             burnin = 3000, nitt = 13000, thin = 10)
  plot(fit.mcmc)
  est <- apply(fit.mcmc$Sol, 2, median)
  plot(est, coef(fit.glmmadmb.c)) # pretty close
  abline(0, 1)
  #ci95 <- HPDinterval(fit.mcmc$Sol)
  #data.frame(ss = i, b.name = names(unlist(B))[-2], b = unlist(B)[-2], est = est, ci95)
  
  var(as.matrix(fit.mcmc$X) %*% est)
  var(as.matrix(fit.mcmc$X)[ , -11:-12] %*% est[-11:-12])
  var(as.matrix(fit.mcmc$X)[ , c(-5:-7, -11:-12)] %*% est[c(-5:-7, -11:-12)])
}

# slow, so skip for now, but might be useful for partial R^2
if(FALSE) {
  # fit full model and all drop1 models
  start.time <- Sys.time()
  fit.glmmadmb.c.list <-
    lapply(form.c, function(form) {
      print(form)
      glmmadmb(formula(paste(form, "+ (1 | fr.ward) + (1 | to.ward)")),
               family = "truncnbinom", data = datw1)
    })#, mc.cores = detectCores())
  Sys.time() - start.time
  
  summary(fit.glmmadmb.c.list$full)
  
  # compare models to full model
  mod.comp.c <-
    sapply(fit.glmmadmb.c.list[fixed.eff], function(fit) {
      n <- fit.glmmadmb.c.list$full$n
      mss.full <- 
        var(model.matrix(fit.glmmadmb.c.list$full) %*% fixef(fit.glmmadmb.c.list$full)) * (n - 1)
      mss.red <- var(model.matrix(fit) %*% fixef(fit)) * (n - 1)
      c(p = anova(fit, fit.glmmadmb.c.list$full)[2, "Pr(>Chi)"],
        eff.size = 1 - mss.red/mss.full)
    })
}





plot(density(datw$dist.km[datw$caprine > 0], adjust = adjust), lty = 2, xlim = c(0, 350),
     main = "Distribution of movement distances", xlab = "kilometres")
lines(density(datw$dist.km[datw$cattle > 0], adjust = adjust), lty = 1)
legend("topright", legend = c("Cattle", "Sheep & goats"), lty = 1:2)

