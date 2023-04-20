
# parameters ####
# patch 1 par s culex rural
kv1 <- 80000 # carrying capacity
Reggs1 <- 0.8 # rate eggs are laid 
Heggs1 <- 0.4 # hatching rate eggs
meggs1 <- 0.00001 # death rate eggs 
Peggs1 <- 0.94 # probability of being viable eggs
mpup1 <- 0.135 # death rate pupae
pupmaturation1 <- 0.3 # pupal maturation rate (Jen - 3-7 days)
pupsurvival1 <- 0.3 # pupal probability of surviving to maturation
mv1 <- 0.063 # mortality adult vector (Jen - 7-30 days)
lv1 <- 0.83 # latent period adult vector 
probVinfected1 <- 0.78 # probability of vector becoming infected from a bite
biterate1 <- 0.06 # rate that vectors bite host (Jen - 3-6 days)

y1  <- 0.2 # host clearance
lh1 <- 0.5 # host latent period remove latent period in hosts? 
mh1 <- 0.00043 # host mortality 
xh1 <- 0.0037 # host reproduction
kh1 <- 2000 # host carrying capacity 
probHinfected1 <- 0.1 # prob infected from a bite (0-1)
ah1 <- 0.00125 # additional parasite induced mortality

# patch 2 pars aedes peri
kv2 <- 80000
Reggs2 <- 0.6
Heggs2 <- 0.23
Peggs2 <- 0.84 #https://www.scielo.br/j/bjb/a/QvsZdPBrqsYRZVGxCdM7jZB/?format=pdf&lang=en
meggs2 <- 0.001###
mpup2 <- 0.12
pupmaturation2 <- 0.3
pupsurvival2 <- 0.6
mv2 <- 0.07
lv2 <- 0.83
probVinfected2 <- 0.78
biterate2 <- 0.06

y2  <- 0.2
lh2 <- 0.5
mh2 <- 0.00043
xh2 <- 0.0037
kh2 <- 800
probHinfected2 <- 0.022
ah2 <- 0.01


pH12 <- 0.0001 # host movement from rural to PU 
pH21 <- 0.001 # host movement from PU to rural

pars = c(kv1, Reggs1, Heggs1, meggs1,
         Peggs1, mpup1, pupmaturation1,
         pupsurvival1, mv1, lv1, probVinfected1, 
         biterate1, y1, lh1, mh1, xh1, kh1, probHinfected1,
         ah1, kv2, Reggs2, Heggs2, Peggs2, meggs2, mpup2,
         pupmaturation2, pupsurvival2, mv2, lv2, 
         probVinfected2, biterate2, y2, lh2, mh2,
         xh2, kh2, probHinfected2, ah2, pH12,pH21)

I_V01 = 1000
I_V02 = 1000

S_V01 = 25000
S_V02 = 25000

init <- c(kh1, 0, 0, 0, 
          1000000, 0, S_V01, 0, I_V01,
          kh2, 0, 0, 0, 
          100, 0, S_V02, 0, I_V02)

start = 0
stop = 365*100
t = seq(start, stop, by = 1)