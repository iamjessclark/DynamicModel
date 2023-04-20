rm(list = ls())
## plotting results table as a chart
sim.res.tab <- (read.csv("interventions.csv"))
sim.res.tab$Reduction <- sim.res.tab$red.prop * 100
sim.res.tab$se <- sim.res.tab$stand.error * 100

library(ggplot2)
str(sim.res.tab)

library("RColorBrewer")

display.brewer.all()

#sim.res.tab$Reduction <- sim.res.tab$Reduction * 100

#sim.res.tab$se <- sim.res.tab$se * 100

levels(sim.res.tab$Intervention)

x <- ggplot(sim.res.tab, aes(x = Intervention, y = Reduction, fill=Ro))+
  geom_bar(position = position_dodge(), stat = "identity") +
  geom_errorbar(aes(ymin= Reduction-se, ymax = Reduction + se),
                width=.2,
                position = position_dodge(.9))+
  xlab("Intervention type") +
  ylab("% Reduction in population cumulative incidence ")+
  scale_fill_hue(name="Pathogen",
                 breaks=c("1", "2"),
                 labels=c("Fast", "Slow")) + #, palette = "rainbow")+
  scale_fill_brewer(palette = "Set1", labels=c("Fast", "Slow"), name="Pathogen")+
  #ggtitle("Reduction in population cumulative incidence of disease after 1 year with simulated \n fast and slow pathogens and different interventions applied at a ward level ")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 50
                                   ,hjust = 1))
  

ggsave("PhilTransFig3.pdf", x, width = 7, height = 5)


