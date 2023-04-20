
SolPatch %>%
  filter(!is.na(Status), Agent == "vector") %>%
  ggplot()+
  geom_line(aes(x=time/365, y=value, colour=Status))+
  facet_wrap(Patch ~ Status, scales = "free", nrow = 2, ncol = 6) +
  theme_bw()+
  theme(axis.text.x = element_text(size = 10), 
        axis.title = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16), 
        strip.text = element_text(size = 12))+
  xlab("years")+
  scale_x_continuous(breaks=seq(0,100,10))+
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.1)))

SolPatch %>%
  filter(!is.na(Status), Agent == "host") %>%
  ggplot()+
  geom_line(aes(x=time/365, y=value, colour=Status))+
  facet_wrap(Patch ~ Status, scales = "free", nrow = 2, ncol = 5) +
  scale_x_continuous(breaks=seq(0,100,10))+
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.1)))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10), 
        axis.title = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16), 
        strip.text = element_text(size = 12))+
  xlab("years")



SolPatch %>%
  filter(!is.na(Status), Patch == "rural", Agent == "host") %>%
  ggplot()+
  geom_line(aes(x=time/365, y=value, colour=Status))+
  facet_wrap(. ~ Status, scales = "free", ncol = 5) +
  theme_bw()+
  theme(axis.text.x = element_text(size = 10), 
        axis.title = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16), 
        strip.text = element_text(size = 12))+
  xlab("years")+
  scale_x_continuous(breaks=seq(0,100,10))+
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.1)))

SolPatch %>%
  filter(!is.na(Status), Patch == "rural", Agent == "vector") %>%
  ggplot()+
  geom_line(aes(x=time/365, y=value, colour=Status))+
  facet_wrap(. ~ Status, scales = "free", ncol = 6) +
  theme_bw()+
  theme(axis.text.x = element_text(size = 10), 
        axis.title = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16), 
        strip.text = element_text(size = 12))+
  xlab("years")+
  scale_x_continuous(breaks=seq(0,100,10))+
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.1)))

SolPatch %>%
  filter(!is.na(Status), Patch == "peri-urban", Agent == "host") %>%
  ggplot()+
  geom_line(aes(x=time/365, y=value, colour=Status))+
  facet_wrap(. ~ Status, scales = "free", ncol = 5) +
  theme_bw()+
  theme(axis.text.x = element_text(size = 10), 
        axis.title = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16), 
        strip.text = element_text(size = 12))+
  xlab("years")+
  scale_x_continuous(breaks=seq(0,100,10))+
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.1)))

SolPatch %>%
  filter(!is.na(Status), Patch == "peri-urban", Agent == "vector") %>%
  ggplot()+
  geom_line(aes(x=time/365, y=value, colour=Status))+
  facet_wrap(. ~ Status, scales = "free", ncol = 6) +
  theme_bw()+
  theme(axis.text.x = element_text(size = 10), 
        axis.title = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16), 
        strip.text = element_text(size = 12))+
  xlab("years")+
  scale_x_continuous(breaks=seq(0,100,10))+
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.1)))
