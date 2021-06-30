library(tidyverse)

# Figure Seed morphology

## Plot figure A

read.csv("data/traits.csv") %>%
  select(Taxon, Family,
         Seed.mass, Seed.number) %>%
  mutate(Family = ifelse(Family %in% c("Leguminosae", "Poaceae"), Family, "Others"),
         Family = ifelse(Family %in% "Leguminosae", "Fabaceae", Family)) %>%
  gather(Trait, Value, Seed.mass:Seed.number) %>%
  mutate(Family = fct_relevel(Family, c("Poaceae", "Fabaceae", "Others"))) %>%
  mutate(Trait = fct_relevel(Trait, c("Seed.mass", "Seed.number"))) %>%
  mutate(Trait = fct_recode(Trait, 
                            "Seed mass (log)" = "Seed.mass",
                            "Seed number (log)" = "Seed.number")) %>%
  mutate(Value = ifelse(! Trait %in% c("Seed roundness"), log(Value), Value)) %>%
  ggplot(aes(Family, Value, color = Family, fill = Family)) + 
  geom_violin(alpha = 0.5, draw_quantiles = c(0.25, 0.5, 0.75)) +
  facet_wrap( ~ Trait, strip.position = "left", scales = "free_y", ncol = 1) +
  scale_y_continuous(labels = scales::number_format(accuracy = .1)) +
  ggthemes::theme_tufte() +
  theme(panel.background = element_rect(color = "grey96", fill = "grey96"),
        strip.text = element_text(size = 12, color = "black"),
        axis.text = element_text(size = 11, color = "black"),
        axis.title = element_blank(),
        legend.position = "none", strip.placement = "outside") +
  geom_jitter(shape = 16, position = position_jitter(0.05)) +
  scale_color_manual(values = c("yellowgreen",  "gold", "darkorchid")) +
  scale_fill_manual(values = c("yellowgreen",  "gold", "darkorchid")) -> f2a; f2a

## Plot figure B

read.csv("data/traits.csv") %>%
  select(Taxon, Family) %>%
  merge(read.csv("data/morphometrics.csv"), by = "Taxon") %>%
  select(Taxon, Family) %>%
  unique %>%
  group_by(Family) %>%
  tally() %>%
  pull(Family) -> fams

read.csv("data/traits.csv") %>%
  select(Taxon, Family) %>%
  merge(read.csv("data/morphometrics.csv"), by = "Taxon") %>% 
  filter(Family %in% fams) %>%
  mutate(Family = ifelse(Family %in% "Compositae", "Asteraceae", Family),
         Family = ifelse(Family %in% "Leguminosae", "Fabaceae", Family),
         Family = ifelse(Family %in% c("Poaceae", "Fabaceae"), Family, "Others")) -> morpho

cent <- aggregate(cbind(Width, Length) ~ Family, data = morpho, FUN = mean)
segs <- merge(morpho, setNames(cent, c("Family", "oDCA1", "oDCA2")), by = "Family", sort = FALSE)

morpho %>%
  ggplot(aes(Width, Length, color = Family)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_segment(aes(x = 0, xend = 15, y = 0, yend = 15), linetype = "dashed", color = "black") +
  geom_point(aes(color = Family), alpha = 0.2, size = 3) +
  geom_label(data = cent, size = 4, aes(label = Family, fill = Family), color = "white", show.legend = F) + 
  ggthemes::theme_tufte() + 
  coord_cartesian(xlim = c(0, 9), ylim = c(0, 14)) +
  theme(legend.position = "none",
        panel.background = element_rect(color = "grey96", fill = "grey96"),
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black")) +
  scale_color_manual(values = c("gold", "darkorchid", "yellowgreen")) +
  scale_fill_manual(values = c("gold", "darkorchid", "yellowgreen")) +
  annotate("text", x = 1.4, y = 13.5, label = "Elongated\nseeds", color = "grey50") +
  annotate("text", x = 7.7, y = 4, label = "Rounded\nseeds", color = "grey50") +
  scale_x_continuous(name = "Width (mm)") + 
  scale_y_continuous(name = "Length (mm)") -> f2b

## Merge

gridExtra::arrangeGrob(f2a, f2b, nrow = 1) -> f2

## Export

ggsave(f2, file = "results/figures/Fig2.png", 
       path = NULL, scale = 1, width = 173, height = 110, units = "mm", dpi = 600)

# Figure MCMC

## Load models

load("results/models/allspp.RData")
load("results/models/leguminosae.RData")
load("results/models/poaceae.RData")
load("results/models/rest.RData")

## Bind effect sizes data frame

rbind(`All species:Scarification` = summary(allspp)$solutions[2, ],
      `All species:Stratification` = summary(allspp)$solutions[3, ],
      `All species:GA3` = summary(allspp)$solutions[4, ],
      `All species:Temperature` = summary(allspp)$solutions[5, ],
      `All species:Alternating` = summary(allspp)$solutions[6, ],
      `All species:Light` = summary(allspp)$solutions[7, ],
      
      `Poaceae:Scarification` = summary(poaceae)$solutions[2, ],
      `Poaceae:Stratification` = summary(poaceae)$solutions[3, ],
      `Poaceae:GA3` = summary(poaceae)$solutions[4, ],
      `Poaceae:Temperature` = summary(poaceae)$solutions[5, ],
      `Poaceae:Alternating` = summary(poaceae)$solutions[6, ],
      `Poaceae:Light` = summary(poaceae)$solutions[7, ],
      
      `Fabaceae:Scarification` = summary(leguminosae)$solutions[2, ],
      `Fabaceae:Stratification` = summary(leguminosae)$solutions[3, ],
      `Fabaceae:GA3` = summary(leguminosae)$solutions[4, ],
      `Fabaceae:Temperature` = summary(leguminosae)$solutions[5, ],
      `Fabaceae:Alternating` = summary(leguminosae)$solutions[6, ],
      `Fabaceae:Light` = summary(leguminosae)$solutions[7, ],
      
      `Others:Scarification` = summary(rest)$solutions[2, ],
      `Others:Stratification` = summary(rest)$solutions[3, ],
      `Others:GA3` = summary(rest)$solutions[4, ],
      `Others:Temperature` = summary(rest)$solutions[5, ],
      `Others:Alternating` = summary(rest)$solutions[6, ],
      `Others:Light` = summary(rest)$solutions[7, ]) %>%
  data.frame %>%
  rownames_to_column(var = "Model") %>%
  separate(Model, into = c("Group", "Effect"), sep = ":") %>%
  mutate(Effect = fct_relevel(Effect, c("Light", "Alternating", "Temperature", 
                                        "GA3", "Stratification", "Scarification")),
         Group = fct_relevel(Group, c("All species", "Poaceae", "Fabaceae", "Others"))) -> efectos

## Plot figure

efectos %>%
  ggplot(aes(y = Effect, x = post.mean, 
             xmin = l.95..CI, xmax = u.95..CI,
             color = Effect)) +
  facet_wrap(~ Group, scales = "free_x", nrow = 1) +
  geom_point(size = 2) +
  labs(x = "Effect size") +
  geom_errorbarh(height = .3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ggthemes::theme_tufte() +
  theme(legend.position = "none", axis.title.y = element_blank(),
        panel.background = element_rect(color = "grey96", fill = "grey96"),
        axis.text.y = element_text(size = 12,
                                   color = c("gold",
                                             "purple", 
                                             "turquoise4", 
                                             "yellowgreen", 
                                             "steelblue1", 
                                             "olivedrab")),
        axis.text.x = element_text(size = 7.5, color = "black"),
        strip.text.x = element_text(size = 12)) +
  scale_color_manual(values = c("gold",
                                "purple", 
                                "turquoise4", 
                                "yellowgreen", 
                                "steelblue1", 
                                "olivedrab")) -> f3

## Export

ggsave(f3, file = "results/figures/Fig3.png", 
       path = NULL, scale = 1, width = 173, height = 65, units = "mm", dpi = 600)

# Figure Ordination with environment

## Prepare data frame
 
read.csv("data/germination.csv") %>%
  mutate(Germination = Germinated / Germinable) %>%
  group_by(Taxon) %>%
  summarise(Temperature = weighted.mean(Tmean, w = Germination),
            Alternating = weighted.mean(Alternating, w = Germination),
            Light = weighted.mean(Light, w = Germination),
            Scarification = weighted.mean(Scarification, w = Germination),
            Stratification = weighted.mean(Stratification, w = Germination),
            GA3 = weighted.mean(GA3, w = Germination)) %>%
  merge(read.csv("data/traits.csv"), by = "Taxon") %>%
  mutate(Family = ifelse(Family %in% c("Leguminosae", "Poaceae"), Family, "Others"),
         Family = ifelse(Family %in% "Leguminosae", "Fabaceae", Family)) %>%
  mutate(Family = fct_relevel(Family, c("Poaceae", "Fabaceae", "Others"))) %>%
  select(Taxon, Family, Temperature:Stratification, Seed.mass, Seed.number, bio06:pH) %>%
  na.omit() -> traits

## Do PCA

traits %>%
  gather(Trait, Value, Temperature:Seed.mass) %>%
  mutate(Value = log(0.01+ Value)) %>%
  spread(Trait, Value) %>%
  select(-c(Taxon, Family)) %>%
  FactoMineR::PCA(graph = FALSE) -> pca

traits %>%
  select(Taxon, Family) %>%
  cbind(pca$ind$coord[, 1:2]) %>%
  filter(Dim.2 > -3) -> pcaInds

pca$var$coord[, 1:2] %>%
  data.frame %>%
  rownames_to_column(var = "Variable") %>%
  mutate(Variable = fct_recode(Variable, `Seed mass` = "Seed.mass", `Seed number` = "Seed.number",
                               `Embryo:seed` = "Embryo")) %>%
  mutate(Type = ifelse(Variable %in% c("bio06", "bio14", "pH", "Sand"), "Environment", "Seed")) -> pcaVars

## Plot PCA

ggplot(pcaInds, aes(x = Dim.1, y = Dim.2)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(aes(color = Family), alpha = 0.8, size = 8) +
  ggthemes::theme_tufte() + 
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 12, color = "black"),
        panel.background = element_rect(color = "grey96", fill = "grey96"),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black")) +
  scale_x_continuous(name = paste("Axis 1 (", round(pca$eig[1, 2], 0),
                                  "% variance explained)", sep = "")) + 
  scale_y_continuous(name = paste("Axis 2 (", round(pca$eig[2, 2], 0), 
                                  "% variance explained)", sep = "")) +
  coord_cartesian(xlim = c(-4.60, 4), ylim = c(-3.5, 3)) +
  geom_segment(data = pcaVars, aes(x = 0, y = 0, xend = 4*Dim.1, yend = 4*Dim.2)) +
  geom_label(data = filter(pcaVars, ! Variable %in% c("Light", "Temperature", "pH", "bio14", "Stratification")), aes(x = 4*Dim.1, y = 4*Dim.2, label = Variable, fill = Type),  show.legend = FALSE, size = 4) +
  geom_label(data = filter(pcaVars, Variable %in% c("Light")), aes(x = .45 + (4*Dim.1), y = 4*Dim.2, label = Variable), size = 4) +
  geom_label(data = filter(pcaVars, Variable %in% c("Temperature")), aes(x = -.60 + (4*Dim.1), y = -.15 + 4*Dim.2, label = Variable), size = 4) +
  geom_label(data = filter(pcaVars, Variable %in% c("pH")), aes(x = .05 + (4*Dim.1), y = -0.1 + (4*Dim.2), label = Variable, fill = Type),  show.legend = FALSE, size = 4) +
  geom_label(data = filter(pcaVars, Variable %in% c("Stratification")), aes(x = -.4 + (4*Dim.1), y =(4*Dim.2), label = Variable, fill = Type),  show.legend = FALSE, size = 4) +
  geom_label(data = filter(pcaVars, Variable %in% c("bio14")), aes(x = -.24 + (4*Dim.1), y = 0.25 + (4*Dim.2), label = Variable, fill = Type),  show.legend = FALSE, size = 4) +
  geom_segment(
    x = -4.2, y = 1.5,
    xend = -4.2, yend = -1.5,
    lineend = "round", # See available arrow types in example above
    linejoin = "mitre",
    size = 1, 
    arrow = arrow(length = unit(0.1, "inches")),
    colour = "grey50" # Also accepts "red", "blue' etc
  ) +
  geom_segment(
    x = -4.2, y = -1.5,
    xend = -4.2, yend = 1.5,
    lineend = "round", # See available arrow types in example above
    linejoin = "mitre",
    size = 1, 
    arrow = arrow(length = unit(0.1, "inches")),
    colour = "grey50" # Also accepts "red", "blue' etc
  ) +
  geom_segment(
    x = 1.5, y = -3,
    xend = -1.5, yend = -3,
    lineend = "round", # See available arrow types in example above
    linejoin = "mitre",
    size = 1, 
    arrow = arrow(length = unit(0.1, "inches")),
    colour = "grey50" # Also accepts "red", "blue' etc
  ) +
  geom_segment(
    x = -1.5, y = -3,
    xend = 1.5, yend = -3,
    lineend = "round", # See available arrow types in example above
    linejoin = "mitre",
    size = 1, 
    arrow = arrow(length = unit(0.1, "inches")),
    colour = "grey50" # Also accepts "red", "blue' etc
  ) +
  annotate("text", x = -4.2, y = -2.1, label = "Gap-sensing\nMany seeds\nPoaceae", color = "grey50") +
  annotate("text", x = -4.2, y = 2.1, label = "Hardseedness\nHeavy seeds\nFabaceae", color = "grey50") +
  annotate("text", x = -2.4, y = -3, label = "Warm-dry climes\nAcidic soils", color = "grey50") +
  annotate("text", x = 2.4, y = -3, label = "Cool-wet climes\nAlkaline soils", color = "grey50") +
  
  annotate("text", x = 4, y = 2.9, label = "Poaceae", color = "yellowgreen", hjust = 1) +
  annotate("text", x = 4, y = 2.6, label = "Fabaceae", color = "gold", hjust = 1) +
  annotate("text", x = 4, y = 2.3, label = "Other families", color = "darkorchid", hjust = 1) +
  
  scale_color_manual(values = c("yellowgreen",  "gold", "darkorchid")) +
  scale_fill_manual(values = c("grey65",  "white"))-> f4;f4

## Export

ggsave(f4, file = "results/figures/Fig4.png", 
       path = NULL, scale = 1, width = 173, height = 140, units = "mm", dpi = 600)

# Figure other habitats

read.csv("data/germination.csv") %>%
  filter(Reference == "SOSPRADERAS") %>%
  filter(GA3 == 0) %>%
  merge(read.csv("data/traits.csv")) %>%
  filter((Family == "Leguminosae" & Scarification == 1) | Family != "Leguminosae") %>%
  group_by(Taxon, Seed.mass, Tmean) %>%
  summarise(G = sum(Germinated) / sum(Germinable)) %>%
  spread(Tmean, G) %>%
  group_by(Taxon) %>%
  mutate(Habitat = "Mesic meadows") %>%
  rename(F14 = `9`, F22 = `17`, F30 = `25`) %>%
  select(Taxon, Habitat, F14:F30) %>%
  rbind(read.csv("data/others.csv")) %>%
  group_by() -> df1

## Do PCA

df1 %>%
  select(-c(Taxon, Habitat)) %>%
  FactoMineR::PCA(graph = FALSE, scale = FALSE) -> pca

df1 %>%
  select(c(Taxon, Habitat)) %>%
  cbind(pca$ind$coord[, 1:2]) -> pcaInds

pca$var$coord[, 1:2] %>%
  data.frame %>%
  rownames_to_column(var = "Variable") %>%
  mutate(Variable = fct_recode(Variable, `22/12 ºC` = "F22", 
                               `30/20 ºC` = "F30", `14/4 ºC` = "F14")) -> pcaVars

cent <- aggregate(cbind(Dim.1, Dim.2) ~ Habitat, data = pcaInds, FUN = mean)
segs <- merge(pcaInds, setNames(cent, c("Habitat", "oDCA1", "oDCA2")), by = "Habitat", sort = FALSE)

## Plot PCA

ggplot(pcaInds, aes(x = Dim.1, y = Dim.2)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_segment(data = pcaVars, aes(x = 0, y = 0, xend = 3*Dim.1, yend = 3*Dim.2), size = 0.3) +
  geom_label(data = pcaVars, aes(x = 3*Dim.1, y = 3*Dim.2, label = Variable), size = 4) +
  geom_segment(data = segs, mapping = aes(xend = oDCA1, yend = oDCA2, color = Habitat), show.legend = F, alpha = 0.5) +
  geom_point(data = cent, size = 6, aes(color = Habitat), show.legend = T) + 
  ggthemes::theme_tufte() + 
  theme(legend.position = "top",
        panel.background = element_rect(color = "grey96", fill = "grey96"),
        legend.text = element_text(size = 11, color = "black"),
        legend.title = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black")) +
  scale_x_continuous(name = paste("Axis 1 (", round(pca$eig[1, 2], 0),
                                  "% variance explained)", sep = "")) + 
  scale_y_continuous(name = paste("Axis 2 (", round(pca$eig[2, 2], 0), 
                                  "% variance explained)", sep = "")) +
  coord_cartesian(xlim = c(-.9, 1.2), ylim = c(-.65, .65)) +
  geom_segment(
    x = -0.4, y = -0.6,
    xend = 0.4, yend = -0.6,
    lineend = "round", # See available arrow types in example above
    linejoin = "mitre",
    size = 1, 
    arrow = arrow(length = unit(0.1, "inches")),
    colour = "grey50" # Also accepts "red", "blue' etc
  ) +
  geom_segment(
    x = 0.4, y = -0.6,
    xend = -0.4, yend = -0.6,
    lineend = "round", # See available arrow types in example above
    linejoin = "mitre",
    size = 1, 
    arrow = arrow(length = unit(0.1, "inches")),
    colour = "grey50" # Also accepts "red", "blue' etc
  ) +
  geom_segment(
    x = -0.75, y = -0.3,
    xend = -0.75, yend = 0.3,
    lineend = "round", # See available arrow types in example above
    linejoin = "mitre",
    size = 1, 
    arrow = arrow(length = unit(0.1, "inches")),
    colour = "grey50" # Also accepts "red", "blue' etc
  ) +
  geom_segment(
    x = -0.75, y = 0.3,
    xend = -0.75, yend = -0.3,
    lineend = "round", # See available arrow types in example above
    linejoin = "mitre",
    size = 1, 
    arrow = arrow(length = unit(0.1, "inches")),
    colour = "grey50" # Also accepts "red", "blue' etc
  ) +
  annotate("text", x = -.75, y = -0.4, label = "Cold-cued\ngermination", color = "grey50") +
  annotate("text", x = -.75, y = 0.4, label = "Warm-cued\ngermination", color = "grey50") +
  annotate("text", x = .55, y = -0.6, label = "More\ngermination", color = "grey50") +
  annotate("text", x = -.55, y = -0.6, label = "Less\ngermination", color = "grey50") +
  scale_color_manual(values = c("#fdda62",  "saddlebrown", "#3dc6fa", "#5fc765")) -> f5;f5


ggsave(f5, file = "results/figures/Fig5.png", 
       path = NULL, scale = 1, width = 173, height = 140, units = "mm", dpi = 600)