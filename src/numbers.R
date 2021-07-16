# Read data

read.csv("../data/germination.csv") -> germination
read.csv("../data/traits.csv") -> traits
load("../results/models/allspp.RData")

# Calculate numbers

germination %>% pull(Taxon) %>% length -> MSrecords
germination %>% pull(Taxon) %>% unique %>% length -> MSspeciesN
germination %>% pull(Country) %>% unique %>% length -> MScountries
sum(germination$Germinable) -> MSseedsN
min(germination$Tmean) -> MSminT
max(germination$Tmean) -> MSmaxT
germination %>% filter(Alternating == 1) %>% pull(Taxon) %>% length -> MSaltY
germination %>% filter(Alternating == 0) %>% pull(Taxon) %>% length -> MSaltN
germination %>% filter(Light == 1) %>% pull(Taxon) %>% length -> MSlightY
germination %>% filter(Light == 0) %>% pull(Taxon) %>% length -> MSlightN
germination %>% filter(Stratification == 0) %>% pull(Taxon) %>% length -> MSstratN
germination %>% filter(Stratification == 1) %>% pull(Taxon) %>% length -> MSstratY
germination %>% filter(Scarification == 0) %>% pull(Taxon) %>% length -> MSscarN
germination %>% filter(Scarification == 1) %>% pull(Taxon) %>% length -> MSscarY
germination %>% filter(GA3 == 0) %>% pull(Taxon) %>% length -> MSgaN
germination %>% filter(GA3 == 1) %>% pull(Taxon) %>% length -> MSgaY

# Calculate lambda http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11/chapter-11-1-simple-model-mcmcglmm

lambda <- allspp$VCV[,"animal"]/(allspp$VCV[,"animal"] + allspp$VCV[,"units"])

mean(lambda) %>% round(2) -> MSlambda
coda::HPDinterval(lambda)[, 1] %>% round(2) -> MSlambdal
coda::HPDinterval(lambda)[, 2] %>% round(2) -> MSlambdah

# Random effects

summary(allspp)$Gcovariances[1, 1] %>% round(2) -> MSphylo
summary(allspp)$Gcovariances[1, 2] %>% round(2) -> MSphylol
summary(allspp)$Gcovariances[1, 3] %>% round(2) -> MSphyloh

summary(allspp)$Gcovariances[3, 1] %>% round(2) -> MSpopulation
summary(allspp)$Gcovariances[3, 2] %>% round(2) -> MSpopulationl
summary(allspp)$Gcovariances[3, 3] %>% round(2) -> MSpopulationh

summary(allspp)$Gcovariances[2, 1] %>% round(2) -> MSsource
summary(allspp)$Gcovariances[2, 2] %>% round(2) -> MSsourcel
summary(allspp)$Gcovariances[2, 3] %>% round(2) -> MSsourceh

# PCA numbers

read.csv("../data/germination.csv") %>%
  mutate(Germination = Germinated / Germinable) %>%
  group_by(Taxon) %>%
  summarise(Temperature = weighted.mean(Tmean, w = Germination),
            Alternating = weighted.mean(Alternating, w = Germination),
            Light = weighted.mean(Light, w = Germination),
            Scarification = weighted.mean(Scarification, w = Germination),
            Stratification = weighted.mean(Stratification, w = Germination),
            GA3 = weighted.mean(GA3, w = Germination)) %>%
  merge(read.csv("../data/traits.csv"), by = "Taxon") %>%
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

pca$var$contrib


