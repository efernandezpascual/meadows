library(tidyverse)

# Read tree

phangorn::nnls.tree(cophenetic(ape::read.tree("data/meadowstree.tree")), 
                    ape::read.tree("data/meadowstree.tree"), rooted = TRUE) -> 
  nnls_orig

# Set number of iterations

nite = 500000
nthi = 50
nbur = 50000

# Set priors from https://doi.org/10.1017/S0960258517000332

list(R = list(V = 1, nu = 50), 
     G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500), 
              G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500),
              G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 500)
     )) -> 
  priors

# Read data

read.csv("data/germination.csv") %>%
  merge(read.csv("data/traits.csv"), by = "Taxon") %>%
  mutate(ID = gsub(" ", "_", Taxon), animal = ID) %>% 
  select(ID, animal, Family,
         Reference:Population, Scarification:GA3,
         Tmean:Light, Germinated:Germinable) -> 
  data

# All species

MCMCglmm::MCMCglmm(cbind(Germinated, Germinable - Germinated) ~ scale(Scarification) + 
                     scale(Stratification) + scale(GA3) + scale(Tmean) + 
                     scale(Alternating) + scale(Light), 
                   random = ~ animal + Reference + Population:Reference, 
                   family = "multinomial2", pedigree = nnls_orig, prior = priors, data = data,
                   nitt = nite, thin = nthi, burnin = nbur, 
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                   pr = FALSE, pl = FALSE) -> allspp

summary(allspp)
save(allspp, file = "results/models/allspp.Rdata")
rm(allspp)

# Poaceae

data %>% filter(Family == "Poaceae") -> poaceaedf

MCMCglmm::MCMCglmm(cbind(Germinated, Germinable - Germinated) ~ scale(Scarification) + 
                     scale(Stratification) + scale(GA3) + scale(Tmean) + 
                     scale(Alternating) + scale(Light), 
                   random = ~ animal + Reference + Population:Reference, 
                   family = "multinomial2", pedigree = nnls_orig, prior = priors, data = poaceaedf,
                   nitt = nite, thin = nthi, burnin = nbur, 
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                   pr = FALSE, pl = FALSE) -> poaceae

summary(poaceae)
save(poaceae, file = "results/models/poaceae.Rdata")
rm(poaceaedf, poaceae)

# Fabaceae

data %>% filter(Family == "Leguminosae") -> leguminosaedf

MCMCglmm::MCMCglmm(cbind(Germinated, Germinable - Germinated) ~ scale(Scarification) + 
                     scale(Stratification) + scale(GA3) + scale(Tmean) + 
                     scale(Alternating) + scale(Light), 
                   random = ~ animal + Reference + Population:Reference, 
                   family = "multinomial2", pedigree = nnls_orig, prior = priors, data = leguminosaedf,
                   nitt = nite, thin = nthi, burnin = nbur, 
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                   pr = FALSE, pl = FALSE) -> leguminosae

summary(leguminosae)
save(leguminosae, file = "results/models/leguminosae.Rdata")
rm(leguminosaedf, leguminosae)

# Rest

data %>% filter(! Family %in% c("Poaceae", "Leguminosae")) -> restdf

MCMCglmm::MCMCglmm(cbind(Germinated, Germinable - Germinated) ~ scale(Scarification) + 
                     scale(Stratification) + scale(GA3) + scale(Tmean) + 
                     scale(Alternating) + scale(Light), 
                   random = ~ animal + Reference + Population:Reference, 
                   family = "multinomial2", pedigree = nnls_orig, prior = priors, data = restdf,
                   nitt = nite, thin = nthi, burnin = nbur, 
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, 
                   pr = FALSE, pl = FALSE) -> rest

summary(rest)
save(rest, file = "results/models/rest.Rdata")
rm(restdf, rest)