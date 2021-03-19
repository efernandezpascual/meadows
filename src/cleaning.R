library(tidyverse); library(openxlsx)

# Clean species data

read.xlsx("../#data/plots/data/sospraderas/Inventarios SOS PRADERAS privacidad.xlsx", sheet = 1) %>% # These are the relevés in the DB
  merge(read.xlsx("../#data/plots/data/sospraderas/Inventarios SOS PRADERAS privacidad.xlsx", sheet = 2), # These are the species names in the DB
        by.x = "codigo_especie", by.y = "ID") %>%
  merge(read.csv("../#data/plots/data/sospraderas/Especies revisadas 20210108.csv")) %>% # Revised spp names
  mutate(valor = fct_recode(valor, # Braun-Blanquet to cover scale
                            "87.5" = "5",
                            "62.5" = "4",
                            "37.5" = "3",
                            "15" = "2",
                            "5" = "1",
                            "1" = "+", 
                            "0.1" = "r"),
         valor = as.character(valor), 
         valor = as.numeric(valor)) %>%
  rename(Plot = Cod_Sitio,
         Cover = valor) %>%
  dplyr::select(Plot, Species, Cover) %>% 
  arrange(Plot) %>%
  group_by(Plot, Species) %>%
  summarise(Cover = max(Cover)) -> # remove species that may have been merged by taxonomy
  species # clean species tibble

# Clean header data

read.xlsx("../#data/plots/data/sospraderas/Inventarios SOS PRADERAS privacidad.xlsx", sheet = 3) %>% # This is the header data in the database
  merge(read.csv("../#data/plots/data/sospraderas/Coordenadas ETRS.csv")) %>%
  mutate(Date = convertToDate(Data)) %>%
  filter(Date > 2001-01-01 & !is.na(Latitude)) %>% # Only 21st Century, with coordinates
  mutate(`SIC.(LIC)` = as.character(`SIC.(LIC)`), 
         Region = ifelse(`SIC.(LIC)` %in% "PNPE", "Cantabrian Mountains", Pais),
         Region = fct_recode(Region, 
                             "French Pyrenees" = "França", 
                             "Spanish Pyrenees" = "Espanha",
                             "Northern Portugal" = "Portugal")) %>%
  dplyr::select(Cod_Sitio, `Área.Mínima`, Date, xcoord, ycoord, `Nivel.de.Precisão`, Region) %>%
  rename(Plot = Cod_Sitio,
         Longitude = xcoord, Latitude = ycoord,
         Area = `Área.Mínima`, Precision = `Nivel.de.Precisão`) %>%
  filter(! Region %in% "French Pyrenees") %>%
  filter(! Plot %in% "419") -> # Plot with only presence (+) Malva moschata, transcription error
  header0

# Clean soil data

read.xlsx("../#data/plots/data/sospraderas/Analisis_Suelos_SOSPraderas_16-02-2018.xlsx") %>% # Extra data not in DB, this file contains the same information as the DB but also additional parameters
  dplyr::select(-c(`B.(mg/kg)`, `ACCE.(%)`, `CCE.(%)`, X3)) %>% # Mostly empty cells
  na.omit %>% # Drop the rows used to label groups
  merge(read.csv("../#data/plots/data/sospraderas/plotIDs.csv"), by = "X2") %>% # This is a file I had to do manually to fix the mess of IDs in the extra soil data file
  dplyr::select_if(~ any(is.numeric(.))) %>% # Remove non-numeric variables
  rename(Plot = Cod_Sitio,
         pH = `pH.H2O.(1:2.5)`,
         Fertility = `CICE.(cmol+/kg)`,
         Sand = `Arena.(%)`) %>%
  select(Plot, pH, Sand) -> soils

# Get CHELSA bioclimatic variables
# CHELSA rasters are stored outside project folder, must download from https://chelsa-climate.org/bioclim/

bio06 <- raster::raster("../#data/maps/CHELSA/CHELSA_bio10_06.tif") # Min temperature coldest month
bio14 <- raster::raster("../#data/maps/CHELSA/CHELSA_bio10_14.tif") # Precipitation driest month

pts <- header0
sp::coordinates(pts) <- ~Longitude+Latitude # Make points a spatial object
sp::proj4string(pts) <- sp::CRS(SRS_string = "OGC:CRS84") # Asign projection
Chelsa <- data.frame(sp::coordinates(pts)) # Create data.frame with the Enscobase points
Chelsa$bio06 <- raster::extract(bio06, pts)/10 # Extract from CHELSA layer, correct units
Chelsa$bio14 <- raster::extract(bio14, pts) # Extract from CHELSA layer

# Update header0 with environment

header0 %>%
  merge(Chelsa, by = c("Longitude", "Latitude")) %>%
  merge(soils) ->
  header

rm(header0, Chelsa, soils, bio06, bio14, pts)

# Calculate SNCs

header %>%
  merge(species) %>%
  select(Species, bio06:Sand, Cover) %>%
  gather(Trait, Value, bio06:Sand) %>%
  group_by(Species, Trait) %>%
  summarise(Value = weighted.mean(Value, Cover)) %>%
  spread(Trait, Value) %>%
  rename(Taxon = Species) -> sncs

# Calculate species dominance

species %>%
  group_by(Plot) %>%
  summarise(TCover = sum(Cover)) %>% # Max cover in each plot
  merge(species) %>%
  mutate(RCover = Cover/TCover) %>% # Relative species cover by plot's total cover
  merge(header) %>%
  group_by(Species) %>%
  summarise(Frequency = length(Species), Abundance = sum(Cover)) %>% # Relative abundance
  mutate(Abundance = scales::rescale(Abundance, to = c(1, 100))) %>% # Cumulative relative abundance (rescaled 1-100)
  mutate(Frequency = 100 * Frequency / length(unique(species$Plot))) %>%
  mutate(Dominance = ifelse(Abundance >= 12, "Dominant", "Subordinate"),
         Dominance = ifelse(Abundance < 2, "Transient", Dominance)) %>%
  select(Species, Abundance, Dominance) %>%
  arrange(-Abundance) %>%
  rename(Taxon = Species) -> dominances

# Clean seedlot data

read.csv("../#data/germination/data/sospraderas/Muestras de semillas.csv") %>% # Seedlot data from project's final report
  merge(read.csv("../#data/germination/data/sospraderas/Localidades de recolección.csv"), by = "ID.Localidad") %>% # Collecting site data from project's final report
  mutate(Date = as.Date(Fecha.de.recolección, "%d/%m/%Y")) %>%
  rename(Longitude = xcoord, 
         Latitude = ycoord,
         Seedlot = ID.Recolección,
         Tabela.de.Registo.de.Inventários.de.Espécies = Taxón,
         Collectors = Recolectores,
         Site = Descripción,
         Municipality = Municipio,
         Province = Provincia,
         Region = Región,
         Country = País) %>%
  merge(read.csv("../#data/germination/data/sospraderas/Especies revisadas 20210108.csv")) %>% # Revised spp names
  dplyr::select(Seedlot,
                Species,
                Collectors,
                Date,
                Longitude,
                Latitude,
                Site,
                Municipality,
                Province,
                Region,
                Country) ->
  seedlots

# Clean morphometric analysis data

read.csv("../#data/germination/data/sospraderas/Accesiones morfo.csv") %>% # Seedlot information for the morphological analysis
  filter(Tipo.de.espécies == "Espécies silvestres") %>% # Commercial seed lots of the project, are of no interest here
  dplyr::select(ID, Label) %>% 
  merge(read.csv("../#data/germination/data/sospraderas/Análise morfométrica.csv"), by = "Label") %>% # Results of the morphometric analysis as provided by Madalena Vaz
  dplyr::select(-c(Label, X.1)) %>%
  rename(Seedlot = ID) %>% # This is the seedlot ID code as used in the final report of the project
  merge(seedlots, by = "Seedlot") %>%
  rename(Taxon = Species) %>%
  group_by(Taxon) %>%
  summarise(Length = mean(Height), Width = mean(Width)) -> morphometrics

# Clean seed mass

read.csv("../#data/germination/data/enscobase meadows/Species_domimance_in_ENSCOBASE.csv") %>%
  mutate(Taxon = gsub("_", " ", animal)) %>%
  select(Taxon, mass) %>% 
  rename(Seed.mass = mass) %>%
  unique -> mass1

read.csv("../#data/germination/data/sospraderas/SID Seed Mass.csv") %>%
  rename(Taxon = Species,
         Seed.mass = Seed.mass.mg) -> mass2

read.csv("../#data/seedmass/results/seedmass.csv") %>%
  rename(Taxon = TPLName) %>%
  filter(Taxon %in% species$Species) -> mass3

rbind(mass1, mass2, mass3) %>%
  group_by(Taxon) %>%
  summarise(Seed.mass = mean(Seed.mass)) -> mass

rm(mass1, mass2, mass3)

# Get dormancy classes

read.csv("../#data/baskin/results/dormancy.csv") %>%
  rename(Taxon = TPLName) %>%
  filter(Taxon %in% species$Species) -> dormancy1

read.csv("../#data/germination/data/enscobase meadows/Species_domimance_in_ENSCOBASE.csv") %>%
  mutate(Taxon = gsub("_", " ", animal)) %>%
  select(Taxon, dormancy) %>% 
  rename(Dormancy = dormancy) %>%
  filter(! Taxon %in% dormancy1$Taxon) %>%
  unique -> dormancy2

rbind(dormancy1, dormancy2) -> dormancy

rm(dormancy1, dormancy2)

# Merge species traits

dominances %>%
  filter(Dominance != "Transient") %>%
  merge(sncs, all.x = TRUE) %>%
  merge(mass, all.x = TRUE) %>%
  merge(morphometrics, all.x = TRUE) %>%
  merge(dormancy, all.x = TRUE) %>%
  merge(read.csv("../#data/tpl/results/TPLNames.csv")) %>%
  select(Taxon, Family, Abundance, Dominance, bio06, bio14, pH, Sand, Seed.mass, Length, Width, Dormancy) -> traits0

rm(dominances, dormancy, mass, morphometrics, sncs, species, header)

# Extract higher taxonomical ranks

traits0$Family %>%
  unique %>%
  as.character() -> families

taxize::classification(families, db = "gbif") -> 
  taxonomy

as.data.frame(t(sapply(names(taxonomy), function (x) taxonomy[[x]] [, 1])[c(5, 4), ])) %>%
  remove_rownames() %>%
  rename(Family = V1, order = V2) -> orders

traits0 %>%
  merge(orders, all.x = TRUE) %>%
  merge(read.csv("../#data/germination/data/sospraderas/clades.csv"), all.x = TRUE) %>%
  rename(Order = order, Clade = clade) %>%
  select(Taxon, Family, Order, Clade,
         Abundance, Dominance, bio06, bio14, pH, Sand, Seed.mass, Length, Width, Dormancy) -> traits

rm(orders, families, taxonomy, traits0)

# Clean germination data from SOS PRADERAS

read.csv("../#data/germination/data/sospraderas/IDs.csv") %>% # File with the different IDs used during the project
  dplyr::select(ID.Recolección, ID.Biesca) %>%
  merge(read.csv("../#data/germination/data/sospraderas/Germinación.csv"), by = "ID.Biesca") %>% # Germination dataset as provided with the project's final report
  merge(seedlots, by.x = "ID.Recolección", by.y = "Seedlot") %>%
  mutate(Scarification = ifelse(Pretratamiento == "Escarificado", 1, 0), # Binary scarification
         GA3 = ifelse(Pretratamiento == "GA3", 1, 0), # Binary GA3
         Stratification = 0, # All fresh
         Light = 1, # All in light 22/12h
         Alternating = 1, # All alternating 10C
         Germinated = Germinadas.S1 + Germinadas.S2 + Germinadas.S3 + Germinadas.S4,
         Germinable = Germinated + Normales,
         Temperature = gsub("C", "", Temperatura),
         Reference = "SOSPRADERAS") %>%
  separate(Temperature, into = c("T1", "T2"), sep = "/") %>%
  mutate(Database = "SOSPRADERAS") %>%
  mutate(Tmean = (as.numeric(T1) + as.numeric(T2)) / 2) %>% # Average germination temperature
  rename(Taxon = Species,
         Population = ID.Recolección,
         Pretreatment = Pretratamiento,
         Temperature = Temperatura,
         Dish = Placa,
         Germinated.w1 = Germinadas.S1,
         Germinated.w2 = Germinadas.S2,
         Germinated.w3 = Germinadas.S3,
         Germinated.w4 = Germinadas.S4,
         Mouldy = Contaminadas,
         Empty = Vacías) %>%
  select(Taxon, Database, Reference, Country, Population,
         Scarification, Stratification, GA3,
         Tmean, Alternating, Light,
         Germinated, Germinable) ->
  germination1

# Germination data from Enscobase (Angelino Carta)

read.csv("../#data/germination/data/enscobase meadows/Species_domimance_in_ENSCOBASE.csv") %>% 
  filter(aged == 0) %>% # Remove ageing experiments
  filter(chemical != "KNO3") %>% # remove, rare treatment
  mutate(Taxon = gsub("_", " ", animal)) %>% # Original taxon name follows EuroMed
  mutate(Stratification = ifelse(cold2 == 1 | warm2 == 1, 1, 0)) %>% # Binary stratification yes/no
  mutate(GA3 = ifelse(chemical == "GA3", 1, 0)) %>% # Binary GA3 yes/no
  mutate(Database = "Enscobase") %>%
  rename(Population = acc.orig, # Seedlot ID
         Reference = germ_inst, # Reference or Lab ID
         Country = country_collection, # Country of collection
         Tmean = T2, # Mean germination temperature
         Light = L2,
         Alternating = Alt2,
         Scarification = scarification,
         Germinated = y,
         Germinable = n) %>%
  select(Taxon, Database, Reference, Country, Population,
         Scarification, Stratification, GA3,
         Tmean, Alternating, Light,
         Germinated, Germinable) ->
  germination2

# Germination from SylvanSeeds

read.csv("../#data/germination/results/TBMF_Database.csv") %>%
  filter(! Country %in% c("USA", "Argentina", "Uruguay", "Canada", "Australia", 
                          "Brazil", "Pakistan", "South Korea", "Chile", "New Zealand",
                          "Syria", "Japan", "Iran")) %>% # Remove non-European records
  mutate(GA3 = 0) %>% # All non-GA3
  mutate(Database = "SylvanSeeds") %>%
  mutate(Stratification = ifelse(Stratification == "N", 0, 1)) %>%
  mutate(Light = ifelse(Light == "Y" | is.na(Light), 1, 0)) %>%
  mutate(Alternating = ifelse(Alternating == "N", 0, 1)) %>%
  mutate(Scarification = ifelse(Scarification == "N", 0, 1)) %>%
  rename(Taxon = ï..Taxon,
         Germinable = Number_seeds) %>%
  select(Taxon, Database, Reference, Country, Population,
         Scarification, Stratification, GA3,
         Tmean, Alternating, Light,
         Germinated, Germinable) ->
  germination3

# Merge databases

rbind(germination1, germination2, germination3) %>%
  mutate(Country = ifelse(Country == "United-Kingdom", "UK", Country)) %>%
  mutate(Country = ifelse(Country == "España", "Spain", Country)) %>%
  mutate(Country = ifelse(Country == "Francia", "France", Country)) %>%
  filter(Germinated <= Germinable) %>% # Check no cases of more germination than seeds
  filter(Germinable > 0) %>% # At least one germinable seed!
  group_by(Taxon, Database, Reference, Country, Population,
           Scarification, Stratification, GA3,
           Tmean, Alternating, Light) %>%
  summarise(Germinated = sum(Germinated), Germinable = sum(Germinable)) %>% # Merge dishes
  filter(Taxon %in% traits$Taxon) %>%
  group_by() -> germination

rm(germination1, germination2, germination3)

# Phylogenetic tree

traits %>%
  filter(Taxon %in% germination$Taxon) %>%
  separate(Taxon, into = c("Genus", "Species"), sep = " ") %>%
  mutate(species = paste(Genus, Species),
         genus = Genus,
         family = Family) %>%
  select(species, genus, family) %>%
  unique %>%
  mutate(family = fct_recode(family, 
                             "Asteraceae" = "Compositae",
                             "Fabaceae" = "Leguminosae",
                             "Asphodelaceae" = "Xanthorrhoeaceae")) %>%
  arrange(species) %>%
  na.omit -> 
  ranks1

library(V.PhyloMaker)

phylo.maker(sp.list = ranks1, 
            tree = GBOTB.extended, 
            nodes = nodes.info.1, 
            scenarios = "S3") ->
  tree

rm(ranks1)

# Save

write.tree(tree$scenario.3, file = "data/meadowstree.tree")
write.csv(germination, file = "data/germination.csv", row.names = FALSE)
traits %>%
  filter(Taxon %in% germination$Taxon) %>%
  write.csv(file = "data/traits.csv", row.names = FALSE)
