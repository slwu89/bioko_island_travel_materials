library(here)
library(raster)
library(sp)
library(rgdal)

here()

# Travel frequency and duration ####
# Load in 2018 off-island travel data; examine
data.off <- fread(here("data/raw/2018_travel_data/trips_to_EG.csv"))
View(data.off)

# very crude multinomial probability of going to each of these districts on the mainland
# also a very crude model of duration of travel to each district
dt.off <- data.off[, .(t.off.prob = .N, duration = mean(nights, na.rm = TRUE)), by = "dest_district"]
dt.off$t.off.prob <- dt.off$t.off.prob/(sum(dt.off$t.off.prob))

# Find mean PR in the different locations #### 
dt.off$PR <- 0
# Load in PR map from MAP - 2017 median surface
EG.map <- raster::raster(here("data/clean/country_profiles_data/2019_Global_PfPR_GNQ_2017.tiff"))
EG.map.extent <- projectExtent(EG.map, crs(EG.shp))
EG.map.projected <- projectRaster(EG.map, EG.map.extent)
plot(EG.map.projected)
EG.mainland <- crop(EG.map.projected, raster::extent(5e05, 759708.5, .5e05, 2.7e05))
plot(EG.mainland)
plot(EG.shp, add = TRUE)

plot(EG.mainland)
plot(EG.shp[EG.shp$admin2 == "Micomeseng",], add = TRUE)
holder <- values(mask(EG.mainland, EG.shp[EG.shp$admin2 == "Micomeseng",]))
dt.off[dest_district == "micomiseng"]$PR <- mean(holder[!is.na(holder)])

plot(EG.mainland)
plot(EG.shp[EG.shp$admin2 == "Bata",], add = TRUE)
holder <- values(mask(EG.mainland, EG.shp[EG.shp$admin2 == "Bata",]))
dt.off[dest_district == "bata"]$PR <- mean(holder[!is.na(holder)])

plot(EG.mainland)
plot(EG.shp[EG.shp$admin2 == "Niefang",], add = TRUE)
holder <- values(mask(EG.mainland, EG.shp[EG.shp$admin2 == "Niefang",]))
dt.off[dest_district == "niefang"]$PR <- mean(holder[!is.na(holder)])

plot(EG.mainland)
plot(EG.shp[EG.shp$admin2 == "Mongomo",], add = TRUE)
holder <- values(mask(EG.mainland, EG.shp[EG.shp$admin2 == "Mongomo",]))
dt.off[dest_district == "mongomo"]$PR <- mean(holder[!is.na(holder)])

plot(EG.mainland)
plot(EG.shp[EG.shp$admin2 == "Ebibeyín",], add = TRUE)
holder <- values(mask(EG.mainland, EG.shp[EG.shp$admin2 == "Ebibeyín",]))
dt.off[dest_district == "ebebiyin"]$PR <- mean(holder[!is.na(holder)])

plot(EG.mainland)
plot(EG.shp[EG.shp$admin2 == "Evinayong",], add = TRUE)
holder <- values(mask(EG.mainland, EG.shp[EG.shp$admin2 == "Evinayong",]))
dt.off[dest_district == "evinayong"]$PR <- mean(holder[!is.na(holder)])

plot(EG.mainland)
plot(EG.shp[EG.shp$admin2 == "Nsok Nsomo",], add = TRUE)
holder <- values(mask(EG.mainland, EG.shp[EG.shp$admin2 == "Nsok Nsomo",]))
dt.off[dest_district == "nsok nsomo"]$PR <- mean(holder[!is.na(holder)])

plot(EG.mainland)
plot(EG.shp[EG.shp$admin2 == "Nsok Nsomo",], add = TRUE)
holder <- values(mask(EG.mainland, EG.shp[EG.shp$admin2 == "Nsok Nsomo",]))
dt.off[dest_district == "nsok nsomo"]$PR <- mean(holder[!is.na(holder)])

plot(EG.mainland)
plot(EG.shp[EG.shp$admin2 == "Aconibe",], add = TRUE)
holder <- values(mask(EG.mainland, EG.shp[EG.shp$admin2 == "Aconibe",]))
dt.off[dest_district == "akonibe"]$PR <- mean(holder[!is.na(holder)])

plot(EG.mainland)
plot(EG.shp[EG.shp$admin2 == "Añisok",], add = TRUE)
holder <- values(mask(EG.mainland, EG.shp[EG.shp$admin2 == "Añisok",]))
dt.off[dest_district == "anisok"]$PR <- mean(holder[!is.na(holder)])

plot(EG.mainland)
plot(EG.shp[EG.shp$admin2 == "Nsork",], add = TRUE)
holder <- values(mask(EG.mainland, EG.shp[EG.shp$admin2 == "Nsork",]))
dt.off[dest_district == "nsork"]$PR <- mean(holder[!is.na(holder)])

plot(EG.mainland)
plot(EG.shp[EG.shp$admin2 == "Cogo",], add = TRUE)
holder <- values(mask(EG.mainland, EG.shp[EG.shp$admin2 == "Cogo",]))
dt.off[dest_district == "cogo"]$PR <- mean(holder[!is.na(holder)])

plot(EG.mainland)
plot(EG.shp[EG.shp$admin2 == "Acurenam",], add = TRUE)
holder <- values(mask(EG.mainland, EG.shp[EG.shp$admin2 == "Acurenam",]))
dt.off[dest_district == "acurenam"]$PR <- mean(holder[!is.na(holder)])

plot(EG.mainland)
plot(EG.shp[EG.shp$admin2 == "Mbini",], add = TRUE)
holder <- values(mask(EG.mainland, EG.shp[EG.shp$admin2 == "Mbini",]))
dt.off[dest_district == "mbini"]$PR <- mean(holder[!is.na(holder)])

# Need just 2 more:
#dt.off[dest_district == "oyala"]$PR <- mean(holder[!is.na(holder)])
#dt.off[dest_district == "mbere"]$PR <- mean(holder[!is.na(holder)])


