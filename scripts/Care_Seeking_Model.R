####
# Raw data - Care Seeking Behavior and Treatment
# 
# We model the treatment cascade first as the probability of symptoms appearing followed by the probability of following through with seeking treatment.
#
# Both of these quantities are modeled using the symptoms and treatment reported in the MIS data.
#
# October 2, 2019
#
####

library(here)

# Load in the raw data, where people report their treatment seeking behavior and symptoms
raw.2015.2017.data <- fread(here("data/raw/2015-2017_survey_data/summaries.csv"))
care.fever <- raw.2015.2017.data[,.(areaId, n, pf, pffv, pftto)]

# Load in cleaned data, which has more accurate ad2 sorting
clean.data <- fread(here("data/clean/aggregated_2015_2018_travel_data.csv"))
# travel distance data set, we'll use this for the distance to malabo:
reg.travel.dist <- fread(here("data/clean/travel_dist_by_region.csv"))

# merge data before fitting
care.fever <- merge(clean.data[year == 2015, .(areaId, ad2, pop)], care.fever, by = "areaId")
care.fever <- merge(care.fever, reg.travel.dist[,.(areaId, dist.mal)])
care.fever[is.na(n)]$n <- 0
care.fever[is.na(pf)]$pf <- 0
care.fever[is.na(pffv)]$pffv <- 0
care.fever[is.na(pftto)]$pffto <- 0

sum(h$pffv)/sum(h$pf)

# probability of reporting fever, given that one had reported an infection
# perform binomial model fit, similar to how we fit the travel frequency model (Travel_frequency_model.R)
h <- glm(cbind(pffv, pf - pffv) ~ pop + dist.mal,
         data = care.fever, 
         family = binomial(link = logit))
# aic = 483.4
# the model that used ad2 as an indicator covariate made it so that Ureka was a very weird outlier, with probability of fever = 0
# taking distance from malabo into account, it appears that the farther away one is the more likely one will have symptoms
care.fever$prob.fever <- predict(h, data = care.fever, type = "response")

## some summary statistics
# hist(care.fever$prob.fever)
# summary(care.fever$prob.fever)
# mean = 0.1170
## compare to the results from
# sum(care.fever$pffv)/sum(care.fever$pf)
# 0.1116336

# probability of reporting seeking treatment, given that one had a reported symptoms
h <- glm(cbind(pftto, pffv - pftto) ~ pop + dist.mal,
         data = care.fever, 
         family = binomial(link = logit))

care.fever$prob.treatment <- predict(h, data = care.fever, type = "response")

## summary statistics
hist(care.fever$prob.treatment)
summary(care.fever$prob.treatment)
## mean: 0.5744
## compare to:
#sum(care.fever$pftto)/sum(care.fever$pffv)
## 0.602

fwrite(care.fever, here("data/clean/Care_Seeking_Model_estimates.csv"))


# Of course, the problem is actually that for this particular simulation application
# we will have to use a single number for fever + care seeking behavior.
# Yes, these numbers very likely vary across the island, but in this case
# the macro.pfsi model only accommodates a single value.