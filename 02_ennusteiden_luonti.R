getwd()
library(bayesPop)
library(wpp2024)
library(scales)
library(geofacet)
library(tidyverse)
library(geofi)
library(writexl)
library(bayesMig)
set.seed(123)
#tehdään ennusteet vaihteittan. 
#1. syntyvyys
#2. elinajanodote
#3. muuttoliike
###########
#fertility
###########
my.tfr.file <- "data/tfr.txt"
head(read.delim(my.tfr.file, check.names = FALSE))
nat.tfr.dir <- "data/wpp2024_projections/TFR1unc/sim20241101"
nat.tfr.dir <- "data/wpp2024_projections/TFR1unc/sim20241101"

subnat.tfr.dir <- "results/tfr"
tfr.preds <- tfr.predict.subnat(246, my.tfr.file = my.tfr.file,
                                sim.dir = nat.tfr.dir, output.dir = subnat.tfr.dir,
                                annual = TRUE, start.year = 2024, end.year = 2040,
                                nr.traj = 10000)

tfr.preds <- get.regtfr.prediction(subnat.tfr.dir)
tfr.subnat <-tfr.preds[["246"]]
nat.tfr.pred <- get.tfr.prediction(nat.tfr.dir)

###########
#elinajanodote
###########
my.e0M.file <- "data/e0M.txt"
my.e0F.file <- "data/e0F.txt"
head(read.delim(my.e0F.file, check.names = FALSE))
nat.e0.dir <- "data/wpp2024_projections/e01/sim20241101"
subnat.e0.dir <- "results/e0"

e0.preds <- e0.predict.subnat(246, my.e0.file = my.e0F.file,
                              sim.dir = nat.e0.dir, output.dir = subnat.e0.dir,
                              annual = TRUE, start.year = 2024, end.year = 2040,
                              predict.jmale = TRUE, my.e0M.file = my.e0M.file,
                              nr.traj = 10000)

e0.subnat <- get.rege0.prediction(subnat.e0.dir, 246)
e0.trajectories.plot(e0.subnat, "Uusimaa", both.sexes = TRUE)
e0.trajectories.plot(e0.subnat, 100, both.sexes = FALSE)

#############
# Migration
#############
my.mig.file <- "data/mig_rates.txt"
head(read.delim(my.mig.file, check.names = FALSE))

subnat.mig.dir <- "results/mig"
mc<- run.mig.mcmc(nr.chains = 3, iter = 5000,
                  output.dir = subnat.mig.dir,
                  present.year = 2024, start.year = 1990,
                  my.mig.file = my.mig.file, annual = TRUE,
                  replace.output = TRUE, verbose.iter = 1000)

mig.subnat <- mig.predict(sim.dir = subnat.mig.dir, nr.traj = 10000,
                          burnin = 1000, end.year = 2050, save.as.ascii = 1000)

mig.subnat <- get.mig.prediction(sim.dir = subnat.mig.dir)
mig.trajectories.plot(mig.subnat, "Uusimaa", nr.traj = 30)
mig.trajectories.plot(mig.subnat, "Koko maa", nr.traj = 30)

#############################
# Population projections
###############################

location.file <- "data/wafips.txt"
popM0.file<- "data/popM.txt"
popF0.file<- "data/popF.txt"
mxM.file <- "data/mxM.txt"
mxF.file<- "data/mxF.txt"
mig.file <- "data/mig_counts.txt"
mig.traj.file<- file.path(subnat.mig.dir, 
                          "predictions/ascii_trajectories.csv")
subnat.tfr.results<- file.path(subnat.tfr.dir, "subnat/c246")
subnat.e0.results <- file.path(subnat.e0.dir,"subnat_ar1/c246")
subnat.pop.dir <- "results/pop"
pop.subnat <- pop.predict.subnat(output.dir = subnat.pop.dir,
                                 locations = location.file, default.country = 246,
                                 verbose = TRUE, annual = TRUE, wpp.year = 2024,
                                 present.year = 2023, end.year = 2040,
                                 nr.traj = 10000, replace.output = TRUE,
                                 inputs = list(
                                   popM = popM0.file, popF = popF0.file,
                                   mxM = mxM.file,mxF = mxF.file,
                                   mig = mig.file, migtraj = mig.traj.file,
                                   tfr.sim.dir = subnat.tfr.results,
                                   e0F.sim.dir = subnat.e0.results, 
                                   e0M.sim.dir = "joint_"  ),
                                 mig.age.method = "rc",
                                 mig.is.rate = c(FALSE, TRUE),
                                 keep.vital.events = TRUE,
                                 pasfr.ignore.phase2 = TRUE
                                 )

pop.aggr <- pop.aggregate.subnat(pop.subnat, regions = 100,
                                locations = location.file)

pop.subnat <- get.pop.prediction(subnat.pop.dir)
pop.aggr <- get.pop.aggregation(sim.dir = subnat.pop.dir)

#############################
# results
###############################

# Tallennukset Exceliin
pop.trajectories.table(pop.subnat,  expression = "E100{0}") |> as.data.frame() |> write_xlsx("E0.xlsx")
pop.trajectories.table(pop.subnat, expression ="P100") |> as.data.frame() |> write_xlsx("Pop.xlsx")
pop.trajectories.table(pop.subnat,expression = "F100") |> as.data.frame() |> write_xlsx("trf.xlsx")
pop.trajectories.table(pop.subnat, expression ="G100") |> as.data.frame() |> write_xlsx("migrants.xlsx")

# Aluekohtaiset ennusteet
data <- map_dfr(setdiff(1:21, c(3, 20)), ~ {
  pop.trajectories.table(pop.subnat, expression =paste0("P", .x), pi = c(80, 95)) |> 
    as.data.frame() |> 
    rownames_to_column("year") |> 
    mutate(reg_code = .x)
}) |> 
  filter(year > 2010) |> 
  left_join(maakunta_koodit_export) |> 
  mutate(ennuste = "bayesPop")

combined <- data |> 
  bind_rows(readRDS("data/tk_ennuste.rds")) |> 
  filter(year < 2041, !is.na(name))

write_xlsx(combined, "alueelliset väestöennusteet koko väestö.xlsx")

# Kuva koko väestöstä
pop_all<-ggplot(combined, aes(as.numeric(year), median, color = ennuste, fill = ennuste)) +
  geom_ribbon(aes(ymin = `0.025`, ymax = `0.975`), alpha = 0.2) +
  geom_ribbon(aes(ymin = `0.1`, ymax = `0.9`), alpha = 0.4) +
  geom_line(size = 1) +
  geom_vline(xintercept = 2024, linetype = "dashed") +
  facet_geo(~name, grid = grid_maakunta, scale = "free_y") +
  scale_y_continuous(labels = label_number(scale = 1e-3, suffix = "t")) +
  theme_minimal()  +
  labs(
    y = "",
    x = "",
    title = ""
  )
ggsave("figs/all_pop.svg",pop_all, width = 8, height = 10)

# Ikäprofiili Pohjois-Savolle
ika_profile <- map_dfr(0:100, ~ {
  pop.trajectories.table(pop.subnat, expression =paste0("P11[", .x, "]"), pi = c(80,95)) |> 
    as.data.frame() |> 
    rownames_to_column("year") |> 
    mutate(age = .x, reg_code = 11)
}) |> 
  left_join(maakunta_koodit_export) |> 
  filter(year %in% c(2023, 2040), name == "Pohjois-Savo")

age_pattern<-ggplot(ika_profile, aes(age, median, color = year, fill = year)) +
  geom_ribbon(aes(ymin = `0.025`, ymax = `0.975`), alpha = 0.2) +
  geom_ribbon(aes(ymin = `0.1`, ymax = `0.9`), alpha = 0.4) +
  geom_line(size = 1) +
  scale_y_continuous(labels = label_number(scale = 1e-3, suffix = "t")) +
  theme_minimal()  +
  labs(
    y = "",
    x = "",
    title = ""
  )
ggsave("figs/age_pattern.svg",age_pattern, width = 8, height = 10)

# 80+ väestö
data80 <- map_dfr(setdiff(1:21, c(3,20)), ~ {
  pop.trajectories.table(pop.subnat, expression =paste0("P", .x, "[80:130]"), pi = c(80,95)) |> 
    as.data.frame() |> 
    rownames_to_column("year") |> 
    mutate(reg_code = .x)
}) |> filter(year > 2010) |> left_join(maakunta_koodit_export) |> mutate(ennuste="bayesPop")

write_xlsx(data80, "alueelliset väestöennusteet 80+ väestö.xlsx")

older<-ggplot(data80, aes(as.numeric(year), median, color = ennuste, fill = ennuste)) +
  geom_ribbon(aes(ymin = `0.025`, ymax = `0.975`), alpha = 0.2) +
  geom_ribbon(aes(ymin = `0.1`, ymax = `0.9`), alpha = 0.4) +
  geom_line(size = 1) +
  geom_vline(xintercept = 2024, linetype = "dashed") +
  facet_geo(~name, grid = grid_maakunta, scale = "free_y") +
  scale_y_continuous(labels = label_number(scale = 1e-3, suffix = "t")) +
  theme_minimal() + theme(legend.position = "none") +
  labs(
    y = "",
    x = "",
    title = ""
  )
ggsave("figs/older.svg",older, width = 8, height = 10)

# Huoltosuhde
data_support <- map_dfr(setdiff(1:21, c(3,20)), ~ {
  pop.trajectories.table(pop.subnat, expression =paste0("(P", .x, "[0:18]+P", .x, "[65:130])/P", .x, "[19:64]"), pi=c(80,95)) |> 
    as.data.frame() |> 
    rownames_to_column("year") |> 
    mutate(reg_code=.x)
}) |> filter(year > 2010) |> left_join(maakunta_koodit_export) |> mutate(ennuste="bayesPop")

write_xlsx(data_support, "alueelliset väestöennusteet huoltosuhde.xlsx")

huolto<-ggplot(data_support, aes(as.numeric(year), median, fill=ennuste, color=ennuste)) +
  geom_ribbon(aes(ymin=`0.025`, ymax=`0.975`), alpha=0.2) +
  geom_ribbon(aes(ymin=`0.1`, ymax=`0.9`), alpha=0.4) +
  geom_line(size=1) +
  facet_geo(~name, grid=grid_maakunta) +
  geom_vline(xintercept=2024, linetype="dashed") +
  theme_minimal()  
  ggsave("figs/support.svg", huolto, width=8, height=10)

# Miksi Lapin väestö kasvaa vuonna 2040
test <- get.pop.ex("P19", pop.subnat)
sum(test["2040", ] > test["2024", ]) / ncol(test)

# Hedelmällisyys ja muuttoliike
observed <- get.pop.ex("G100", pop.subnat, as.dt=TRUE, observed=TRUE) |> rename(Muuttaneita=indicator) |> 
  left_join(get.pop.ex("P100", pop.subnat, as.dt=TRUE, observed=TRUE)) |> rename(pop=indicator) |> 
  left_join(get.pop.ex("F100", pop.subnat, as.dt=TRUE, observed=TRUE)) |> rename(Syntyvyys=indicator)

base <- get.pop.ex("G100", pop.subnat, as.dt=TRUE) |> rename(Muuttaneita=indicator) |> 
  left_join(get.pop.ex("P100", pop.subnat, as.dt=TRUE)) |> rename(pop=indicator) |> 
  left_join(get.pop.ex("F100", pop.subnat, as.dt=TRUE)) |> rename(Syntyvyys=indicator)

migration <- base |> group_by(trajectory) |> filter(any(year == "2040")) |> ungroup()

pop_change_df <- migration |> filter(year %in% c("2023", "2040")) |> 
  select(trajectory, year, pop) |> pivot_wider(names_from=year, values_from=pop, names_prefix="pop_") |> 
  mutate(pop_change=if_else(pop_2040 < pop_2023, "Väkiluku laskee", "Väkiluku nousee")) |> 
  select(trajectory, pop_change)

df_long <- migration |> left_join(pop_change_df, by="trajectory") |> 
  pivot_longer(c(Muuttaneita, Syntyvyys), names_to="type", values_to="value")

median_df <- df_long |> group_by(year, type, pop_change) |> summarize(value=median(value, na.rm=TRUE), .groups="drop") |> mutate(trajectory="median")

sampled_trajectories <- df_long |> filter(trajectory != "median") |> group_by(type, pop_change) |> distinct(trajectory) |> slice_sample(n=10)

df_sampled <- df_long |> semi_join(sampled_trajectories, by=c("trajectory", "type", "pop_change")) |> mutate(trajectory=as.character(trajectory))

df_obs_long <- observed |> pivot_longer(c(Muuttaneita, Syntyvyys), names_to="type", values_to="value") |> mutate(trajectory="observed") |> crossing(pop_change=c("Väkiluku laskee", "Väkiluku nousee")) |> filter(year>2010)

df_plot <- bind_rows(df_sampled, median_df, df_obs_long) |> filter(year>2010)

ggplot() +
  geom_line(data=df_sampled, aes(as.numeric(year), value, group=trajectory), color="black", alpha=0.5) +
  geom_line(data=median_df, aes(as.numeric(year), value), color="red", size=1.2) +
  geom_line(data=df_obs_long, aes(as.numeric(year), value), color="blue", size=1) +
  facet_grid(type ~ pop_change, scales="free_y") +
  theme_minimal() +
  labs(title="Hedelmällisyys ja muuttoliike, Lapin väkiluvun muutosarviot", x="", y="")

# 80+ väestön osuus
data_older_share <- map_dfr(setdiff(1:21, c(3, 20)), ~ {
  pop.trajectories.table(pop.subnat, expression =paste0("P", .x, "[80:130]/P", .x), nr.traj=0, pi=c(80,95)) |> 
    as.data.frame() |> rownames_to_column("year") |> mutate(reg_code=.x)
}) |> filter(year>2010) |> left_join(maakunta_koodit_export) |> mutate(ennuste="bayesPop")

write_xlsx(data_older_share, "alueelliset väestöennusteet 80+ väestön osuus.xlsx")

older_share<- ggplot(data_older_share, aes(as.numeric(year), median, color=ennuste)) +
  geom_ribbon(aes(ymin=`0.025`, ymax=`0.975`), alpha=0.2) +
  geom_ribbon(aes(ymin=`0.1`, ymax=`0.9`), alpha=0.4) +
  geom_line(size=1) +
  facet_geo(~name, grid=grid_maakunta) +
  geom_vline(xintercept=2024, linetype="dashed") +
  theme_minimal() +
  theme(legend.position="bottom") 

ggsave("older_share.svg",older_share, width=8, height=10)


