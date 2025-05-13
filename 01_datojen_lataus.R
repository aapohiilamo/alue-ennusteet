####
library(pxweb)
library(readxl)
library(dplyr)
library(readr)
library(sotkanet)
library(tidyverse)
#Ladataan datat suoraan TK:lta. 

#seuraavat datat
#1. metatiedot (kuntien ja maakuntien linkit)
#2. väestö
#3. syntyvyys
#4. elinajanodote (tämä sotkanetistä) (muista ajaa my;s koodi jolla jaetaan nama ikakohtaiseksi)
#5. muuttoliike
#6. tk:n ennuste


#####################################
#1. metatiedot (kuntien ja maakuntien linkit)
###################################
####maakunta koodit https://stat.fi/fi/luokitukset/tupa https://media.stat.fi/A7H6ohk0S8qafyCM4bfDaz/DICwBPn5Q8uifsM6dwow
maakunta_koodit <- read_excel(
  "meta_data/kunnat_ja_kuntapohjaiset_alueet_2025_suomeksi_ruotsiksi_englanniksi (1).xlsx", 
  range = "A2:AJ310"
) |>
  select(
    reg_code = 'Maakunnan koodi',
    name = Maakunta
  ) |>
  distinct() |>
  bind_rows(
    tibble(reg_code = "100", name = "Koko maa")
  )

maakunta_koodit_export <- maakunta_koodit |>
  mutate(
    country_code = 100,
    reg_code = as.numeric(reg_code),
    location_type = case_when(
      reg_code == 100 ~ 0,
      TRUE ~ 4
    )  ) |>
  select(country_code, reg_code, name, location_type) |>
  arrange(-reg_code)

write_delim(maakunta_koodit_export, "data/wafips.txt", delim = "\t", na = "", quote_escape = "none")

#####meta data 

meta <- read_excel("meta_data/kunnat_ja_kuntapohjaiset_alueet_2025_suomeksi_ruotsiksi_englanniksi (1).xlsx", 
                   range = "A2:AJ310") |> 
  distinct()


#####################################
#2.väestö
###################################
#population counts
pxweb_query_list <- 
  list("Alue"=c("SSS","KU020","KU005","KU009","KU010","KU016","KU018","KU019","KU035","KU043","KU046","KU047","KU049","KU050","KU051","KU052","KU060","KU061","KU062","KU065","KU069","KU071","KU072","KU074","KU075","KU076","KU077","KU078","KU079","KU081","KU082","KU086","KU111","KU090","KU091","KU097","KU098","KU102","KU103","KU105","KU106","KU108","KU109","KU139","KU140","KU142","KU143","KU145","KU146","KU153","KU148","KU149","KU151","KU152","KU165","KU167","KU169","KU170","KU171","KU172","KU176","KU177","KU178","KU179","KU181","KU182","KU186","KU202","KU204","KU205","KU208","KU211","KU213","KU214","KU216","KU217","KU218","KU224","KU226","KU230","KU231","KU232","KU233","KU235","KU236","KU239","KU240","KU320","KU241","KU322","KU244","KU245","KU249","KU250","KU256","KU257","KU260","KU261","KU263","KU265","KU271","KU272","KU273","KU275","KU276","KU280","KU284","KU285","KU286","KU287","KU288","KU290","KU291","KU295","KU297","KU300","KU301","KU304","KU305","KU312","KU316","KU317","KU318","KU398","KU399","KU400","KU407","KU402","KU403","KU405","KU408","KU410","KU416","KU417","KU418","KU420","KU421","KU422","KU423","KU425","KU426","KU444","KU430","KU433","KU434","KU435","KU436","KU438","KU440","KU441","KU475","KU478","KU480","KU481","KU483","KU484","KU489","KU491","KU494","KU495","KU498","KU499","KU500","KU503","KU504","KU505","KU508","KU507","KU529","KU531","KU535","KU536","KU538","KU541","KU543","KU545","KU560","KU561","KU562","KU563","KU564","KU309","KU576","KU577","KU578","KU445","KU580","KU581","KU599","KU583","KU854","KU584","KU592","KU593","KU595","KU598","KU601","KU604","KU607","KU608","KU609","KU611","KU638","KU614","KU615","KU616","KU619","KU620","KU623","KU624","KU625","KU626","KU630","KU631","KU635","KU636","KU678","KU710","KU680","KU681","KU683","KU684","KU686","KU687","KU689","KU691","KU694","KU697","KU698","KU700","KU702","KU704","KU707","KU729","KU732","KU734","KU736","KU790","KU738","KU739","KU740","KU742","KU743","KU746","KU747","KU748","KU791","KU749","KU751","KU753","KU755","KU758","KU759","KU761","KU762","KU765","KU766","KU768","KU771","KU777","KU778","KU781","KU783","KU831","KU832","KU833","KU834","KU837","KU844","KU845","KU846","KU848","KU849","KU850","KU851","KU853","KU857","KU858","KU859","KU886","KU887","KU889","KU890","KU892","KU893","KU895","KU785","KU905","KU908","KU092","KU915","KU918","KU921","KU922","KU924","KU925","KU927","KU931","KU934","KU935","KU936","KU941","KU946","KU976","KU977","KU980","KU981","KU989","KU992"),
       "Ikä"=c("SSS","000","001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042","043","044","045","046","047","048","049","050","051","052","053","054","055","056","057","058","059","060","061","062","063","064","065","066","067","068","069","070","071","072","073","074","075","076","077","078","079","080","081","082","083","084","085","086","087","088","089","090","091","092","093","094","095","096","097","098","099","100-"),
       "Sukupuoli"=c("SSS","1","2"),
       "Vuosi"=c("1972","1973","1974","1975","1976","1977","1978","1979","1980","1981","1982","1983","1984","1985","1986","1987","1988","1989","1990","1991","1992","1993","1994","1995","1996","1997","1998","1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017","2018","2019","2020","2021","2022","2023","2024"),
       "Tiedot"=c("vaesto"))
px_data <- 
  pxweb_get(url = "https://statfin.stat.fi/PXWeb/api/v1/fi/StatFin/vaerak/statfin_vaerak_pxt_11re.px",
            query = pxweb_query_list)

px_data_frame <- as.data.frame(px_data, column.name.type = "text", variable.value.type = "text")

vaesto <- px_data_frame|> 
  rename(Kunta=Alue,
         pop='Väestö 31.12.',
         age=Ikä) |>
  mutate(age=case_when(age=="100 -"~100, TRUE~as.numeric(age))) |>
  left_join(meta) |>
  group_by(Maakunta ,Vuosi,age, Sukupuoli) |>
  summarise(counts=sum(pop)) |>
  ungroup() |>
  rename(name=Maakunta) |>
  mutate(name=case_when(is.na(name)~"Koko maa", TRUE~name))|>
  left_join(maakunta_koodit) |>
  mutate(reg_code=as.numeric(reg_code))

popF <- vaesto |>
  filter(Sukupuoli == "Naiset") |>
  select(reg_code, name, age, year = Vuosi, value = counts) |>
  pivot_wider(names_from = year, values_from = value) |>
  arrange(reg_code, age) |> filter(!is.na(age)) |>
  arrange(-reg_code)

write_delim(popF, "data/popF.txt", delim = "\t", na = "", quote_escape = "none")

popM <- vaesto |>
  filter(Sukupuoli == "Miehet") |>
  select(reg_code, name, age, year = Vuosi, value = counts) |>
  pivot_wider(names_from = year, values_from = value) |>
  arrange(reg_code, age) |> filter(!is.na(age)) |>
  arrange(-reg_code)


write_delim(popM, "data/popM.txt", delim = "\t", na = "", quote_escape = "none")


#####################################
#3. syntyvyys
###################################

pxweb_query_list <- 
  list("Maakunta"=c("SSS","MA1","MK01","MK02","MK04","MK05","MK06","MK07","MK08","MK09","MK10","MK11","MK12","MK13","MK14","MK15","MK16","MK17","MK18","MK19","MA2","MK21"),
       "Tiedot"=c("tfr"),
       "Vuosi"=c("1990","1991","1992","1993","1994","1995","1996","1997","1998","1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017","2018","2019","2020","2021","2022","2023"))

# Download data 
px_data <- 
  pxweb_get(url = "https://statfin.stat.fi/PXWeb/api/v1/fi/StatFin/synt/statfin_synt_pxt_12du.px",
            query = pxweb_query_list)

# Convert to data.frame 
px_data_frame <- as.data.frame(px_data, column.name.type = "text", variable.value.type = "text")
tfr <- px_data_frame |>
  filter(!str_starts(Maakunta, "MA")) |>  # 1. Drop regions starting with "MA"
  mutate(
    Maakunta = str_remove(Maakunta, "^MK\\d+\\s*"),  # 2. Remove "MK01 ", "MK02 ", etc.
    country_code = 246,
    name = case_when(Maakunta=="KOKO MAA" ~ "Koko maa", TRUE~ Maakunta),  # 2. Remove "MK01 ", "MK02 ", etc.
  )  |>
  left_join(maakunta_koodit, by = c("name")) |>  # Join to get region codes
  select(country_code, reg_code, name, year = Vuosi, value = Kokonaishedelmällisyysluku) |>
  pivot_wider(names_from = year, values_from = value) |>
  mutate(reg_code=as.numeric(reg_code)) |>
  arrange(-reg_code)

# Write to .txt file
write_delim(tfr, "data/tfr.txt", delim = "\t", na = "", quote_escape = "none")


#####
#4. elinajanodote
#############
dat <- get_sotkanet(indicators = 4011)

region_codes <- dat |> select(reg_code = region, name = region.title) |>  distinct()

e0F <- dat|>
  filter(gender == "female") |>
  arrange(year) |>
  mutate(country_code = 246) |>
  select(country_code, name = region.title, year, value = primary.value) |>
  left_join(maakunta_koodit) |>
  select(country_code, reg_code, name, year, value) |>
  pivot_wider(names_from = year, values_from = value) |>
  mutate(reg_code=as.numeric(reg_code)) |>
  arrange(-reg_code)

# Write to tab-delimited txt file (without row names, and with raw column names)
write_delim(e0F, "data/e0F.txt", delim = "\t", na = "", quote_escape = "none")

e0M <-  dat|>
  filter(gender == "male") |>
  arrange(year) |>
  mutate(country_code = 246) |>
  select(country_code, name = region.title, year, value = primary.value) |>
  left_join(maakunta_koodit) |>
  select(country_code, reg_code, name, year, value) |>
  pivot_wider(names_from = year, values_from = value) |>
  mutate(reg_code=as.numeric(reg_code)) |>
  arrange(-reg_code)

write_delim(e0M, "data/e0M.txt", delim = "\t", na = "", quote_escape = "none")


####### MUUTTOLIIKE
pxweb_query_list <- 
  list("Vuosi"=c("1990","1991","1992","1993","1994","1995","1996","1997","1998","1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017","2018","2019","2020","2021","2022","2023"),
       "Alue"=c("SSS","KU020","KU005","KU009","KU010","KU016","KU018","KU019","KU035","KU043","KU046","KU047","KU049","KU050","KU051","KU052","KU060","KU061","KU062","KU065","KU069","KU071","KU072","KU074","KU075","KU076","KU077","KU078","KU079","KU081","KU082","KU086","KU111","KU090","KU091","KU097","KU098","KU102","KU103","KU105","KU106","KU108","KU109","KU139","KU140","KU142","KU143","KU145","KU146","KU153","KU148","KU149","KU151","KU152","KU165","KU167","KU169","KU170","KU171","KU172","KU176","KU177","KU178","KU179","KU181","KU182","KU186","KU202","KU204","KU205","KU208","KU211","KU213","KU214","KU216","KU217","KU218","KU224","KU226","KU230","KU231","KU232","KU233","KU235","KU236","KU239","KU240","KU320","KU241","KU322","KU244","KU245","KU249","KU250","KU256","KU257","KU260","KU261","KU263","KU265","KU271","KU272","KU273","KU275","KU276","KU280","KU284","KU285","KU286","KU287","KU288","KU290","KU291","KU295","KU297","KU300","KU301","KU304","KU305","KU312","KU316","KU317","KU318","KU398","KU399","KU400","KU407","KU402","KU403","KU405","KU408","KU410","KU416","KU417","KU418","KU420","KU421","KU422","KU423","KU425","KU426","KU444","KU430","KU433","KU434","KU435","KU436","KU438","KU440","KU441","KU475","KU478","KU480","KU481","KU483","KU484","KU489","KU491","KU494","KU495","KU498","KU499","KU500","KU503","KU504","KU505","KU508","KU507","KU529","KU531","KU535","KU536","KU538","KU541","KU543","KU545","KU560","KU561","KU562","KU563","KU564","KU309","KU576","KU577","KU578","KU445","KU580","KU581","KU599","KU583","KU854","KU584","KU588","KU592","KU593","KU595","KU598","KU601","KU604","KU607","KU608","KU609","KU611","KU638","KU614","KU615","KU616","KU619","KU620","KU623","KU624","KU625","KU626","KU630","KU631","KU635","KU636","KU678","KU710","KU680","KU681","KU683","KU684","KU686","KU687","KU689","KU691","KU694","KU697","KU698","KU700","KU702","KU704","KU707","KU729","KU732","KU734","KU736","KU790","KU738","KU739","KU740","KU742","KU743","KU746","KU747","KU748","KU791","KU749","KU751","KU753","KU755","KU758","KU759","KU761","KU762","KU765","KU766","KU768","KU771","KU777","KU778","KU781","KU783","KU831","KU832","KU833","KU834","KU837","KU844","KU845","KU846","KU848","KU849","KU850","KU851","KU853","KU857","KU858","KU859","KU886","KU887","KU889","KU890","KU892","KU893","KU895","KU785","KU905","KU908","KU092","KU915","KU918","KU921","KU922","KU924","KU925","KU927","KU931","KU934","KU935","KU936","KU941","KU946","KU976","KU977","KU980","KU981","KU989","KU992"),
       "Tiedot"=c("koknetmuutto"))

# Download data 
px_data <- 
  pxweb_get(url = "https://statfin.stat.fi/PXWeb/api/v1/fi/StatFin/muutl/statfin_muutl_pxt_11ae.px",
            query = pxweb_query_list)

# Convert to data.frame 
mig_data_raw <- as.data.frame(px_data, column.name.type = "text", variable.value.type = "text")
muutto<- mig_data_raw |> rename(Kunta=Alue) |>
  left_join(meta) |>
  group_by(Maakunta ,Vuosi) |>
  summarise(counts=sum(Kokonaisnettomuutto)) |>
  ungroup() |>
  rename(name=Maakunta) |>
  left_join(maakunta_koodit) |>
  mutate(reg_code=as.numeric(reg_code))

###crude immigration rate
pop <- vaesto |>
  filter(Sukupuoli == "Yhteensä") |>
  select(reg_code, name, age, Vuosi, value = counts) |>
  arrange(reg_code, age) |> filter(is.na(age))

muutto_crude<- mig_data_raw |> rename(Kunta=Alue) |>
  left_join(meta) |>
  mutate(Maakunta=case_when(is.na(Maakunta)~"Koko maa", TRUE~Maakunta))|>
  group_by(Maakunta ,Vuosi) |>
  summarise(counts=sum(Kokonaisnettomuutto)) |>
  ungroup() |>
  rename(name=Maakunta) |>
  left_join(maakunta_koodit) |>
  mutate(reg_code=as.numeric(reg_code)) |>
  left_join(pop)|>
  mutate(crude_migration=counts/value)

###väestön ennakkotiedot crude immigration

pxweb_query_list <- 
  list("Alue"=c("SSS","KU020","KU005","KU009","KU010","KU016","KU018","KU019","KU035","KU043","KU046","KU047","KU049","KU050","KU051","KU052","KU060","KU061","KU062","KU065","KU069","KU071","KU072","KU074","KU075","KU076","KU077","KU078","KU079","KU081","KU082","KU086","KU111","KU090","KU091","KU097","KU098","KU102","KU103","KU105","KU106","KU108","KU109","KU139","KU140","KU142","KU143","KU145","KU146","KU153","KU148","KU149","KU151","KU152","KU165","KU167","KU169","KU170","KU171","KU172","KU176","KU177","KU178","KU179","KU181","KU182","KU186","KU202","KU204","KU205","KU208","KU211","KU213","KU214","KU216","KU217","KU218","KU224","KU226","KU230","KU231","KU232","KU233","KU235","KU236","KU239","KU240","KU320","KU241","KU322","KU244","KU245","KU249","KU250","KU256","KU257","KU260","KU261","KU263","KU265","KU271","KU272","KU273","KU275","KU276","KU280","KU284","KU285","KU286","KU287","KU288","KU290","KU291","KU295","KU297","KU300","KU301","KU304","KU305","KU312","KU316","KU317","KU318","KU398","KU399","KU400","KU407","KU402","KU403","KU405","KU408","KU410","KU416","KU417","KU418","KU420","KU421","KU422","KU423","KU425","KU426","KU444","KU430","KU433","KU434","KU435","KU436","KU438","KU440","KU441","KU475","KU478","KU480","KU481","KU483","KU484","KU489","KU491","KU494","KU495","KU498","KU499","KU500","KU503","KU504","KU505","KU508","KU507","KU529","KU531","KU535","KU536","KU538","KU541","KU543","KU545","KU560","KU561","KU562","KU563","KU564","KU309","KU576","KU577","KU578","KU445","KU580","KU581","KU599","KU583","KU854","KU584","KU592","KU593","KU595","KU598","KU601","KU604","KU607","KU608","KU609","KU611","KU638","KU614","KU615","KU616","KU619","KU620","KU623","KU624","KU625","KU626","KU630","KU631","KU635","KU636","KU678","KU710","KU680","KU681","KU683","KU684","KU686","KU687","KU689","KU691","KU694","KU697","KU698","KU700","KU702","KU704","KU707","KU729","KU732","KU734","KU736","KU790","KU738","KU739","KU740","KU742","KU743","KU746","KU747","KU748","KU791","KU749","KU751","KU753","KU755","KU758","KU759","KU761","KU762","KU765","KU766","KU768","KU771","KU777","KU778","KU781","KU783","KU831","KU832","KU833","KU834","KU837","KU844","KU845","KU846","KU848","KU849","KU850","KU851","KU853","KU857","KU858","KU859","KU886","KU887","KU889","KU890","KU892","KU893","KU895","KU785","KU905","KU908","KU092","KU915","KU918","KU921","KU922","KU924","KU925","KU927","KU931","KU934","KU935","KU936","KU941","KU946","KU976","KU977","KU980","KU981","KU989","KU992"),
       "Vuosineljännes"=c("2024Q1","2024Q2","2024Q3","2024Q4"),
       "Tiedot"=c("koknetmuutto","vaesto"))

# Download data 
px_data <- 
  pxweb_get(url = "https://statfin.stat.fi/PXWeb/api/v1/fi/StatFin/vamuu/statfin_vamuu_pxt_11lk.px",
            query = pxweb_query_list)

# Convert to data.frame 
px_data_frame <- as.data.frame(px_data, column.name.type = "text", variable.value.type = "text")

Väkiluku2024 <- px_data_frame |> rename(Kunta=Alue) |>
  filter(Vuosineljännes=="2024Q4*") |>
  left_join(meta) |>
  mutate(Maakunta=case_when(is.na(Maakunta)~"Koko maa", TRUE~Maakunta))|>
  group_by(Maakunta) |>
  summarise(pop=sum(Väkiluku),)|>
  rename(name=Maakunta) |>
  left_join(maakunta_koodit) |>
  mutate(reg_code=as.numeric(reg_code))

muutto_crude2024<- px_data_frame |> rename(Kunta=Alue) |>
  left_join(meta) |>
  mutate(Maakunta=case_when(is.na(Maakunta)~"Koko maa", TRUE~Maakunta))|>
  group_by(Maakunta) |>
  summarise(counts=sum(Kokonaisnettomuutto),
            ) |>
  ungroup() |>
  rename(name=Maakunta) |>
  left_join(maakunta_koodit) |>
  mutate(reg_code=as.numeric(reg_code)) |>
  left_join(Väkiluku2024)|>
  mutate(crude_migration=counts/pop) |> 
  mutate(Vuosi="2024")
 
 mig_rates <- muutto_crude |> 
  bind_rows(muutto_crude2024)  |>
    select(country_code=reg_code, name,  year = Vuosi, value = crude_migration) |>
    pivot_wider(names_from = year, values_from = value) |>
    arrange(country_code) |>
    mutate(country_code=as.numeric(country_code)) |>
    arrange(-country_code)
  
  write_delim(mig_rates, "data/mig_rates.txt", delim = "\t", na = "", quote_escape = "none")
 
  ####### MUUTTOLIIKE  counts

  muutto<- muutto   |>
    mutate(reg_code=case_when(is.na(reg_code)~100, TRUE~reg_code))   |>
    mutate(name=case_when(is.na(name)~"Koko maa", TRUE~name))   |>
    bind_rows(muutto_crude2024)
  
  mig <- muutto |>
    mutate(
      country_code = 100,
      location_type = 4,  # or whatever constant value makes sense
      Vuosi = as.character(Vuosi)  # ensure years are character before pivot
    ) |>
    select(name, country_code, reg_code, location_type, year = Vuosi, value = counts) |>
    pivot_wider(names_from = year, values_from = value) |>
    arrange(reg_code) |>
    mutate(reg_code=as.numeric(reg_code)) 
  
   # Write to file
  write_delim(mig, "data/mig_counts.txt", delim = "\t", na = "", quote_escape = "none")
  
  
###tk;n ennuste

pxweb_query_list <- 
  list("Alue"=c("SSS","KU020","KU005","KU009","KU010","KU016","KU018","KU019","KU035","KU043","KU046","KU047","KU049","KU050","KU051","KU052","KU060","KU061","KU062","KU065","KU069","KU071","KU072","KU074","KU075","KU076","KU077","KU078","KU079","KU081","KU082","KU086","KU111","KU090","KU091","KU097","KU098","KU102","KU103","KU105","KU106","KU108","KU109","KU139","KU140","KU142","KU143","KU145","KU146","KU153","KU148","KU149","KU151","KU152","KU165","KU167","KU169","KU170","KU171","KU172","KU176","KU177","KU178","KU179","KU181","KU182","KU186","KU202","KU204","KU205","KU208","KU211","KU213","KU214","KU216","KU217","KU218","KU224","KU226","KU230","KU231","KU232","KU233","KU235","KU236","KU239","KU240","KU320","KU241","KU322","KU244","KU245","KU249","KU250","KU256","KU257","KU260","KU261","KU263","KU265","KU271","KU272","KU273","KU275","KU276","KU280","KU284","KU285","KU286","KU287","KU288","KU290","KU291","KU295","KU297","KU300","KU301","KU304","KU305","KU312","KU316","KU317","KU318","KU398","KU399","KU400","KU407","KU402","KU403","KU405","KU408","KU410","KU416","KU417","KU418","KU420","KU421","KU422","KU423","KU425","KU426","KU444","KU430","KU433","KU434","KU435","KU436","KU438","KU440","KU441","KU475","KU478","KU480","KU481","KU483","KU484","KU489","KU491","KU494","KU495","KU498","KU499","KU500","KU503","KU504","KU505","KU508","KU507","KU529","KU531","KU535","KU536","KU538","KU541","KU543","KU545","KU560","KU561","KU562","KU563","KU564","KU309","KU576","KU577","KU578","KU445","KU580","KU581","KU599","KU583","KU854","KU584","KU588","KU592","KU593","KU595","KU598","KU601","KU604","KU607","KU608","KU609","KU611","KU638","KU614","KU615","KU616","KU619","KU620","KU623","KU624","KU625","KU626","KU630","KU631","KU635","KU636","KU678","KU710","KU680","KU681","KU683","KU684","KU686","KU687","KU689","KU691","KU694","KU697","KU698","KU700","KU702","KU704","KU707","KU729","KU732","KU734","KU736","KU790","KU738","KU739","KU740","KU742","KU743","KU746","KU747","KU748","KU791","KU749","KU751","KU753","KU755","KU758","KU759","KU761","KU762","KU765","KU766","KU768","KU771","KU777","KU778","KU781","KU783","KU831","KU832","KU833","KU834","KU837","KU844","KU845","KU846","KU848","KU849","KU850","KU851","KU853","KU857","KU858","KU859","KU886","KU887","KU889","KU890","KU892","KU893","KU895","KU785","KU905","KU908","KU092","KU915","KU918","KU921","KU922","KU924","KU925","KU927","KU931","KU934","KU935","KU936","KU941","KU946","KU976","KU977","KU980","KU981","KU989","KU992"),
       "Vuosi"=c("2024","2025","2026","2027","2028","2029","2030","2031","2032","2033","2034","2035","2036","2037","2038","2039","2040","2041","2042","2043","2044","2045"),
       "Sukupuoli"=c("SSS","1","2"),
       "Ikä"=c("SSS","000","001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042","043","044","045","046","047","048","049","050","051","052","053","054","055","056","057","058","059","060","061","062","063","064","065","066","067","068","069","070","071","072","073","074","075","076","077","078","079","080","081","082","083","084","085","086","087","088","089","090","091","092","093","094","095","096","097","098","099","100-"),
       "Tiedot"=c("vaesto_e24"))

px_data <- 
  pxweb_get(url = "https://statfin.stat.fi/PXWeb/api/v1/fi/StatFin/vaenn/statfin_vaenn_pxt_14wx.px",
            query = pxweb_query_list)

ennuste_data <- as.data.frame(px_data, column.name.type = "text", variable.value.type = "text")

ennuste_data_ <- ennuste_data |>
  rename(Kunta=Alue) |>
  mutate(         count=as.numeric(`Väestö 31.12. (ennuste 2024)`)) |>
  left_join(meta) |>
  filter(Sukupuoli=="Yhteensä", Ikä =="Yhteensä") |>
  group_by(Vuosi, Maakunta) |>
  summarise(median=sum(count)) |>
  ungroup() |>
  select(year=Vuosi, name=Maakunta,median) |>
  mutate(ennuste="Tilastokeskus")

saveRDS(ennuste_data_, file = "data/tk_ennuste.rds")

