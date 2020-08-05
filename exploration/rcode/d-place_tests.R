# Run in terminal: git clone https://github.com/D-PLACE/dplace-data

d <- read.csv(file = "dplace-data/datasets/Binford/data.csv")

str(d)

#construct full dataframe with all variables

a  <- which(d$var_id %in% c("B033", "B034", "B008", "B009", "B010", "B012", "B029", "B031")) # specificy variables of interest here. 
soc_id <- unique(d$soc_id[a])
df <- data.frame(soc_id,
                 money = rep(NA, length(soc_id)),
                 spec = rep(NA, length(soc_id)),
                 mob = rep(NA, length(soc_id)),
                 dens = rep(NA, length(soc_id)),
                 coop_size = rep(NA, length(soc_id)),
                 agg_size = rep(NA, length(soc_id)),
                 pol_comp = rep(NA, length(soc_id)),
                 class = rep(NA, length(soc_id))
                 )
for(i in 1:nrow(d)){
  for(j in 1:nrow(df)){
    if(d$soc_id[i] == df$soc_id[j]){
      if(d$var_id[i] == "B033") df$money[j] <-d$code[i] 
      if(d$var_id[i] == "B034") df$spec[j] <-d$code[i] 
      if(d$var_id[i] == "B008") df$dens[j] <-d$code[i] 
      if(d$var_id[i] == "B009") df$mob[j] <-d$code[i] 
      if(d$var_id[i] == "B010") df$coop_size[j] <-d$code[i] 
      if(d$var_id[i] == "B012") df$agg_size[j] <- d$code[i]
      if(d$var_id[i] == "B029") df$pol_comp[j] <- d$code[i]
      if(d$var_id[i] == "B031") df$class[j] <- d$code[i]
    }
  }
}


summary(lm(df$spec ~  df$money + df$pol_comp + df$dens + df$class))

spec <- d$code[d$var_id == "B034"]
density <- d$code[d$var_id == "B008"]
mobility <- d$code[d$var_id == "B009"]
sub_coop <- d$code[d$var_id == "B010"]

