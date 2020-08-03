########################################################
############# Intial tests #############################

spec<- apply(df_income[,1:13], 1, theil)
plot(spec ~ df_income$market_int)
spec_pp <- spec/demographics$hh_adult_eq
plot((spec_pp ~ df_income$market_int))
plot((spec ~ demographics$hh_adult_eq))
plot((spec ~ demographics$hh_over_18))
plot((df_income$market_int ~ demographics$hh_over_18))
plot((spec ~ demographics$hh_over_18 ))
plot((df_income$market_int ~ demographics$females_under18))
plot((spec ~ demographics$females_under18))
plot((df_income$market_int ~ demographics$hhh_sex))
plot((spec ~ demographics$hhh_sex))
plot((spec ~ demographics$hhh_sex))
plot((df_income$market_int ~ demographics$hhh_sex))
plot((df_income$market_int ~ demographics$female_over18))
plot((df_income$market_int ~ demographics$male_over18))
plot((spec ~ demographics$female_over18))
plot(spec ~ demographics$male_over18)
plot(spec ~ as.numeric(demographics$hhh_age))
plot(log(df_income$total_inc) ~ spec)
plot(log(df_income$total_inc) ~ df_income$market_int)
plot(log(wealth$total_wealth) ~ df_income$market_int)
plot(log(wealth$total_wealth) ~ spec)
plot(spec ~ demographics$hhh_edu)
plot(df_income$market_int ~ demographics$hhh_edu)
plot(df_income$market_int ~ demographics$father_education)
plot(spec ~ demographics$father_education)
plot(spec ~ demographics$mother_education)
plot(df_income$market_int ~ demographics$mother_education)
plot(df_income$market_int ~ demographics$hh_u_18_a_5)
plot(spec ~ demographics$hh_u_18_a_5)



plot((spec_pp ~ demographics$hh_adult_eq))
plot((spec_pp ~ demographics$hh_over_18))
plot((df_income$market_int ~ demographics$hh_over_18))
plot((spec_pp ~ demographics$hh_over_18 ))
plot((df_income$market_int ~ demographics$females_under18))
plot((spec_pp ~ demographics$females_under18))
plot((df_income$market_int ~ demographics$hhh_sex))
plot((spec_pp ~ demographics$hhh_sex))
plot((spec_pp ~ demographics$hhh_sex))
plot((df_income$market_int ~ demographics$hhh_sex))
plot((df_income$market_int ~ demographics$female_over18))
plot((df_income$market_int ~ demographics$male_over18))
plot((spec_pp ~ demographics$female_over18))
plot(spec_pp ~ demographics$male_over18)
plot(spec_pp ~ as.numeric(demographics$hhh_age))
plot(log(df_income$total_inc) ~ spec_pp)
plot(log(df_income$total_inc) ~ df_income$market_int)
plot(log(wealth$total_wealth) ~ df_income$market_int)
plot(log(wealth$total_wealth) ~ spec_pp)
plot(spec_pp ~ demographics$hhh_edu)
plot(df_income$market_int ~ demographics$hhh_edu)
plot(df_income$market_int ~ demographics$father_education)
plot(spec_pp ~ demographics$father_education)
plot(spec_pp ~ demographics$mother_education)
plot(df_income$market_int ~ demographics$mother_education)
plot(df_income$market_int ~ demographics$hh_u_18_a_5)
plot(spec_pp ~ demographics$hh_u_18_a_5)



#THE GOAL IS TO QUANTIFY HOW MANY DIFFRENT INCOME SOURCES ARE IN THE DATA SET 
# You gotta add up each tiny source of income and obtain a theil measure. 

income_list <- select(income, ends_with("quant"))
prod_list <- data.frame(select(farming , ends_with("prod")))


