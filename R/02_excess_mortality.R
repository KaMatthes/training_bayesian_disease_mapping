# load your packages

rm(list=ls())
source("R/00_setup.R")


# define plot parameters
axis_text_size <- 20
plot_title_size <- 25
lwd_size <- 1.5
point_size <- 3
pd <- position_dodge(width=0.4) # space between the line in the coefficient plot
fatten_size <- 8 # make the lines thicker

# load your data
dtd <- read.csv("data/deaths_copen.csv", sep = ";") 
  
# load your shapefiles
dts <- st_read("data/shapefiles/Polygonbasis_183_eli.shp") %>%
  rename (District = BEZNR)

# create your neighboorhood list, QUEEN is the default
nbk <- poly2nb(dts, dts$District)  # build neighbours list
nb2INLA("data/Districts", nbk) # build spatial neighbours for INLA

# link data and spatial data

dt <- poly2nb(dts, dts$District) %>%
  attr("region.id") %>%
  as.data.frame() %>%
  rename(District = ".") %>%
  mutate(Region = 1:183)  %>%
  left_join(dtd) %>%
  mutate(
    mx = deaths/pop,
    death_mod=ifelse(Year ==1918, NA, deaths) # to exclude 1918 in the baseline, 1918 will be predicted given 1913 - 1917
    ) 

#####################
# Statistical model #
#####################

# Penalized Complexity prior

hyper.bym <- list(
  # Overall variability
  # P(SD > 1) = 0.01
  # Only a 1% prior probability that the standard deviation of the spatial effect is greater than 1.
  theta1 = list(
    prior = "pc.prec",
    param = c(1, 0.01)
  ),
  
  # Spatial mixing
  # P(phi < 0.5) = 0.5
  #Only a 1% prior probability that the standard deviation of the spatial effect is greater than 1.
  theta2 = list(
    prior = "pc",
    param = c(0.5, 0.5)
  )
)

# Priors for  Year -> secular trend

hyper.year <- list(
  # 
  theta = list(prior="pc.prec", 
               param=c(1, 0.01)))


# formular for share_industry

formula <-   death_mod ~ 1 + offset(log(pop)) + 
  f(Year,model = "rw2", hyper = hyper.year, constr = TRUE, scale.model=TRUE) + 
  f(Region, model="bym2", graph="data/Districts", hyper = hyper.bym ) # random effect, spatial 

set.seed(20260818) # to get always the sanme results

# INLA 
inla.mod <- inla(formula,
                 data=dt,
                 family="nbinomial", # in case of overdispersion
                 #family = "poisson",
                 control.predictor = list(compute = TRUE), # Compute posterior marginals of linear predictors
                 control.compute = list(
                   config = TRUE,  # to use inla.posterior.sample later
                   dic = TRUE,
                   waic = TRUE,
                   cpo = TRUE)) # for model comparison

summary(inla.mod)


# get the posterior distribution with all information
# It generates random draws (here 1000 )from INLA's approximation to get the posterior distribution
post.samples <- inla.posterior.sample(n = 1000, result = inla.mod, seed=20261808)

# get the predicted values from 1000 samples and transform the prediction from the 
# log scale back to the mortality-rate scale

predlist <- do.call(cbind, lapply(post.samples, function(X)
  exp(X$latent[startsWith(rownames(X$latent), "Pred")]))) 

# unlist the predicted values
dtr <-array(unlist( predlist), dim=c(dim(dt)[1], 1000)) %>%
  as.data.frame() %>%
  cbind(dt)

# create a matrix of the 1000 samples to summarize them in the next step
sample_matrix <- as.matrix(dtr %>%  select(starts_with("V")))

# get the data, estimate excess mortality
results <- dtr %>%
  select(Year, Region, District) %>%
  mutate(
    fit = rowMedians(sample_matrix),
    LL  = rowQuantiles(sample_matrix, probs = 0.025),
    UL  = rowQuantiles(sample_matrix, probs = 0.975)
  ) %>% 
  left_join(dt) %>%
  mutate(
  # excess mortality
    exc = deaths - fit,
    exc_ll = deaths-UL,
    exc_ul = deaths-LL,
  # calculate the p-score in percentage
    p_score = exc/fit *100,
  # define with observed death is above UL and defined as excess
    col_sig = ifelse(deaths > UL, "yes", "no"))

# Plot 4 districts as examples

dty <- results %>%
  filter(District %in% c(111,2500,2207,1901)) # filter which districts

ggplot(dty) +
  geom_ribbon(aes(ymin=LL, ymax=UL,x=Year),fill="grey80",linetype=1, alpha=1) +
  geom_line(aes(x=Year, y=fit),col="grey20",lwd=1) +
  geom_point(aes(x=Year, y=deaths, col=col_sig),size=point_size ) +
  facet_wrap(~District , ncol=2, scales = "free_y") +
  xlab("Year") +
  ylab("Deaths") +
  scale_color_manual("For a 95% CrI:",
                     breaks=c("yes","no"),
                     # labels=c("Observed deaths", "Expected deaths" ),
                     values=c("red", "black"))+
  theme_bw() +
  theme(
    strip.text = element_text(size = plot_title_size),
    axis.text = element_text(size=axis_text_size),
    axis.title  = element_text(size=axis_text_size),
    legend.position = c(0.2,0.9),
    legend.text=element_text(size=axis_text_size),
    legend.title =element_text(size=axis_text_size),
    plot.title = element_text(size=plot_title_size),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank()) 

ggsave("excess_mortality.png",h=10,w=22)

####################
# Plot the results #
####################

# link results and shape files
dt2 <- dts %>%
  left_join(results) %>%
  filter(Year == 1918) %>%
  mutate(
   # calculate quantiles for the maps
    p_score_quar = cut( p_score,
                        breaks = c(min(p_score, na.rm = TRUE),0,quantile(p_score[p_score > 0], probs = seq(0.2, 1, 0.2), na.rm = TRUE)),
                        include.lowest = TRUE,
                        dig.lab = 3)
						)

# plot map with p-scores

ggplot(dt2) +
  geom_sf(aes(fill = p_score_quar),color = "black", linewidth = 0.1) +
  scale_fill_brewer(
    palette = "YlOrRd",
    name = "P-score") +
  labs(title = "P-score") +
  theme_bw() +
  theme(
    axis.text = element_blank(),
    axis.title  = element_text(size=axis_text_size),
    plot.title = element_text(size=plot_title_size),
    legend.title  = element_text(size=axis_text_size),
    legend.text  = element_text(size=axis_text_size),
    legend.position = "bottom"
  )

ggsave("maps_p_score.png",h=10,w=20)

# Which demographic and socioeconomic factors were associated with higher mortality during the 1918–1919 influenza pandemic?

dt3 <- dt2 %>% # scale data so that it is better comparable 
  mutate(
      dens_doc= scale(dens_doc),
      share_male= scale(share_male),
      share_5_14= scale(share_5_14),
      share_20_39= scale(share_20_39),
      share_60= scale(share_60),
      share_industry= scale(share_industry),
      gdp = scale(gdp),
      dens_pop= scale(dens_pop),
      houshold_house =scale(houshold_house),
      household_size= scale(household_size),
  )

# robust linear regression
Mod1 <- coef(summary(rlm(p_score ~ gdp, data=dt3)))
Mod2 <- coef(summary(rlm(p_score ~ share_industry, data=dt3)))
Mod3 <- coef(summary(rlm(p_score ~ dens_doc, data=dt3)))
Mod4 <- coef(summary(rlm(p_score ~ houshold_house, data=dt3)))
Mod5 <- coef(summary(rlm(p_score ~ household_size, data=dt3)))
Mod6 <- coef(summary(rlm(p_score ~ share_male, data=dt3)))
Mod7 <- coef(summary(rlm(p_score ~ share_5_14, data=dt3)))
Mod8 <- coef(summary(rlm(p_score ~ share_20_39, data=dt3)))
Mod9 <-  coef(summary(rlm(p_score ~ share_60, data=dt3)))
Mod10 <- coef(summary(rlm(p_score ~ dens_pop, data=dt3)))

# combine the results of the model
res_uni <- rbind(Mod1, Mod2, Mod3, Mod4, Mod5, Mod6, Mod7, Mod8, Mod9,Mod10) %>%
  data.frame() %>%
  mutate(Cofactor=row.names(.)) %>%
  filter( Cofactor=="gdp" | Cofactor=="share_industry" |   Cofactor=="dens_doc"|  Cofactor=="houshold_house"| Cofactor=="household_size"|
            Cofactor=="share_male" | Cofactor=="share_5_14" |   Cofactor=="share_20_39"| 
            Cofactor=="share_60" |   Cofactor=="dens_pop") %>%
  mutate(
  # get coefficient and 95% confidence intervals
  est= round(Value,2),
  Cl = round(Value - 1.96* Std..Error,2),
  Cu = round(Value + 1.96* Std..Error,2),
  Univariate = paste0(est," (",Cl,"-",Cu, ")"),
  # change names of the factors
  Cofactor = recode(Cofactor, 
                              "dens_doc" = "Private physicians per km2",
                              "share_male"  = "Share of men",
                              "share_5_14"  = "Share of 5-14 years old",
                              "share_20_39"    =  "Share of 20-39 years old",
                              "share_60"  = "Share of >= 60 years old",
                              "share_industry"     = "Share of industry",
                              "gdp" = "GDP per capita",
                              "dens_pop"    = "Population density",
                              "houshold_house"  = "Households per house",
                              "household_size" = "Household size"),
# reorder factors in the order they should appear in the plot
	Cofactor = factor(Cofactor, 
                           levels = c("Population density",
                                      "GDP per capita",
                                      "Share of industry",
                                      "Private physicians per km2",
                                      "Household size",
                                      "Households per house",
                                      "Share of men",
                                      "Share of 5-14 years old",
                                      "Share of 20-39 years old" ,
                                      "Share of >= 60 years old"))
  )


# Plot the coefficient

ggplot(res_uni , aes(x=Cofactor,ymin=Cl, ymax=Cu,y= est), position=pd) + 
  geom_hline(yintercept=0, colour="grey", lwd=lwd_size) + 
  geom_pointrange(position=pd, fatten=fatten_size,lwd=lwd_size) +
  scale_x_discrete(limits = rev)+
  xlab("") +
  ylab("Regression coefficients and 95% CI") +
  theme_bw()+
  theme(
    axis.text= element_text(size=axis_text_size),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size=axis_text_size),
    title =element_text(size=axis_text_size))+
  coord_flip()

ggsave("coefficient_pscore.png",h=10,w=15)


