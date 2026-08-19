# load your packages

rm(list=ls())
source("R/00_setup.R")

# define plot parameters
axis_text_size <- 20
plot_title_size <- 25
lwd_size <- 1.5


# load your data

dtc <- read.csv("data/data_cases_copen.csv", sep = ";") %>%
  mutate(gdp_z = scale(gdp))

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
  left_join(dtc)

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

# formular for share_industry

formula <- cases ~ 1 + share_industry+ offset(log(pop)) +
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

# plot fix effect

# Get marginal fixed effect of GPP
ind_marg <- inla.mod$marginals.fixed[["share_industry"]]

# Transform whole posterior to RR scale
ind_RR <- inla.tmarginal(exp, ind_marg)

# Median and 95% CrI
ind_median <- inla.qmarginal(0.5, ind_RR )
ind_LL     <- inla.qmarginal(0.025, ind_RR )
ind_UL     <- inla.qmarginal(0.975, ind_RR)

# Plot distribution

ggplot(ind_RR, aes(x = x, y = y)) +
  geom_line(linewidth = 1) +
  geom_area(alpha = 0.1) +
  geom_vline(xintercept = ind_median, lwd=lwd_size) +
  geom_vline(xintercept = c(ind_LL, ind_UL),
             linetype = "dashed") +
  labs(
    title = "Posterior distribution of share of industry",
    x = "Risk Ratio",
    y = "Posterior density"
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size=axis_text_size),
    axis.title  = element_text(size=axis_text_size),
    plot.title = element_text(size=plot_title_size)
  )

ggsave("share_industry.png",h=10,w=20)

# get the posterior distribution, not only the fitted value

# We draw 1,000 plausible realizations from the fitted posterior distribution
post.samples <- inla.posterior.sample(n = 1000, result = inla.mod, seed=20261808)

# get the predicted values from 1000 samples and transform the prediction from the log scale back to the mortality-rate scale

predlist <- do.call(cbind, lapply(post.samples, function(X)
  exp(X$latent[startsWith(rownames(X$latent), "Pred")]))) 

# unlist the predicted values
dtr <-array(unlist( predlist), dim=c(dim(dt)[1], 1000)) %>%
  as.data.frame() %>%
  cbind(dt)
  
# create a matrix of the 1000 samples to summarize in the next step
sample_matrix <- as.matrix(dtr %>%  select(starts_with("V")))

# summarize the 1000 samples per district

results <- dtr %>%
  select(Region, District) %>%
  mutate(
    fit = rowMedians(sample_matrix),
    LL  = rowQuantiles(sample_matrix, probs = 0.025),
    UL  = rowQuantiles(sample_matrix, probs = 0.975)
  ) %>% 
  left_join(dt) %>%
  mutate(
  # calculate quantiles for the maps
    fit_quar = ntile(fit, 5),
    fit_quar= factor(fit_quar,levels = 1:5,labels = c("Q1", "Q2", "Q3", "Q4", "Q5")),
    case_quar = ntile(cases, 5),
    case_quar= factor(case_quar,levels = 1:5,labels = c("Q1", "Q2", "Q3", "Q4", "Q5"))
    
  )

####################
# Plot the results #
####################

# link results with shapefiles
dt2 <- dts %>%
  left_join(results)  

# plot observed cases in quantile
ggplot(dt2) +
  geom_sf(aes(fill = case_quar),color = "black", linewidth = 0.1) +
  scale_fill_brewer(
    palette = "YlOrRd",
    name = "Observed Cases"
  ) +
  labs(title = "Observed Cases") +
  theme_bw() +
  theme(
    axis.text = element_blank(),
    axis.title  = element_text(size=axis_text_size),
    plot.title = element_text(size=plot_title_size),
    legend.title  = element_text(size=axis_text_size),
    legend.text  = element_text(size=axis_text_size),
    legend.position = "bottom"
  )

ggsave("maps_observed.png",h=10,w=20)

# plot fitted values in quantile
ggplot(dt2) +
  geom_sf(aes(fill = fit_quar),color = "black", linewidth = 0.1) +
  scale_fill_brewer(
    palette = "YlOrRd",
    name = "Fitted values"
  ) +
  labs(title = "Fitted values") +
  theme_bw() +
  theme(
    axis.text = element_blank(),
    axis.title  = element_text(size=axis_text_size),
    plot.title = element_text(size=plot_title_size),
    legend.title  = element_text(size=axis_text_size),
    legend.text  = element_text(size=axis_text_size),
    legend.position = "bottom"
  )

ggsave("maps_fitted.png",h=10,w=20)

# bivariate map
# create the data for a bivariate map

map_biv <- dt2 %>%
  bi_class(
    x = fit,
    y = share_industry,
    style = "quantile",
    dim = 3
  )

# plot a bivariate map

biplot <- ggplot(map_biv) +
  geom_sf(aes(fill = bi_class),color = "white",linewidth = 0.15 ) +
  bi_scale_fill(pal = "DkBlue",dim = 3) +
  bi_theme() +
  theme(legend.position = "none")

# legende for bivariate map

legend <- bi_legend(
  pal = "DkBlue",
  dim = 3,
  xlab = "Incidence",
  ylab = "Share Industry",
  size = 12) +
  theme(
    axis.title = element_text(size = axis_text_size),
    axis.text = element_text(size = axis_text_size)
  )

# combine plot and legend
ggdraw() +
  draw_plot(biplot, 0, 0.05, 1, 1) +
  draw_plot(legend, 0.05, 0.7, 0.25, 0.25)

ggsave("maps_bivariate.png",h=10,w=20)