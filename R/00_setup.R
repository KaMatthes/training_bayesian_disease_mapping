

if (!require(tidyverse)) install.packages("tidyverse")
if (!require(sf)) install.packages("sf")
if (!require(spdep)) install.packages("spdep")
if (!require(matrixStats)) install.packages("matrixStats")
if (!require(biscale)) install.packages("biscale")
if (!require(cowplot)) install.packages("cowplot")
if (!require(MASS)) install.packages("MASS")
if (!require(conflicted)) install.packages("conflicted")



library(tidyverse)
library(sf)
library(spdep)
library(INLA)
library(matrixStats)
library(biscale)
library(cowplot)
library(MASS)
library(conflicted)

# conflicted functions
conflict_prefer("select", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("recode", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("summarise", "dplyr")
conflict_prefer("arrange", "dplyr")