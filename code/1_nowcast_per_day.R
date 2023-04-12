# Load packages
library(tidyverse)
library(cmdstanr)
library(lubridate)
library(RcppRoll)
library(parallel)
library(posterior)
library(bayesplot)
theme_set(theme_bw())

source("./daily_nowcast_per_strat.R")
source("./daily_nowcast_per_strat_postproc.R")

