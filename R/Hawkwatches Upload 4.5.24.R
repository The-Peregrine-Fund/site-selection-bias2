# load packages
library(Rfast)
library(tidyr)
library(lme4)
library(ggplot2)
library(tidyverse)
library(viridis)
library(ggpubr)
library(MASS)
library(cowplot)

# Analyze hawkwatch data to determine shape of distributions and correlation bt mean and SD
t=read.csv=read.csv("data//HW 3.27.24.csv",header=T,se=',',na.strings=c("","NA"))

t=subset(t,t$species=='Pandion haliaetus')
t$rate=t$count/t$hrs
am=data.frame(aggregate(t$rate,list(site=t$site),mean,na.rm=T))
as=data.frame(aggregate(t$rate,list(site=t$site),sd,na.rm=T))
d=data.frame(site=am$site,m=am$x,s=as$x)

fitdistr(am$x,densfun='gamma')
 summary(lm(s~m,data=d))


# run simulations 
# simulation settings
nsites = 100 ; nsurveys = 1 ; nyears = 20 ; mean.beta = 0 ; sd.lam = c(0, 0)
shape = c(1.2535605,1.8666828,0.6058704,0.44195343,0.3917880,0.253360307,1.3542088,0.4658104,0.9922213,0.3687364,0.7664877)
rate =  c(0.6363898,0.2439588,3.0894622,0.09669700,1.2133352,0.001807080,2.0306067,0.5118202,0.4418974,0.4886916,0.6909364)
my_dep = c(0.31061,0.25797,0.262666,0.23625,0.388404,0.24113,0.42263,0.52013,0.3660,0.25706,0.41207) # coefficient between temporal variance and abundance
nsim = 1000
year <- (1:nyears) - 1
df_simul_results <- data.frame(matrix(NA, nrow = nsim, ncol = length(my_dep)))
UCL <- data.frame(matrix(NA, nrow = nsim, ncol = length(my_dep)))
LCL <- data.frame(matrix(NA, nrow = nsim, ncol = length(my_dep)))
colnames(df_simul_results) <- paste0("case_", 1:length(my_dep))

# store start time
(start.time_total <- Sys.time())

# set random seed
set.seed(112823)

# run the simulations
for (j in 1:length(my_dep)) {
  dep_ratio <- my_dep[j]
  shape1 <- shape[j]
  rate1 <- rate[j]
  for (i in 1:nsim) {
    
    alpha <- rgamma(n = nsites, shape = shape1, rate = rate1) # rate is the parameter to adjust the distance of outliers
    beta <- rnorm(n = nsites, mean = mean.beta, sd = sd.lam[2]) # population trend, set to zero with no variability
    lam <- N <- p <- ep <- array(NA, dim = c(nsites, nyears)) # create a matrix with nsites rows and nyears columns for all parameters
    
    # fill the matrix with the expected abundance
    for (t in (0:nyears)) { 
      lam[, t] <- exp(log(alpha) + beta * year[t])
    }
    
    # fill the matrix with the realised abundance
    for (t in (1:nyears)) { 
      N[, t] <- rnorm(n = nsites, mean = lam[, t], sd = dep_ratio*lam[, t]) # temporal SD proportional to expected site abundance
    }
    
    # set negative abundances to zero (as did the authors)
    N <- ifelse(N<0, 0, N)
    
    # select the 10 sites with the greatest initial abundance
    vect_10_selected_sites <- c(Rfast::nth(x = N[, 1], k = 10, num.of.nths = 10, descending = TRUE, index.return = TRUE))
    N_selected <- N[vect_10_selected_sites, ]
    
    # convert simulated data from matrix to dataframe format
    df_counts <- as.data.frame.array(N_selected)
    df_counts <- cbind(site = as.factor(c(1:10)), df_counts) # add column with site id
    df_counts <- pivot_longer(df_counts, !site, names_to = "year", values_to = "count", names_prefix = "V") # convert the dataframe to long format
    df_counts$year <- as.integer(df_counts$year)
    
    # fit the linear model
    my_model <- lmer(log(df_counts$count+1) ~ year + (1|site), data = df_counts)
    
    # compute the apparent_pop_change
    estimated_beta <- coef(my_model)$site$year[1]
    apparent_pop_change <- exp(estimated_beta*19)
    df_simul_results[i, j] <- apparent_pop_change
    CI=confint(my_model,year,level=0.90)
    UCL[i,j]=exp(CI[4,2]*19)
    LCL[i,j]=exp(CI[4,1]*19)
    
    print(c(i,j))
    
  } # i
} # j

# print end time and total execution time
(end.time_total <- Sys.time())
(time.taken_total <- difftime(end.time_total, start.time_total, units = "mins"))
colMedians(df_simul_results) # mean estimated population change across all simulations
apply(df_simul_results,2,function (df_simul_results) mean(df_simul_results<1))
apply(UCL,2,function (UCL) mean(UCL<0.81))


# convert the table to long format 
df_all_results_var_dependence_gradient <- rbind(data.frame(Species = spec[1], estim_pop_change = df_simul_results$case_1),
                                                  data.frame(Species = spec[2], estim_pop_change = df_simul_results$case_2),
                                                  data.frame(Species = spec[3], estim_pop_change = df_simul_results$case_3),
                                                  data.frame(Species = spec[4], estim_pop_change = df_simul_results$case_4),
                                                  data.frame(Species = spec[5], estim_pop_change = df_simul_results$case_5),
                                                  data.frame(Species = spec[6], estim_pop_change = df_simul_results$case_6),
                                                  data.frame(Species = spec[7], estim_pop_change = df_simul_results$case_7),
                                                  data.frame(Species = spec[8], estim_pop_change = df_simul_results$case_8),
                                                  data.frame(Species = spec[9], estim_pop_change = df_simul_results$case_9),
                                                  data.frame(Species = spec[10], estim_pop_change = df_simul_results$case_10),
                                                  data.frame(Species = spec[11], estim_pop_change = df_simul_results$case_11))


ucl1 <- rbind(data.frame(Species = spec[1], UCL1 = UCL[,1]),
                                                data.frame(Species = spec[2], UCL1 = UCL[,2]),
                                                data.frame(Species = spec[3], UCL1 = UCL[,3]),
                                                data.frame(Species = spec[4], UCL1 = UCL[,4]),
                                                data.frame(Species = spec[5], UCL1 = UCL[,5]),
                                                data.frame(Species = spec[6], UCL1 = UCL[,6]),
                                                data.frame(Species = spec[7], UCL1 = UCL[,7]),
                                                data.frame(Species = spec[8], UCL1 = UCL[,8]),
                                                data.frame(Species = spec[9], UCL1 = UCL[,9]),
                                                data.frame(Species = spec[10], UCL1 = UCL[,10]),
                                                data.frame(Species = spec[11], UCL1 = UCL[,11]))
# make figures
ucl=data.frame(ucl=ucl, Species=spec)
m=merge(ucl,df_all_results_var_dependence_gradient,by=c("Species"))

m$Species <- 
  fct_reorder(m$Species, -m$estim_pop_change)

# plot the results
p3 <- ggplot(m, aes(x = Species, y = estim_pop_change)) + 
  geom_rect(data = NULL, aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0.81), fill = "lightgrey") +
  geom_rect(data = NULL, aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0.5), fill = "darkgrey") +
  geom_boxplot(fill = NA) +
  theme_classic() +
  theme(text = element_text(size = 15),
       axis.text.x=element_blank())+
  geom_hline(yintercept = 1, linetype = "dashed") +
  xlab("") +
  ylab("Apparent\nPopulation\nChange")+
  ylim(0,1.65)

p3

ucl1$Species <- 
  fct_reorder(ucl1$Species, -m$estim_pop_change)

p4 <- ggplot(ucl1, aes(x = Species, y = UCL1)) + 
  geom_rect(data = NULL, aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0.81), fill = "lightgrey") +
  geom_rect(data = NULL, aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0.5), fill = "darkgrey") +
  geom_boxplot(fill = NA) +
  theme_classic() +
  theme(text = element_text(size = 15),
         axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  xlab("Simulated Species") +
  ylab("90% UCL of Apparent\nPopulation\nChange")+
  ylim(0,1.65)

p4
