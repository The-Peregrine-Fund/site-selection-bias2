# load packages
library(Rfast)
library(tidyr)
library(lme4)
library(ggplot2)
library(tidyverse)
library(viridis)
library(ggpubr)
library(cowplot)

theme_x=(function (base_size = 16, base_family = "") 
{
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          strip.background = element_rect(colour = "black", 
                                          size = 0.5), legend.key = element_blank(),
          legend.position=c(0.1, .9))
})


#### Figure 1 ####
nsites = 100 ; nsurveys = 1 ; nyears = 10 ; mean.beta = 0 ; sd.lam = c(0, 0)
my_rate = c(0.0005, 0.005, 0.05)
nsim = 1000
year <- (1:nyears) - 1
saved_N <- list(NA, NA, NA)

# set random seed
set.seed(123)

# run the simulations
for (j in 1:length(my_rate)) {
  rate <- my_rate[j]
  
  alpha <- rgamma(n = nsites, shape = 0.1, rate = rate) # rate is the parameter to adjust the distance of outliers
  beta <- rnorm(n = nsites, mean = mean.beta, sd = sd.lam[2]) # population trend, set to zero with no variability
  lam <- N <- p <- ep <- array(NA, dim = c(nsites, nyears)) # create a matrix with nsites rows and nyears columns for all parameters
  
  # fill the matrix with the expected abundance
  for (t in (0:nyears)) { 
    lam[, t] <- exp(log(alpha) + beta * year[t])
  }
  
  # fill the matrix with the realised abundance
  for (t in (1:nyears)) { 
    N[, t] <- rnorm(n = nsites, mean = lam[, t], sd = lam[, t]) # temporal SD equal to expected site abundance
  }
  
  # set negative abundances to zero (as did the authors)
  N <- ifelse(N<0, 0, N)
  
  # save the matrix with the realised abundance for plotting
  saved_N[[j]] <- N
  
} # j


# format data for plotting
format_my_data <- function(N_matrix) {
  df_N <- as.data.frame.matrix(N_matrix)
  df_N <- df_N[order(df_N$V1, decreasing = TRUE), ]
  colnames(df_N) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
  df_N <- cbind(site_ID = c(1:100), site_status = c(rep("selected", times = 10), rep("not_selected", times = 90)), df_N)
  
  # convert to long format
  df_N_long <- df_N %>%
    pivot_longer(cols = -c(site_ID, site_status), names_to = "year", values_to = "abundance")
  
  # set right column types
  df_N_long$site_ID <- as.factor(df_N_long$site_ID)
  df_N_long$site_status <- as.factor(df_N_long$site_status)
  df_N_long$year <- as.integer(df_N_long$year)
  df_N_long$abundance <- as.double(df_N_long$abundance)
  
  return(df_N_long)
}

df_abundances_0005 <- format_my_data(N_matrix = saved_N[[1]])
df_abundances_005 <- format_my_data(N_matrix = saved_N[[2]])
df_abundances_05 <- format_my_data(N_matrix = saved_N[[3]])

# make the plots
my_palette <- viridis(n = 6)

plot_rate_0005 <- ggplot(df_abundances_0005, aes(x = year, y = abundance, group = site_ID)) +
  geom_line(data = subset(df_abundances_0005, df_abundances_0005$site_ID == "1"), color = my_palette[1], size = 1) +
  geom_line(data = subset(df_abundances_0005, df_abundances_0005$site_ID == "2"), color = my_palette[2], size = 1) +
  geom_line(data = subset(df_abundances_0005, df_abundances_0005$site_ID == "3"), color = my_palette[3], size = 1) +
  geom_line(data = subset(df_abundances_0005, df_abundances_0005$site_ID == "4"), color = my_palette[4], size = 1) +
  geom_line(data = subset(df_abundances_0005, df_abundances_0005$site_ID == "5"), color = my_palette[5], size = 1) +
  geom_line(data = subset(df_abundances_0005, df_abundances_0005$site_ID == "6"), color = my_palette[6], size = 1) +
  theme_classic() +
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5)) +
  xlab("Year") +
  ylab("Abundance") +
  ggtitle("Rate = 0.0005")

plot_rate_005 <- ggplot(df_abundances_005, aes(x = year, y = abundance, group = site_ID)) +
  geom_line(data = subset(df_abundances_005, df_abundances_005$site_ID == "1"), color = my_palette[1], size = 1) +
  geom_line(data = subset(df_abundances_005, df_abundances_005$site_ID == "2"), color = my_palette[2], size = 1) +
  geom_line(data = subset(df_abundances_005, df_abundances_005$site_ID == "3"), color = my_palette[3], size = 1) +
  geom_line(data = subset(df_abundances_005, df_abundances_005$site_ID == "4"), color = my_palette[4], size = 1) +
  geom_line(data = subset(df_abundances_005, df_abundances_005$site_ID == "5"), color = my_palette[5], size = 1) +
  geom_line(data = subset(df_abundances_005, df_abundances_005$site_ID == "6"), color = my_palette[6], size = 1) +
  theme_classic() +
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5)) +
  xlab("Year") +
  ylab("Abundance") +
  ggtitle("Rate = 0.005")

plot_rate_05 <- ggplot(df_abundances_05, aes(x = year, y = abundance, group = site_ID)) +
  geom_line(data = subset(df_abundances_05, df_abundances_05$site_ID == "1"), color = my_palette[1], size = 1) +
  geom_line(data = subset(df_abundances_05, df_abundances_05$site_ID == "2"), color = my_palette[2], size = 1) +
  geom_line(data = subset(df_abundances_05, df_abundances_05$site_ID == "3"), color = my_palette[3], size = 1) +
  geom_line(data = subset(df_abundances_05, df_abundances_05$site_ID == "4"), color = my_palette[4], size = 1) +
  geom_line(data = subset(df_abundances_05, df_abundances_05$site_ID == "5"), color = my_palette[5], size = 1) +
  geom_line(data = subset(df_abundances_05, df_abundances_05$site_ID == "6"), color = my_palette[6], size = 1) +
  theme_classic() +
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5)) +
  xlab("Year") +
  ylab("Abundance") +
  ggtitle("Rate = 0.05")


cavu=read.csv("data//CAVU 10 years long 11.28.23.csv",header=T,se=',',na.strings=c("","NA"))

plot_real=ggplot(cavu)+
  geom_line(aes(x=year,y=count,color=fct_reorder(site, -count)),size=1)+
  theme_classic() +
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5)) +
  xlab("Year") +
  ylab("Abundance") +
  ggtitle("Cape Vulture")+
  scale_color_viridis_d()+
  theme(legend.position=c('none'))
plot_real



# arrange the plots together to make the figure
Figure1 <- ggarrange(plot_rate_0005, plot_rate_005, plot_rate_05,plot_real,
                     labels = c("(a)", "(b)", "(c)", "(d)"),
                     ncol = 2, nrow = 2)

Figure1

# Regress mean vs SD of Cape Vulture data
avc=data.frame(aggregate(cavu$count,list(site=cavu$site),mean,na.rm=t))
avc$av=avc$x
avc$x=NULL

stdc=data.frame(aggregate(cavu$count,list(site=cavu$site),sd,na.rm=t))
avc$up=avc$av+stdc$x
avc$down=avc$av-stdc$x
avc$sd=stdc$x
avc

res=lm(sd~av,data=avc)
summary(res)
pred=data.frame(predict(res,se=T))
pred$up=pred$fit+1.96*pred$se.fit
pred$down=pred$fit-1.96*pred$se.fit
pred$av=avc$av

p1=ggplot()+theme_x()+xlab('Abundance')+ylab('Standard Deviation')+
  geom_point(aes(x=av,y=sd), data=avc,size=3)+theme(text = element_text(size = 19))+
  geom_line(aes(x=av,y=fit),data=pred)+ 
  geom_ribbon(aes(x=av,ymax=up,ymin=down),alpha=0.5,data=pred)
p1

#### Figure 2 ####
# 1 simulated plot

# simulation settings
nsites = 100 ; nsurveys = 1 ; nyears = 20 ; mean.beta = 0; sd.lam = c(0,0)
rate = 0.0005
my_dep = c(0.06920)
nsim = 1000
year <- (1:nyears) - 1
saved_N_dependence <- list(NA, NA, NA)

# set random seed
set.seed(112823)

# run the simulations
for (j in 1:length(my_dep)) {
  dep_ratio <- my_dep[j]
  
  alpha <- rgamma(n = nsites, shape = 0.1, rate = rate) # rate is the parameter to adjust the distance of outliers
  beta <- rnorm(n = nsites, mean = mean.beta, sd = sd.lam[2]) # population trend, set to zero with no variability
  lam <- N <- p <- ep <- array(NA, dim = c(nsites, nyears)) # create a matrix with nsites rows and nyears columns for all parameters
  
  # fill the matrix with the expected abundance
  for (t in (0:nyears)) { 
    lam[, t] <- exp(log(alpha) + beta * year[t])
  }
  
  # fill the matrix with the realized abundance
  for (t in (1:nyears)) { 
    N[, t] <- rnorm(n = nsites, mean = lam[, t], sd = dep_ratio*lam[, t]) # temporal SD proportional to expected site abundance
  }
  
  # set negative abundances to zero (as did the authors)
  N <- ifelse(N<0, 0, N)
  
  # save the matrix with the realized abundance for plotting
  saved_N_dependence[[j]] <- N
  
} # j

format_my_data <- function(N_matrix) {
  df_N <- as.data.frame.matrix(N_matrix)
  df_N <- df_N[order(df_N$V1, decreasing = TRUE), ]
  colnames(df_N) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                      "11","12","13","14", "15", "16", "17", "18", "19", "20")
  df_N <- cbind(site_ID = c(1:100), site_status = c(rep("selected", times = 10), rep("not_selected", times = 90)), df_N)
  
  # convert to long format
  df_N_long <- df_N %>%
    pivot_longer(cols = -c(site_ID, site_status), names_to = "year", values_to = "abundance")
  
  # set right column types
  df_N_long$site_ID <- as.factor(df_N_long$site_ID)
  df_N_long$site_status <- as.factor(df_N_long$site_status)
  df_N_long$year <- as.integer(df_N_long$year)
  df_N_long$abundance <- as.double(df_N_long$abundance)
  
  return(df_N_long)
}

# format data for plotting
df_abundances_dep02 <- format_my_data(N_matrix = saved_N_dependence[[1]])

# make the plots
my_palette <- viridis(n = 10)

p2 <- ggplot(df_abundances_dep02, aes(x = year, y = abundance, group = site_ID)) +
  geom_line(color = "lightgrey", size = 1) +
  geom_line(data = subset(df_abundances_dep02, df_abundances_dep02$site_ID == "1"), color = my_palette[1], size = 1) +
  geom_line(data = subset(df_abundances_dep02, df_abundances_dep02$site_ID == "2"), color = my_palette[2], size = 1) +
  geom_line(data = subset(df_abundances_dep02, df_abundances_dep02$site_ID == "3"), color = my_palette[3], size = 1) +
  geom_line(data = subset(df_abundances_dep02, df_abundances_dep02$site_ID == "4"), color = my_palette[4], size = 1) +
  geom_line(data = subset(df_abundances_dep02, df_abundances_dep02$site_ID == "5"), color = my_palette[5], size = 1) +
  geom_line(data = subset(df_abundances_dep02, df_abundances_dep02$site_ID == "6"), color = my_palette[6], size = 1) +
  geom_line(data = subset(df_abundances_dep02, df_abundances_dep02$site_ID == "7"), color = my_palette[7], size = 1) +
  geom_line(data = subset(df_abundances_dep02, df_abundances_dep02$site_ID == "8"), color = my_palette[8], size = 1) +
  geom_line(data = subset(df_abundances_dep02, df_abundances_dep02$site_ID == "9"), color = my_palette[9], size = 1) +
  geom_line(data = subset(df_abundances_dep02, df_abundances_dep02$site_ID == "10"), color = my_palette[10], size = 1) +
  theme_x() +
  theme(text = element_text(size = 19))+
  xlab("Year") +
  ylab("Abundance") 

p2

# run many simulations for analysis

# simulation settings
nsites = 100 ; nsurveys = 1 ; nyears = 20 ; mean.beta = 0 ; sd.lam = c(0, 0)
rate = 0.0005
my_dep = c(0.07) # coefficient between temporal variance and abundance
nsim = 1000
year <- (1:nyears) - 1
df_simul_results <- data.frame(matrix(NA, nrow = nsim, ncol = length(my_dep)))
colnames(df_simul_results) <- paste0("case_", 1:length(my_dep))

# store start time
(start.time_total <- Sys.time())

# set random seed
set.seed(112823)

# run the simulations
for (j in 1:length(my_dep)) {
  dep_ratio <- my_dep[j]
  
  for (i in 1:nsim) {
    
    alpha <- rgamma(n = nsites, shape = 0.1, rate = rate) # rate is the parameter to adjust the distance of outliers
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
    
    print(c(i,j))
    
  } # i
} # j

# print end time and total execution time
(end.time_total <- Sys.time())
(time.taken_total <- difftime(end.time_total, start.time_total, units = "mins"))
colMeans(df_simul_results) # mean estimated population change across all simulations


# convert the table to long format 
df_all_results_var_dependence_gradient <- rbind(data.frame(var_dependence = my_dep[1], estim_pop_change = df_simul_results$case_1))
df_all_results_var_dependence_gradient$var_dependence <- as.factor(df_all_results_var_dependence_gradient$var_dependence)

# save simulation output
mean(df_all_results_var_dependence_gradient$estim_pop_change<1)

# plot the results
p3 <- ggplot(df_all_results_var_dependence_gradient, aes(x = var_dependence, y = estim_pop_change)) + 
  geom_rect(data = NULL, aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 1), fill = "lightgrey") +
  geom_boxplot(fill = NA) +
  theme_classic() +
  theme(text = element_text(size = 19)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  xlab("Coefficient") +
  ylab("Apparent\nPopulation Change")


p3

# arrange the plots together to make the figure
Figure2=plot_grid(p1,p2,p3,labels=c("(a)","(b)","(c)"),
                  ncol = 3,label_size = 19,
                  align = 'hv')

Figure2
