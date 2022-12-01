#*******************************#
#  Integrated Population Model  #
#       Varanger Red Fox        #
#*******************************#

library(coda)
library(ggplot2)
library(reshape2)
library(viridis)
library(plyr)

# NOTE: 
# This has to be adjusted for every set of models you want to compare

# Loading & reformatting posteriors #
#-----------------------------------#

## Model A1
load('210507_RedFox_IPM_A1.RData')
out.mat <- as.matrix(RF.IPM.test)
#out.mat <- as.matrix(RF.IPM.test)[5001:15000,] # Drop first chain

for(i in 1:4){
	RF.IPM.test[[i]][which(is.na(RF.IPM.test[[i]]))] <- 0
}
#plot(RF.IPM.test, ask = T)
# --> Convergence is okay, but mixing is not so good for the majority of parameters
# --> All non-immigration parameters are well estimated
# --> SD of immigrant numbers is at the upper boundary, mean is at the lower boundary -> normal approximation may not work so well

sam.A1 <- melt(out.mat)
sam.A1$Model <- 'A1 (Bristol)'  


## Model A2
load('210507_RedFox_IPM_A2.RData')
out.mat <- as.matrix(RF.IPM.test)
#out.mat <- as.matrix(RF.IPM.test)[5001:15000,]

for(i in 1:4){
	RF.IPM.test[[i]][which(is.na(RF.IPM.test[[i]]))] <- 0
}
#plot(RF.IPM.test, ask = T)
# --> Same as for A1

sam.A2 <- melt(out.mat)
sam.A2$Model <- 'A2 (N-Sweden)'  

## Model A3
load('210507_RedFox_IPM_A3.RData')
out.mat <- as.matrix(RF.IPM.test)

for(i in 1:4){
	RF.IPM.test[[i]][which(is.na(RF.IPM.test[[i]]))] <- 0
}
#plot(RF.IPM.test, ask = T)
# --> Same as A1, but chain 4 diverged, then converged again in the middle

sam.A3 <- melt(out.mat)
sam.A3$Model <- 'A3 (Meta-all)'  


## Model A4
load('210507_RedFox_IPM_A4.RData')
out.mat <- as.matrix(RF.IPM.test)
#out.mat <- as.matrix(RF.IPM.test)[5001:15000,]

for(i in 1:4){
	RF.IPM.test[[i]][which(is.na(RF.IPM.test[[i]]))] <- 0
}
#plot(RF.IPM.test, ask = T)
# --> Generally like A1, but chain 3 only coverged with others at the very end

sam.A4 <- melt(out.mat)
sam.A4$Model <- 'A4 (Meta-sub)'  


## Model B1
load('210507_RedFox_IPM_B1.RData')
out.mat <- as.matrix(RF.IPM.test)
#out.mat <- as.matrix(RF.IPM.test)[5001:15000,]

for(i in 1:4){
	RF.IPM.test[[i]][which(is.na(RF.IPM.test[[i]]))] <- 0
}
#plot(RF.IPM.test, ask = T)
# --> Chains 1 & 4 converged to the same, unrealistic solution. Chains 2 & 3 converged to a more plausible solution, but only in the last 3rd of posterior samples



sam.B1 <- melt(out.mat)
sam.B1$Model <- 'B1 (Hoening-fixR)'  


## Model B1
load('210507_RedFox_IPM_B2.RData')
out.mat <- as.matrix(RF.IPM.test)
#out.mat <- as.matrix(RF.IPM.test)[5001:15000,]

for(i in 1:4){
	RF.IPM.test[[i]][which(is.na(RF.IPM.test[[i]]))] <- 0
}
#plot(RF.IPM.test, ask = T)
# --> Similar to B1, but much worse. Chain 1 maybe on the way to converge with 2 & 3, chain 4 complete BS stray. A disaster.

sam.B2 <- melt(out.mat)
sam.B2$Model <- 'B2 (Hoening-priorR)'  


# Combining data #
#----------------#

## Bind data together
sam.data <- rbind(sam.A1, sam.A2, sam.A3, sam.A4, sam.B1, sam.B2)

## Re-name columns
colnames(sam.data) <- c('Sample', 'Parameter', 'Estimate', 'Model')


# Plotting - Vital rates #
#------------------------#

## Define VR parameters 
VR.params <- c(paste0('Mu.mH[', 1:5, ']'), paste0('Mu.mO[', 1:5, ']'), 'JuvAdRatio',
               paste0('Mu.Psi[', 2:5, ']'), paste0('Mu.rho[', 2:5, ']'), 
               'mean.S0', 'Mu.Imm',
               'sigma.mH', 'sigma.Psi', 'sigma.rho', 'sigma.Imm')

## Plot VRs
pdf('PosteriorsComp_ModelsAB_VRs.pdf', width = 11, height = 8)

ggplot(subset(sam.data, Parameter %in% VR.params), aes(x = Estimate, group = Model)) +
        geom_density(aes(color = Model), fill = NA) +
        facet_wrap(~Parameter, scales = 'free') +
        ggtitle('Vital rate parameters') +
        scale_color_viridis(discrete = T) +
        theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), plot.title = element_text(face = 'bold'))

dev.off()


# Plotting - Total population sizes #
#-----------------------------------#

## Set parameters to plot
N.params <- c(paste0('N.tot[', 1:16, ']'))

## Plot 
pdf('PosteriorsComp_ModelsAB_Ntot.pdf', width = 11, height = 8)

ggplot(subset(sam.data, Parameter %in% N.params), aes(x = Estimate, group = Model)) +
  geom_density(aes(color = Model), fill = NA) +
  facet_wrap(~Parameter, scales = 'free') +
  ggtitle('Population sizes') +
  scale_color_viridis(discrete = T) +
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), plot.title = element_text(face = 'bold'))

dev.off()


# Plotting - Breeding population sizes #
#--------------------------------------#

## Set parameters to plot
B.params <- c(paste0('B.tot[', 1:16, ']'))

## Plot 
pdf('PosteriorsComp_ModelsAB_Btot.pdf', width = 11, height = 8)

ggplot(subset(sam.data, Parameter %in% B.params), aes(x = Estimate, group = Model)) +
  geom_density(aes(color = Model), fill = NA) +
  facet_wrap(~Parameter, scales = 'free') +
  ggtitle('Breeding population sizes') +
  scale_color_viridis(discrete = T) +
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), plot.title = element_text(face = 'bold'))

dev.off()

# Plotting - Recruit numbers #
#----------------------------#

## Set parameters to plot
R.params <- c(paste0('R.tot[', 1:16, ']'))

## Plot 
pdf('PosteriorsComp_ModelsAB_Rtot.pdf', width = 11, height = 8)

ggplot(subset(sam.data, Parameter %in% R.params), aes(x = Estimate, group = Model)) +
  geom_density(aes(color = Model), fill = NA) +
  facet_wrap(~Parameter, scales = 'free') +
  ggtitle('Recruit numbers') +
  scale_color_viridis(discrete = T) +
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), plot.title = element_text(face = 'bold'))

dev.off()


# Plotting - Immigrant numbers #
#------------------------------#

## Set parameters to plot
I.params <- c(paste0('Imm[', 2:16, ']'))

## Plot 
pdf('PosteriorsComp_ModelsAB_Imm.pdf', width = 11, height = 8)

ggplot(subset(sam.data, Parameter %in% I.params), aes(x = Estimate, group = Model)) +
  geom_density(aes(color = Model), fill = NA) +
  facet_wrap(~Parameter, scales = 'free') +
  ggtitle('Immigrant numbers') +
  scale_color_viridis(discrete = T) +
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), plot.title = element_text(face = 'bold'))

dev.off()

# Plotting - Population size over time #
#--------------------------------------#

## Summarize posterior into median, 90% and 95% CIs
data.sum <- ddply(sam.data, .(Parameter, Model), summarise, 
                  median = median(Estimate, na.rm = T), 
                  lCI_90 = quantile(Estimate, probs = 0.05, na.rm = T), uCI_90 = quantile(Estimate, probs = 0.95, na.rm = T), 
                  lCI_50 = quantile(Estimate, probs = 0.25, na.rm = T), uCI_50 = quantile(Estimate, probs = 0.75, na.rm = T))

## Subset data 
data.Ntot <- subset(data.sum, Parameter%in%N.params)

data.Ntot$indexT <- rep(c(1:15), each = 6)
data.Ntot$Year <- data.Ntot$indexT+2003
data.Ntot <- data.Ntot[order(data.Ntot$Year),]

pdf("Comp_ModelsAB_Ntot_time.pdf", width = 8, height = 4)

ggplot(data.Ntot, aes(x = Year, y = median, group = Model)) + 
  geom_line(aes(color = Model)) + 
  geom_ribbon(aes(ymin = lCI_90, ymax = uCI_90, fill = Model), alpha = 0.1) +
  ylab('Total population size') + 
  scale_x_continuous(breaks = c(2003:2018)) + 
  scale_color_viridis(discrete = T) +
  scale_fill_viridis(discrete = T) +
  theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))

ggplot(data.Ntot, aes(x = Year, y = median, group = Model)) + 
  geom_line(aes(color = Model)) + 
  geom_line(aes(y = lCI_90, color = Model), linetype = 'dotted', size = 0.25) + 
  geom_line(aes(y = uCI_90, color = Model), linetype = 'dotted', size = 0.25) + 
  ylab('Total population size') + 
  scale_x_continuous(breaks = c(2003:2018)) + 
  scale_color_viridis(discrete = T) +
  theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))

dev.off()



