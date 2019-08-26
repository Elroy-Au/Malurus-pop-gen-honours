## Population Genomics of M. assimilis 
## Elroy Au

### ================================== ###
### Samples Distribution Map           ###
### ================================== ###

# Generate a map of the distribution of M. assimilis tissue samples across Australia 
# coded with collection date 

# Load libraries
library(oz)       # map of Australia for plotting 
library(readxl)   # for reading excel files 
                  # read in excel file: X29_7_Samples_Long_Lat_Date
                  # located in: "Documents/University/Honours 2019/Australian Birds
                  # /tissues/tissue samples/samples-plot/

oz()              #plot a map of Australia with state lines 

lat1 <- X29_7_Samples_Long_Lat_Date$`Latitude-1`        #set latitude, longitude
lat2 <- X29_7_Samples_Long_Lat_Date$`Latitude-2`       
lat3 <- X29_7_Samples_Long_Lat_Date$`Latitude-3`        
lat4 <- X29_7_Samples_Long_Lat_Date$`Latitude-4`        
lat5 <- X29_7_Samples_Long_Lat_Date$`Latitude-5`        
lat6 <- X29_7_Samples_Long_Lat_Date$`Latitude-6`       
lat7 <- X29_7_Samples_Long_Lat_Date$`Latitude-7`
long1 <- X29_7_Samples_Long_Lat_Date$`Longitude-1`
long2 <- X29_7_Samples_Long_Lat_Date$`Longitude-2`
long3 <- X29_7_Samples_Long_Lat_Date$`Longitude-3`
long4 <- X29_7_Samples_Long_Lat_Date$`Longitude-4`
long5 <- X29_7_Samples_Long_Lat_Date$`Longitude-5`
long6 <- X29_7_Samples_Long_Lat_Date$`Longitude-6`
long7 <- X29_7_Samples_Long_Lat_Date$`Longitude-7`

points(long1,lat1, col = "#85C1E9", bg = "#85C1E9", pch = 24)   #plot points 
points(long2,lat2, col = "#5DADE2", bg = "#5DADE2", pch = 24)
points(long3,lat3, col = "#3498DB", bg = "#3498DB", pch = 24)
points(long4,lat4, col = "#2E86C1", bg = "#2E86C1", pch = 24)
points(long5,lat5, col = "#2874A6", bg = "#2874A6", pch = 24)
points(long6,lat6, col = "#21618C", bg = "#21618C", pch = 24)
points(long7,lat7, col = "#1B4F72", bg = "#1B4F72", pch = 24)

### ================================== ###
### SNP Parameters                     ###
### ================================== ###

# Plotting parameters calculated for SNP filtering including:
# Coverage Depth, Individual Missingness, Average Site Depth
# Site Missingness

### ================================== ###
### Coverage Depth                     ###
### ================================== ###

# Load libraries

library(ggplot2)    # for plotting
library(readxl)     # for reading excel files
                    # read in excel file: coverage
                    # located in: Documents/University/Honours 2019/Everything/
                    # /coverage/coverage.xlsx

ggplot (coverage, aes(coverage$`Coverage Depth`)) + 
  geom_histogram(breaks=seq(0,8, by =1), col="grey1", aes(fill=..count..)) 
+ scale_fill_gradient("Count", high = "paleturquoise4", low = "paleturquoise1") 
+ labs(x="Coverage Depth", y = "Frequency") + theme_bw()  # plot a histogram
                                                          # spanning from 0 - 8
                                                          # with bin size = 1

### ================================== ###
### Average Site Depth                 ###
### ================================== ###

# Load libraries 

library(ggplot2)    # for plotting
library(readxl)     # for reading excel files
                    # read in excel file: coverage
                    # located in: Documents/University/Honours 2019/Everything/
                    # /snp-filtering/mean site depth.xlsx

depth <- mean_site_depth_trimmed$`Mean Site Depth`
hist(depth, breaks=100, main = "Average Site Depth",    
     xlab = "Average Site Depth")                          # plot a histogram
                                                           # with bin size = 100

mean <- mean(depth)     # calculate the average depth
SD <- sd(depth)         # calculate the standard deviation
sdevs <- sd * 3         # calculate 3 standard deviations
cutoff <- sdevs + mean  # the cuttof i used was 3 standard deviations 
                        # from the mean. Everything above this has a  
                        # high chance of being errors, e.g. paralogs

abline(v = mean, col = "red", lwd = 2)                     # mark the mean 
abline(v = cutoff, col = "blue", lwd = 2)                  # mark 3 standard 
                                                           # from the mean

### ================================== ###
### Individual Missingness             ###
### ================================== ###

# Load libraries

library(ggplot2)    # for plotting
library(readxl)     # for reading excel files
                    # read in excel file: coverage
                    # located in: Documents/University/Honours 2019/Everything/
                    # /snp-filtering/missingness.xlsx

ggplot (missingness, aes(missingness)) + 
  geom_histogram(breaks=seq(0,0.8, by =0.1), col="grey1", aes(fill=..count..)) 
+ scale_fill_gradient("Count", high = "red4", low = "palevioletred2") 
+ labs(x="Individual Missingness", y = "Frequency") + theme_bw()

### ================================== ###
### Site Missingness                   ###
### ================================== ###

# Load libraries

library(ggplot2)    # for plotting
library(readxl)     # for reading excel files
                    # read in excel file: coverage
                    # located in: Documents/University/Honours 2019/Everything/
                    # /snp-filtering/site-missingness.xlsx
                    # this is a subset of 100,000 estimates of site-missingness 
                    # taken from the site-missingness estimates across all sites

ggplot (site_missingness, aes(site_missingness$`site-missingness`)) 
+ geom_histogram(breaks=seq(0,1, by =0.1), col="grey1", aes(fill=..count..)) 
+ scale_fill_gradient("Count", high = "goldenrod2", low = "lightgoldenrod2") 
+ labs(x="Site Missingness", y = "Frequency") + theme_bw()

### ================================== ###
### Data Analysis                      ###
### ================================== ###

### ================================== ###
### PCA                                ###
### ================================== ###

# Load libraries

library(ggplot2)    # for plotting
library(readxl)     # for reading excel files 
library(RcppCNPy)   # for reading numpy files 

covmatrix <- npyLoad("/Users/Elroy/downloads/PCA.cov.npy") 
Eigenvalues <- eigen(covmatrix)$values          # calculate eigenvalues
Eigenvectors <- eigen(covmatrix)$vectors        # calculate eigenvectors
PC <- as.matrix(covmatrix) %*% Eigenvectors     # calculate the principal components
                                                # by multiplying the covariance matrix 
                                                # against the eigenvectors 

print(round(Eigenvalues/sum(Eigenvalues) * 100, digits = 2))
variance <- (round(Eigenvalues/sum(Eigenvalues) * 100, digits = 2)) 
barplot(variance) # make shift scree plot. bar plot of the proportion of variance
                  # explained by each of the principal components 
round(cumsum(Eigenvalues)/sum(Eigenvalues) * 100, digits = 2)   # display the 
                                                                # cumulative 
                                                                # variance explained
                                                                # by the principal
                                                                # components 

PC1 <- PC[,1]
PC2 <- PC[,2]
PC3 <- PC[,3]

PCA <- data.frame(PC1, PC2, PC3)          # create a dataframe of the principal
                                          # components. move this into an excel 
                                          # file containing the sample metadata
                                          # for ease of plotting 

ggplot(no_outliers_PCA_metadata, aes(x=PC1, y=PC2)) 
+ geom_point(aes(shape = Period, color = STATE), position = "jitter", size = 2.5) 
+ scale_shape_manual(values = c(3, 4, 8, 15, 17, 18, 19, 25)) 
+ theme_bw()                              # plot principal components with 
                                          # the collection date and geographic
                                          # location (state) metdata 

ggplot (no_outliers_PCA_metadata, aes(x=PC1, y=PC2)) 
+ geom_point(aes(shape = Period, color = MISSINGNESS), position = "jitter", size = 2.5) 
+ scale_shape_manual(values = c(3, 4, 8, 15, 17, 18, 19, 25)) 
+ scale_color_continuous(high = "black", low = "turquoise2") 
+ theme_bw()                              # plot the principal components with
                                          # the collection date and individual 
                                          # sample missingness metadata 

ggplot(no_outliers_PCA_metadata, aes(x=LONGITUDE, y=LATITUDE, color = PC1)) 
+ geom_point( size = 2.5) + scale_color_gradient2() # plot the latitude and 
                                                    # longitude against variation 
                                                    # explained by PC1. Reveals 
                                                    # geographic structure 

a <- ggplot(no_outliers_PCA_metadata, aes(PC1,PC2)) 
+ geom_point(color = 'cadetblue3') + theme_bw()

a + geom_text(aes(label = rownames(no_outliers_PCA_metadata)),
              size = 3)                       # identify outlier points using 
                                              # rownames

require("ggrepel")

a + geom_text_repel(aes(label = rownames(no_outliers_PCA_metadata)), 
                    size = 3)             # use "ggrepel" to make labels clearer
                                          # NOTE: this package is so far 
                                          # unavailable using R 3.6

### ================================== ###
### Linear Regression                  ###
### ================================== ###

library(ggplot2)    # for plotting
library(readxl)     # for reading excel files 
library(tidyverse)  # for correlation matrix

regression <- lm(PC1~MISSINGNESS, 
                 data = no_outliers_PCA_metadata)  # linearegression model
                                                   # of PC1 and individual 
                                                   # missingness 
summary(regression)                                # print values

cPC1 <- no_outliers_PCA_metadata$PC1
cMiss <- no_outliers_PCA_metadata$MISSINGNESS
cor.test(cMiss, cPC1, method=c("pearson"))         # calculate a correlation
                                                   # co-efficient. other methods
                                                   # include "spearman" and 
                                                   # "kendall" 

cor(chr27_sorted %>% select(PC1, PC2, Year, LATITUDE, LONGITUDE), 
    use="pairwise.complete.obs")                   # perform a correlation 
                                                   # matrix on the variables PC1,
                                                   # PC2, Year, Lat, Long


ggplot(no_outliers_PCA_metadata, aes(y=PC1, x=Missingness)) # plot linear regression
+ geom_point() + geom_smooth(method="lm") + theme_bw()

### ================================== ###
### Inbreeding                         ###
### ================================== ###

ggplot (inbreedvalues, aes(inbreed)) 
+ geom_histogram(col="grey1", fill="lightblue1") + theme_bw()

ggplot(no_outliers_PCA_metadata, aes(MISSINGNESS,INBREEDING)) 
+ geom_point(color = 'cadetblue3') + theme_bw()

ggplot(no_outliers_PCA_metadata, aes(Year,INBREEDING)) 
+ geom_point(color = 'firebrick2') + theme_bw()


### ================================== ###
### NGSadmix - ADMIXTURE               ###
### ================================== ###

pop <- read.table("/Users/Elroy/Documents/University
                  /Honours 2019/Everything/Data Analysis/
                  ADMIX/Pop-Data.txt", 
                  fill = TRUE, header = FALSE)          # fill in blanks = TRUE
                                                        # there is no header

barplot(t(q2)[,ord],col=2:10, space=0, border=NA,
        xlab="Individuals",                             # create a barplot of   
        ylab="Admixture Proportions (K=2)")             # the individual 
                                                        # admixture proportions

# need to fix this. currently, the labels are crowded
# and unreadable 

text(tapply(1:nrow(pop), pop[ord,1], mean),-0.05,       # add individual sample
     unique(pop[ord,1]),xpd=T)                          # labels on x axis


abline(v=cumsum(sapply(unique(pop[ord,1]), 
                       function(x){sum(pop[ord,1]==x)})),
       col=1,lwd=1.2)                                   # add lines between each
                                                        # individual







