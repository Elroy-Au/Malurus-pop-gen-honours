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

ggplot(no_outliers_metadata, aes(x=PC1, y=PC2)) 
+ geom_point(aes(shape = Period, color = STATE), position = "jitter", size = 2.5) 
+ scale_shape_manual(values = c(3, 4, 8, 15, 17, 18, 19, 25)) 
+ theme_bw()                              # plot principal components with 
                                          # the collection date and geographic
                                          # location (state) metdata 

ggplot (no_outliers_metadata, aes(x=PC1, y=PC2)) 
+ geom_point(aes(shape = Period, color = MISSINGNESS), position = "jitter", 
             size = 2.5) + scale_shape_manual(values = 
                                                c(3, 4, 8, 15, 17, 18, 19, 25)) 
+ scale_color_continuous(high = "black", low = "turquoise2") + theme_bw()                              
                                          # plot the principal components with
                                          # the collection date and individual 
                                          # sample missingness metadata 

ggplot(no_outliers_metadata, aes(x=LONGITUDE, y=LATITUDE, color = PC1)) 
+ geom_point(size = 2.5) + scale_color_gradient2()  # plot the latitude and 
                                                    # longitude against variation 
                                                    # explained by PC1. Reveals 
                                                    # geographic structure 

# Basic PCA / scatterplot

a <- ggplot(no_outliers_metadata, aes(PC1,PC2)) 
+ geom_point(color = 'cadetblue3') + theme_bw()

# Identify 

a + geom_text(aes(label = rownames(no_outliers_metadata)),
              size = 3)                       # identify outlier points using 
                                              # rownames

library(ggrepel)

a + geom_text_repel(aes(label = rownames(no_outliers_metadata)), 
                    size = 3)             # use "ggrepel" to make labels clearer

### ================================== ###
### Linear Regression                  ###
### ================================== ###

library(ggplot2)    # for plotting
library(readxl)     # for reading excel files 
library(tidyverse)  # for correlation matrix

regression <- lm(PC1~MISSINGNESS, 
                 data = no_outliers_metadata)  # linearegression model
                                                   # of PC1 and individual 
                                                   # missingness 
summary(regression)                                # print values

cPC1 <- no_outliers_metadata$PC1
cMiss <- no_outliers_metadata$MISSINGNESS
cor.test(cMiss, cPC1, method=c("pearson"))         # calculate a correlation
                                                   # co-efficient. other methods
                                                   # include "spearman" and 
                                                   # "kendall" 

cor(chr27_sorted %>% select(PC1, PC2, Year, LATITUDE, LONGITUDE), 
    use="pairwise.complete.obs")                   # perform a correlation 
                                                   # matrix on the variables PC1,
                                                   # PC2, Year, Lat, Long


ggplot(no_outliers_metadata, aes(y=PC1, x=Missingness)) # plot linear regression
+ geom_point() + geom_smooth(method="lm") + theme_bw()

### ================================== ###
### Inbreeding                         ###
### ================================== ###

ggplot (inbreedvalues, aes(inbreed)) 
+ geom_histogram(col="grey1", fill="lightblue1") + theme_bw()

ggplot(no_outliers_metadata, aes(MISSINGNESS,INBREEDING)) 
+ geom_point(color = 'cadetblue3') + theme_bw()

ggplot(no_outliers_metadata, aes(Year,INBREEDING)) 
+ geom_point(color = 'firebrick2') + theme_bw()


### ================================== ###
### NGSadmix - ADMIXTURE               ###
### ================================== ###

# load files

pop <- read.table("/Users/Elroy/Documents/University
                  /Honours 2019/Everything/Data Analysis/
                  ADMIX/Pop-Data.txt", 
                  fill = TRUE, header = FALSE)          # fill in blanks = TRUE
                                                        # there is no header

q1 <- read.table("/Users/Elroy/Documents/University/
                 Honours 2019/Everything/Data Analysis/
                 ADMIX/admix.1.txt", fill = TRUE, header = FALSE)

q2 <- read.table("/Users/Elroy/Documents/University/
                 Honours 2019/Everything/Data Analysis/
                 ADMIX/admix.2.txt", fill = TRUE, header = FALSE)

q3 <- read.table("/Users/Elroy/Documents/University/
                 Honours 2019/Everything/Data Analysis/
                 ADMIX/admix.3.txt", fill = TRUE, header = FALSE)

q4 <- read.table("/Users/Elroy/Documents/University/
                 Honours 2019/Everything/Data Analysis/
                 ADMIX/admix.4.txt", fill = TRUE, header = FALSE)

q5 <- read.table("/Users/Elroy/Documents/University/
                 Honours 2019/Everything/Data Analysis/
                 ADMIX/admix.5.txt", fill = TRUE, header = FALSE)

q6 <- read.table("/Users/Elroy/Documents/University/
                 Honours 2019/Everything/Data Analysis/
                 ADMIX/admix.6.txt", fill = TRUE, header = FALSE)

ord <- order (pop[,2])  # order according to population 

# plot admixture

barplot(t(q2)[,ord],col=2:10, space=0, border=NA,
        xlab="Individuals",                             # create a barplot of   
        ylab="Admixture Proportions (K=2)")             # the individual 
                                                        # admixture proportions

barplot(t(q3)[,ord],col=brewer.pal(n=3, name="RdBu"),   # choose colours using
        space=0, border=NA, xlab="Individuals",         # RColourBrewer
        ylab="Admixture Proportions (K=3)")

text(tapply(1:nrow(pop), pop[ord,2], mean),-0.05,       # add individual 
     unique(pop[ord,2]),xpd=T)                          # population labels


abline(v=cumsum(sapply(unique(pop[ord,1]), 
                       function(x){sum(pop[ord,1]==x)})),
       col="white",lwd=0.5)                             # add lines between each
                                                        # individual

abline(v=cumsum(sapply(unique(pop[ord,2]), 
                       function(x){sum(pop[ord,2]==x)})),
       col=1,lwd=2)                                     # add lines between each
                                                        # population

### ================================== ###
### ADMIXTURE SCATTER PLOTS            ###
### ================================== ###

# load libraries

library(ggplot2)        # plotting
library(ggrepel)        # ggplot add-on: text repel
library(scatterpie)     # ggplot add-on: pie charts
library(sf)             # read shapefile

# load shapefile of Australia

australia <- st_read("/Users/Elroy/Documents/University/Honours 2019/
                     Everything/Data Analysis/Aus Shapefile/australia.shp")

# create a vector plot of Australia 

a <- ggplot() + geom_sf(data = australia, color = "black", fill = "white") 

### K=2

PC1 <- admix_geodata$PC1
PC2 <- admix_geodata$PC2
Lat <- admix_geodata$LATITUDE
Long <- admix_geodata$LONGITUDE
K1 <- admix_geodata$K2.1
K2 <- admix_geodata$K2.2

data <- data.frame(PC1,PC2,Lat,Long,K1,K2)

ggplot(data, aes(PC1, PC2)) 
+ geom_scatterpie(aes(x = PC1, y = PC2), data = data, cols = c("K1", "K2")) 
+ scale_fill_manual(values=c("#F39C12", "#C39BD3", 
                             "#82E0AA", "#F7DC6F", "#EC7063", "#85C1E9")) 
+ theme_light()

a + geom_scatterpie(aes(x = Long, y = Lat), data = data, 
                    cols = c("K1", "K2")) 
+ scale_fill_manual(values=c("#F39C12", "#C39BD3", "#82E0AA", 
                             "#F7DC6F", "#EC7063", "#85C1E9")) 
+ theme_light()

### K=3

PC1 <- admix_geodata$PC1
PC2 <- admix_geodata$PC2
Lat <- admix_geodata$LATITUDE
Long <- admix_geodata$LONGITUDE
K1 <- admix_geodata$K3.1
K2 <- admix_geodata$K3.2
K3 <- admix_geodata$K3.3

data <- data.frame(Lat,Long,PC1,PC2,K1,K2,K3)

ggplot(data, aes(PC1, PC2)) 
+ geom_scatterpie(aes(x = PC1, y = PC2), 
                  data = data, cols = c("K1", "K2", "K3")) 
+ scale_fill_manual(values=c("#F39C12", "#C39BD3", "#82E0AA", 
                             "#F7DC6F", "#EC7063", "#85C1E9")) + theme_light()

a + geom_scatterpie(aes(x = Long, y = Lat), 
                    data = data, cols = c("K1", "K2", "K3")) 
+ scale_fill_manual(values=c("#F39C12", "#C39BD3", "#82E0AA", 
                             "#F7DC6F", "#EC7063", "#85C1E9")) + theme_light()

### K=4

PC1 <- admix_geodata$PC1
PC2 <- admix_geodata$PC2
Lat <- admix_geodata$LATITUDE
Long <- admix_geodata$LONGITUDE
K1 <- admix_geodata$K4.1
K2 <- admix_geodata$K4.2
K3 <- admix_geodata$K4.3
K4 <- admix_geodata$K4.4

data <- data.frame(Lat,Long,PC1,PC2,K1,K2,K3,K4)

ggplot(data, aes(PC1, PC2)) 
+ geom_scatterpie(aes(x = PC1, y = PC2), data = data, 
                  cols = c("K1", "K2", "K3", "K4")) 
+ scale_fill_manual(values=c("#F39C12", "#C39BD3", "#82E0AA", 
                             "#F7DC6F", "#EC7063", "#85C1E9")) + theme_light()

a + geom_scatterpie(aes(x = Long, y = Lat), data = data, 
                    cols = c("K1", "K2", "K3", "K4")) 
+ scale_fill_manual(values=c("#F39C12", "#C39BD3", "#82E0AA", 
                             "#F7DC6F", "#EC7063", "#85C1E9")) + theme_light()

### K=5

PC1 <- admix_geodata$PC1
PC2 <- admix_geodata$PC2
Lat <- admix_geodata$LATITUDE
Long <- admix_geodata$LONGITUDE
K1 <- admix_geodata$K5.1
K2 <- admix_geodata$K5.2
K3 <- admix_geodata$K5.3
K4 <- admix_geodata$K5.4
K5 <- admix_geodata$K5.5

data <- data.frame(Lat,Long,PC1,PC2,K1,K2,K3,K4,K5)

ggplot(data, aes(PC1, PC2)) 
+ geom_scatterpie(aes(x = PC1, y = PC2), data = data, 
                  cols = c("K1", "K2", "K3", "K4", "K5")) 
+ scale_fill_manual(values=c("#F39C12", "#C39BD3", "#82E0AA", 
                             "#F7DC6F", "#EC7063", "#85C1E9")) + theme_light()

a + geom_scatterpie(aes(x = Long, y = Lat), data = data, 
                    cols = c("K1", "K2", "K3", "K4", "K5")) 
+ scale_fill_manual(values=c("#F39C12", "#C39BD3", "#82E0AA", 
                             "#F7DC6F", "#EC7063", "#85C1E9")) + theme_light()

### K=6

PC1 <- admix_geodata$PC1
PC2 <- admix_geodata$PC2
Lat <- admix_geodata$LATITUDE
Long <- admix_geodata$LONGITUDE
K1 <- admix_geodata$K1
K2 <- admix_geodata$K2
K3 <- admix_geodata$K3
K4 <- admix_geodata$K4
K5 <- admix_geodata$K5
K6 <- admix_geodata$K6

data <- data.frame(Lat,Long,PC1,PC2,K1,K2,K3,K4,K5,K6)

ggplot(data, aes(PC1, PC2)) 
+ geom_scatterpie(aes(x = PC1, y = PC2), data = data, 
                  cols = c("K1", "K2", "K3", "K4", "K5", "K6")) 
+ scale_fill_manual(values=c("#F39C12", "#C39BD3", "#82E0AA", 
                             "#F7DC6F", "#EC7063", "#85C1E9")) + theme_light()

a + geom_scatterpie(aes(x = Long, y = Lat), data = data, 
                    cols = c("K1", "K2", "K3", "K4", "K5", "K6")) 
+ scale_fill_manual(values=c("#F39C12", "#C39BD3", "#82E0AA", 
                             "#F7DC6F", "#EC7063", "#85C1E9")) + theme_light()

# add labels 
# need to define ID (sample ID) or cID (collection ID) as a value 
# in the dataframe
# too many samples leading to excessive cluttering

ID <- admix_geodata$`COLLECTION ID`

ggplot(data, aes(PC1, PC2, label = ID)) 
+ geom_scatterpie(aes(x = PC1, y = PC2), data = data, 
                  cols = c("K1", "K2", "K3", "K4", "K5", "K6")) 
+ scale_fill_manual(values=c("#F39C12", "#C39BD3", "#82E0AA", 
                             "#F7DC6F", "#EC7063", "#85C1E9")) 
+ geom_label_repel() + theme_light()


### ================================== ###
### MIXED MODELS IN SOMMER             ###
### ================================== ###

# load the covariance matrix 
# generated using PCangsd 

genmatrix <- read.table("/Users/Elroy/Documents/University/Honours 2019/
                Everything/Data Analysis/Mixed Model/covariance-matrix.txt")

# convert into a matrix and check data structure 

M <- as.matrix(genmatrix)
M[1:4,1:4]

###           A37114     A36182     A35583     A35582
### A37114 1.07916844 0.02958575 0.03612941 0.03744407
### A36182 0.02958575 1.10712433 0.03258741 0.02230000
### A35583 0.03612941 0.03258741 1.05975997 0.03946498
### A35582 0.03744407 0.02230000 0.03946498 1.08377147

# specify a single fixed effect 
# random effects are specified using the ~vs(x) function
# which denotes a covariance matrix applied at the random effect x 
# a known covariance matrix is denoted and provided using the function Gu()
# unknown covariance matrices are created using functions like vs(ds())
# rcov is a formula specifying the name of the error term, i.e. rcov= ~ units
# data provides the dataset

mix1 <- mmer(fixed=INBREEDING~Year,
            random=~vs(ID,Gu=M),
            rcov=~units,
            data=inbreedtestdata)

# specify multiple fixed effects
# the operater "+" can be used to include multiple effects (random or fixed)

mix2 <- mmer(fixed=INBREEDING~Year + LATITUDE + LONGITUDE,
            random=~vs(ID,Gu=M),
            rcov=~units,
            data=inbreedtestdata)

# provide a summary of the model results

summary(mix)
