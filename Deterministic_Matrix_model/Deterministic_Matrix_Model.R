
setwd("/Users/hannajackson/Dropbox/BISC-838-pop-models/HW1_Deterministic")
rm(list=ls())

require(popbio)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~ Set up the population matrix ~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
stages <- c("Egg/hatchling",
            "Small Juvenile",
            "Large Juvenile",
            "Subadult",
            "Early Breeders",
            "Remigrants",
            "Mature")

myvector <- c(0     ,0     ,0     ,0     ,127   ,4     ,80,
              0.6747,0.7370,0     ,0     ,0     ,0     ,0,
              0     ,0.0486,0.6610,0     ,0     ,0     ,0,
              0     ,0     ,0.0147,0.6907,0     ,0     ,0,
              0     ,0     ,0     ,0.0518,0     ,0     ,0,
              0     ,0     ,0     ,0     ,0.8091,0     ,0,
              0     ,0     ,0     ,0     ,0     ,0.8091,0.8089)

## Make sure that's correct
length(myvector) == 7*7

## Construct our matrix 
M <- matrix(myvector,
            byrow = TRUE,
            nrow  = 7,
            dimnames = list(stages,stages))

## Look at the eigenvalues and eigenvectors 
eigens <- eigen(M)

## The leading eigenvalue is associated with the first class - 0.945 
print(lambda <- Re(eigens$values)[which.max(Re(eigens$values))])





## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~ Get and plot the eigen analysis: ~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Run the analysis two ways, we'll need both for plotting later 
out1 <- eigen.analysis(M, zero = FALSE)
out2 <- eigen.analysis(M, zero = TRUE)

## The lambda (0.945031) we just got matches the output from
##    eigen.analysis in the previous section 
lambda
out1$lambda1


## ~~~~~~~~~~ Plot matricies to PDF via image2: ~~~~~~~~~~
pdf(file = "Matrix.fig.pdf", width = 8, height = 8)
layout(matrix(1:4, nrow = 2,byrow = TRUE))
## Set up the plotting options 
par(oma = c(0.1,0.1,0.1, 0.1), 
    mgp = c(2,0.2,0),
    tcl = 0, 
    cex.axis = 1.3, 
    pty = 's')

## Getting the margins to line up ended up being tricky for this, I
##    ended up setting this parameter 'ma' that then gets used in the mar
##    argument of each pannel: 
ma <- 8

## Plot each of the matricies
## ~~~ Pannel 1: ~~~
image2(M,cex=0.9, border="black", round = 2,
       mar = c(0,ma,ma,0))
mtext("Projection matrix", line=2, side=3, cex=1.5)

## ~~~ Pannel 2: ~~~
image2(out1$elasticities,
       cex=0.9, border="black", round = 2,
       labels=3,mar = c(0,0,ma,0))
mtext("Elasticity matrix", line=2, side=3, cex=1.5)

## ~~~ Pannel 3: ~~~
image2(out1$sensitivities,
       cex=0.9, border="black", round = 2,
       labels=2,mar = c(0,ma,0,0))
mtext("Sensitivity matrix all elements", line=1, side=3, cex=1.5)

## ~~~ Pannel 4: ~~~
image2(out2$sensitivities,
       cex=0.9, border="black", round = 2, labels=NA,mar = c(ma/2,0,ma/2,0))
mtext("Sensitivity matrix", line=1,side=3,cex=1.5)

dev.off() 



## ~~~~~~~~~~ Plot Elasticities and Sensitivities vs Stage ~~~~~~~~~~

## Get the data read for plotting - extract only the parameters of
## interest for both sensitivities and elasticities

## Set up a data structure to hold only the parameters we need and
## then fill it
sens.plot <- matrix(NA, nrow=7,ncol=3, dimnames=list(stages,c("P","G","F")))

sens.plot[1:7,"P"] <- out1$sensitivities[cbind(1:7,1:7)]
sens.plot[1:6,"G"] <- out1$sensitivities[cbind(2:7,1:6)]
sens.plot[2:7,"F"] <- out1$sensitivities[cbind(rep(1,6),2:7)]

elast.plot <- matrix(NA, nrow=7,ncol=3, dimnames=list(stages,c("P","G","F")))

elast.plot[1:7,"P"] <- out1$elasticities[cbind(1:7,1:7)]
elast.plot[1:6,"G"] <- out1$elasticities[cbind(2:7,1:6)]
elast.plot[2:7,"F"] <- out1$elasticities[cbind(rep(1,6),2:7)]


## Plot those to PDF:  
pdf(file="Elast.sens.pdf", width=6, height=8)

## Setup plotting options 
layout(matrix(1:2,ncol=1))
cols <- c("blue","purple","forestgreen")
par(oma = c(7,5,2, 2),
    mar = c(c(0, 0, 0.2, 0)),
    mgp = c(2,0.2,0),
    tcl = 0, 
    cex.axis = 1.3,
    lwd=2.5)


## ~~~ Pannel 1: Sensitivities through the stages ~~~
plot(NA,
     xlim = c(1,7),
     ylim = c(0,max(sens.plot,na.rm=TRUE)),
     xlab = "",
     ylab = "",
     xaxt = 'n',
     las  = 1)
mtext(text="Sensitivities",side=2,line=3,cex=1.3)
lines(sens.plot[,"P"], col=cols[2]) ## purple 
lines(sens.plot[,"G"], col=cols[1]) ## blue
lines(sens.plot[,"F"], col=cols[3]) ## green 
legend("topright",
       legend = c("Transitioning to next stage",
                  "Remaining in same stage",
                  "Reproductive output"), 
       col = cols,
       cex = 1.2,
       pch = 15, 
       bty = "n",
       ncol = 1)


## ~~~ Pannel 2: Elasticities through the stages ~~~ 
plot(NA,
     xlim = c(1,7),
     ylim = c(0,max(elast.plot,na.rm=TRUE)),
     xlab = "",
     ylab = "",
     xaxt = 'n',
     las = 1)
mtext(text = "Stage",side=1,line=5.5,cex=1.3)
mtext(text = "Elasticities",side=2,line=3,cex=1.3)
text(x = 1:7+0.1,
     y = par('usr')[3],
     labels = stages,
     srt = 35,
     adj = c(1.1,1.1),
     xpd = NA,
     cex = 1.2)
lines(elast.plot[,"P"], col = cols[2])
lines(elast.plot[,"G"], col = cols[1])
lines(elast.plot[,"F"], col = cols[3])

## Save the PDF 
dev.off()





## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~ Make and plot the population projection ~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Plot projection of the population through time
nsim <- 40
## my.simulation <- pop.projection(M, rep(1,7),nsim)
my.simulation <- pop.projection(M, rep(20,7),nsim)

## Specify colours using RColorBrewer

## install.packages("RColorBrewer")
library(RColorBrewer)
cols <- brewer.pal(7,"Dark2")

## Extract just the population sizes 
sizes <- my.simulation$stage.vectors


## Actually plot it 
pdf(file = "Reaching_Stable_age.pdf", width=7, height=8.3)
layout(matrix(c(1,1,2,2,3,4),nrow=3,byrow=TRUE))

par(oma = c(5.5,2,0.1, 2),
    mar = c(0.2, 3, 1, 3),
    mgp = c(2,0.2,0),
    tcl = 0, 
    cex.axis = 1.3)


## ~~~ Pannel 1: Plot the number of individuals vs time ~~~ 

## Make an empty plot of the correct size 
plot(NA,
     xlim = c(0,nsim), ylim = c(0,max(sizes)),
     xlab = "", ylab = "",
     las = 1, xaxt = "n")
mtext("Number of individuals", side=2, line=3)

## If I plotted it as it is, all the classes with low values would
##   overlap and you would only be able to see the last one plotted. So
##   I'm going to make those classes dotted lines, and I'll jitter
##   their x axes systematically so that the dotted lines don't also overlap
lines  <- c(1,1,1,3,3,3,3) ## line types 
x.vals <- matrix(NA, ncol=40, nrow=7)
x.vals[1,] <- 0:39
x.vals[2,] <- 0:39
x.vals[3,] <- 0:39
x.vals[4,] <- 0:39+(1/4)
x.vals[5,] <- 0:39+(2/4)
x.vals[6,] <- 0:39+(3/4)
x.vals[7,] <- 0:39+(4/4)

## plot each row of the matrix we created as a line on that empty plot
for(ii in 1:7){
  lines(x=x.vals[ii,], y=sizes[ii,], col=cols[ii],lwd=3, lty=lines[ii])
}
legend("topright",
       legend = stages, 
       col    = cols,
       cex    = 1.2,
       lwd    = 6, 
       bty    = "n",
       ncol   = 1)

## Adding a 3rd axis on the right side for the total population
## Need to scale it to the scale of the current plot: 
scaled.sizes <- colSums(sizes) * max(sizes) / max(colSums(sizes))
lines(scaled.sizes, col="grey50",lwd=4,lty=2)

## Then also need to scale the axis in the same way: 
desired.axis <- c(0,100,200,300,400)
axis(side     = 4,
     at       = desired.axis * max(sizes)/max(desired.axis),
     labels   = desired.axis,
     col.axis = "grey40",
     las      =1 )
## And add the label 
mtext("Total number of individuals", side=4, line=3,col='grey40')


## ~~~ Pannel 2: Plot relative pop sizes ~~~

## Making x.vals for this plot too 
x.vals <- matrix(NA, ncol=40, nrow=7)
x.vals[1,] <- 0:39
x.vals[2,] <- 0:39
x.vals[3,] <- 0:39
x.vals[4,] <- 0:39+(1/8)
x.vals[5,] <- 0:39+(2/8)
x.vals[6,] <- 0:39+(3/8)
x.vals[7,] <- 0:39+(4/8)

## First make the object that tells us what proportion of the
##    individuals at any one time are in each class

## Start by just duplicating 'sizes' and then modify iteratively 
prop.sizes <- sizes 
for (ii in 1:nsim){
  prop.sizes[,ii] <- sizes[,ii] / sum(sizes[,ii])
}

## Plot that 
plot(NA,
     xlim = c(0,nsim),
     ylim = c(0,max(prop.sizes)),
     xlab = "",
     ylab = "",
     las  = 1)
mtext("Proportion of individuals", side=2, line=3)
mtext("Time step", side=1, line=1.5)
for(ii in 1:7){
  lines(x=x.vals[ii,], y=prop.sizes[ii,], col=cols[ii],lwd=3,lty=lines[ii])
}


## ~~~ Pannel 3: Plot the Stable Age Distribution ~~~
barplot(out1$stable.stage+0.01,
        col=cols, ylim=c(0,1),las=1, xaxt='n')
mtext("Stable stage proportion", side=2, line=3)
text(x = c(0.8,2,3.2,4.4,5.6,6.7,8),
     y = par('usr')[3],
     labels = stages,
     srt = 35,
     adj = c(1.1,1.1),
     xpd = NA,
     cex = 1.2)


## ~~~ Pannel 4: Plot the Reproductive Values 
barplot(out1$repro.value+0.01,
        col=cols,las=1, xaxt='n', yaxt='n',
        ylim=c(0,max(out1$repro.value+30)))
axis(side=2,
     labels=c(0,100,200,300,400,500),
     at=c(0,100,200,300,400,500))
mtext("Reproductive Value", side=2, line=3)
text(x = c(0.8,2,3.2,4.4,5.6,6.7,8),
     y = par('usr')[3],
     labels = stages,
     srt = 35,
     adj = c(1.1,1.1),
     xpd = NA,
     cex = 1.2)

dev.off()
