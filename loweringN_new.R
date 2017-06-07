# Combination of the means plots (Fig 6)

# first load in the first bit of loweringN.R, then instead of the plot commands there, do these:
# for now, quick and ugly without ggplot 
M <- max(meanES1,meanES55,meanES25,meanES79)
plot(c(10000,100000),c(0,M),type="n",xlab="cut-off (x 1000)",ylab="mean p-value",xaxt="n")
axis(1,c(0,10000,25000,50000,75000,100000),c("0","10","25","50","75","100"))
lines(n,meanES1,type="l",col="dark blue",lwd=2)
points(n,meanES1,type="p",col="dark blue",lwd=2,pch=19)

lines(n,meanES55,type="l", col="dark green",lwd=2)
points(n,meanES55,type="p",col="dark green",lwd=2,pch=19)

lines(n,meanES25,type="l",col="light blue",lwd=2)
points(n,meanES25,type="p",col="light blue",lwd=2,pch=17)

lines(n,meanES79,type="l",col="light green",lwd=2)
points(n,meanES79,type="p",col="light green",lwd=2,pch=17)


M <- max(nrNA1,nrNA55,nrNA25,nrNA79)
plot(c(10000,100000),c(0,M),type="n",xlab="cut-off (x 1000)",ylab="nr excluded values",xaxt="n")
axis(1,c(0,10000,25000,50000,75000,100000),c("0","10","25","50","75","100"))
lines(n,nrNA1,type="l",col="dark blue",lwd=2)
points(n,nrNA1,type="p",col="dark blue",lwd=2,pch=19)

lines(n,nrNA55,type="l", col="dark green",lwd=2)
points(n,nrNA55,type="p",col="dark green",lwd=2,pch=19)

lines(n,nrNA25,type="l",col="light blue",lwd=2)
points(n,nrNA25,type="p",col="light blue",lwd=2,pch=17)

lines(n,nrNA79,type="l",col="light green",lwd=2)
points(n,nrNA79,type="p",col="light green",lwd=2,pch=17)
