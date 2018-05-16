Objective:  Create plots of the gamma distributions that we are using
for the nearly neutral paper.

The following succeeds at this, but raises questions about the gamma
distribution parameters that we are using.  These are:

# dfe_adv_prob=0.01
# dfe_adv_alpha=1.0
# dfe_adv_theta=0.01
# dfe_disadv_alpha=0.2
# dfe_disadv_theta=0.5

with alternatives:
# dfe_adv_theta=0.02
# dfe_adv_theta=0.005

Reference website:
http://www.statmethods.net/advgraphs/probability.html

# The following lines can be pasted directly into R to create the gamma PDFs plot.

dev.off()
posx = seq(0.00,.04,length=200)
negx = seq(-0.04,.00,length=200)
nnegx = seq(0.04,0.0,length=200)
colors= c("red","darkblue","coral4","cyan","orange","brown")
shapes = c(1.0,1.0,1.0,0.2,0.2)
scales = c(0.01,0.02,0.005,0.5,1.0)
linetypes = c(2,4,3,5,6)
labels=c(
    expression(paste(alpha,"=0.2, ",beta,"=0.5 Deleterious")),
    expression(paste(alpha,"=0.2, ",beta,"=1.0 Deleterious")),
    expression(paste(alpha,"=1, ",beta,"=0.01 Mixed")),
    expression(paste(alpha,"=1, ",beta,"=0.02 Mixed")),
    expression(paste(alpha,"=1, ",beta,"=0.005 Mixed")))
hx = dgamma(posx,shape=shapes[1],scale=scales[1])
plot(posx,hx,type="l",lty=linetypes[1],log="",col=colors[1],xlab="Selection Coefficient",ylab="Density",ylim=c(0.1,100.0),xlim=c(-0.04,0.04),yaxt="n",main="gamma probability density functions")
axis(2,at=c(),pos=0.0)
for (i in 2:3){lines(posx, dgamma(posx,shape=shapes[i],scale=scales[i]), lwd=2, col=colors[i],lty=linetypes[i])}
for (i in 4:5){lines(negx, dgamma(nnegx,shape=shapes[i],scale=scales[i]), lwd=2, col=colors[i], lty=linetypes[i])}
legend(-0.043,80,inset=0.5,lty=linetypes,labels,lwd=2,col=colors,cex=0.85,bty="n")
dev.copy(png,'gamma_dists_2sided.png')    # copies what is on the screen to a file.
dev.off()
dev.copy(tiff,'gamma_dists_2sided.tiff')    # copies what is on the screen to a file.
dev.off()
dev.copy(pdf,'gamma_dists_2sided.pdf')    # copies what is on the screen to a file.
dev.off()


