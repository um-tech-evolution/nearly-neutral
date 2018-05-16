# R script to produce plot the mean for N from 25 to 6400 and for various N_mut values.
# Author:  Alden Wright  2/28/2018  
#setwd("../../continuous_variable_evolution/data/1_17_18")
#cv = read.csv("N6400_na1_N_mut_bi3_0_fs1_ntrials1000.csv",comment.char="#")  # this constructs the "old plot"
cv = read.csv("N6400_na1_N_mut_bi3_0_fs1_ntrials1000combined.csv",comment.char="#")
sub_cv=subset(cv,N>=25,select=c(N,N_mut,fitness_mean))
agg_cv=aggregate(sub_cv,by=list(sub_cv$N,sub_cv$N_mut),FUN=mean)
colors= c("green","red","darkblue","coral4","cyan","orange","brown")
pchars = c(15,16,17,18,19,20,14)
labels=c(
"Nmut = 0.125",
"Nmut = 0.25",
"Nmut = 0.5",
"Nmut = 1.0"
)
linetypes = c(1,2,3,4,5,6,7,8,9)
Nmut = c(0.125,0.25,0.5,1.0)
#dev.off()  # uncomment if a graphics window is already open
i=1
plot(agg_cv$N[agg_cv$N_mut==Nmut[i]],agg_cv$fitness_mean[agg_cv$N_mut==Nmut[i]],
	log="x",xaxt="n",pch=pchars[1],col=colors[1],xlab="N",ylab="fitness",ylim=c(0.9,1.0),xlim=c(25.0,6400.0))
axis(1,at=c(25,50,100,200,400,800,1600,3200,6400))
lines(agg_cv$N[agg_cv$N_mut==Nmut[i]],agg_cv$fitness_mean[agg_cv$N_mut==Nmut[i]],col=colors[1],lty=linetypes[1])
for( i in 2:4){
  points(agg_cv$N[agg_cv$N_mut==Nmut[i]],agg_cv$fitness_mean[agg_cv$N_mut==Nmut[i]],col=colors[i],pch=pchars[i])
  lines(agg_cv$N[agg_cv$N_mut==Nmut[i]],agg_cv$fitness_mean[agg_cv$N_mut==Nmut[i]],col=colors[i],lty=linetypes[i])
}
title("Non-neutral Fitness Mean with burn-in=3.0")
legend(1200,0.94,labels,col=colors,pch=pchars,lty=linetypes,bty="n")
dev.copy(pdf,"non_neutral_fit_mean_N_mut_bi3_0.pdf")    # copies what is on the screen to a file.
dev.off()
dev.copy(png,"non_neutral_fit_mean_N_mut_bi3_0.png")    # copies what is on the screen to a file.
dev.off()
dev.copy(tiff,"non_neutral_fit_mean_N_mut_bi3_0.tiff")    # copies what is on the screen to a file.
dev.off()
