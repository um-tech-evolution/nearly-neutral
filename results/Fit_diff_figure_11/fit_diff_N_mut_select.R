# R script to produce plot fit diffs for N from 25 to 6400 and for various N_mut values.
# Author:  Alden Wright  2/28/2018  
#setwd("../../continuous_variable_evolution/data/1_17_18")
#cv = read.csv("N6400_na1_N_mut_bi3_0_fs1_ntrials1000.csv",comment.char="#")  # these are the "old plots"
cv = read.csv("N6400_na1_N_mut_bi3_0_fs1_ntrials1000combined.csv",comment.char="#")
sub_cv=subset(cv,N>=25,select=c(N,N_mut,fit_diff_neg_fract,fit_diff_neg_neutral,fit_diff_pos_neutral,fit_diff_pos_fract))
agg_cv=aggregate(sub_cv,by=list(sub_cv$N,sub_cv$N_mut),FUN=mean)
colors= c("green","red","darkblue","coral4","cyan","orange","brown")
pchars = c(15,16,17,18,19,20,14)
labels=c(
  "fit diff neg neutral",
  "fit diff pos neutral",
  "fit diff neg fract",
  "fit diff pos fract"
)
linetypes = c(1,2,3,4,5,6,7,8,9)
Nmut = c( 1.0, 0.5, 0.25, 0.125)
#dev.off()  # uncomment if a graphics window is already open
plot(agg_cv$N[agg_cv$N_mut==Nmut[1]],agg_cv$fit_diff_neg_neutral[agg_cv$N_mut==Nmut[1]],
	log="x",xaxt="n",pch=pchars[1],col=colors[1],xlab="N",ylab="fraction",ylim=c(0.0,0.5),xlim=c(25.0,6400.0))
axis(1,at=c(25,50,100,200,400,800,1600,3200,6400))
lines(agg_cv$N[agg_cv$N_mut==Nmut[1]],agg_cv$fit_diff_neg_neutral[agg_cv$N_mut==Nmut[1]],col=colors[1],lty=linetypes[1])
points(agg_cv$N[agg_cv$N_mut==Nmut[1]],agg_cv$fit_diff_pos_neutral[agg_cv$N_mut==Nmut[1]],col=colors[2],pch=pchars[2])
lines(agg_cv$N[agg_cv$N_mut==Nmut[1]],agg_cv$fit_diff_pos_neutral[agg_cv$N_mut==Nmut[1]],col=colors[2],lty=linetypes[2])
points(agg_cv$N[agg_cv$N_mut==Nmut[1]],agg_cv$fit_diff_neg_fract[agg_cv$N_mut==Nmut[1]],col=colors[4],pch=pchars[4])
lines(agg_cv$N[agg_cv$N_mut==Nmut[1]],agg_cv$fit_diff_neg_fract[agg_cv$N_mut==Nmut[1]],col=colors[4],lty=linetypes[4])
points(agg_cv$N[agg_cv$N_mut==Nmut[1]],agg_cv$fit_diff_pos_fract[agg_cv$N_mut==Nmut[1]],col=colors[3],pch=pchars[3])
lines(agg_cv$N[agg_cv$N_mut==Nmut[1]],agg_cv$fit_diff_pos_fract[agg_cv$N_mut==Nmut[1]],col=colors[3],lty=linetypes[3])
#title("Non-neutral Fitness Differences NMut=1")
legend(800,0.12,labels,col=colors,pch=pchars,lty=linetypes,bty="n")
dev.copy(pdf,"fit_diff_NMut=1.pdf")    # copies what is on the screen to a file.
dev.off()
dev.copy(png,"fit_diff_NMut=1.png")    # copies what is on the screen to a file.
dev.off()
dev.copy(tiff,"fit_diff_NMut=1.tiff")    # copies what is on the screen to a file.
dev.off()
