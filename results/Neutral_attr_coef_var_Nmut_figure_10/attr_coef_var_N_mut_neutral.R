# R script to produce plot the attribute coefficient of variation for N from 25 to 6400 and for various N_mut values.
# Author:  Alden Wright  1/18/2018  
#setwd("../../continuous_variable_evolution/data/12_5_17")
cv = read.csv("N6400_na1_N_mut_bi3_0_neutral_ntrials1000.csv",comment.char="#")
sub_cv=subset(cv,N>=25,select=c(N,N_mut,attribute_mean,attribute_coef_var))
agg_cv=aggregate(sub_cv,by=list(sub_cv$N,sub_cv$N_mut),FUN=mean)
colors= c("green","red","darkblue","coral4","cyan","orange","brown")
pchars = c(15,16,17,18,19,20,14)
labels=c(
"N_mut = 0.125",
"N_mut = 0.25",
"N_mut = 0.5",
"N_mut = 1.0"
)
linetypes = c(1,2,3,4,5,6,7,8,9)
Nmut = c(0.125,0.25,0.5,1.0)
#dev.off()  # uncomment if a graphics window is already open
i=1
plot(agg_cv$N[agg_cv$N_mut==Nmut[i]],agg_cv$attribute_coef_var[agg_cv$N_mut==Nmut[i]],
	log="x",xaxt="n",pch=pchars[1],col=colors[1],xlab="N",ylab="attribute coeffcient of variation",ylim=c(0.0,0.2),xlim=c(25.0,6400.0))
axis(1,at=c(25,50,100,200,400,800,1600,3200,6400))
lines(agg_cv$N[agg_cv$N_mut==Nmut[i]],agg_cv$attribute_coef_var[agg_cv$N_mut==Nmut[i]],col=colors[1],lty=linetypes[1])
for( i in 2:4){
  points(agg_cv$N[agg_cv$N_mut==Nmut[i]],agg_cv$attribute_coef_var[agg_cv$N_mut==Nmut[i]],col=colors[i],pch=pchars[i])
  lines(agg_cv$N[agg_cv$N_mut==Nmut[i]],agg_cv$attribute_coef_var[agg_cv$N_mut==Nmut[i]],col=colors[i],lty=linetypes[i])
}
#title("Neutral Attribute Coefficient of Variation with burn-in=3.0")
legend(1200,0.15,labels,col=colors,pch=pchars,lty=linetypes,bty="n")
dev.copy(pdf,"neutral_attr_coef_var_N_mut_bi3_0.pdf")    # copies what is on the screen to a file.
dev.off()
dev.copy(png,"neutral_attr_coef_var_N_mut_bi3_0.png")    # copies what is on the screen to a file.
dev.off()
dev.copy(tiff,"neutral_attr_coef_var_N_mut_bi3_0.tiff")    # copies what is on the screen to a file.
dev.off()
