# R script to produce plot the attribute coefficient of variation for N from 25 to 6400 and for various mutation standard deviations.
# Author:  Alden Wright  12/14/2017  (verified 1/18/18, 7/2/18)
#setwd("../../continuous_variable_evolution/data/9_30_17")
cv = read.csv("N6400_2_na1_mutstddev_bi400_multnerr_neutral.csv",comment.char="#")
sub_cv=subset(cv,N>=25,select=c(N,mutation_stddev,attribute_mean,attribute_coef_var))
agg_cv=aggregate(sub_cv,by=list(sub_cv$N,sub_cv$mutation_stddev),FUN=mean)
colors= c("green","red","darkblue","coral4","cyan","orange","brown")
pchars = c(15,16,17,18,19,20,14)
labels=c(
"mutation stddev = 0.002",
"mutation stddev = 0.005",
"mutation stddev = 0.01",
"mutation stddev = 0.02",
"mutation stddev = 0.04"
)
linetypes = c(1,2,3,4,5,6,7,8,9)
mtsd = c(0.002,0.005,0.01,0.02,0.04)
#dev.off()  # uncomment if a graphics window is already open
i=1
plot(agg_cv$N[agg_cv$mutation_stddev==mtsd[i]],agg_cv$attribute_coef_var[agg_cv$mutation_stddev==mtsd[i]],
	log="x",xaxt="n",pch=pchars[1],col=colors[1],xlab="N",ylab="attribute coeffcient of variation",ylim=c(0.0,1.0),xlim=c(25.0,6400.0))
axis(1,at=c(25,50,100,200,400,800,1600,3200,6400))
lines(agg_cv$N[agg_cv$mutation_stddev==mtsd[i]],agg_cv$attribute_coef_var[agg_cv$mutation_stddev==mtsd[i]],col=colors[1],lty=linetypes[1])
for( i in 2:5){
  points(agg_cv$N[agg_cv$mutation_stddev==mtsd[i]],agg_cv$attribute_coef_var[agg_cv$mutation_stddev==mtsd[i]],col=colors[i],pch=pchars[i])
  lines(agg_cv$N[agg_cv$mutation_stddev==mtsd[i]],agg_cv$attribute_coef_var[agg_cv$mutation_stddev==mtsd[i]],col=colors[i],lty=linetypes[i])
}
#title("Neutral Attribute Coefficient of Variation with burn-in=400")
legend(25,0.9,labels,col=colors,pch=pchars,lty=linetypes,bty="n")
dev.copy(png,"neutral_attr_coef_var_bi400.png")    # copies what is on the screen to a file.
dev.off()
dev.copy(tiff,"neutral_attr_coef_var_bi400.tiff")    # copies what is on the screen to a file.
dev.off()
dev.copy(pdf,"neutral_attr_coef_var_bi400.pdf")    # copies what is on the screen to a file.
dev.off()
