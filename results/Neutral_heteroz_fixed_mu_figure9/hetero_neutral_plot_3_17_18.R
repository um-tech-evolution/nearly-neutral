setwd("../../../nearly_neutral/data/12_11_17")

neut = read.csv("N25_6400_nn_neut_mutstddev.csv",comment.char="#")

neut025 = subset(neut,mu==0.025,select=c(N,w_heteroz,expected_w_heteroz,mu))
neut01 = subset(neut,mu==0.01,select=c(N,w_heteroz,expected_w_heteroz,mu))
neut005 = subset(neut,mu==0.005,select=c(N,w_heteroz,expected_w_heteroz,mu))
neut0025 = subset(neut,mu==0.0025,select=c(N,w_heteroz,expected_w_heteroz,mu))

mu_list= c(0.025,0.01,0.005,0.0025)
ymin=0.00;  ymax=1.0; legend_max = 0.6; legend_right=1000 
# settings for the other nmu values have not been determined
  
neutstats=subset(neut,mu==nmu,select=c(N,mu,expected_w_heteroz,w_heteroz))

colors= c("red","green","darkblue","coral4","cyan","orange","brown")
pchars = c(15,16,17,18,19,20,14)
labels=c(
paste0("Neutral µ = ",mu_list[2]),
paste0("Neutral µ = ",mu_list[1]),
paste0("Neutral µ = ",mu_list[3]),
paste0("Neutral µ = ",mu_list[4])
)
linetypes = c(1,2,3,4,5,6,7)
dev.off()
plot(neut01$N,neut01$expected_w_heteroz,xaxt="n",log="x",pch=pchars[1],col=colors[1],xlab="N",ylab="Heterozygosity",ylim=c(ymin,ymax),xlim=c(25.0,6400.0))
axis(1,at=c(25,50,100,200,400,800,1600,3200,6400))
lines(neut01$N,neut01$expected_w_heteroz,col=colors[1],lty=linetypes[1])
points(neut025$N,neut025$expected_w_heteroz,col=colors[2],pch=pchars[2])
lines(neut025$N,neut025$expected_w_heteroz,col=colors[2],lty=linetypes[2])
points(neut005$N,neut005$expected_w_heteroz,col=colors[3],pch=pchars[3])
lines(neut005$N,neut005$expected_w_heteroz,col=colors[3],lty=linetypes[3])
points(neut0025$N,neut0025$expected_w_heteroz,col=colors[4],pch=pchars[4])
lines(neut0025$N,neut0025$expected_w_heteroz,col=colors[4],lty=linetypes[4])
title("Neutral Heterozygosity for fixed µ")
legend(legend_right,legend_max,labels,col=colors,pch=pchars,lty=linetypes,bty="n")
dev.copy(png,'Neutral_Heterozygosity_by_mu.png')    # copies what is on the screen to a file.
dev.off()
dev.copy(tiff,'Neutral_Heterozygosity_by_mu.tiff')    # copies what is on the screen to a file.
dev.off()
dev.copy(pdf,'Neutral_Heterozygosity_by_mu.pdf')    # copies what is on the screen to a file.
dev.off()

