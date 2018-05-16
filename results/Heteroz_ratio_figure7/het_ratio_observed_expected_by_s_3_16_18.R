# Create plot of H1/H0 2*(4*N*s-1+exp(-4*N*s))/(4*N*s(1-exp(-4*N*s)))
#   which is equation 5.25, page 250 of Hartl & Clark 4th edition.
# This version of 3/16/18 does not include the counts, only neutral and selected.
# H1 is heterozygosity with selection, H0 is neutral heterozygosity
# This version (6/1/17) relies on the spreadsheet to calculate the H1/H0 ratio and the observed ratios.
# The spreadsheet is built by concatenating the results for different selection coefficient.
# Remember that the formulas must be "filled" for each section of the spreadsheet.
setwd("../../data/6_1_17/2_22_17/")   # Doesn't work, but shows the right directory
setwd("nearly-neutral/data/9_13_17fixed/")   # Doesn't work, but shows the right directory

# Remember that when saving from xlsx to csv, you need to check that the N_list and mu_list parameter lines are preceded by a #
s_all = read.csv("in_fixed_mulist_N25_6400_s_all_ratios.csv",comment.char="#")
s_neut=subset(s_all,s==0.0,select=c(s,N,ave_heteroz,observed_heteroz_ratio,expected_heteroz_ratio,ave_innov_count,count_ratio))
s0005=subset(s_all,s==-0.0005,select=c(s,N,ave_heteroz,observed_heteroz_ratio,expected_heteroz_ratio,ave_innov_count,count_ratio))
s001=subset(s_all,s==-0.001,select=c(s,N,ave_heteroz,observed_heteroz_ratio,expected_heteroz_ratio,ave_innov_count,count_ratio))
s002=subset(s_all,s==-0.002,select=c(s,N,ave_heteroz,observed_heteroz_ratio,expected_heteroz_ratio,ave_innov_count,count_ratio))
s005=subset(s_all,s==-0.005,select=c(s,N,ave_heteroz,observed_heteroz_ratio,expected_heteroz_ratio,ave_innov_count,count_ratio))
s01=subset(s_all,s==-0.01,select=c(s,N,ave_heteroz,observed_heteroz_ratio,expected_heteroz_ratio,ave_innov_count,count_ratio))
s02=subset(s_all,s==-0.02,select=c(s,N,ave_heteroz,observed_heteroz_ratio,expected_heteroz_ratio,ave_innov_count,count_ratio))

colors= c("green","red","darkblue","coral4","cyan","orange","brown")
pchars = c(15,16,17,18,19,20,14)
linetypes = c(1,2,3,4,5,6,7)
labels=c("expected heteterozygosity","counts")
#labels=c("observed heterozygosity","expected heteterozygosity")
#labels=c("observed heterozygosity","expected heteterozygosity","counts")

# Uncomment the line corresponding to the selection coefficient s that you want.
#sss = setNames(s0005, names(s0005)); legend_right=50; legend_max = 0.3;
#sss = setNames(s001, names(s001)); legend_right=50; legend_max = 0.3;
#sss = setNames(s002, names(s002)); legend_right=50; legend_max = 0.3;
#sss = setNames(s005, names(s005)); legend_right=50; legend_max = 0.3;
#sss = setNames(s01, names(s01)); legend_right=270; legend_max = 0.95;
sss = setNames(s02, names(s02)); legend_right=270; legend_max = 0.95;
dev.off()
plot(sss$N,sss$expected_heteroz_ratio,type="l",log="x",ylim=c(0.05,1.00),col=colors[1],lty=1,xlab="N",ylab="Ratio",xaxt="n",yaxt="n")
axis(1,at=sss$N)
axis(2,at=c(0.1,0.3,0.5,0.7,0.9))
points(sss$N,sss$expected_heteroz_ratio,col=colors[1],pch=pchars[1])
#lines(sss$N,sss$observed_heteroz_ratio,col=colors[2],lty=2,pch=pchars[2])
#points(sss$N,sss$observed_heteroz_ratio,col=colors[2],pch=pchars[2])
lines(sss$N,sss$count_ratio,col=colors[3],lty=3,pch=pchars[3])
points(sss$N,sss$count_ratio,col=colors[3],pch=pchars[3])
legend(legend_right,legend_max,labels,col=colors,pch=pchars,lty=linetypes)
title(paste0("Ratio of Neutral to Expected for s = ",sss$s[1]))
dev.copy(png,paste0("het_ratio_neut_to_exp_s_",sss$s[1],".png"))
dev.off()
dev.copy(tiff,paste0("het_ratio_neut_to_exp_s_",sss$s[1],".tiff"))
dev.off()
dev.copy(pdf,paste0("het_ratio_neut_to_exp_s_",sss$s[1],".pdf"))
dev.off()



