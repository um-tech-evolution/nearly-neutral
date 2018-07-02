#setwd("../9_9_17")

neut = read.csv("nn_neut_Nmulist_fm45_N25_6400.csv",comment.char="#")
del0_5 = read.csv("nn_del_Nmulist_fm45_N25_6400_theta0.5.csv",comment.char="#")
del1 = read.csv("nn_del_Nmulist_fm45_N25_6400_theta1.csv",comment.char="#")
mixed05 = read.csv("nn_mixed_Nmulist_fm45_prob0.01_N25_6400_theta0.5_0.005.csv",comment.char="#")
mixed1 = read.csv("nn_mixed_Nmulist_fm45_prob0.01_N25_6400_theta0.5_0.01.csv",comment.char="#")
mixed2 = read.csv("nn_mixed_Nmulist_fm45_prob0.01_N25_6400_theta0.5_0.02.csv",comment.char="#")

# neut = read.csv("nn_neut_Nmulist_fm45_N50_6400.csv",comment.char="#")
# del0_5 = read.csv("nn_del_Nmulist_fm45_N50_6400_theta0.5.csv",,comment.char="#")
# del1 = read.csv("nn_del_Nmulist_fm45_N50_6400_theta1.csv",comment.char="#")
# mixed05=read.csv("nn_mixed_Nmulist_fm45_prob0.01_N50_6400_theta0.5_0.005.csv",comment.char="#")
# mixed1=read.csv("nn_mixed_Nmulist_fm45_prob0.01_N50_6400_theta0.5_0.01.csv",comment.char="#")
# mixed2=read.csv("nn_mixed_Nmulist_fm45_prob0.01_N50_6400_theta0.5_0.02.csv",comment.char="#")
#nmu=0.5      # theta = 1.0
nmu=1.0      # theta = 2.0
#nmu=2.0      # theta = 4.0
#nmu=4.0      # theta = 8.0
if (nmu==1.0) { ymin=0.0;  ymax=0.7; legend_max = 0.21; legend_right=22 }
# settings for the other nmu values have not been determined
  
neutstats=subset(neut,N_mu==nmu,select=c(N,N_mu,expected_w_heteroz,w_heteroz,expected_richness,average_richness,IQV))
del0_5stats=subset(del0_5,N_mu==nmu,select=c(N,N_mu,expected_w_heteroz,w_heteroz,expected_richness,average_richness,IQV))
del1stats=subset(del1,N_mu==nmu,select=c(N,N_mu,expected_w_heteroz,w_heteroz,expected_richness,average_richness,IQV))
mixed05stats=subset(mixed05,N_mu==nmu,select=c(N,N_mu,expected_w_heteroz,w_heteroz,expected_richness,average_richness,IQV))
mixed1stats=subset(mixed1,N_mu==nmu,select=c(N,N_mu,expected_w_heteroz,w_heteroz,expected_richness,average_richness,IQV))
mixed2stats=subset(mixed2,N_mu==nmu,select=c(N,N_mu,expected_w_heteroz,w_heteroz,expected_richness,average_richness,IQV))

# Add column "type" to dataframes
neutstats$type = 1
del0_5stats$type=2
del1stats$type=3
mixed05stats$type=4
mixed1stats$type=5
mixed2stats$type=6
# concatenate dataframes
all_stats = rbind(neutstats, del0_5stats,del1stats,mixed05stats,mixed1stats,mixed2stats)
colors= c("green","red","darkblue","coral4","cyan","orange","brown")
pchars = c(15,16,17,18,19,20,14)
labels=c(
"Neutral",
expression(paste("Deleterious ",beta,"=0.5")),
expression(paste("Deleterious ",beta,"=1.0")),
expression(paste("Mixed ",beta,"=0.005")),
expression(paste("Mixed ",beta,"=0.1")),
expression(paste("Mixed ",beta,"=0.2")),
"Neutral expected")
linetypes = c(1,2,3,4,5,6,7)
dev.off()
i=1
plot(all_stats$N[all_stats$type==i],all_stats$expected_w_heteroz[all_stats$type==i],xaxt="n",log="x",pch=pchars[7],col=colors[7],xlab="N",ylab="Heterozygosity",ylim=c(0.0,0.7),xlim=c(25.0,6400.0))
axis(1,at=c(25,50,100,200,400,800,1600,3200,6400))
lines(all_stats$N[all_stats$type==i],all_stats$expected_w_heteroz[all_stats$type==i],col=colors[7],lty=linetypes[7])
for( i in 1:6 ){
   points(all_stats$N[all_stats$type==i],all_stats$w_heteroz[all_stats$type==i],pch=pchars[i],col=colors[i])
  lines(all_stats$N[all_stats$type==i],all_stats$w_heteroz[all_stats$type==i],col=colors[i],lty=i)
}
#title("Heterozygosity with  θ = 2")
legend(legend_right,legend_max,labels,col=colors,pch=pchars,lty=linetypes,bty="n")
dev.copy(png,'Heterozygosity_theta_2.png')    # copies what is on the screen to a file.
dev.off()
dev.copy(tiff,'Heterozygosity_theta_2.tiff')    # copies what is on the screen to a file.
dev.off()
dev.copy(pdf,'Heterozygosity_theta_2.pdf')    # copies what is on the screen to a file.
dev.off()

