# R code which produces the richness plot or plots for the infinite alleles model for the Nearly Neutral paper
# This code can be pasted into R or read into R using the R source command.
#if( flag2_27 ){ setwd("../02_27_17") }
setwd("nearly-neutral/data/9_9_17")   # May not work, but shows the folder containing the data

del05 = read.csv("nn_del_Nmulist_fm45_N25_6400_theta0.5.csv",comment.char="#")
del10 = read.csv("nn_del_Nmulist_fm45_N25_6400_theta1.csv",comment.char="#")
mix05 = read.csv("nn_mixed_Nmulist_fm45_prob0.01_N25_6400_theta0.5_0.005.csv",comment.char="#")
mix1 = read.csv("nn_mixed_Nmulist_fm45_prob0.01_N25_6400_theta0.5_0.01.csv",comment.char="#")
mix2 = read.csv("nn_mixed_Nmulist_fm45_prob0.01_N25_6400_theta0.5_0.02.csv",comment.char="#")
neut = read.csv("nn_neut_Nmulist_fm45_N25_6400.csv",comment.char="#")

# these must come before subset statements
#nmu=0.5
#nmu=1.0
#nmu=2.0
nmu=4.0
neutstats=subset(neut,N_mu==nmu,select=c(N,N_mu,expected_w_heteroz,w_heteroz,expected_richness,average_richness,IQV))
del05stats=subset(del05,N_mu==nmu,select=c(N,N_mu,expected_w_heteroz,w_heteroz,expected_richness,average_richness,IQV))
del10stats=subset(del10,N_mu==nmu,select=c(N,N_mu,expected_w_heteroz,w_heteroz,expected_richness,average_richness,IQV))
mix05stats=subset(mix05,N_mu==nmu,select=c(N,N_mu,expected_w_heteroz,w_heteroz,expected_richness,average_richness,IQV))
mix1stats=subset(mix2,N_mu==nmu,select=c(N,N_mu,expected_w_heteroz,w_heteroz,expected_richness,average_richness,IQV))
mix2stats=subset(mix1,N_mu==nmu,select=c(N,N_mu,expected_w_heteroz,w_heteroz,expected_richness,average_richness,IQV))
# Add column "type" to dataframes
neutstats$type = 1
del05stats$type=2
del10stats$type=3
mix05stats$type=4
mix1stats$type=5
mix2stats$type=6
# concatenate dataframes
all_stats = rbind(neutstats, del05stats, del10stats,mix05stats,mix1stats,mix2stats)

# all nmu for richness
exp = all_stats$expected_richness; avg = all_stats$average_richness
ylabel = "Richness"
if( nmu == 0.5 ){ ymax = 10.; ymin = 2.40; legend_max =ymax; legend_right = 22}
if( nmu == 1.0 ){ ymax = 17.; ymin = 3.30; legend_max =ymax; legend_right = 22 }
if( nmu == 2.0 ){ ymax = 32.; ymin = 5.00; legend_max =32.; legend_right = 22 }
if( nmu == 4.0 ){ ymax = 54.; ymin = 10.0; legend_max =54.; legend_right = 22 }

colors= c("green","red","darkblue","coral4","cyan","orange","brown")
pchars = c(15,16,17,18,19,20,14)
linetypes = c(1,2,3,4,5,6,7)
labels=c(
  expression(paste("Neutral ")),
  expression(paste("Deleterious ",beta,"=0.5")),
  expression(paste("Deleterious ",beta,"=1.0")),
  expression(paste("Mixed ",beta,"=0.005")),
  expression(paste("Mixed ",beta,"=0.01")),
  expression(paste("Mixed ",beta,"=0.02")),
  expression(paste("Neutral expected "))   # comment out for IQV but not heterozygosity or richness
)
dev.off()
i=1
plot(all_stats$N[all_stats$type==i],avg[all_stats$type==i],log="x",xaxt="n",pch=pchars[1],col=colors[1],xlab="N",ylab=ylabel,ylim=c(ymin,ymax),xlim=c(25.0,6400.0))
axis(1,at=c(25,50,100,200,400,800,1600,3200,6400))
lines(all_stats$N[all_stats$type==i],avg[all_stats$type==i],col=colors[1],lty=linetypes[1])
for( i in 2:6 ){
  points(all_stats$N[all_stats$type==i],avg[all_stats$type==i],pch=pchars[i],col=colors[i])
  lines(all_stats$N[all_stats$type==i],avg[all_stats$type==i],col=colors[i],lty=i)
}
i=7  # Neutral expected.  Comment out for IQV since no IQV expected
points(all_stats$N[all_stats$type==1],exp[all_stats$type==1],pch=pchars[i],col=colors[i]) # comment out for IQV
lines(all_stats$N[all_stats$type==1],exp[all_stats$type==1],col=colors[i],lty=i)  # comment out for IQV
title(paste0(ylabel," with θ = ",(2.0*nmu)))
legend(legend_right,legend_max,labels,col=colors,pch=pchars,lty=linetypes,bty="n")
dev.copy(pdf,paste0(ylabel,"_theta_",(2.0*nmu),".pdf"))    # copies what is on the screen to a file.
dev.off()
dev.copy(png,paste0(ylabel,"_theta_",(2.0*nmu),".png"))    # copies what is on the screen to a file.
dev.off()
dev.copy(tiff,paste0(ylabel,"_theta_",(2.0*nmu),".tiff"))    # copies what is on the screen to a file.
dev.off()
