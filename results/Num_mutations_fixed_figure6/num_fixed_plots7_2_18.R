# R script to produce plot of number of fixations for a fixed mu using the infinite sites model
#setwd("../../../nearly-neutral/data/9_13_17")
neut = read.csv("in_neut_mulist_N25_6400.csv",comment.char="#")
del5 = read.csv("in_del_mulist_N25_6400_theta0.5.csv",comment.char="#")
del1 = read.csv("in_del_mulist_N25_6400_theta1.csv",comment.char="#")
mixed05 = read.csv("in_mixed_mulist_N25_6400_theta0.5_0.005.csv",comment.char="#")
mixed1 = read.csv("in_mixed_mulist_N25_6400_theta0.5_0.01.csv",comment.char="#")
mixed2 = read.csv("in_mixed_mulist_N25_6400_theta0.5_0.02.csv",comment.char="#")

# uncomment the line with the nmu (mu) value that you want to use.
#nmu=0.0250; ymax = 17000; ymin=2000; legend_max = 12000 
#nmu=0.0100; ymax = 6900; ymin=900; legend_max = 4800
nmu=0.0050; ymax = 3400; ymin=400; legend_max = 2300
#nmu=0.0025; ymax = 1800; ymin = 200; legend_max = 1200 

neutstats=subset(neut,mu==nmu,select=c(N,mu,num_extinct,num_fixed,fraction_fixed,ave_fixed_time))
del5stats=subset(del5,mu==nmu,select=c(N,mu,num_extinct,num_fixed,fraction_fixed,ave_fixed_time))
del1stats=subset(del1,mu==nmu,select=c(N,mu,num_extinct,num_fixed,fraction_fixed,ave_fixed_time))
mixed05stats=subset(mixed05,mu==nmu,select=c(N,mu,num_extinct,num_fixed,fraction_fixed,ave_fixed_time))
mixed1stats=subset(mixed1,mu==nmu,select=c(N,mu,num_extinct,num_fixed,fraction_fixed,ave_fixed_time))
mixed2stats=subset(mixed2,mu==nmu,select=c(N,mu,num_extinct,num_fixed,fraction_fixed,ave_fixed_time))

# Add column "type" to dataframes
neutstats$type = 1
del5stats$type=2
del1stats$type=3
mixed05stats$type=4
mixed1stats$type=5
mixed2stats$type=6
# concatenate dataframes
all_stats = rbind(neutstats, 
   #del04stats, del08stats, del2stats, 
  del5stats, del1stats, mixed05stats,mixed1stats,mixed2stats) 
colors= c("green","red","darkblue","coral4","cyan","orange","brown")
pchars = c(15,16,17,18,19,20,14)
labels=c(
"Neutral",
expression(paste("Deleterious disadv ",beta,"=0.5")),
expression(paste("Deleterious disadv ",beta,"=1.0")),
expression(paste("Mixed adv ",beta,"=0.005")),
expression(paste("Mixed adv ",beta,"=0.01")),
expression(paste("Mixed adv ",beta,"=0.02")))
linetypes = c(1,2,3,4,5,6,7,8,9)
dev.off()
i=1
plot(all_stats$N[all_stats$type==i],all_stats$num_fixed[all_stats$type==i],log="x",xaxt="n",pch=pchars[1],col=colors[1],xlab="N",ylab="Number fixed",ylim=c(ymin,ymax),xlim=c(25.0,6400.0))
axis(1,at=c(25,50,100,200,400,800,1600,3200,6400))
lines(all_stats$N[all_stats$type==i],all_stats$num_fixed[all_stats$type==i],col=colors[1],lty=linetypes[1])
for( i in 2:6 ){
  points(all_stats$N[all_stats$type==i],all_stats$num_fixed[all_stats$type==i],pch=pchars[i],col=colors[i])
  lines(all_stats$N[all_stats$type==i],all_stats$num_fixed[all_stats$type==i],col=colors[i],lty=i)
}
nmu_string = paste0(" = ",nmu)
#title(expression(paste("Number mutations fixed with ",mu,nmu_string)))
legend(25,legend_max,labels,col=colors,pch=pchars,lty=linetypes,bty="n")
dev.copy(png,paste0("Num_mutations_fixed_mu_",nmu,".png"))    # copies what is on the screen to a file.
dev.off()
dev.copy(tiff,paste0("Num_mutations_fixed_mu_",nmu,".tiff"))    # copies what is on the screen to a file.
dev.off()
dev.copy(pdf,paste0("Num_mutations_fixed_mu_",nmu,".pdf"))    # copies what is on the screen to a file.
dev.off()

