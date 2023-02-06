# author: Li Juan

library(wesanderson)
args=commandArgs(trailingOnly = TRUE)

popSize=5000
add=c(0.01)
epi=c(0, 0.01, 0.1)

pdf("demo.pdf", width=5, height=2)

par(mfrow=c(1,1), mai=c(0.4, 0.42, 0.2,0.08), mgp=c(2,0.5,0), cex = 0.7,  las=1, tck=-0.02)

plot(NULL, xlim=c(0, 1600), ylim=c(-8, -1), ylab="", xlab="", xaxt="n",  yaxt="n")
#title(main=paste("Ruggedness=", e/add, sep=""), cex.main=0.8, line=0.2)
mtext("Generation", side=1, line=1.5, cex=0.7)
mtext("Population fitness variance", side=2, line=2, cex=0.7, las=0)

generation=c(0, seq(200, 1600, 200))
axis(1, at=generation, labels = generation)
axis(2, at = seq(-9, -1, 1), labels = expression(10^-9, 10^-8, 10^-7, 10^-6, 10^-5, 10^-4, 10^-3, 10^-2, 10^-1))

colors=wes_palette("Zissou1", 10, type="continuous")[c(1,4,9)]

text(750, -1,  substitute(bold("Ruggedness:")), adj=c(1,1.85), cex=1)
legend(750, -1,  legend = c(0, 1, 10), col=colors, pch=rep(19, 3),bty ="n", cex=1, horiz=TRUE)

i=1

for(e in epi){

  data=read.table(paste(args[1],"population_details/RMF_L15_mu0.01_stda0_stdb", e, ".dat", sep=""),
                  header=TRUE, comment.char = "")

  for(FLid in 0:9){
    tmp=subset(data, data[,1]==FLid)
    print(head(tmp))
    points(tmp$t, log10(tmp$var_fitness), type="l", col=colors[i])

  }

  i=i+1
}


dev.off()
