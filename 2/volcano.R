

#install.packages("ggplot2")


library(ggplot2)           
logFCfilter=0.5             
adj.P.Val.Filter=0.05      
inputFile="all.txt"        
setwd("C:\\Users\\lyj19\\Desktop\\vtegse19151logfc0.5\\04")       


rt=read.table(inputFile, header=T, sep="\t", check.names=F)

Sig=ifelse((rt$adj.P.Val<adj.P.Val.Filter) & (abs(rt$logFC)>logFCfilter), ifelse(rt$logFC>logFCfilter,"Up","Down"), "Not")
rt=cbind(rt, Sig=Sig)


p=ggplot(rt, aes(logFC, -log10(adj.P.Val)))+
    geom_point(aes(col=Sig))+
    scale_color_manual(values=c("blue", "black", "red"))+
    xlim(-2,2)+
    labs(title = " ")+
    geom_vline(xintercept=c(-logFCfilter,logFCfilter), col="green", cex=1, linetype=2)+
    geom_hline(yintercept= -log10(adj.P.Val.Filter), col="green", cex=1, linetype=2)+
    theme(plot.title=element_text(size=16, hjust=0.5, face="bold"))
p=p+theme_bw()


pdf(file="volcano.pdf", width=6, height=5.1)
print(p)
dev.off()