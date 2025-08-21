##########SCatter plot
dat = as.data.frame(allDiff)

library(data.table)
dat2= reshape2::melt(dat,id=3)
dat2$value = as.numeric(dat2$value)

pcreg = function(ds1, ds2) {
  #Principle components regression to calculate slope 
  r = prcomp(~ds1+ds2)
  slope <- r$rotation[2,1] / r$rotation[1,1]
  intercept <- r$center[2] - slope*r$center[1]
  return(list(slope,intercept))
}

fit_1= pcreg(dat$ESCC, dat$LUSC)
fit_2 = pcreg(dat$ESCC, dat$CSCC)
fit_3 = pcreg(dat$ESCC, dat$HNSC)
fit_4 = pcreg(dat$ESCC, dat$CESC)

dat2$variable = as.character(dat2$variable)

dat2$variable = gsub("LUSC", paste("LUSC, slope=", signif(fit_1[[1]],2), sep=""), dat2$variable)
dat2$variable = gsub("CSCC", paste("CSCC, slope=", signif(fit_2[[1]],2), sep=""), dat2$variable)
dat2$variable = gsub("HNSC", paste("HNSC, slope=", signif(fit_3[[1]],2), sep=""), dat2$variable)
dat2$variable = gsub("CESC", paste("CESC, slope=", signif(fit_4[[1]],2), sep=""), dat2$variable)


library(ggplot2)

SCCs_DEGs=ggplot(dat2,aes(x=ESCC,y=value,color=variable)) + 
  scale_color_manual(values = c("#DBE0ED","#CC88B0", "#998DB7","#F4CEB4"))+ 
  geom_point(alpha=.9, size = 0.7) + 
  geom_abline(slope=1, lty=2) + xlim(-4,4) + ylim(-4,4) + 
  geom_abline(slope=fit_1[[1]], intercept = fit_1[[2]], color="#F4CEB4",size = 0.7) + 
  geom_abline(slope=fit_2[[1]], intercept = fit_2[[2]], color="#CC88B0",size = 0.7) +
  geom_abline(slope=fit_3[[1]], intercept = fit_3[[2]], color="#998DB7",size = 0.7) + 
  geom_abline(slope=fit_4[[1]], intercept = fit_4[[2]], color="#DBE0ED",size = 0.7) +
  xlab("ESCC (log2FC)") + ylab("SCCs(log2FC)") +
  coord_fixed(ratio=1)+
  theme(panel.grid=element_blank(), 
        panel.background=element_rect(color="black", fill="transparent"))##theme

SCCs_DEGs