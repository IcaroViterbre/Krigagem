rm(list=ls())


library(geoR)
library(sp) 
library(ggplot2)
library(gstat)
library(automap)
library(RColorBrewer)
library(graphics)
library(raster)
library(ggspatial)
library(vroom)


DadosTotAjus <- as.data.frame(apply(vroom::vroom("Aleatoria.txt"), 2, as.numeric))

Borders <- as.data.frame(apply(vroom::vroom("BordasYL.txt"), 2, as.numeric))

DadosPlot <- DadosTotAjus

n <- length(DadosPlot$z)

limsup <- (sort(DadosPlot$z,partial=n)[n]+sort(DadosPlot$z,partial=n-1)[n-1]+
             sort(DadosPlot$z,partial=n-2)[n-2]+  sort(DadosPlot$z,partial=n-3)[n-3]+
             sort(DadosPlot$z,partial=n-4)[n-4])/5 

liminf <- (sort(DadosPlot$z,partial=1)[1]+sort(DadosPlot$z,partial=2)[2]+
             sort(DadosPlot$z,partial=3)[3]+  sort(DadosPlot$z,partial=4)[4]+
             sort(DadosPlot$z,partial=5)[5])/5


Dados <- geoR::as.geodata(DadosTotAjus,head = T, coords.col = 1:2, data.col = 3) 
Dados$borders <- Borders 
distmax <- 1*(variog(Dados, option ="bin", estimator.type ="modulus")$max.dist) 
semivar <- variog(Dados, estimator = "modulus",max.dist = distmax) 
 
sp::coordinates(DadosTotAjus) <- ~ x+y
 
vario.fit <- automap::autofitVariogram(z~1,
                              DadosTotAjus,
                              model =c("Exp","Sph","Gau"),
                              kappa=c(0.05,sq(0.2,2,0.1),5,10))


model <- as.character(vario.fit$var_model[[2,1]])
cov.model <- if(model=="Exp"){"exponential"}else{if(model=="Sph"){"spherical"}else{"gaussian"}}
semivarols <- variofit(semivar, cov.model = cov.model , ini=c((vario.fit$var_model[[1,2]]+vario.fit$var_model[[2,2]]),vario.fit$var_model[[2,3]]), nugget=vario.fit$var_model[[1,2]], weights="equal")
 
grid <- as.matrix(expand.grid(seq(min(Dados$borders[,1]),max(Dados$borders[,1]),l=200), seq(min(Dados$borders[,2]),max(Dados$borders[,2]),l=200)))

ind <- point.in.polygon(grid[,1],grid[,2],Dados$borders[,1],Dados$border[,2])

grid <- grid[which(ind==1),]
krg <- krige.conv(Dados,location=grid,krige=krige.control(cov.model=semivarols$cov.model,cov.pars=semivarols$cov.pars,nugget=semivarols$nugget))

r3 <- SpatialPointsDataFrame(grid, data.frame(predict = (krg$predict)))
gridded(r3) <- TRUE
r3 <- as(r3,'RasterLayer')
test_spdf3 <- as(r3, "SpatialPixelsDataFrame")
test_df3 <- as.data.frame(test_spdf3)
colnames(test_df3) <- c("value", "x", "y")
Final3 <- data.frame(long = test_df3$x, lat = test_df3$y, z = test_df3$value) 
Final3$z[Final3$z<liminf] <- liminf
Final3$z[Final3$z>limsup] <- limsup

palett <- rev(brewer.pal(11, "Spectral"))


ggplot(Final3, aes(x=long, y=lat, fill=z))+ 
  geom_raster()+ 
  geom_point(data = DadosPlot,aes(x=x, y=y, fill=z), 
             col = "black", alpha=0.7)+
  scale_fill_gradientn(limits=c(liminf,limsup),colours = palett )+ 
  scale_color_gradientn(limits=c(liminf,limsup),colors = palett )+
  theme_classic()+
  labs(title = " ",
       subtitle= " ",
       x=" ",
       y= "",
       fill= " ")+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        panel.grid= element_blank(),
        axis.ticks = element_blank(),
        panel.border=element_blank(),
        line=element_blank())

