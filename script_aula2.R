##
## Probit e Logit espacial
##


# Biliotecas
library(spdep)
library(spatialprobit)

# Abra o arquivo 'gm10.shp'
gm10.shp <- readShapePoly(file.choose(), 
                          IDvar = "CODMUN6", 
                          verbose=TRUE, 
                          proj4string=CRS("+proj=longlat +datum=WGS84"))

# Vari�vel dicot�mica
gm10.shp@data$y <- gm10.shp@data$TCVPA10 > mean(gm10.shp@data$TCVPA10)
table(gm10.shp@data$y)

# Matriz de vizinhos espaciais
coords <- coordinates(gm10.shp)
w1 <- kNearestNeighbors(coords[,1],coords[,2], k = 1)

# Especifica��o
esp <- y ~ GM + POLMPC09 + RICPOB10 + RENDA10 + DESOCU10 + CHEFA10 + DPOP10

# SAR Probit
mod1 <- sarprobit(esp, w1, gm10.shp@data)
summary(mod1)
impacts(mod1)
logLik(mod1)
plot(mod1)

# SAR Tobit
mod2 <- sartobit(esp, w1, gm10.shp@data)
summary(mod2)
impacts(mod2)
plot(mod2)

# SEM Probit
mod3 <- semprobit(esp, w1, gm10.shp@data)
summary(mod3)
logLik(mod3)
plot(mod3)

# SEM Tobit
mod4 <- sartobit(esp, w1, gm10.shp@data)
summary(mod4)
plot(mod4)







##
## Painel de Dados Espaciais
##

# Biblioteca
library(plm)
library(splm)

# Abra o arquivo 'guarda.shp'
guarda.shp <- readShapePoly(file.choose(), 
                          IDvar = "CODMUN6", 
                          verbose=TRUE, 
                          proj4string=CRS("+proj=longlat +datum=WGS84"))

# Matriz de vizinhos espaciais
IDs <- row.names(as.data.frame(guarda.shp))
coords <- coordinates(guarda.shp)
w1 <- matrix <- nb2listw(knn2nb(knearneigh(coords, k=1),row.names=IDs),style="W")

# Vari�veis para o modelo
dados <- guarda.shp@data
dados <- subset(dados, select=c("CODMUN6", "TCVPA00", "TCVPA10","RENDA00", "RENDA10", "THEIL00","THEIL10"))


# Vari�veis defasadas
dados$LAGTCVPA00 <- lag.listw(w1, dados$TCVPA00)
dados$LAGTCVPA10 <- lag.listw(w1, dados$TCVPA10)
dados$LAGRENDA00 <- lag.listw(w1, dados$RENDA00)
dados$LAGRENDA10 <- lag.listw(w1, dados$RENDA10)
dados$LAGTHEIL00 <- lag.listw(w1, dados$THEIL00)
dados$LAGTHEIL10 <- lag.listw(w1, dados$THEIL10)

# Empilhar dados
dadosp <- data.frame(codmun=as.integer(), ano=as.integer(), tcvpa=as.numeric(), renda=as.numeric(), theil=as.numeric(), lagtcvpa=as.numeric(), lagrenda=as.numeric(), lagtheil=as.numeric())
for(i in 1:nrow(dados)){
  for(a in 1:2){
    dadosp <- rbind(dadosp, data.frame(codmun=dados[i, "CODMUN6"],
                                       ano=ifelse(a==1,"2000","2010"),
                                       tcvpa=dados[i,1+a],
                                       renda=dados[i,3+a],
                                       theil=dados[i,5+a],
                                       lagtcvpa=dados[i,7+a],
                                       lagrenda=dados[i,9+a],
                                       lagtheil=dados[i,11+a]))
  }
}

# Especifica��o
esp2 <- tcvpa ~ renda + theil

# Painel n�o espacial, efeitos fixos
fe <- plm(esp2, data=dadosp)
summary(fe)

# Painel n�o espacial, efeitos aleat�rios
re <- plm(esp2, data=dadosp, model = "random")
summary(re)

# Teste de Hausman
ph <- phtest(fe, re)
print(ph)

# Teste de Pesaran CD
cd <- pcdtest(esp2, data=dadosp)
print(cd)

# Modelo tradicional (n�o espacial)
modTrad <- plm(esp2, data = dadosp)
summary(modTrad)

# Modelos SAR, SEM, SAC

# Modelo SAR
modSAR <- spml(esp2, data=dadosp, listw=w1, lag=TRUE, model="within", effect="individual", spatial.error="none")
summary(modSAR)
impSAR <- impacts(modSAR, listw=w1, time=2)
summary(impSAR)

# Modelo SEM
modSEM <- spml(esp2, data=dadosp, listw=w1, lag=FALSE, model="within", effect="individual", spatial.error="b")
summary(modSEM)

# Modelo SAC
modSAC <- spml(esp2, data=dadosp, listw=w1, lag=TRUE, model="within", effect="individual", spatial.error="b")
summary(modSAC)
impSAC <- impacts(modSAC, listw=w1, time=2)
summary(impSAC)

# Modelos SDM, SDEM, SLX

esp3 <- tcvpa ~ lagrenda + lagtheil

# Modelo SDM
modSDM <- spml(esp3, data=dadosp, listw=w1, lag=TRUE, model="within", effect="individual", spatial.error="none")
summary(modSDM)
impSDM <- impacts(modSDM, listw=w1, time=2)
summary(impSDM)

# Modelo SDEM
modSDEM <- spml(esp3, data=dadosp, listw=w1, lag=FALSE, model="within", effect="individual", spatial.error="b")
summary(modSDEM)

# Modelo SLX
modSLX <- plm(esp3, data=dadosp)
summary(modSLX)
