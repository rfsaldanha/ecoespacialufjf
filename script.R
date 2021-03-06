# Pacotes
library(maptools)
library(rgeos)

# Abra o arquivo 'gm10.shp'
gm10.shp <- readShapePoly(file.choose(), 
                          IDvar = "CODMUN6", 
                          verbose=TRUE, 
                          proj4string=CRS("+proj=longlat +datum=WGS84"))

# Plotar o mapa
plot(gm10.shp)

# Pacotes
library(spdep)

# Matriz queen
w1 <- nb2listw(poly2nb(gm10.shp, queen = TRUE))
summary(w1)

# Matrix queen padronizada na linha
w1.w <- nb2listw(poly2nb(gm10.shp, queen=TRUE), style="W")
summary(w1.w)

# Matriz rook
w2 <- nb2listw(poly2nb(gm10.shp, queen = FALSE))
summary(w2)

# Matriz rook padronizada globalmente
w2.c <- nb2listw(poly2nb(gm10.shp, queen = FALSE), style = "C")
summary(w2.c)

# Dist�ncia inversa
coords <- coordinates(gm10.shp)
nb <- dnearneigh(coords, 0, 1000)
dlist <- nbdists(nb, coords)
dlist <- lapply(dlist, function(x) 1/x)
w3 <- nb2listw(nb, glist=dlist)
summary(w3)

# Dist�ncia inversa padronizada pelo n�mero de vizinhos
w3.u <- nb2listw(nb, glist=dlist, style="U")
summary(w3.u)

# K vizinhos espaciais

# N�mero de permuta��es
per <- 999

# N�mero m�ximo de k vizinhos testados
kv <- 20

# Nome dos registros
IDs <- row.names(gm10.shp@data)

# Cria��o da tabela que ir� receber a estat�stica I de Moran e signific�ncia
# para cada k n�mero de vizinhos testado
res.pesos <- data.frame()

# In�cio do loop
for(k in 1:kv)
{
  # Armazenando n�mero k atual
  res.pesos[k,1] <- k
  # Calculando o I e signific�ncia para o k atual
  moran.k <- moran.mc(gm10.shp@data$TCVPA10,
                      listw=nb2listw(knn2nb(knearneigh(coords, k=k),
                                            row.names=IDs),style="W"),
                      nsim=per)
  # Armazenando o valor I para o k atual
  res.pesos[k,2] <- moran.k$statistic
  # Armazenando o p-value para o k atual
  res.pesos[k,3] <- moran.k$p.value
}

# Ver a tabela de k vizinhos, I de Moran e signific�ncia
res.pesos

# Sendo todos significativos, iremos usar o k que retornou o maior valor I
maxi <- which.max(res.pesos[,2])

# Cria��o da matriz usando o k escolhido
w5 <- nb2listw(knn2nb(knearneigh(coords, k=maxi),row.names=IDs),style="W")


# Importando um arquivo GAL do Geoda
w6 <- nb2listw(read.gal(file.choose(), override.id = TRUE))

# Modelo n�o espacial
esp <- TCVPA10 ~ GM + POLMPC09 + RICPOB10 + RENDA10 + DESOCU10 + CHEFA10 + DPOP10
mod1 <- lm(formula = esp, data = gm10.shp@data)
summary(mod1)

# Teste dos res�duos: Moran I
mod1.moran <- lm.morantest(model = mod1, listw = w5)
mod1.moran

# Teste dos res�duos: Multiplicador de Lagrange
mod1.lagrange <- lm.LMtests(model = mod1, listw = w5,
                            test = c("LMerr","RLMerr","LMlag","RLMlag","SARMA"))
mod1.lagrange


# SAR
mod1.sar <- lagsarlm(formula = esp, data = gm10.shp@data, listw = w5)
summary(mod1.sar)
impacts(mod1.sar, listw = w5)

# CAR
mod1.car <- spautolm(formula = esp, data = gm10.shp@data, listw = w5, family = "CAR")
summary(mod1.car)

# SEM
mod1.sem <- errorsarlm(formula = esp, data = gm10.shp@data, listw = w5)
summary(mod1.sem)

# SEM por GMM
mod1.semGMM <- GMerrorsar(formula = esp, data = gm10.shp@data, listw = w5)
summary(mod1.semGMM)

# SAC
mod1.sac <- sacsarlm(formula = esp, data = gm10.shp@data, listw = w5)
summary(mod1.sac)
impacts(mod1.sac, listw = w5)

# SAC por GMM
mod1.sacGMM <- gstsls(formula = esp, data = gm10.shp@data, listw = w5)
summary(mod1.sacGMM)
impacts(mod1.sacGMM, listw = w5)

# SDM
mod1.sdm <- lagsarlm(formula = esp, data = gm10.shp@data, listw = w5, type = "mixed")
summary(mod1.sdm)
impacts(mod1.sdm, listw = w5)

# SDEM
mod1.sdem <- errorsarlm(formula = esp, data = gm10.shp@data, listw = w5, etype = "emixed")
summary(mod1.sdem)

# SMA
mod1.sma <- spautolm(formula = esp, data = gm10.shp@data, listw = w5, family = "SMA")
summary(mod1.sma)





# NOT RUN

# Modelos com heterocedasticidade
library(sphet)
# HAC
gm10.shp.utm <- spTransform(gm10.shp, CRS("+init=epsg:32723"))
coords.utm <- coordinates(gm10.shp.utm)
coords.utm <- as.data.frame(coords.utm)
coords.utm$id  <- row.names(coords.utm)
coords.utm <- coords.utm[,c(3,1,2)]
names(coords.utm) <- c("id","lat","long")

coords.dist <- distance(coord = coords.utm, 
                        region.id = coords.utm$id, 
                        output = TRUE, 
                        file = "dist.gwt", 
                        type = "NN", 
                        nn=6, 
                        firstline = FALSE,
                        shape.name = "gm10.shp.utm")
coords.dist <- read.gwt2dist(file = "dist.gwt")
mod1.hac <- stslshac(formula = esp, data = gm10.shp@data, listw = w5, distance = coords.dist)
mod1.hec <- gstslshet(formula = esp, data = gm10.shp@data, listw = w5)
mod1.spreg <- spreg(formula = esp, data = gm10.shp@data, listw = w5, het = TRUE, verbose = FALSE)

# END NOT RUN