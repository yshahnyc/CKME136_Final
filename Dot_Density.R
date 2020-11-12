#Install (if needed) libraries
install.packages("dplyr")
install.packages("maptools")
install.packages("rgeos")
install.packages("tidyverse")
install.packages("rgdal")
install.packages("raster")
install.packages("ggmap")
install.packages("ggplot2")

#Load Libraries
library(dplyr)
library(sp)
library(maptools)
library(rgeos)
library(tidyverse)
library(rgdal)

#Download population data from Statistics Canada Census Website and save in working directory
BasePopulation = read.csv("98-400-X2016003_English_CSV_data.csv")

#For Referencing location codes
PlaceNames = read.csv("Geo_starting_row_CSV.csv")

#Clean table for total population and create city table. Find the correct code for GEO_CSD in PlaceNames table. Using "Orangeville"" as an example.
TotalPopulation = subset(BasePopulation,`DIM..Age..in.single.years..and.average.age..127.` == "Total - Age", select = c("GEO_CODE..POR.","ALT_GEO_CODE","Dim..Sex..3...Member.ID...1...Total...Sex"))
DAPopulation = subset(TotalPopulation, nchar(TotalPopulation$GEO_CODE..POR.) == max(nchar(TotalPopulation$GEO_CODE..POR.)))
colnames(DAPopulation) = c("DAUID","GEO_ID","POP")
DAPopulation$GEO_CSD = as.numeric(substr(DAPopulation$GEO_ID,0,7))
DAPopulation$POP = as.numeric(as.character(DAPopulation$POP))
CityPopulation = subset(DAPopulation, GEO_CSD =='3522014', select = c("DAUID", "POP"))

#Load Shapefiles from correct working directory and merge dissemination area data (same city as in previous step: here Orangeville)
DACanada = spTransform(readOGR(dsn = "DA_Boundaries/lda_000b16a_e.shp"), CRS("+proj=longlat +datum=WGS84"))
DACity = subset(DACanada, CSDNAME == "Orangeville")
Merge_CityPopulation = subset(merge(DACity@data, CityPopulation, by.x = "DAUID", all = TRUE, sort = FALSE), select = c("DAUID","POP"))

#Remove large files to preserve working memory (if needed)
remove(BasePopulation, TotalPopulation, DAPopulation, DACanada)

#Test sample without clustering
CityTest = select(Merge_CityPopulation, POP)
CityTest$POP[is.na(CityTest$POP)] = 0
CityTest.pts = dotsInPolys(DACity, as.integer(CityTest$POP), f="regular")
projection(CityTest.pts) = CRS("+proj=longlat +datum=WGS84")

#Base plot of DA Boundaries
png("DACity_Boundaries.png", units = "in", width = 10, height = 8, res = 1000)
plot(DACity, lwd = 1, border = "black")
dev.off()

#Base plot of CityTest points
png("CityTest_Points.png", units = "in", width = 10, height = 8, res = 1000)
plot(DACity, lwd = 0.1, border = "black")
plot(CityTest.pts, add = T, pch = 16, cex = 0.1, col = "light blue")
dev.off()

#Load additional libraries
library(raster)
library(ggmap)
library(ggplot2)

#Get Google map for location. Here, we're using Orangeville.
GoogleCity = get_googlemap(center = "Orangeville, Ontario, Canada", zoom = 12, scale = 2, format = "png8" ,maptype = "satellite", filename = "orangevilleGoogle", color = "color")

#ggmap to raster function (copied from R. Lovelace on Github)
ggmap_rast = function(map){
  map_bbox = attr(map, 'bb') 
  .extent = extent(as.numeric(map_bbox[c(2,4,1,3)]))
  my_map = raster(.extent, nrow= nrow(map), ncol = ncol(map))
  rgb_cols = setNames(as.data.frame(t(col2rgb(map))), c('red','green','blue'))
  red = my_map
  values(red) = rgb_cols[['red']]
  green = my_map
  values(green) = rgb_cols[['green']]
  blue = my_map
  values(blue) = rgb_cols[['blue']]
  stack(red,green,blue)
}

#convert Google Map capture to raster file. Crop extent to dissemination area boundaries
CityRaster = ggmap_rast(GoogleCity)
CityRaster_cropped = crop(CityRaster, DACity)
DFCity = as.data.frame(CityRaster_cropped)
plotRGB(CityRaster_cropped, 3,2,1)

#Find k for cluster analysis and plot elbow diagram
k_est = (nrow(DFCity)-1)*sum(apply(DFCity,2,var))
  for (i in 2:10){
    set.seed(1234)
    k_est[i] <- sum(kmeans(DFCity, centers=i)$withinss)}
png("Elbow_Graph.png", units = "in", width = 5, height = 5, res = 500)
plot(1:10, k_est, type="b", xlab = "Number of Clusters", ylab = "Within groups sum of squares")
dev.off()

#Land use classification and check for NAs
City_k <- kmeans(DFCity, centers = 3, iter.max = 100, nstart = 10)
City_cluster = as.matrix(City_k$cluster)
CityRaster_cropped2 = as.data.frame(CityRaster_cropped)
CityRaster_cropped2$class = City_k$cluster
CityRaster2 = setValues(CityRaster_cropped, CityRaster_cropped2$class)
plot(CityRaster2, legend = FALSE, col = c("orange", "red", "blue"), main = FALSE, 3)

#Identify residential areas and select the correct value class
mask1 = setValues(CityRaster2, NA)
mask1[CityRaster2 == 1] = 1
plot(mask1, col="orange", legend = FALSE, main = FALSE, 3)

#Convert to shp and export to GIS program for further geoprocessing (R's spatial functions do not work here)
CityPolygon_res = rasterToPolygons(mask1, fun = NULL, n = 4, na.rm = TRUE, dissolve = FALSE)
projection(CityPolygon_res) = CRS("+proj=longlat +datum=WGS84")
CityPolygon_res2 = crop(CityPolygon_res, DACity)
CityPolygon_res3 = gSimplify(CityPolygon_res2, tol = 0.005, topologyPreserve = TRUE)
projection(CityPolygon_res3) = CRS("+init=epsg:3347")
CityPolygon_res3 = gBuffer(CityPolygon_res3, byid=TRUE, width=0)
projection(CityPolygon_res3) = CRS("+proj=longlat +datum=WGS84")
writeOGR(DACity, dsn = "DA_Boundaries", layer = "DACity", driver = "ESRI Shapefile", overwrite_layer = TRUE)
writeOGR(CityPolygon_res2, dsn = "DA_Boundaries", layer = "CityPolygon_res", driver = "ESRI Shapefile", overwrite_layer = TRUE)

#Once residential polygon has been intersected with dissemination areas in GIS program
CityPolygon_intersect = spTransform(readOGR(dsn = "DA_Boundaries/CityPolygon_intersect.shp"),CRS("+proj=longlat +datum=WGS84"))
Merge_CityPolygon = merge(CityPolygon_intersect@data, subset(Merge_CityPopulation, select = c("DAUID","POP")), by.x = "DAUID", all = TRUE, sort = FALSE)

#Link to population by dissemination area and plot points
CityTest.pts2 = dotsInPolys(CityPolygon_intersect, as.integer(Merge_CityPolygon$POP), f="regular")
png("CityTest_Points2.png", units = "in", width = 10, height = 8, res = 1000)
plot(DACity, lwd = 0.1, border = "gray")
plot(CityTest.pts2, add = T, pch = 16, cex = 0.25, col = "light blue")
dev.off()

#Check on Google Maps overlay
png("PointsOverlay.png", units = "in", width = 10, height = 8, res = 1000)
plotRGB(CityRaster_cropped, 3,2,1)
plot(CityTest.pts2, add = T, pch = 16, cex = 0.25, col = "light blue")
dev.off()

#Export for GIS use
projection(CityTest.pts2) = CRS("+proj=longlat +datum=WGS84")
writeOGR(CityTest.pts2, dsn = "DA_Boundaries", layer = "CityTest_pts2", driver = "ESRI Shapefile", overwrite_layer = TRUE)
