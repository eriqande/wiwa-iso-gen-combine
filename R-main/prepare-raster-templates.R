

## Quick script to prepare the rasters that give the locations of the 
## different clusters found by Ruegg et al.  Project things accordingly
## and sample to correct resolution, etc.


#### LOAD LIBRARIES  ####
stopifnot(
  library(raster, logical.return = TRUE),
  library(rgdal, logical.return = TRUE),
  library(RCurl, logical.return = TRUE),
  library(ncdf4, logical.return = TRUE),
  library(digest, logical.return = TRUE),
  library(dplyr, logical.return = TRUE),
  library(tidyr, logical.return = TRUE)
)



# this is an example raster from Kristina's isotope probability surfaces.  Note that I just
# gave it a projection of +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0
# because it didn't have one.  
iso_raster_template <- readRDS("data/wiwa101-raster.rds")



#### here we read in the range map as spatial polygons, then we turn it ####
# into spatial lines and then we rasterize it at the same resolution
# as the isotope stuff

# read it
wiwa.SPDF <- readOGR(dsn="./private/Wilsonia_pusilla_9151_NS", layer="Wilsonia_pusilla_9151_NS")

# once we have done that, project it, same as iso_raster_template
wiwa.SPDF.proj <- spTransform(wiwa.SPDF,  projection(iso_raster_template))   # here it is projected same as wiwa101 iso raster

# then whittle it down to just the breeding range 
wiwa.breed.SPDF <- wiwa.SPDF.proj[wiwa.SPDF.proj$SEASONAL==2,]

# make spatial lines of it 
wibreed_lines <- as(wiwa.breed.SPDF, "SpatialLinesDataFrame")

# rasterize from the polygons
wi_poly_rast <- rasterize(wiwa.breed.SPDF, iso_raster_template)

# rasterize via lines too.  But don't actually do it because it is a huge bother
# because the shape file has so many lines in it, I think.
# wi_lines_rast <- rasterize(wibreed_lines, iso_raster_template)

# now, wi_poly_rast can be used as a mask, like this
plot(mask(iso_raster_template, wi_poly_rast))



#### Now we want to get the different cluster areas from the Geneland analysis ####
ProbAFile <- "./data/proba.pop.membership.txt"

pp <- read.table(ProbAFile, header=T) %>%
  tbl_df

# here are the names of the regions.  add them on there:
regions <- c("CalSierra", "Basin.Rockies", "Eastern",  "AK.EastBC.AB",  "Wa.To.NorCalCoast",  "CentCalCoast")

names(pp)[1:8] <- c("long", "lat", regions)



# then try rasterizing.  Get a whole list of them and turn them into a rasterBrick
ppb <- brick(lapply(3:8, function(x)  rasterFromXYZ(pp[,c(1,2,x)], crs=CRS("+proj=longlat"))))


# then we have to project it to wiwa.crs:
ppb1 <- projectRaster(ppb, crs=projection(iso_raster_template))

# now get it on the same grid as the wi_poly_rast  (we have to do this before doing the
# masking part because resampling each layer overlays the cells.)
ppb_resam <- resample(ppb1, wi_poly_rast)

# look at em
plot(ppb_resam)

# now, put a 1 in the layer that has the highest posterior prob and NAs elsewhere
ppb_ones <- calc(ppb_resam, fun = function(x) {y<-rep(0, length(x)); y[which.max(x)] <- 1; y })
names(ppb_ones) <- names(ppb_resam)

ppb_masked <- mask(ppb_ones, wi_poly_rast)

# now have a look at that:
plot(ppb_masked)



# stick them all together on the same level
plot(calc(ppb_masked, fun = sum, na.rm = TRUE))

# look! we can sum up the values in any raster layer like this:
cellStats(ppb_masked, sum)

# and we can multiply these things by values like this:
ppb_masked * iso_raster_template


#### Now do a few things with the actual genetic assignments ######
load("data/WIWA-main-carryover-variables.Rda")
GenAss <- WM.gr %>% tbl_df

# now, just get WIWA101 out of there.  Actually, let's just put it all into a good long format
LongAss <- GenAss %>% 
  select(Sample_Name.1, AK.EastBC.AB:Eastern) %>%
  tidyr::gather(., key = region, value = posterior, -Sample_Name.1)

LongAss$region <- factor(LongAss$region, levels = regions)  # to get these in the right order

LongAss2 <- LongAss %>%
  arrange(Sample_Name.1, region)

# now we can get a named vector of posteriors for WIWA101
tmp <- LongAss2 %>%
  filter(Sample_Name.1 == "WIWA101")

posts <- tmp$posterior

# this is not the way to do it, but we will do it this way:
# apply the posterior probs "smeared out within each cluster"
x <- ppb_masked
for(i in names(posts)) {
  x[[i]] <- x[[i]]  * posts[i] / cellStats(x[[i]], sum)
}
# then we just sum all the layers in x into one big layer of genetic
# posterior probs.  These then will sum to 1
xx <- calc(x, fun = sum)


# then get the isotope assignment probs masked by the breeding range
# and then normalized to sum to one.  
iso_probs <- mask(iso_raster_template, wi_poly_rast) 
iso_probs <- iso_probs / cellStats(iso_probs, sum)

# and now we take the product of those with the isotopes, and normalize so they
# sum to 1
final <- xx * iso_probs
final <- final / cellStats(final, sum)


# compare the different data types.
par(mfrow = c(3,1))
plot(iso_probs)
plot(xx)
plot(final)


# so, what we need to send to the Rstudio server for Kristina is the 
# rasterBrick with the regions, and the data frame WM.gr, but we will want
# to name these different things.
writeRaster(ppb_masked, filename = "RegionBrick.nc", bylayer = FALSE)   # can read this with brick("RegionBrick.grd")
saveRDS(LongAss2, file = "genetic-posteriors.rds", compress = "xz")
