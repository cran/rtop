library(rtop)
library(sf)
options(error = recover)
rpath = system.file("extdata",package="rtop")
observations = st_read(rpath,"observations")
# Create a column with the specific runoff:
observations$obs = observations$QSUMMER_OB/observations$AREASQKM
predictionLocations = st_read(rpath,"predictionLocations")

params = list(gDist = TRUE, cloud = FALSE)
# Create a column with the specific runoff:
observations$obs = observations$QSUMMER_OB/observations$AREASQKM
# Build an object
rtopObj = createRtopObject(observations, predictionLocations, 
                           params = params)


rtopObj = rtopFitVariogram(rtopObj)

# Predicting at prediction locations
rtopObj = rtopKrige(rtopObj)

# Cross-validation
rtopObj = rtopKrige(rtopObj,cv=TRUE)
cor(rtopObj$predictions$observed,rtopObj$predictions$var1.pred)
