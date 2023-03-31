## IGNORE_RDIFF_BEGIN
library(rtop)
library(sf)
set.seed(1)
options(error = recover)
rpath = system.file("extdata",package="rtop")
observations = st_read(rpath,"observations")
# Create a column with the specific runoff:
observations$obs = observations$QSUMMER_OB/observations$AREASQKM
predictionLocations = st_read(rpath,"predictionLocations")
## IGNORE_RDIFF_END

params = list(gDist = TRUE, cloud = FALSE)
# Create a column with the specific runoff:
observations$obs = observations$QSUMMER_OB/observations$AREASQKM
# Build an object
rtopObj = createRtopObject(observations, predictionLocations, 
                           params = params, formulaString = "obs ~1")


rtopObj = rtopFitVariogram(rtopObj, iprint = -1)

# Predicting at prediction locations
rtopObj = rtopKrige(rtopObj)

# Cross-validation
rtopObj = rtopKrige(rtopObj,cv=TRUE)
print(cor(rtopObj$predictions$observed,rtopObj$predictions$var1.pred), 4)
