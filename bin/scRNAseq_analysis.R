custom_functions = "./bin/custom_functions/"

library(devtools)
library(Antler)

# load custom functions
sapply(list.files(custom_functions, full.names = T), source)


