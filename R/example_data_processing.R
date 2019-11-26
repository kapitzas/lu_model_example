require(raster)
mask <- readRDS("~/OneDrive - The University of Melbourne/land use model/MDAP examples/data/mask_eur.rds")
covs <- readRDS("~/OneDrive - The University of Melbourne/land use model/MDAP examples/data/cov.rds")
lu <- readRDS("~/OneDrive - The University of Melbourne/land use model/MDAP examples/data/lu.rds")

sizes <- table(mask[])
inds <- which(!is.na(mask[]))
countries <- mask[inds]
example_country <- which(countries%in%c(32)) #Germany

mask[which(!mask[]%in%c(32))] <- NA
mask <- trim(mask)

covs <- covs[example_country,]
lu <- lu[example_country,]
saveRDS(covs, "./data/cov.rds")
saveRDS(lu, "./data/lu.rds")
saveRDS(mask, "./data/mask.rds")
