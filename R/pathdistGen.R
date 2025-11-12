#' @name pathdistGen
#'
#' @title Generate a stack of path distance raster objects
#' @description Generate a stack of path accumulated distance raster objects
#'
#' @param sf_ob sf object with point geometries
#' @param costras SpatRaster cost raster
#' @param range numeric. Range of interpolation neighborhood
#' @param yearmon character. String specifying the name of the sf_ob
#' @param progressbar logical show progressbar during processing?
#'
#' @return Multi-layer SpatRaster object of path distances
#'
#' @importFrom terra res classify writeRaster c hist
#' @importFrom gdistance transition accCost
#' @importFrom utils setTxtProgressBar
#' @importFrom sf st_coordinates st_crs
#' @export
#'
#' @examples
#' library(sf)
#' sf_ob <- data.frame(rnorm(2))
#' xy   <- data.frame(x = c(4, 2), y = c(8, 4))
#' sf_ob <- st_as_sf(cbind(sf_ob, xy), coords = c("x", "y"))
#'
#' m <- matrix(NA, 10, 10)
#' costras <- rast(m)
#' costras[] <- runif(ncell(costras), min = 1, max = 10)
#' # introduce spatial gradient
#' for (i in 1:nrow(costras)) {
#'   costras[i, ] <- costras[i, ] + i
#'   costras[, i] <- costras[, i] + i
#' }
#'
#' rstack <- pathdistGen(sf_ob, costras, 100, progressbar = FALSE)
pathdistGen <- function(sf_ob, costras, range, yearmon = "default",
                        progressbar = TRUE) {

  if (!inherits(sf_ob, "sf")) {
    stop("sf_ob object must be of class sf")
  }

  if (!identical(sf::st_crs(sf_ob), sf::st_crs(costras))) {
    stop("Point data projection and cost raster projections do not match,
    		 see sf::st_transform")
  }

  ipdw_range <- range / terra::res(costras)[1] / 2 # this is a per cell distance

  # start interpolation
  # calculate conductances hence 1/max(x)
  trans <- gdistance::transition(raster::raster(costras), function(x) 1 / max(x),
    directions = 16)
  i <- 1
  coord <- sf_ob[i, ]
  costsurf <- gdistance::accCost(trans, t(matrix(st_coordinates(coord))))
  costsurf <- terra::rast(costsurf)

  ipdw_dist <- terra::hist(costsurf, plot = FALSE)$breaks[2]
  if (ipdw_dist < ipdw_range) {
    ipdw_dist <- ipdw_range
  }

  if (progressbar == TRUE) {
    pb <- utils::txtProgressBar(max = nrow(sf_ob),
      style = 3)
  }

  for (i in seq_len(nrow(sf_ob))) {
    coord <- sf_ob[i, ]
    costsurf <- gdistance::accCost(trans, t(matrix(st_coordinates(coord))))
    costsurf_reclass <- raster::reclassify(costsurf, c(ipdw_dist, +Inf, NA,
      ipdw_range, ipdw_dist, ipdw_range))
    # raster cells are 1 map unit
    costsurf_scaled <- ((ipdw_range / costsurf_reclass)^2)
    costsurf_scaled_reclass <- raster::reclassify(costsurf_scaled,
      c(-Inf, 1, 0))

    # check output with - zoom(A4,breaks=seq(from=0,to=5,by=1),col=rainbow(5))
    # showTmpFiles()

    raster::writeRaster(costsurf_scaled_reclass,
      filename = file.path(tempdir(), paste(yearmon, "A4ras", i,
        ".grd", sep = "")), overwrite = TRUE, NAflag = -99999)
    if (progressbar == TRUE) {
      setTxtProgressBar(pb, i)
    }
  }
  if (progressbar == TRUE) {
    close(pb)
  }

  # create raster stack
  raster_flist <- list.files(path = file.path(tempdir()),
    pattern = paste(yearmon, "A4ras*", sep = ""),
    full.names = TRUE)
  raster_flist <- raster_flist[grep(".grd", raster_flist, fixed = TRUE)]
  as.numeric(gsub(".*A4ras([0123456789]*)\\.grd$",
    "\\1", raster_flist)) -> fileNum
  raster_flist <- raster_flist[order(fileNum)]
  rstack <- raster::stack(raster_flist)
  rstack <- raster::reclassify(rstack, cbind(-99999, NA))
  file.remove(list.files(path = file.path(tempdir()),
    pattern = paste(yearmon, "A4ras*", sep = ""), full.names = TRUE))

  rstack <- new("pathDist", rstack,
    range = range)

  return(rstack)
}

pathDist   <- setClass("pathDist",
  slots = c(range = "numeric"), contains = "RasterBrick")
