# 09_ConstrainCO2Flux_getSAUdata.R
# created August 25, 2022 by Jamie Collins, jcollins@edf.org

# a helper script to retrieve Sea Around Us data using the SAU API, now that the
# rsearoundus package seems to have disappeared from CRAN

# the functions were copied from the source files still available on Aug. 25
# at https://github.com/robsalasco/rseaaroundus

getapibaseurl <- function() "http://api.seaaroundus.org"

# call API with GET and return data
callapi <- function(path, args = list(), ...) {
  conn <- crul::HttpClient$new(
    url = getapibaseurl(),
    headers = list(`X-Request-Source` = "r"),
    opts = list(followredirects = TRUE, ...)
  )
  resp <- conn$get(path = path, query = args)
  resp$raise_for_status()
  jsonlite::fromJSON(resp$parse("UTF-8"))$data
}

catchdata("eez", 840, measure="tonnage", dimension="gear")

catchdata <- function(region, id, measure="tonnage", dimension="taxon",
                      limit=10, chart=FALSE, ...) {
  
  # create url path and query parameters
  path <- paste("api/v1", region, measure, dimension, "", sep="/")
  args <- list(region_id = id, limit = limit)
  
  # call API
  data <- callapi(path, args, ...)
  
  if (is.null(data)) return(data.frame(NULL))
  
  # extract data from response
  values <- data$values
  years <- values[[1]][,1]
  cols <- lapply(values, function(v) { v[,2] })
  
  # create dataframe
  df <- data.frame(years, cols, stringsAsFactors = FALSE)
  df[is.na(df)] <- 0
  colnames(df) <- tolower(c("years", data$key))
  
  # return dataframe
  if (!chart) {
    return(df)
    
    # return chart
  } else {
    ylabel <- ifelse(measure == "tonnage", "Catch (t x 1000)",
                     "Real 2005 value (million US$)")
    charttitle <- toupper(paste(region, id, measure, "by",
                                dimension, sep=" "))
    
    df <- Reduce(function(...) rbind(...), lapply(colnames(df), function(name) {
      coldata <- df[,name] / ifelse(measure=="tonnage", 10^3, 10^6)
      data.frame(year=years, data=coldata, dim=rep(name, nrow(df)))
    }))
    
    spectral <- c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08b",
                  "#ffffbf", "#e6f598", "#abdda4", "#66c2a5",
                  "#3288bd", "#5e4fa2", '#666', '#f88', '#88f', '#8f8', '#800', '#080',
                  '#008', '#888', '#333')
    chartcolors <- rep(spectral, 10)
    
    ggplot(df, aes(year, data)) +
      geom_area(aes(fill=dim), position="stack") +
      theme(legend.position="top", legend.key.size=ggplot2::unit(0.5, "cm")) +
      scale_x_continuous(breaks=seq(min(years), max(years), 10),
                         expand=c(0, 0)) +
      scale_y_continuous(breaks=pretty_breaks(n=10), expand=c(0, 0)) +
      guides(fill=guide_legend(title=NULL, direction="horizontal",
                               byrow=TRUE, ncol=5)) +
      scale_fill_manual(values = chartcolors) +
      labs(y = ylabel, x = "Year") +
      ggtitle(charttitle)
  }
}

listregions <- function(region) {
  
  # create url
  path <- paste("api/v1", region, "?nospatial=true", sep="/")
  
  # call API
  data <- callapi(path)
  
  # # create url
  # baseurl <- getapibaseurl()
  # url <- paste(baseurl, region, "?nospatial=true", sep="/")
  # 
  # # call API
  # data <- callapi(url)
  
  # extract data from response
  df <- data.frame(data, row.names='id')
  return(df)
}

# retrieve the data we are looking for
# the codes used for the SAU regions appear to have *some* overlap with the UN M49
# country codes, but there are some weird idiosyncrasies

SAUregions <- listregions("eez")

SAUcatchdata <- as.data.frame(matrix(data = NA,
                                     nrow = 13,
                                     ncol = nrow(SAUregions)))
rownames(SAUcatchdata) <- c("Country","SAU_code",c(2009:2018),"Sum_2009-2018")
SAUcatchdata[1,] <- SAUregions$title
SAUcatchdata[2,] <- as.numeric(rownames(SAUregions))

for (i in 1:ncol(SAUcatchdata)) {
  
  these.catchdata <- catchdata("eez", SAUcatchdata[2,i], measure="tonnage", dimension="gear")
  catchdata.Yearswanted.ind <- which(these.catchdata$years %in% c(2009:2018))
  bottomtrawlData.yearsWanted <- these.catchdata[catchdata.Yearswanted.ind,
                                                 c("bottom trawl")]
  
  if (!is.null(bottomtrawlData.yearsWanted)) {
    
    SAUcatchdata[3:12,i] <- bottomtrawlData.yearsWanted
    SAUcatchdata[13,i] <- mean(bottomtrawlData.yearsWanted, na.rm = T)
  
  }
  
}

# save a copy of this data now that we've collected it

SAU_benthicCatchdata <- SAUcatchdata
SAU_benthicCatchdata <- as.data.frame(t(SAU_benthicCatchdata))
  
save(SAU_benthicCatchdata, file = "data/global_trawling/derived/SAU_catchdata/SAU_benthicCatchdata_2009_2018.RData")

SAU_benthicCatchdata[,3:13] <- as.numeric(unlist(SAU_benthicCatchdata[,3:13]))

SAU_benthicCatchdata.sorted <- SAU_benthicCatchdata[sort(SAU_benthicCatchdata$`Sum_2009-2018`, decreasing = TRUE,
     index.return = TRUE, na.last = TRUE)$ix,c(1,13)]

# need to do some combining, manually will be easiest ...


