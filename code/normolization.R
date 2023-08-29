centering <- function(data) {
  col_names <- colnames(data)
  row_names <- rownames(data)
  data <- as.matrix(data)
  df <- apply(data, 2, function(x) x - mean(x, na.rm=TRUE))
  colnames(df) <- col_names
  rownames(df) <- row_names
  return(df)
}

auto_scale <- function(data) {
  col_names <- colnames(data)
  row_names <- rownames(data)
  data <- as.matrix(data)
  df <- apply(data, 2, function(x) (x - mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE))
  colnames(df) <- col_names
  rownames(df) <- row_names
  return(df)
}

range_scale <- function(data) {
  col_names <- colnames(data)
  row_names <- rownames(data)
  data <- as.matrix(data)
  df <- apply(data, 2, function(x) (x - mean(x, na.rm=TRUE))/(max(x) - min(x)))
  colnames(df) <- col_names
  rownames(df) <- row_names
  return(df)
}

minmax_scale <- function(data) {
  col_names <- colnames(data)
  row_names <- rownames(data)
  data <- as.matrix(data)
  df <- apply(data, 2, function(x) (x - min(x))/(max(x) - min(x)))
  colnames(df) <- col_names
  rownames(df) <- row_names
  return(df)
}

pareto_scale <- function(data) {
  col_names <- colnames(data)
  row_names <- rownames(data)
  data <- as.matrix(data)
  df <- apply(data, 2, function(x) (x - mean(x, na.rm=TRUE))/sqrt(sd(x, na.rm=TRUE)))
  colnames(df) <- col_names
  rownames(df) <- row_names
  return(df)
}

vast_scale <- function(data) {
  col_names <- colnames(data)
  row_names <- rownames(data)
  data <- as.matrix(data)
  df <- apply(data, 2, function(x) (x - mean(x, na.rm=TRUE))*mean(x, na.rm=TRUE)/var(x, na.rm=TRUE))
  colnames(df) <- col_names
  rownames(df) <- row_names
  return(df)
}

level_scale <- function(data) {
  col_names <- colnames(data)
  row_names <- rownames(data)
  data <- as.matrix(data)
  df <- apply(data, 2, function(x) (x - mean(x, na.rm=TRUE))/mean(x, na.rm=TRUE))
  colnames(df) <- col_names
  rownames(df) <- row_names
  return(df)
}

log_transform <- function(data, base=2) {
  col_names <- colnames(data)
  row_names <- rownames(data)
  data <- as.matrix(data)
  if (min(data) < 0) {
    stop("Log transformation could not process negative data")
  }
  if (base == "e") {
    base = exp(1)
  } else {
    base = as.numeric(base) 
  }
  df <- log(1+data, base = base)
  colnames(df) <- col_names
  rownames(df) <- row_names
  return(df)
}

power_transform <- function(data) {
  col_names <- colnames(data)
  row_names <- rownames(data)
  data <- as.matrix(data)
  if (min(data) < 0) {
    stop("Power transformation could not process negative data")
  }
  df <- sqrt(data)
  colnames(df) <- col_names
  rownames(df) <- row_names
  return(df)
}


#average
object = paste0(
  '~0+',
  paste0(
    "data[,",
    1:length(x = group.by),
    "]",
    collapse = ":"
  )
)

tmp <- apply(combine_all, 1,function(x){
  
})