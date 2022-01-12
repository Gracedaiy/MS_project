
evaluation <- function(k){ 
    rv = rnorm(10000, 0, 1)
    write.csv(rv, file = paste('data', k, '.csv', sep = ""), row.names = FALSE)
}


args = commandArgs(trailingOnly = TRUE)
k = as.numeric(args[1])
evaluation(k)
