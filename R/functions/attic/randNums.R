# function that generates n random integers with fixed sum s within boundaries lb and hb

randNums <- function(n,minVal,maxVal,s){
  hit <- FALSE
  while(hit==FALSE){
    total <- 0
    count <- 1
    nums <- rep(NA, n)
    
    while(total < s && count < n){
      r <- round(runif(1, min=minVal, max=maxVal))
      total <- total + r
      
      nums[count] <- r
      count <- count + 1
    }
    #nums[count:n] <- 0
    if(total == s){
      hit <- TRUE
    }
  }
  return(nums)
}
