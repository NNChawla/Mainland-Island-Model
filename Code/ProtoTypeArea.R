#####################
sample <- c(1, 2, 4, 5, 7, 8)

r <- y[[1]][sample]
s <- y[[2]][sample]
g <- y[[3]][sample]
a <- y[[4]][sample, sample]
init.x <- y[[5]][sample]

######################

mougi_model <- function(t,x,parms){
  dx <- x * (r - s*x + g * (a %*% x) - (t(a) %*% x))
  list(dx)
}

n.integrate <- function(time=time, init.x= init.x, model=model){
  t.out <- seq(time$start,time$end,length=time$steps)
  as.data.frame(lsoda(init.x, t.out, model, parms = parms))
}

# Integration window
time <- list(start = 0, end = 100, steps = 100)
# dummy variable for lvm() function defined above
parms <- c(0) ### dummy variable (can have any numerical value)


out <- n.integrate(time, init.x, model = mougi_model)

#list(L, xxx, n, r, s, g, a, init.x)