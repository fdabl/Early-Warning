library('MESS')
library('lemon')
library('readr')
library('dplyr')
library('qgraph')
library('ggplot2')
library('deSolve')
library('nleqslv')
library('numDeriv')
library('rootSolve')
library('latex2exp')
library('gridExtra')
library('rootSolve')
library('RColorBrewer')
source('Code/helpers.R')


ptheme <- theme(
  plot.title = element_text(hjust = .5),
  text = element_text(size = 14),
  legend.key = element_blank(),
  legend.title = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  panel.border = element_blank(),
  legend.position = 'top',
  axis.line = element_blank(),
  strip.background = element_blank()
)

cm <- 1.25
cl <- 1.25
ca <- 1.5

# Creates the subplots in Figure 1
# I combined them in Adobe Illustrator afterwards
####################################################################
###### Figure 1: Critical Slowing Down in Unidimensional Systems ###
####################################################################
bifurcation <- function(equation, cs, x = seq(0, 10, .001), r = 1, K = 10, b = 1) {
  n <- length(cs)
  roots <- matrix(NA, nrow = n, ncol = 4)
  
  for (i in seq(n)) {
    root <- uniroot.all(function(x) equation(x, cs[i], r = r, K = K), c(0, 10))
    roots[i, ] <- c(rep(NA, 4 - length(root)), root)
  }
  
  roots
}

logistic_sn <- function(x, cparam, r = 1, K = 10, b = 1) r*x * (1 - x/K) - cparam * x^2 / (b + x^2)

pdf('Figures/Fig-1-Bifurcation.pdf', width = 7, height = 6)
cs <- seq(1, 3, length.out = 500)
roots <- bifurcation(logistic_sn, cs)

u <- roots[, 3]
csel <- which(cs > 1.786 & cs < 2.6)
ix <- which(roots[, 4] < 3)[1] - 1

plot(
  cs[seq(ix)], roots[seq(ix), 4], type = 'l', ylim = c(0, 10), lwd = 2, xlim = c(1, 3),
  axes = FALSE, xlab = 'Predation Rate', ylab = expression(x^'*'), main = 'Bifurcation Diagram', font.main = 1,
  cex.main = 1.75, cex.lab = 1.5, cex.axis = 2
)
lines(cs[seq(ix + 1, length(cs))], roots[seq(ix+1, length(cs)), 4], lwd = 2)
lines(cs, roots[, 2], lwd = 2)
lines(cs[csel], roots[, 3][csel], lty = 'dashed', lwd = 2)
lines(c(1, 3), c(0, 0), lty = 2, lwd = 2)
axis(1)
axis(2, las = 2)
points(1.786, roots[csel[1], 2], pch = 20, col = 'gray76', cex = 2)
points(1.786, roots[csel[1], 2], cex = 1.5)
points(2.6, roots[csel[length(csel)], 4] - .25, pch = 20, col = 'gray76', cex = 2)
points(2.6, roots[csel[length(csel)], 4] - .25, cex = 1.5)

# arrows(1.3, 3, 1.3, 8, length = .1, lwd = 1.5)
# arrows(1.7, 3, 1.7, 7.5, length = .1, lwd = 1.5)
# arrows(2.1, 3, 2.1, 6.5, length = .1, lwd = 1.5)
# arrows(2.2, 2.2, 2.2, 0.8, length = .1, lwd = 1.5)
# arrows(2.5, 3, 2.5, 0.8, length = .1, lwd = 1.5)
# arrows(2.8, 5, 2.8, 0.8, length = .1, lwd = 1.5)
dev.off()


pdf('Figures/Fig-1-Potentials.pdf', width = 9, height = 4)
par(mfrow = c(1, 2))
cols <- RColorBrewer::brewer.pal(3, 'Set1')

potential_grazing <- function(x, cparam = 1, r = 1, K = 10) {
  cparam*x - cparam*atan(x) + r*x^2*(2*x - 3*K) / (6*K)
}

x <- seq(0, 10, .001)
ballcex <- 5.575

# Far from the bifurcation
plot(
  x, potential_grazing(x, cparam = 1.9), type = 'l',
  main = 'Higher Resilience', xaxs = 'i', yaxs = 'i', xlab = 'x',
  axes = FALSE, ylab = 'Potential', col = cols[2], lwd = 1.5, ylim = c(-4, 4),
  font.main = 1, cex.main = cm, cex.lab = cl, cex.axis = ca
)

y <- potential_grazing(x, cparam = 1.9)

polygon(c(min(x), x, max(x)), c(-10, y, -10), col = grDevices::adjustcolor(cols[2], alpha.f = .60))
axis(1)
axis(2, las = 2, at = seq(-4, 4, 2))

p <- uniroot.all(function(x) logistic_sn(x, c = 1.9), interval = c(0, 10))[-1]
points(p[3], -1.98, pch = 20, cex = ballcex, col = 'gray70')

# Close to the bifurcation
plot(
  x, potential_grazing(x, cparam = 2.4), type = 'l',
  main = 'Lower Resilience', xaxs = 'i', yaxs = 'i', xlab = 'x',
  axes = FALSE, ylab = 'Potential', col = cols[2], lwd = 1.5, ylim = c(-4, 4),
  font.main = 1, cex.main = cm, cex.lab = cl, cex.axis = ca
)

y <- potential_grazing(x, cparam = 2.4)

polygon(c(min(x), x, max(x)), c(-4, y, -5), col = grDevices::adjustcolor(cols[2], alpha.f = .60))
axis(1)
axis(2, las = 2, at = seq(-4, 4, 2))

p1 <- uniroot.all(function(x) logistic_sn(x, c = 2.4), interval = c(0, 10))[-1]
points(p1[3], 0.77, pch = 20, cex = ballcex, col = 'gray70')
dev.off()

LV_saddlenode <- function(t, X, params) {
  r <- params[['r']]
  c <- params[['c']]
  b <- params[['b']]
  K <- params[['K']]
  esd <- params[['esd']]
  e <- rnorm(1, 0, esd)
  dX <- r * X * (1 - X / K) - (c)*X^2 / (b^2 + X^2) + e
  list(dX)
}

get_params_saddle_node <- function(c, esd = 0.40) {
  list('r' = 1, 'c' = c, 'b' = 1, 'K' = 10, 'esd' = esd)
}

# Similar to van Nes & Scheffer (2007) plot
create_impulse <- function(y0, times1 = seq(0, 30, .001), times2 = seq(30, 50, .001), cparam = 2) {
  out1 <- ode(y = y0, times = times1, func = LV_saddlenode, parms = get_params_saddle_node(cparam))
  out2 <- ode(y = out1[nrow(out1), 2] * 0.50, times = times2, func = LV_saddlenode, parms = get_params_saddle_node(cparam))
  out <- rbind(out1, out2)
}

cbind(
  uniroot.all(function(x) logistic_sn(x, c = 1.9), interval = c(0, 10))[-1],
  uniroot.all(function(x) logistic_sn(x, c = 2.4), interval = c(0, 10))[-1]
)

set.seed(1)
y01 <- 7.5160459
y02 <- 6.2635324
out1 <- create_impulse(y0 = y01, cparam = 1.90)
out2 <- create_impulse(y0 = y02, cparam = 2.40)

set.seed(2)
x1 <- ode(
  y = y01, times = seq(0, 1000, .01), func = LV_saddlenode,
  parms = get_params_saddle_node(1.90, 0.1)
)

x2 <- ode(
  y = y02, times = seq(0, 1000, .01), func = LV_saddlenode,
  parms = get_params_saddle_node(2.40, 0.1)
)

x1 <- x1[seq(1, nrow(x1), 100), 2]
x2 <- x2[seq(1, nrow(x2), 100), 2]

f <- function(x) {
  cbind(x[-1], x[-length(x)]) / mean(x)
}


cols <- brewer.pal('Set1', n = 4)
pdf('Figures/Fig-1-Higher-Resilience-EWS.pdf', width = 8, height = 4)
par(mfrow = c(1, 2))
plot(
  out1[, 1], out1[, 2] / y01, type = 'l', axes = FALSE,
  xlab = 'Time t', ylab = expression(x[t] / bar(x)), ylim = c(0, 1.2), col = cols[2], lty = 1,
  main = 'Recovery Rate', lwd = 2,
  font.main = 1, cex.main = cm, cex.lab = cl, cex.axis = ca
)
lines(c(0, 50), c(1, 1), lty = 2)
lines(out1[, 1], out1[, 2] / y01, lwd = 2, col = cols[2])
axis(1)
axis(2, las = 1)

plot(
  f(x1), pch = 16, cex = 0.75, axes = FALSE,
  ylim = c(0.80, 1.20), xlim = c(0.80, 1.20), main = 'Autocorrelation',
  xlab = expression(x[t] / bar(x)), ylab = expression(x[t+1] / bar(x)),
  col = scales::alpha(cols[2], 0.4),
  font.main = 1, cex.main = cm, cex.lab = cl, cex.axis = ca
)
axis(1)
axis(2, las = 2)
dev.off()



pdf('Figures/Fig-1-Lower-Resilience-EWS.pdf', width = 8, height = 4)
par(mfrow = c(1, 2))
plot(
  out2[, 1], out2[, 2] / y02, type = 'l', axes = FALSE,
  xlab = 'Time t', ylab = expression(x[t] / bar(x)), ylim = c(0, 1.2), col = cols[2], lty = 1,
  main = 'Recovery Rate', lwd = 2,
  font.main = 1, cex.main = cm, cex.lab = cl, cex.axis = ca
)
lines(c(0, 50), c(1, 1), lty = 2)
lines(out2[, 1], out2[, 2] / y02, lwd = 2, col = cols[2])
axis(1)
axis(2, las = 1)

plot(
  f(x2), pch = 16, cex = 0.75, axes = FALSE,
  ylim = c(0.80, 1.20), xlim = c(0.80, 1.20), main = 'Autocorrelation',
  xlab = expression(x[t] / bar(x)), ylab = expression(x[t+1] / bar(x)),
  col = scales::alpha(cols[2], 0.4),
  font.main = 1, cex.main = cm, cex.lab = cl, cex.axis = ca
)
axis(1)
axis(2, las = 2)
dev.off()


#################################################################
###### Figure 2: Lotka-Volterra Bifurcation and Transition ######
#################################################################
x <- c(1, 2, 3, 4)
C <- matrix(c(-.2, .04, -.2, -.2, 
              .04, -.2, -.2, -.2,
              -.2, -.2, -.2, .04, 
              -.2, -.2, .04, -.2), 4, 4, byrow = T)

dxdt <- function(x) {
  dx <- r * x + x * as.vector(C %*% x) + 1
  dx
}

stress_set <- seq(0.5, 1.5, by = 0.01)
searchSpace <- seq(1, 7, 0.5)
grid <- expand.grid(searchSpace, searchSpace )
grid <- cbind(grid[,1], grid[,1], grid[,2], grid[,2])
ansBroy <- list()

for (i in 1:length(stress_set)) {
  r <-  c(1, 1, stress_set[i], stress_set[i]) # loop over r3 & r4
  res <- searchZeros(as.matrix(grid), fn = dxdt, method = 'Broyden', global = 'dbldog')$x
  
  # remove any negative roots
  if (any(res < 0)) {
    ind <- -which( apply(res, 1, function(x) any(x < 0)))
    ansBroy[[i]] <- rbind(res[ind, ])
  } else{
    ansBroy[[i]] <- res
  }
}

## proper counting
mincut <- min(which(sapply(ansBroy, nrow) == 3)) - 1
maxcut <- max(which(sapply(ansBroy, nrow) == 3)) + 1
t1 <- stress_set[1:mincut]
t2 <- stress_set[(mincut+1):(maxcut-1)]
t3 <- stress_set[maxcut:length(stress_set)]
tbelow <- c(t1, t2)
tabove <- c(t2, t3)

## arranging state variables
below1 <- sapply(ansBroy[1:mincut], rbind)
below2 <- sapply(ansBroy[(mincut+1):(maxcut-1)], function(x)x[3,])
below <- cbind(below1, below2)
above2 <- sapply(ansBroy[(mincut+1):(maxcut-1)], function(x)x[1,])
above3 <- sapply(ansBroy[maxcut:length(stress_set)], rbind)
above <- cbind(above2, above3)
unstable <- sapply(ansBroy[(mincut+1):(maxcut-1)], function(x)x[2,])

eigbelow <- array(NA, dim = dim(below))
jacbelow <- Abelow <- list()

for(i in 1:ncol(eigbelow)){
  r <-  c(1, 1, tbelow[i], tbelow[i])
  jacbelow[[i]] <- jacobian(dxdt, x = below[,i])
  eigbelow[,i] <- eigen(jacbelow[[i]])$values
  Abelow[[i]] <- exp(jacbelow[[i]])
  
}

eigabove <- array(NA, dim = dim(above))
jacabove <- Aabove <- list()
for(i in 1:ncol(eigabove)){
  r <-  c(1, 1, tabove[i], tabove[i])
  jacabove[[i]] <- jacobian(dxdt, x = above[,i])
  eigabove[,i] <- eigen(jacabove[[i]])$values
  Aabove[[i]] <- exp(jacabove[[i]])
}

eigunstable <- array(NA, dim = dim(unstable))
for(i in 1:ncol(eigunstable)){
  r <-  c(1, 1, t2[i], t2[i])
  jacunstable <- jacobian(dxdt, x = unstable[,i])
  eigunstable[,i] <- eigen(jacunstable)$values
}


pdf("Figures/Fig-2-Lotka-Volterra-Plot.pdf", width = 12, height = 5)
par(mfrow = c(1, 2))
cols <- RColorBrewer::brewer.pal(n = 4, 'Set1')
plot(NA, xlim = c(.5, 1.5), ylim = c(0, 10), 
     xlab = "r", ylab = "", axes = FALSE,
     main = 'Bifurcation Diagram', font.main = 1,
     cex.main = 1.50, cex.lab = 1.25, cex.axis = 1.5
)
mtext(expression(x^"*"), side = 2, line = 2.5, cex = 1.25)
rect(xleft = stress_set[mincut+1], 0, xright = stress_set[maxcut-1], 10, col = rgb(0.5,0.5,0.5,1/4), border = NA)
lines(tbelow, below[1,], col = cols[1], lwd = 1.5)
lines(tbelow, below[3,], col = cols[2], lwd = 1.5)

lines(tabove, above[1,], col = cols[1], lwd = 1.5)
lines(tabove, above[3,], col = cols[2], lwd = 1.5)

lines(t2, unstable[1,], col = 'grey46', lty = 2, lwd = 1.5)
lines(t2, unstable[3,], col = 'grey46', lty = 2, lwd = 1.5)


legend(
  'topleft', legend = c(expression(x[1]~','~x[2]), expression(x[3]~','~x[4])),
  col = cols[c(1, 2)], lwd = 1.5, bty = 'n', cex = 1.10
)

axis(1, at = seq(.5, 1.5, .1), cex.axis = 1.10)
axis(2, las = 2, cex.axis = 1.10)


set.seed(6)
lotka <- create_lotka_obj(nr_transition_days = 20, nr_baselinedays = 50, noisetype = 'additive', noise = 4)
lotkas <- subsample(lotka, freq = 9)

plot_lotka_rs <- function(
  lotka, xat = c(0, 25, 50, 70, 90) * 100, xlabels = xat / 100, cex.axis = 1, legend.cex = 0.90, ...
) {
  
  cols <- RColorBrewer::brewer.pal(n = 4, 'Set1')
  
  plot(
    lotka$rs, xlab = 'Days',
    ylab = 'r(t)', axes = FALSE, type = 'l', lwd = 2, col = 'black',
    ylim = c(0.6, 1.8), font.main = 1, ...
  )
  
  d <- lotka$dat[, -1]
  d <- (d - min(d)) / (max(d) - min(d)) + 0.60
  
  lines(d[, 1], lwd = 1.5, col = grDevices::adjustcolor(cols[1], alpha.f = .75))
  lines(d[, 2], lwd = 1.5, col = grDevices::adjustcolor(cols[2], alpha.f = .75))
  lines(d[, 3], lwd = 1.5, col = grDevices::adjustcolor(cols[3], alpha.f = .75))
  lines(d[, 4], lwd = 1.5, col = grDevices::adjustcolor(cols[4], alpha.f = .75))
  axis(1, at = xat, labels = xlabels, cex.axis = cex.axis)
  axis(2, las = 2, at = seq(0.6, 1.8, .2), cex.axis = cex.axis)
  
  legend(
    'topleft',
    col = cols, lwd = 1.5, bty = 'n', cex = legend.cex,
    legend = c(expression(x[1]), expression(x[2]), expression(x[3]), expression(x[4]))
  )
}

plot_lotka_rs(
  lotkas, main = 'Example Time Series and r(t)',
  cex.main = 1.50, cex.lab = 1.25, cex.axis = 1.20, legend.cex = 1.10
)
dev.off()



######################################################################
#### Figure 3: Critical Slowing Down in Multidimensional Systems #####
######################################################################
plot_raw <- function(dat, cex = 0.75, alpha = 1) {
  cols <- RColorBrewer::brewer.pal(4, 'Set1')
  plot(
    dat[, 1], dat[, 2], pch = 16, cex = cex, axes = FALSE, type = 'l',
    ylim = c(0, 7), xlab = 'Time', ylab = 'Mood', col = scales::alpha(cols[1], alpha),
    cex.main = 1.75, cex.lab = 1.5, cex.axis = 2
  )
  lines(dat[, 1], dat[, 3], pch = 16, cex = cex, col = scales::alpha(cols[2], alpha))
  lines(dat[, 1], dat[, 4], pch = 16, cex = cex, col = scales::alpha(cols[3], alpha))
  lines(dat[, 1], dat[, 5], pch = 16, cex = cex, col = scales::alpha(cols[4], alpha))
  axis(1)
  axis(2, las = 2)
  
  legend(
    43, 4.5, legend = c(expression(x[1]), expression(x[2]), expression(x[3]), expression(x[4])),
    col = cols, lwd = 1, bty = 'n'
  )
}


plot_autocor <- function(dat, cex = 0.75, alpha = 0.25) {
  cols <- RColorBrewer::brewer.pal(4, 'Set1')
  
  f <- function(dat, sel) {
    d <- dat[, sel]
    cbind(d[-1], d[-1000]) / mean(d)
  }
  
  plot(
    f(dat, 4), pch = 16, cex = cex, axes = FALSE,
    ylim = c(0, 2), xlim = c(0, 2),
    xlab = expression(x[t] / bar(x)), ylab = expression(x[t+1] / bar(x)),
    col = scales::alpha(cols[3], alpha),
    cex.main = 1.75, cex.lab = 1.5, cex.axis = 2
  )
  points(f(dat, 5), pch = 16, cex = cex, col = scales::alpha(cols[4], alpha))
  axis(1)
  axis(2, las = 2)
  
  legend(
    'topleft', legend = c(expression(x[3]), expression(x[4])),
    col = cols[3:4], pch = 16, bty = 'n'
  )
}


plot_graph <- function(mat, directed = FALSE, maximum = .8) {
  layout <- rbind(c(0,1), 
                  c(1,1), 
                  c(0,0), 
                  c(1,0))
  
  qgraph::qgraph(mat, 
                 layout = layout,
                 directed = directed, 
                 # edge.color = "skyblue",
                 edge.labels = TRUE,
                 edge.label.cex = 3,
                 # edge.label.color = "skyblue",
                 labels = c(expression(x[1]), expression(x[2]),
                            expression(x[3]), expression(x[4])), 
                 # lty = m_lty, 
                 vsize = 18, 
                 esize = 16,
                 asize = 12,
                 mar = c(8, 10, 8, 8), 
                 maximum = maximum,
                 fade = TRUE,
                 posCol = cols[2],
                 negCol = cols[1]
  )
}


## I further processe this in Adobe Illustrator
pdf('Figures/Fig-3-Multidimensional-CSD.pdf', width = 12, height = 6)
cols <- RColorBrewer::brewer.pal(n = 4, 'Set1')
par(mfrow = c(2, 3))
set.seed(4)
rs1 <- rep(0.50, 1000)
dat1 <- generate_data(
  time = length(rs1) * .1, rs = rs1,
  timestep = 0.1, noise_sd = 0.50, noisetype = 'additive'
)

set.seed(4)
rs2 <- rep(1.18, 1000)
dat2 <- generate_data(
  time = length(rs2) * .1, rs = rs2,
  timestep = 0.1, noise_sd = 0.50, noisetype = 'additive'
)


lmat <- matrix(1:12, 3, 4, byrow = TRUE)
lmat <- lmat - 1
lo <- layout(lmat, width = c(.1, 1, 1, 1), heights = c(.1, 1, 1, 1))

par(mar = rep(0, 4))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
text(0.5, 0.5, 'Example Time-series', cex = 2.5)

plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
text(0.5, 0.5, 'Autocorrelation', cex = 2.5)

plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
text(0.5, 0.5, 'Cross-correlation', cex = 2.5)

plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
text(0.5, 0.5, 'Higher Resilience', cex = 2, srt = 90)

par(mar = c(5.1, 4, 2, 2))
plot_raw(dat1)
plot_autocor(dat1)
plot_graph(cor(dat1[, -1]))

par(mar = rep(0, 4))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
text(0.5, 0.5, 'Lower Resilience', cex = 2, srt = 90)

par(mar = c(5.1, 4, 2, 2))
plot_raw(dat2)
plot_autocor(dat2)
plot_graph(cor(dat2[, -1]))
dev.off()


################################################################
#### Figure 4: Example Time-series for Simulation Conditions ###
################################################################
set.seed(1)
# lotka6 <- create_lotka_obj(nr_transition_days = 10, nr_baselinedays = 25, noisetype = 'additive', noise = 6)
# lotka8 <- create_lotka_obj(nr_transition_days = 25, nr_baselinedays = 50, noisetype = 'additive', noise = 8)
# lotka10 <- create_lotka_obj(nr_transition_days = 50, nr_baselinedays = 100, noisetype = 'additive', noise = 10)

lotka6 <- create_lotka_obj(nr_transition_days = 25, nr_baselinedays = 50, noisetype = 'additive', noise = 10)
lotka8 <- create_lotka_obj(nr_transition_days = 25, nr_baselinedays = 50, noisetype = 'additive', noise = 8)
lotka10 <- create_lotka_obj(nr_transition_days = 25, nr_baselinedays = 50, noisetype = 'additive', noise = 6)

lotka6s <- subsample(lotka6, freq = 90)
lotka8s <- subsample(lotka8, freq = 180)
lotka10s <- subsample(lotka10, freq = 900)

xlabels1 <- c(0, 25, 35, 55)
xlabels2 <- c(0, 25, 50, 75, 95)
xlabels3 <- c(0, 50, 100, 150, 170)

xlabels1 <- xlabels2
xlabels3 <- xlabels2

pdf('Figures/Fig-4-Example-Simulation-Runs.pdf', width = 12, height = 4)
par(mfrow = c(1, 3))
cm <- 1.75
cl <- 1.75
ca <- 1.50
lc <- 1.35
plot_lotka_rs(
  lotka6s, xat = xlabels1 * 10, cex.main = cm, cex.lab = cl, cex.axis = ca, legend.cex = lc,
  xlabels = xlabels1, main = expression('Sampling 10 x Day and ' ~ sigma[epsilon] ~ ' = 10')
)
plot_lotka_rs(
  lotka8s, xat = xlabels1 * 5, cex.main = cm, cex.lab = cl, cex.axis = ca, legend.cex = lc,
  xlabels = xlabels2, main = expression('Sampling 5 x Day and ' ~ sigma[epsilon] ~ ' = 8')
)
plot_lotka_rs(
  lotka10s, xat = xlabels1 * 1, cex.main = cm, cex.lab = cl, cex.axis = ca, legend.cex = lc,
  xlabels = xlabels3, main = expression('Sampling 1 x Day and ' ~ sigma[epsilon] ~ ' = 6')
)
dev.off()



##################################################
##### Figure 5: ROC Curves Combined Indicator ####
##################################################
# Too big for Github, please send an email if you
# do not want to run all the simulations yourself
datci <- read.csv('Simulation-Results/results-roc-combined-indicator.csv')

# I made them nicer in Adobe Illustrator afterwards
pdf('Figures/Fig-5-ROCs-Combined-Indicator.pdf', width = 10, height = 7)
par(mfrow = c(2, 3))
for (days in c(50, 25)) {
  for (nn in c(4, 6, 8)) {
    
    if (nn == 4) {
      ylab <- 'True Positive Rate'
    } else {
      ylab = ''
    }
    
    if (days == 25) {
      xlab <- 'False Positive Rate'
    } else {
      xlab <- ''
    }
    
    main <- bquote(sigma[epsilon] ~ ' = ' ~ .(nn))
    create_rocplot(
      datci %>% 
        filter(
          noise == nn, nr_transition_days %in% c(NA, days),
          window_size_days == 50, nr_baselinedays == 100
        ), ylab = ylab, xlab = xlab,
      cex.main = 2.5, cex.lab = 1.50, cex.axis = ca, cex.legend = 1.25, main = main
    )
  }
}
dev.off()



###################################
####### Figure 6: AUC Results #####
###################################
datauc <- read_csv('Simulation-Results/results-auc.csv')

dauc <- datauc %>%
  group_by(noise, freq, nr_transition_days, window_size_days, nr_baselinedays, ews) %>% 
  summarize(
    AUC = MESS::auc(c(1, fpr, 0), c(1, tpr, 0), absolutearea = FALSE) # Linearly interpolate
  ) %>% 
  mutate(
    ews = as.character(ews),
    ews = ifelse(ews == 'Combined-Indicator', 'Combined', ews)
  )

# Change for Figure in main text or in appendix
ews_choice <- c('Autocorrelation', 'Variance', 'Cross-correlation', 'CovEigen', 'Combined')
# ews_choice <- c('Mean', 'Skewness', 'Kurtosis', 'Spatial-Skewness', 'Spatial-Kurtosis', 'Spatial-Variance') # Appendix

daucsumtr <- dauc %>% 
  filter(ews %in% ews_choice) %>% 
  group_by(
    ews, noise, freq
  ) %>% 
  summarize(
    meanAUC = mean(AUC),
    sdAUC = sd(AUC)
  ) %>% 
  data.frame() %>% 
  mutate(
    freq = factor(freq, labels = c('Sampling 10x Day', 'Sampling 5x Day', 'Sampling 1x Day'))
  )

daucsumfreq <- dauc %>% 
  filter(ews %in% ews_choice) %>% 
  group_by(
    ews, noise, nr_transition_days
  ) %>% 
  summarize(
    meanAUC = mean(AUC),
    sdAUC = sd(AUC)
  ) %>% 
  data.frame() %>% 
  mutate(
    nr_transition_days = factor(nr_transition_days, labels = 
                 c('Transition Taking 10 Days',
                   'Transition Taking 25 Days',
                   'Transition Taking 50 Days')
    )
  )

daucsumbl <- dauc %>% 
  filter(ews %in% ews_choice) %>% 
  group_by(
    ews, noise, nr_baselinedays
  ) %>% 
  summarize(
    meanAUC = mean(AUC),
    sdAUC = sd(AUC)
  ) %>% 
  data.frame() %>% 
  mutate(
    nr_baselinedays = factor(nr_baselinedays, labels =
                 c('Baseline 25 Days',
                   'Baseline 50 Days',
                   'Baseline 100 Days')
    )
  )

daucsumtr$ews <- reorder(daucsumtr$ews, rep(c(1, 5, 4, 3, 2), each = 12))
daucsumbl$ews <- reorder(daucsumbl$ews, rep(c(1, 5, 4, 3, 2), each = 12))
daucsumbl$nr_baselinedays <- reorder(daucsumbl$nr_baselinedays, rep(c(3, 2, 1), each = 20))
daucsumfreq$ews <- reorder(daucsumfreq$ews, rep(c(1, 5, 4, 3, 2), each = 12))
daucsumfreq$nr_transition_days <- reorder(daucsumfreq$nr_transition_days, rep(c(3, 2, 1), each = 20))

# Appendix
# daucsumtr$ews <- reorder(daucsumtr$ews, rep(c(2, 1, 3, 4, 5, 6), each = 12))
# daucsumfreq$ews <- reorder(daucsumfreq$ews, rep(c(2, 1, 3, 4, 5, 6), each = 12))
# daucsumfreq$nr_transition_days <- factor(
#   daucsumfreq$nr_transition_days,
#   levels = c('Transition Taking 50 Days', 'Transition Taking 25 Days', 'Transition Taking 10 Days')
# )
# daucsumbl$ews <- reorder(daucsumbl$ews, rep(c(2, 1, 3, 4, 5, 6), each = 12))
# daucsumbl$nr_baselinedays <- factor(
#   daucsumbl$nr_baselinedays,
#   levels = c('Baseline 100 Days', 'Baseline 50 Days', 'Baseline 25 Days')
# )
# cols <- rev(brewer.pal(6, 'RdYlBu'))

cols <- c(brewer.pal(4, 'RdYlBu'), 'gray56')
pfreq <- ggplot(
  daucsumtr,
  aes(x = factor(noise), y = meanAUC, color = ews, group = ews)
  ) +
  scale_colour_manual(values = cols) +
  xlab(expression(sigma[epsilon])) +
  ylab('Mean AUC') +
  geom_point(position = position_dodge(0.75), alpha = 1) +
  geom_errorbar(
    aes(ymin = meanAUC - sdAUC, ymax = meanAUC + sdAUC),
    position = position_dodge(0.75), alpha = 1
  ) +
  facet_rep_wrap(~ freq, repeat.tick.labels = TRUE) +
  geom_hline(yintercept = 0.50, color = 'grey76', linetype = 'dotted') +
  scale_x_discrete(breaks = c(2, 4, 6, 8, 10)) +
  scale_y_continuous(breaks = seq(0.30, 1, 0.1), limits = c(0.30, 1)) +
  geom_segment(aes(x = 1, xend = 4, y = -Inf, yend = -Inf), col = 'black')+
  geom_segment(aes(y = 0.30, yend = 1, x = -Inf, xend = -Inf), col = 'black') +
  guides(color = FALSE) +
  ptheme

ptr <- ggplot(
  daucsumfreq,
  aes(x = factor(noise), y = meanAUC, color = ews, group = ews)
) +
  scale_colour_manual(values = cols) +
  xlab(expression(sigma[epsilon])) +
  ylab('Mean AUC') +
  geom_point(position = position_dodge(0.75), alpha = 1) +
  geom_errorbar(
    aes(ymin = meanAUC - sdAUC, ymax = meanAUC + sdAUC),
    position = position_dodge(0.75), alpha = 1
  ) +
  facet_rep_wrap(~ nr_transition_days, repeat.tick.labels = TRUE) +
  geom_hline(yintercept = 0.50, color = 'grey76', linetype = 'dotted') +
  scale_x_discrete(breaks = c(2, 4, 6, 8, 10)) +
  scale_y_continuous(breaks = seq(0.30, 1, 0.1), limits = c(0.30, 1)) +
  geom_segment(aes(x = 1, xend = 4, y = -Inf, yend = -Inf), col = 'black')+
  geom_segment(aes(y = 0.30, yend = 1, x = -Inf, xend = -Inf), col = 'black') +
  guides(color = FALSE) +
  ptheme

pbl <- ggplot(
  daucsumbl,
  aes(x = factor(noise), y = meanAUC, color = ews, group = ews)
) +
  scale_colour_manual(values = cols) +
  ylab('Mean AUC') +
  geom_point(position = position_dodge(0.75), alpha = 1) +
  geom_errorbar(
    aes(ymin = meanAUC - sdAUC, ymax = meanAUC + sdAUC),
    position = position_dodge(0.75), alpha = 1
  ) +
  facet_rep_wrap(~ nr_baselinedays, repeat.tick.labels = TRUE) +
  xlab(expression(sigma[epsilon])) +
  geom_hline(yintercept = 0.50, color = 'grey76', linetype = 'dotted') +
  scale_x_discrete(breaks = c(2, 4, 6, 8, 10)) +
  scale_y_continuous(breaks = seq(0.30, 1, 0.1), limits = c(0.30, 1)) +
  geom_segment(aes(x = 1, xend = 4, y = -Inf, yend = -Inf), col = 'black')+
  geom_segment(aes(y = 0.30, yend = 1, x = -Inf, xend = -Inf), col = 'black') +
  guides(color = FALSE) +
  ptheme

# pdf('Figures/Fig-10-AUC-Appendix.pdf', width = 9, height = 9)
pdf('Figures/Fig-6-AUC-Main.pdf', width = 9, height = 8)
grid.arrange(pfreq, ptr, pbl, nrow = 3)
dev.off()


#######################################
##### Figure 7: Days To Transition ####
#######################################

# Too big for Github, please send an email if you
# do not want to run all the simulations yourself
d <- read.csv('Simulation-Results/results-advance.csv') %>%
  mutate(
    multiplier = ifelse(freq == 90, 1, ifelse(freq == 9, 10, 5)),
    
    freq = factor(freq, labels = c('Sampling 10x Day', 'Sampling 5x Day', 'Sampling 1x Day')),
    
    # Compute theoretical transition
    transition_idx_theoretical = (nr_baselinedays + nr_transition_days) * multiplier,
    
    # Transform first_signal into days (depends on the sampling frequency)
    first_signal = first_signal / multiplier,
    # NB: first_signal gives samples before *actual* transition, not *theoretical* transition
    
    transition_day = (transition_idx - nr_baselinedays * multiplier) / multiplier
  )


dsel <- d %>% 
  filter(
    ews %in% c('CovEigen', 'Cross-correlation', 'Variance', 'Autocorrelation', 'Combined-Indicator')
  ) %>% 
  mutate(
    ews = ifelse(ews == 'Combined-Indicator', 'Combined', ews),
    ews = factor(ews, levels = c('Combined', 'CovEigen', 'Cross-correlation', 'Variance', 'Autocorrelation'))
  )


# Raincloud plots
source('https://raw.githubusercontent.com/RainCloudPlots/RainCloudPlots/master/tutorial_R/R_rainclouds.R')

plot_time2tip <- function(d, y, main) {
  cols <- c('gray', rev(brewer.pal(4, 'RdYlBu')))
  
  ggplot(d, aes(x = ews, y = -(transition_day - first_signal), fill = ews, color = ews)) +
    geom_point(position = position_jitter(width = .15), size = 0.3, alpha = 0.5) +
    xlab('') +
    ylab('Days to Transition') +
    geom_boxplot(width = .2, guides = FALSE, outlier.shape = NA, alpha = 0.5, color = 'black') +
    geom_flat_violin(position = position_nudge(x = 0.2, y = 0), adjust = 1, alpha = 0.5) +
    facet_wrap(~ freq) +
    coord_flip() +
    scale_y_continuous() +
    geom_segment(aes(x = 1, xend = 5, y = -Inf, yend = -Inf), col = 'black') +
    geom_segment(aes(y = y, yend = 0, x = -Inf, xend = -Inf), col = 'black') +
    scale_colour_manual(values = cols) +
    scale_fill_manual(values = cols) +
    ggtitle(main) +
    ptheme +
    theme(
      legend.position = 'none'
    )
}


pdf('Figures/Fig-7-Days-to-Transition.pdf', width = 9, height = 5)
plot_time2tip(
  dsel %>%
    filter(
      nr_transition_days == 50, noise == 6, sigma == 2,
      window_size_days == 50, nr_baselinedays == 100
    ), -60, '')
dev.off()

# True positive rate decreases with sampling frequency
# hence decrease in data from left to right panels
dsel %>%
  filter(
    nr_transition_days == 50, noise == 6,
    window_size_days == 50, nr_baselinedays == 100
  ) %>% 
  group_by(ews, freq) %>% 
  summarize(n = n(), isnotna = sum(!is.na(pred_days)), prop = isnotna / n)


##############################################
##### Table 3: Average Days to Transition ####
##############################################
tab <- d %>% 
  filter(
    ews == 'Combined-Indicator', sigma == 2
  ) %>% 
  group_by(
    noise, freq, nr_transition_days
  ) %>% 
  summarize(
    m = mean(transition_day - first_signal, na.rm = TRUE),
    sd = sd(transition_day - first_signal, na.rm = TRUE)
  )

# Paste into LaTeX
for (i in seq(nrow(tab))) {
  u <- tab[i, ]
  
  if (i %in% c(10, 19, 28)) {
    cat('\n')
  }
  cat(paste0(round(u$m, 2), ' & '))
  # cat(paste0('(', round(u$sd, 2), ') & '))
}


###########################################
#### Figure 9: Boerlijst et al. (2013) ####
###########################################
db <- read.csv('Simulation-Results/results-boerlijst.csv')

pdf('Figures/Fig-9-Boerlijst.pdf', width = 9, height = 3)
par(mfrow = c(1, 3))
cols <- brewer.pal(3, 'Set1')
muP <- db[, ncol(db)]
plot(
  muP, db[, 1], pch = 20, axes = FALSE,
  col = cols[2], ylim = c(.2, 1), ylab = 'Autocorrelation',
  xlab = expression(mu[P]), main = 'Autocorrelation', type = 'l', lwd = 1.5,
  font.main = 1, cex.main = 1.25, cex.lab = 1, cex.axis = 1.5
)
lines(muP, db[, 2], pch = 20, col = cols[3], lwd = 1.5)
lines(muP, db[, 3], pch = 20, col = cols[1], lwd = 1.5)
axis(1, at = c(seq(.45, .53, .02), .552))
axis(2, las = 2, at = seq(.2, 1, .1))
legend(
  'topleft', lty = c(1, 1, 1),
  legend = c('Juveniles', 'Adults', 'Predator'),
  cex = .9, lwd = 1.5, col = c(cols[2], cols[3], cols[1]), box.lty = 0, bty = 'n'
)

plot(
  muP, db[, 1 + 3], pch = 20, axes = FALSE,
  col = cols[2], ylim = c(0, .01), ylab = '',
  xlab = expression(mu[P]), main = 'Standard Deviation', type = 'l', lwd = 1.5,
  font.main = 1, cex.main = 1.25, cex.lab = 1, cex.axis = 1.5
)
mtext('Standard Deviation', side = 2, line = 3.5, cex = 0.63)
lines(muP, db[, 2 + 3], pch = 20, col = cols[3], lwd = 1.5)
lines(muP, db[, 3 + 3], pch = 20, col = cols[1], lwd = 1.5)
axis(1, at = c(seq(.45, .53, .02), .552))
axis(2, las = 2)
legend(
  'topleft', lty = c(1, 1, 1),
  legend = c('Juveniles', 'Adults', 'Predator'),
  cex = .9, lwd = 1.5, col = c(cols[2], cols[3], cols[1]), box.lty = 0, bty = 'n'
)

cols2 <- brewer.pal(3, 'Set2')
zstand <- function(x) (x - mean(x)) / sd(x)
plot(
  muP, zstand(db$Cross_correlation), pch = 20, axes = FALSE,
  col = cols2[1], ylim = c(-4, 4), ylab = 'z-value',
  xlab = expression(mu[P]), main = 'Multivariate Indicators', type = 'l', lwd = 1.5,
  font.main = 1, cex.main = 1.25, cex.lab = 1, cex.axis = 1.5
)
lines(muP, zstand(db$Spatial_Variance), pch = 20, col = cols2[2], lwd = 1.5)
lines(muP, zstand(db$CovEigen), pch = 20, col = cols2[3], lwd = 1.5)
axis(1, at = c(seq(.45, .53, .02), .552))
axis(2, las = 2)

legend(
  'topleft', lty = c(1, 1, 1),
  legend = c('CovEigen', 'Spatial Variance', 'Cross-correlation'),
  cex = .9, lwd = 1.5, col = cols2[c(3, 2, 1)], box.lty = 0, bty = 'n'
)
dev.off()
