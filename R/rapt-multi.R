## Multiple Hit Mark Extraction
multiples <- function(ato, cl.n = 3) {
  mu.ind <- sort(unique(
    c(which(diff(ato$pIndex) == 0), which(diff(ato$pIndex) == 0) + 1)
  ))
  mu.rle <- rle(ato[mu.ind,]$pIndex)
  mu.rle <- data.frame(lengths = mu.rle$lengths, values = mu.rle$values)
  mu.cl <- parallel::makeCluster(cl.n, type = 'FORK')
  mu.ls <- parallel::parLapply(mu.cl, 1:max(mu.rle$lengths), function(x) {
    ind <- which(mu.rle$lengths == x)
    cum <- cumsum(mu.rle$lengths)
    first <- cum[ind-1] + 1
    mu <- mu.ind[first]
    return(mu)
  })
  parallel::stopCluster(mu.cl)
  mu.ls[[1]] <- setdiff(1:dim(ato)[1],mu.ind)
  return(mu.ls)
}


## 2D Double Hit Mass Correlation ##
saxeyPlot <- function(doubles, ms) {
  sax.old <- par(no.readonly = T);
  par(mar = c(1,1,1,1))
  plot(0:1, 0:1, type='n', xlab = '', ylab = '', axes = FALSE)
  par(new = TRUE, plt = c(0.1, 0.85, 0.1, 0.85))
  image(doubles$x, doubles$y, log(doubles$z),
        breaks = seq(floor(min(log(doubles$z[doubles$z > 0]))),
                     floor(max(log(doubles$z))),
                     length.out = 256), col = rainbow(255),
        xlab = '', ylab = '', xaxt = 'n', yaxt = 'n'
  )
  axis(1, tck = -0.015, labels = FALSE)
  axis(1, tick = FALSE, line = -0.5)
  mtext('Mass-to-Charge Ratio (Da)', side = 1, line = 2)
  axis(2, tck = -0.015, labels = FALSE)
  axis(2, tick = F, line = -0.5, las = 1)
  mtext('Mass-to-Charge Ratio (Da)', side = 2, line = 2)
  # Add peak lines
  #abline(h = ms@mass, col = rgb(1,1,1,0.2))
  #abline(v = ms@mass, col = rgb(1,1,1,0.2))
  # Mass spectrum for second hit
  par(new = TRUE, plt = c(0.1, 0.85, 0.85, 0.95))
  plot(ms@intensity~ms@mass, type = 'l',
       xaxs = 'i', yaxs = 'i', bty = 'n', xaxt = 'n', yaxt = 'n',
       xlab = '', ylab = '')
  # Mass spectrum for first hit
  par(new = TRUE, plt = c(0.85, 0.95, 0.1, 0.85))
  plot(ms@mass~ms@intensity, type = 'l',
       xaxs = 'i', yaxs = 'i', bty = 'n', xaxt = 'n', yaxt = 'n',
       xlab = '', ylab = '')
  par(sax.old);
}
