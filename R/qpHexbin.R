qpHexbin <- function (mdata, main = "", ...)
{
  y <- maM(mdata)
  x <- maA(mdata)
  nalim <- !is.na(x) & !is.na(y)
  bin <- hexbin(x[nalim], y[nalim])
  col <- BTY
  maxcnt <- max(bin@count)
##  print(maxcnt)
  colorcut <- seq(0, 1, length = min(17, maxcnt))
  yrange <- c(min(maM(mdata), na.rm = TRUE), max(maM(mdata),
                                na.rm = TRUE) + 1)
  ## move to next traditional graphics plot region
  frame()
  ## Set up grid viewports that correspond to the
  ## current traditional figure region
  vps <- baseViewports()
  pushViewport(vps$inner, vps$figure)
  pushViewport(viewport(gp=gpar(cex=0.6)))
  ## draw a complete hexbin plot with legend
  plot(bin, xlab = "A", ylab = "M",
       ## ylim = yrange,
       main = main, colramp = col, colorcut = colorcut, maxcnt = maxcnt,
       legend = 1,
       lcex = 1,
       ## ysize = 8,
       ## inner = 1)
       newpage=FALSE)
  upViewport(3)
  ## move to next traditional graphics plot region
  ## i.e., skip over region that used to be used for hexbin legend
  frame()
}

qpHexbin <- function (mdata, main = "", ...)
{
  par(mar = c(5, 4, 3, 2))
  y <- maM(mdata)
  x <- maA(mdata)
  nalim <- !is.na(x) & !is.na(y)
  bin <- hexbin(x[nalim], y[nalim])
  col <- BTY
  maxcnt <- max(bin@count)
  ##print(maxcnt)
  colorcut <- seq(0, 1, length = min(17, maxcnt))
  yrange <- c(min(maM(mdata), na.rm = TRUE), max(maM(mdata),
                                na.rm = TRUE) + 1)
  print("1")
  plot(bin@xbnds, bin@ybnds, type = "n", xlab = "A", ylab = "M",
       ylim = yrange, main = main)
  ## set up grid viewports that correspond to the current
  ## traditional plot region
  vps <- baseViewports()
  print("2")
  pushViewport(vps$inner, vps$figure, vps$plot)
  print("3")
  grid.hexagons(bin, colramp = col, colorcut = colorcut, maxcnt = maxcnt)
  print("4")
  upViewport(3)
  print("5")
  par(mar = c(2, 1, 2, 1))
  plot(seq(0, 4, length = 11), -1:9, type = "n", axes = FALSE,
       xlab = "", ylab = "")
  print("6")
  ## set up grid viewports that correspond to the current
  ## traditional plot region
  vps <- baseViewports()
  print("7")
  pushViewport(vps$inner, vps$figure, vps$plot)
  print("9")
  grid.hexlegend(legend = 4, ysize = 8, lcex = 0.8, inner = 1, maxcnt = maxcnt,
                 colorcut = colorcut, colramp = col)
  print("10")
  upViewport(3)
}
