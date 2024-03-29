#' @importFrom ggridges geom_density_ridges
NULL

#' A function for summarize HTS
#' @description A summarizer for HTS
#' @param x a HTS
#' @return A matrix with basic statistics.
#' @examples
#' summaryHTS(subsetHTS(RetHTS, from = 1, to = 10))
#' @export
summaryHTS <- function(x) {
  df <- data.frame(
    Tstamp = numeric(),
    T_period_start = numeric(),
    T_period_end = numeric(),
    mean = numeric(),
    std = numeric(),
    skew = numeric(),
    kurt = numeric(),
    min = numeric(),
    Q1 = numeric(),
    MED = numeric(),
    Q3 = numeric(),
    Max = numeric()
  )
  nr <- length(x@data)
  for (i in 1:nr) {
    tmpo <- x@data[[i]]
    newrow <- data.frame(
      Tstamp = tmpo@tstamp,
      T_period_start = tmpo@period$start,
      T_period_end = tmpo@period$end,
      mean = tmpo@m,
      std = tmpo@s,
      skew = skewH(tmpo),
      kurt = kurtH(tmpo),
      min = tmpo@x[1],
      Q1 = compQ(tmpo, 0.25),
      MED = compQ(tmpo, 0.5),
      Q3 = compQ(tmpo, 0.75),
      Max = tmpo@x[length(tmpo@x)]
    )
    df <- rbind(df, newrow)
  }
  return(df)
}




# OK assigned

plot.gg <- function(x, type = "HISTO", col = "green", border = "black") {
  xlabel <- paste("m= ", format(x@m, nsmall = 4), " std= ", format(x@s, nsmall = 4))

  if (type == "HISTO") {
    lowers <- x@x[1:(length(x@x) - 1)]
    uppers <- x@x[2:length(x@x)]
    maxdens <- 10 / (max(uppers) - lowers[1])
    ampl <- uppers - lowers
    # ampl[which(ampl==0)]=1

    dens <- (x@p[2:length(x@p)] - x@p[1:(length(x@p) - 1)]) / ampl
    dens[which(dens > maxdens)] <- maxdens
    df <- data.frame(xm = lowers, xM = uppers, ym = dens * 0, yM = dens)

    p <- with(df, ggplot(df, aes(xmin = xm, xmax = xM, ymin = ym, ymax = yM)) +
      geom_rect(alpha = 0.7, fill = col, colour = border) +
      geom_vline(xintercept = x@m, colour = col, linetype = "dashed", size = 0.5) +
      xlab(xlabel) +
      ylab("density") +
      ggtitle("Histogram"))
  }
  if (type == "CDF") {
    xs <- x@x
    ps <- x@p
    df <- data.frame(x = xs, cdf = ps)

    p <- with(df, ggplot(df, aes(x = x, y = cdf)) +
      geom_line(colour = border) +
      xlab("x") +
      ylab("p") +
      ggtitle("Cumulative Distribution Function"))
  }
  if (type == "QF") {
    xs <- x@x
    ps <- x@p
    df <- data.frame(x = xs, cdf = ps)

    p <- with(df, ggplot(df, aes(x = cdf, y = x)) +
      geom_line(colour = border) +
      xlab("p") +
      ylab("x") +
      ggtitle("Quantile Function"))
  }
  if (type == "DENS") {
    # generate 200 random points according to the QF
    rn <- 200
    xn <- c(rep(0, rn + 1))
    random_no <- c(0:rn) / rn
    for (i in 0:rn) {
      xn[i + 1] <- compQ(x, random_no[i + 1])
    }

    df <- data.frame(x = xn)
    p <- with(df, ggplot(df, aes(x = x)) +
      geom_density(alpha = 0.7, fill = col, colour = border) +
      xlab(xlabel) +
      ylab("density") +
      ggtitle("Density plot (KDE)"))
    p <- p + geom_vline(xintercept = x@m, colour = "black", linetype = "dashed", size = 0.5)
  }
  if (type == "HBOXPLOT") {
    qua <- c(0, 0.25, 0.5, 0.75, 1)
    df <- data.frame(ymin = 0, lower = 0, middle = 0, upper = 0, ymax = 0)
    for (i in 1:5) {
      df[[i]] <- compQ(x, qua[i])
    }
    df$x <- factor(0)

    p <- with(df, ggplot(df, aes_all(names(df))) +
      geom_boxplot(stat = "identity", fill = col, size = 1, colour = border) +
      ylab(xlabel) +
      xlab("") +
      ggtitle("Horizontal Boxplot") +
      coord_flip())
    # xn=c(0,0,0,0,0)
    # for (i in 1:5){
    #   xn[i]=compQ(x,qua[i])
    # }
    # df=data.frame(x=xn)
    # p=with(df,ggplot(df, aes(x=factor(0), ymin = x[1], lower = x[2], middle =x[3], upper =x[4], ymax =x[5]))+
    #          geom_boxplot(stat = "identity", fill=col, size=1, colour=border)+
    #          ylab(xlabel) +xlab("")+
    #          ggtitle("Horizontal Boxplot")+coord_flip()
    # )
  }
  if (type == "VBOXPLOT") {
    qua <- c(0, 0.25, 0.5, 0.75, 1)
    df <- data.frame(ymin = 0, lower = 0, middle = 0, upper = 0, ymax = 0)
    for (i in 1:5) {
      df[[i]] <- compQ(x, qua[i])
    }
    df$x <- factor(0)

    p <- with(df, ggplot(df, aes_all(names(df))) +
      geom_boxplot(stat = "identity", fill = col, size = 1, colour = border) +
      ylab("x") +
      xlab(xlabel) +
      ggtitle("Vertical Boxplot"))

    # p=with(df,ggplot(df, aes(x=factor(0), ymin = x[1], lower = x[2], middle =x[3], upper =x[4], ymax =x[5]))+
    #          geom_boxplot(stat = "identity", fill=col, size=1, colour=border)+
    #          ylab("x") +xlab(xlabel)+
    #          ggtitle("Vertical Boxplot")
    # )
  }

  return(p)
}

WH.joy <- function(DATA, var) {
  n <- get.MatH.nrows(DATA)
  n_names <- get.MatH.rownames(DATA)
  vp <- c(0:1000) / 1000
  vec1 <- numeric()
  labs <- character()

  for (i in 1:n) {
    tmp <- COMP_Q_VECT(DATA@M[i, var][[1]]@x, DATA@M[i, var][[1]]@p, vp)
    vec1 <- c(vec1, tmp)
    labs <- c(labs, rep(n_names[i], length(tmp)))
  }
  vname <- rep(get.MatH.varnames(DATA)[var], length(labs))

  DF <- data.frame(x = vec1, y = as.factor(labs), z = vname)
  return(DF)
}
WH.joymult <- function(DATA, list_of_vars) { ## AS SUBSTUTE OF PLOT MATH DENS

  DF <- WH.joy(DATA, list_of_vars[1])
  if (length(list_of_vars) > 1) {
    for (i in 2:length(list_of_vars)) {
      tmp <- WH.joy(DATA, list_of_vars[i])
      DF <- rbind(DF, tmp)
    }
  }
  DF$y <- factor(DF$y, levels = rev(get.MatH.rownames(DATA)))
  #  levels(DF$y)=rev(levels(DF$y))
  return(DF)
}

WH.joy.HIST <- function(DATA, var, tr = 1) {
  n <- get.MatH.nrows(DATA)
  n_names <- get.MatH.rownames(DATA)
  vec1 <- numeric()
  vec2 <- numeric()
  vec3 <- numeric()
  labs <- character()
  for (i in 1:n) {
    dom <- DATA@M[i, var][[1]]@x
    widths <- diff(dom)
    widths[which(widths == 0)] <- 1e-16
    p <- diff(DATA@M[i, var][[1]]@p)
    densit <- p / widths
    tmpx <- numeric() # sort(c(dom,dom))
    tmpy <- numeric()
    for (j in 1:length(densit)) {
      tmpx <- c(tmpx, dom[j], dom[j], dom[j + 1], dom[j + 1])
      tmpy <- c(tmpy, 0, densit[j], densit[j], 0)
    }
    #    tmpy=c(tmpy,0)

    vec1 <- c(vec1, tmpx)
    vec2 <- c(vec2, tmpy)
    vec3 <- c(vec3, rep(i, length(tmpx)))
    labs <- c(labs, rep(n_names[i], length(tmpx)))
  }
  vname <- rep(get.MatH.varnames(DATA)[var], length(labs))
  if (tr == 1) {
    treshold <- 2 # quantile(vec2, probs = 0.8)
  } else {
    treshold <- 1.5 * diff(quantile(vec2, probs = c(0.25, 0.75))) + quantile(vec2, probs = c(0.75))
  }
  vec2[vec2 > treshold] <- treshold
  DF <- data.frame(X = vec1, Y = vec3, H = vec2, y = labs, z = as.factor(vname))

  return(DF)
}
WH.joymult.HISTO <- function(DATA, list_of_vars, tr = 1) {
  DF <- WH.joy.HIST(DATA, list_of_vars[1], tr = tr)
  if (length(list_of_vars) > 1) {
    for (i in 2:length(list_of_vars)) {
      tmp <- WH.joy.HIST(DATA, list_of_vars[i], tr = tr)
      DF <- rbind(DF, tmp)
    }
  }
  DF$y <- factor(DF$y, levels = rev(get.MatH.rownames(DATA)))
  return(DF)
}


# OK assigned
plot.M <- function(x, type = "HISTO", border = "black", angL = 330) {
  varsno <- ncol(x@M)
  indno <- nrow(x@M)
  if (varsno == 1 && indno == 1) {
    tmpo <- x@M[1, 1][[1]]
    p <- plot.gg(tmpo, type = type, col = "green", border = border)
  }
  else {
    if (type == "HISTO") {
      df <- WH.joymult.HISTO(DATA = x, list_of_vars = c(1:get.MatH.ncols(x)), tr = 1)
      CC <- get.MatH.ncols(x)
      p <- with(
        df,
        ggplot(df, aes(x = X, y = y, height = H, group = y, fill = z)) +
          geom_density_ridges(stat = "identity", scale = 1.1, alpha = 0.7, color = border) +
          facet_wrap(~z, scales = "free_x", ncol = CC) +
          theme(
            strip.text.x = element_text(face = "bold"),
            axis.text.x = element_text(size = 6),
            axis.text.y = element_text(face = "bold", vjust = 0),
            legend.position = "none",
            axis.title.x = element_blank(), axis.title.y = element_blank()
          )
      )
    }
    if (type == "DENS") {
      df <- WH.joymult(DATA = x, list_of_vars = c(1:get.MatH.ncols(x)))
      CC <- get.MatH.ncols(x)
      p <- with(
        df,
        ggplot(df, aes(x = x, y = y, fill = z)) +
          geom_density_ridges(scale = 1.3, rel_min_height = 0.01, alpha = 0.7) +
          facet_wrap(~z, scales = "free_x", ncol = CC) +
          theme(
            strip.text.x = element_text(face = "bold"),
            axis.text.x = element_text(size = 6),
            axis.text.y = element_text(face = "bold", vjust = 0),
            legend.position = "none",
            axis.title.x = element_blank(), axis.title.y = element_blank()
          )
      )
    }
    if (type == "BOXPLOT") {
      df <- data.frame(
        ymin = numeric(),
        yQ1 = numeric(),
        yME = numeric(),
        yQ3 = numeric(),
        ymax = numeric(),
        indiv = character(), vars = character(), stringsAsFactors = FALSE
      )
      for (i in 1:indno) {
        for (j in 1:varsno) {
          qua <- c(0, 0.25, 0.5, 0.75, 1)
          xn <- c(0, 0, 0, 0, 0)
          tmpo <- x@M[i, j][[1]]
          for (k in 1:5) {
            xn[k] <- compQ(tmpo, qua[k])
          }
          newrow <- data.frame(
            ymin = xn[1],
            yQ1 = xn[2],
            yME = xn[3],
            yQ3 = xn[4],
            ymax = xn[5],
            indiv = rownames(x@M)[i],
            vars = colnames(x@M)[j],
            stringsAsFactors = FALSE
          )
          df <- rbind(df, newrow)
        }
      }
      df$indiv <- as.factor(df$indiv)
      df$indiv <- factor(df$indiv, levels = rownames(x@M))
      df$vars <- as.factor(df$vars)
      df$vars <- factor(df$vars, levels = colnames(x@M))
      p <- with(df, ggplot(df, aes(x = factor(0), ymin = ymin, lower = yQ1, middle = yME, upper = yQ3, ymax = ymax, fill = vars)) +
        geom_boxplot(stat = "identity", size = 0.5, colour = border) +
        facet_grid(vars ~ indiv, scales = "free_y") +
        scale_x_discrete(breaks = NULL) +
        theme(
          legend.position = "none",
          axis.title.x = element_blank(),
          strip.text.x = element_text(size = 12, face = "bold", angle = 30),
          strip.text.y = element_text(size = 12, face = "bold")
        ))
      # print(p)
    }
  }
  return(p)
}

# OK assigned
plot.HTS.1v <- function(x, type = "BOXPLOT", border = "black", maxno.perplot = 15) {
  numofp <- length(x@data)
  numoflines <- ceiling(numofp / maxno.perplot)
  listofP <- list()
  df1 <- summaryHTS(x)
  clint <- c(1:numofp) / numofp
  TT <- max(df1$Tstamp) - min(df1$Tstamp)
  Tini <- min(df1$Tstamp)
  for (plo in 1:numoflines) {
    selected <- c(((plo - 1) * maxno.perplot + 1), min(plo * maxno.perplot, numofp))
    y <- subsetHTS(x, selected[1], selected[2])
    if (type == "BOXPLOT") {
      df <- df1[(selected[1]:selected[2]), ]
      p <- with(df, ggplot(df, aes(
        x = as.factor(Tstamp), ymin = min, lower = Q1, middle = MED,
        upper = Q3, ymax = Max, fill = as.factor(Tstamp)
      )) +
        geom_boxplot(stat = "identity") +
        guides(fill = FALSE) +
        theme(
          legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(vjust = 0.5, size = 8)
        ))
      #   print(p)
      listofP[[plo]] <- p
    }
    if (type == "VIOLIN") {
      nr <- length(y@data)
      df <- data.frame(x = numeric(), t = numeric(), grad = numeric())
      df2 <- data.frame(x = numeric(), t = numeric(), xe = numeric(), te = numeric(), grad = numeric())
      for (i in 1:nr) {
        points <- 100
        xn <- c(rep(0, points))
        random_no <- c(0:points) / points
        for (p in 1:points) {
          xn[p] <- compQ(x@data[[i]], random_no[p])
        }

        newrow <- data.frame(x = xn, t = rep(y@data[[i]]@tstamp), grad = rep((y@data[[i]]@tstamp - Tini) / TT))
        df <- rbind(df, newrow)
        if (i == nr) {
          newrow2 <- data.frame(
            x = y@data[[i]]@m, t = rep(y@data[[i]]@tstamp),
            xe = y@data[[i]]@m, te = rep(y@data[[i]]@tstamp),
            grad = rep((y@data[[i]]@tstamp - Tini) / TT)
          )
        }
        else {
          newrow2 <- data.frame(
            x = y@data[[i]]@m, t = rep(y@data[[i]]@tstamp),
            xe = y@data[[i + 1]]@m, te = rep(y@data[[i + 1]]@tstamp),
            grad = rep((y@data[[i]]@tstamp - Tini) / TT)
          )
        }
        df2 <- rbind(df2, newrow2)
      }
      p <- with(df, ggplot(df, aes(factor(t), x, fill = grad)) +
        geom_violin() +
        geom_point(data = df2, aes(factor(t), x)) +
        geom_line(data = df2, aes(x = factor(t), y = x, group = factor(0)), size = 0.5, alpha = 0.7, linetype = "dashed") +
        scale_fill_gradient2(limits = c(0, 1), low = "red", mid = "white", high = "green", midpoint = 0.5) +
        theme(
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
          legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(angle = 330, vjust = 0.5, size = 10)
        ))

      listofP[[plo]] <- p
    }
  }
  # return(listofP)
  multiplot(listofP)
}


multiplot <- function(plotlist = NULL, ..., file, cols = 1, layout = NULL) {
  # require(grid)
  # require(ggplot2)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots <- length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots / cols)),
      ncol = cols, nrow = ceiling(numPlots / cols)
    )
  }

  if (numPlots == 1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(
        layout.pos.row = matchidx$row,
        layout.pos.col = matchidx$col
      ))
    }
  }
}

#' A function for comparing observed vs predicted histograms
#' @description This function allows the representation of observed vs predicted histograms. It can
#' be used as a tool for interpreting preditive methods (for exampe, the regression of histogrma data)
#' @param PRED a \code{MatH} object with one column, the predicted data
#' @param OBS a \code{MatH} object with one column, the observed data
#' @param type a string. "HISTO" (default), if ones want to compare histograms\cr
#' "CDF", if ones want to compare cumulative distribution functions;\cr
#' "DENS" if ones want to compare approximated densities (using KDE);\cr
#' @param ncolu number of columns in which is arranged the plot, default is 2. If you have a lot of data consider
#' to choose higher values.
#' @return A plot with compared histogram-valued data.
#' @examples
#' ## do a regression
#' pars <- WH.regression.two.components(BLOOD, Yvar = 1, Xvars = c(2:3))
#' ## predict data
#' PRED <- WH.regression.two.components.predict(data = BLOOD[, 2:3], parameters = pars)
#' ## define observed data
#' \dontrun{
#' OBS <- BLOOD[, 1]
#' plotPredVsObs(PRED, OBS, "HISTO")
#' plotPredVsObs(PRED, OBS, "CDF")
#' plotPredVsObs(PRED, OBS, "DENS")
#' }
#' @export
plotPredVsObs <- function(PRED, OBS, type = "HISTO", ncolu = 2) {
  # require ("vioplot")
  nobj <- get.MatH.nrows(PRED)
  maxpercol <- ceiling(nobj / ncolu)
  ListofP <- list()
  M <- data.frame()
  cc <- 0
  for (ind in 1:nobj) {
    cc <- cc + 1
    if (cc > ncolu) cc <- 1
    if (type == "HISTO") {
      x <- OBS@M[ind, 1][[1]]
      lowers <- x@x[1:(length(x@x) - 1)]
      uppers <- x@x[2:length(x@x)]
      maxdens <- 10 / (max(uppers) - lowers[1])
      ampl <- uppers - lowers
      # ampl[which(ampl==0)]=1
      dens <- (x@p[2:length(x@p)] - x@p[1:(length(x@p) - 1)]) / ampl
      dens[which(dens > maxdens)] <- maxdens
      df <- data.frame(
        xm = lowers, xM = uppers, ym = dens * 0, yM = dens,
        Type = rep("OBS", length(lowers)),
        NameO = rep(rownames(OBS@M)[ind], length(lowers)),
        colu = rep(cc, length(lowers))
      )
      x <- PRED@M[ind, 1][[1]]
      lowers <- x@x[1:(length(x@x) - 1)]
      uppers <- x@x[2:length(x@x)]
      maxdens <- 10 / (max(uppers) - lowers[1])
      ampl <- uppers - lowers
      # ampl[which(ampl==0)]=1
      dens <- (x@p[2:length(x@p)] - x@p[1:(length(x@p) - 1)]) / ampl
      dens[which(dens > maxdens)] <- maxdens
      df2 <- data.frame(
        xm = lowers, xM = uppers, ym = dens * 0, yM = dens,
        Type = rep("PRED", length(lowers)),
        NameO = rep(rownames(OBS@M)[ind], length(lowers)),
        colu = rep(cc, length(lowers))
      )
      M <- rbind(M, df, df2)
    }
    if (type == "CDF") {
      x <- PRED@M[ind, 1][[1]]
      xs <- x@x
      ps <- x@p
      df <- data.frame(
        x = xs, cdf = ps,
        Type = rep("PRED", length(xs)),
        NameO = rep(rownames(OBS@M)[ind], length(xs)),
        colu = rep(cc, length(xs))
      )
      x <- OBS@M[ind, 1][[1]]
      xs <- x@x
      ps <- x@p
      df2 <- data.frame(
        x = xs, cdf = ps,
        Type = rep("OBS", length(xs)),
        NameO = rep(rownames(OBS@M)[ind], length(xs)),
        colu = rep(cc, length(xs))
      )
      M <- rbind(M, df, df2)
    }
    if (type == "DENS") {
      # generate 200 random points according to the QF
      rn <- 200
      xn1 <- c(rep(0, rn))
      xn2 <- c(rep(0, rn))
      random_no <- c(0:rn) / rn
      x1 <- PRED@M[ind, 1][[1]]
      x2 <- OBS@M[ind, 1][[1]]
      for (i in 1:rn) {
        xn1[i] <- compQ(x1, random_no[i])
        xn2[i] <- compQ(x2, random_no[i])
      }

      df1 <- data.frame(
        x = xn1,
        Type = rep("PRED", length(xn1)),
        NameO = rep(rownames(OBS@M)[ind], length(xn1)),
        colu = rep(cc, length(xn1))
      )
      df2 <- data.frame(
        x = xn2,
        Type = rep("OBS", length(xn2)),
        NameO = rep(rownames(OBS@M)[ind], length(xn2)),
        colu = rep(cc, length(xn2))
      )
      M <- rbind(M, df1, df2)
    }
  }
  ListofP <- list()
  levels(M$Type) <- c("OBS", "PRED")
  for (CC in 1:ncolu) {
    if (type == "HISTO") {
      p <- with(M, ggplot(subset(M, colu == CC), aes(xmin = xm, xmax = xM, ymin = ym, ymax = yM, fill = Type)) +
        geom_rect(alpha = 0.5) +
        facet_grid(NameO ~ .) +
        theme(legend.position = "bottom", legend.title = element_blank()))
    }
    if (type == "CDF") {
      p <- with(M, ggplot(subset(M, colu == CC), aes(x = x, y = cdf)) +
        geom_line(aes(linetype = Type)) +
        facet_grid(NameO ~ .) +
        theme(
          legend.position = "bottom",
          axis.title.x = element_blank(),
          legend.title = element_blank()
        ))
    }
    if (type == "DENS") {
      p <- with(M, ggplot(subset(M, colu == CC), aes(x = x, fill = Type)) +
        geom_density(alpha = 0.6, colour = "gray") +
        facet_grid(NameO ~ .) +
        theme(
          legend.position = "bottom",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_blank()
        ))
    }
    ListofP <- c(ListofP, list(p))
  }
  multiplot(ListofP, cols = ncolu)
}

#' A function for plotting functions of errors
#' @description This function allows the representation of the difference between observed histograms and
#' the respective predicted ones. It can
#' be used as a tool for interpreting preditive methods (for exampe, the regression of histogrma data)
#' @param PRED a \code{MatH} object with one column, the predicted data
#' @param OBS a \code{MatH} object with one column, the observed data
#' @param type a string. "HISTO_QUA" (default), if ones want to compare histograms quantile differences\cr
#' "HISTO_DEN", if ones want to show the histogram densities differences;\cr
#' "DENS_KDE" if ones want to show the differences between approximated densities (using KDE);\cr
#' @param np number of points considered for density  or quantile computation (default=200).
#' @return A plot with functions of differences between observed and predicted histograms, and a
#' Root Mean Squared value computing by using the L2 Wasserstein distance.
#' @examples
#' ## do a regression
#' pars <- WH.regression.two.components(BLOOD, Yvar = 1, Xvars = c(2:3))
#' ## predict data
#' PRED <- WH.regression.two.components.predict(data = BLOOD[, 2:3], parameters = pars)
#' ## define observed data
#' OBS <- BLOOD[, 1]
#' plot_errors(PRED, OBS, "HISTO_QUA")
#' plot_errors(PRED, OBS, "HISTO_DEN")
#' plot_errors(PRED, OBS, "DENS_KDE")
#' @export
plot_errors <- function(PRED, OBS, type = "HISTO_QUA", np = 200) {
  RMS_W <- 0
  nobj <- get.MatH.nrows(PRED)
  dists <- rep(0, nobj)
  # maxpercol=ceiling(nobj/ncolu)
  ListofP <- list()
  M <- data.frame()
  # M2=data.frame()
  # M3=data.frame()
  minimum <- min(min((get.MatH.stats(OBS, stat = "min"))$mat), min((get.MatH.stats(PRED, stat = "min"))$mat))
  maximum <- max(max((get.MatH.stats(OBS, stat = "max"))$mat), max((get.MatH.stats(PRED, stat = "max"))$mat))
  vals <- minimum + (maximum - minimum) * c(0:np) / np
  range <- maximum - minimum
  NM <- WH.bind.col(OBS, PRED)
  colnames(NM@M) <- c("OBS", "PRED")
  NM <- registerMH(NM)
  for (i in 1:nobj) {
    xo <- NM@M[i, 1][[1]]
    xp <- NM@M[i, 2][[1]]
    dists[i] <- WassSqDistH(xo, xp)
    RMS_W <- RMS_W + 1 / nobj * dists[i]
    if (type == "HISTO_QUA") {
      # tmp=register(xo,xp)
      x <- rep(0, (np + 1))
      p <- rep(0, (np + 1))
      for (j in 0:np) {
        p[j + 1] <- j / np
        x[j + 1] <- compQ(xo, p[j + 1]) - compQ(xp, p[j + 1])
        #       x=xo@x-xp@x
        #       p=xo@p
      }
      df <- data.frame(quantile = x, p = p, NameO = rep(rownames(OBS@M)[i], length(x)))
      M <- rbind(M, df)
    }
    if (type == "HISTO_DEN") {
      fo <- rep(0, np)
      fp <- rep(0, np)
      mp <- 0.5 * (vals[2:np] + vals[1:(np - 1)])
      ran <- (vals[2:np] - vals[1:(np - 1)])
      for (v in 1:np) {
        fo[v] <- compP(xo, vals[v])
        fp[v] <- compP(xp, vals[v])
      }
      do <- (fo[2:np] - fo[1:(np - 1)]) / ran
      dp <- (fp[2:np] - fp[1:(np - 1)]) / ran
      diffden <- do - dp
      df2 <- data.frame(x = mp, density = diffden, NameO = rep(rownames(OBS@M)[i], length(mp)))
      M <- rbind(M, df2)
    }
    if (type == "DENS_KDE") {
      # smooth
      po <- rep(0, np)
      pp <- rep(0, np)
      for (v in 1:np) {
        po[v] <- compQ(xo, (v / np))
        pp[v] <- compQ(xp, (v / np))
      }
      do <- density(po, from = minimum, to = maximum)
      dp <- density(pp, from = minimum, to = maximum)
      diffdkde <- do$y - dp$y
      df3 <- data.frame(x = do$x, density = diffdkde, NameO = rep(rownames(OBS@M)[i], length(do$x)))
      M <- rbind(M, df3)
    }
  }
  if (type == "HISTO_QUA") {
    p <- with(M, ggplot(data = M, aes(x = p, y = quantile, group = NameO, colour = NameO)) +
      geom_line() +
      geom_point() +
      ggtitle("Plot of quantile differences"))
  }
  if (type == "HISTO_DEN") {
    p <- with(M, ggplot(data = M, aes(x = x, y = density, group = NameO, colour = NameO)) +
      geom_line() +
      ggtitle("Plot of density differences"))
  }
  if (type == "DENS_KDE") {
    p <- with(M, ggplot(data = M, aes(x = x, y = density, group = NameO, colour = NameO)) +
      geom_line() +
      ggtitle("Plot of smoothed (KDE) \n density differences"))
  }
  print(p)
  RMS_W <- sqrt(RMS_W)
  print(RMS_W)
  # print(RMS_W/range)
  # print(sqrt(dists))
  # summary(sqrt(dists))
}


