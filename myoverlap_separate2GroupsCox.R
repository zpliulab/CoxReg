my_overlap <- function(x, y){
  # x <- coef_lasso
  # y <- Data_train
  
  coefs.v <- x[,1] %>% { .[. != 0]}
  coefs.v %>% {
    data.frame(gene.name   = names(.),
               coefficient = .,
               stringsAsFactors = FALSE)
  } %>%
    arrange(gene.name) %>%
    knitr::kable()
  
  sele <- rownames(as.matrix(coefs.v))
  gene <- rownames(y[-c(1,2),])
  overlap <- intersect(sele, gene)
  
  lab <- x[,1] %>% { .[. != 0]} %>% as.matrix
  coefs.v <- lab[overlap,]
  
  my <- list(coefs.v, overlap)
  return(my)
}


separate2GroupsCox <- function(chosen.btas, xdata, ydata,
                               probs = c(.5, .5), no.plot = FALSE,
                               # plot.title = 'SurvivalCurves',
                               xlim = NULL, ylim = NULL, expand.yzero = FALSE,
                               legend.outside = FALSE,
                               ...) {
  
  # convert between compatible formats
  if(inherits(chosen.btas, 'numeric')) {
    chosen.btas <- list(chosen.btas)
  }
  
  # convert between compatible formats
  if (inherits(xdata, 'data.frame') ||
      inherits(xdata, 'numeric') ||
      inherits(xdata, 'matrix')) {
    #
    xdata <- as.matrix(xdata)
  }
  
  # checks if main arguments are correct.
  if (!inherits(chosen.btas, 'list')) {
    stop('chosen.btas argument must be a list or vector. ',
         'See documentation ?separate2GroupsCox')
  } else if (!inherits(xdata, 'matrix')) {
    stop('xdata argument must be a matrix, data.frame or vector ',
         'See documentation ?separate2GroupsCox')
  } else if (!inherits(ydata, 'data.frame')) {
    stop('ydata argument must be a data.frame. ',
         'See documentation ?separate2GroupsCox')
  }
  
  if (nrow(xdata) != nrow(ydata)) {
    stop(sprintf('Rows in xdata (%d) and ydata (%d) must be the same',
                 nrow(xdata), nrow(ydata)))
  } else if (!all(ncol(xdata) == vapply(chosen.btas, length, 1))) {
    stop('All or some of the chosen.btas (%s) have different ',
         sprintf('number of variables from xdata (%d)'),
         paste(vapply(chosen.btas, length, 1), collapse = ', '),
         ncol(xdata))
  }
  #
  # creates a matrix from list of chosen.btas
  chosen.btas.mat <- vapply(chosen.btas, function(e){as.vector(e)},
                            rep(1.0, ncol(xdata)))
  # calculate prognostic indexes for each patient and btas
  prognostic.index <- tryCatch(xdata %*% chosen.btas.mat,
                               error = function(err){
                                 stop('xdata is.matrix(.) = ', is.matrix(xdata), '\n',
                                      ' chosen.btas is.matrix(.) = ', is.matrix(chosen.btas.mat), '\n',
                                      '      xdata (nrow)x(ncol) = ', sprintf('%dx%d',
                                                                              nrow(xdata),
                                                                              ncol(xdata)), '\n',
                                      'chosen.btas (nrow)x(ncol) = ', sprintf('%dx%d',
                                                                              nrow(chosen.btas.mat),
                                                                              ncol(chosen.btas.mat)),
                                      '\n', '  error description: ', err)
                               })
  
  colnames(prognostic.index) <- names(chosen.btas)
  futile.logger::flog.debug('')
  futile.logger::flog.debug('prognostic.index', prognostic.index,
                            capture = TRUE)
  prognostic.index.df <- data.frame(time = c(), status = c(), group = c())
  # populate a data.frame with all patients (multiple rows per patients if has
  # multiple btas) already calculate high/low risk groups
  
  for (ix in seq_len(dim(prognostic.index)[2])) {
    # threshold
    #
    #
    temp.group <- array(-1, dim(prognostic.index)[1])
    pi.thres <- stats::quantile(prognostic.index[,ix], probs = c(probs[1],
                                                                 probs[2]))
    
    if (sum(prognostic.index[,ix] <=  pi.thres[1]) == 0 ||
        sum(prognostic.index[,ix] >  pi.thres[2]) == 0) {
      pi.thres[1] <- stats::median(unique(prognostic.index[,ix]))
      pi.thres[2] <- pi.thres[1]
    }
    
    # low risk
    temp.group[prognostic.index[,ix] <=  pi.thres[1]] <- (2 * ix) - 1
    # high risk
    temp.group[prognostic.index[,ix] > pi.thres[2]] <- (2 * ix)
    #
    valid_ix <- temp.group != -1
    #
    prognostic.index.df <- rbind(prognostic.index.df,
                                 data.frame(
                                   pi = prognostic.index[valid_ix, ix],
                                   time   = ydata$time[valid_ix],
                                   status = ydata$status[valid_ix],
                                   group  = temp.group[valid_ix]))
  }
  # factor the group
  prognostic.index.df$group <- factor(prognostic.index.df$group)
  # rename the factor to low / high risk
  new.factor.str <- as.vector(vapply(seq_along(chosen.btas), function(ix) {
    if (!is.null(names(chosen.btas)) && length(names(chosen.btas)) >= ix) {
      e <- names(chosen.btas)[ix]
      as.list(paste0(c('Low risk - ', 'High risk - '), e))
    } else {
      list('Low risk', 'High risk')
    }
  }, list(1,2)))
  
  new.factor.str.l <- as.list(as.character(seq_len(2*length(chosen.btas))))
  names(new.factor.str.l) <- new.factor.str
  
  . <- NULL # Satisfy R CMD CHECK
  prognostic.index.df$group <- prognostic.index.df$group %>%
    list %>%
    c(new.factor.str.l) %>%
    do.call(forcats::fct_collapse, .)
  #
  if (length(levels(prognostic.index.df$group)) == 1) {
    stop('separate2GroupsCox(): There is only one group, cannot create ',
         'kaplan-meir curve with low and high risk groups')
  }
  futile.logger::flog.debug('')
  futile.logger::flog.debug('prognostic.index.df', prognostic.index.df,
                            capture = TRUE)
  #
  # Generate the Kaplan-Meier survival object
  km        <- survival::survfit(survival::Surv(time, status) ~ group,
                                 data = prognostic.index.df)
  km$custom.data <- prognostic.index.df
  futile.logger::flog.debug('')
  futile.logger::flog.debug('kaplan-meier object', km, capture = TRUE)
  # Calculate the logrank test p-value
  surv.prob <- survival::survdiff(survival::Surv(time, status) ~ group,
                                  data = prognostic.index.df)
  futile.logger::flog.debug('')
  futile.logger::flog.debug('surv.prob object', surv.prob, capture = TRUE)
  p_value   <- 1 - stats::pchisq(surv.prob$chisq, df = 1)
  
  futile.logger::flog.debug('')
  futile.logger::flog.debug('pvalue: %g\n', p_value)
  
  if (no.plot) {
    return(list(pvalue = p_value, plot = NULL, km = km))
  }
  
  #
  # Plot survival curve
  #
  # remove group= from legend
  names(km$strata) <- gsub('group=','',names(km$strata))
  # if there are more than 1 btas then lines should have transparency
  if (length(chosen.btas) > 1) {
    my.alpha <- .5
  } else {
    my.alpha <- 1
  }
  
  if (length(chosen.btas) > 1) {
    col.ix <-
      loose.rock::my.colors()[c(1,2,4,3,10,6,12,9,5,7,8,
                                11,13,14,15,16,17)]
  } else {
    col.ix <- c('seagreen', 'indianred2')
  }
  
  p1 <- survminer::ggsurvplot(km,
                              conf.int = F,
                              palette = col.ix,
                              data = prognostic.index.df,
                              ggtheme = ggplot2::theme_bw(),
                              surv.median.line = "hv",  # 自己加
                              risk.table = TRUE, #  自己加
                              ...)
  
  if (expand.yzero)
    p1$plot <- p1$plot + ggplot2::expand_limits(y=.047)
  # limit the x axis if needed
  if (!is.null(xlim))
    p1$plot <- p1$plot + ggplot2::coord_cartesian(xlim=xlim, ylim = ylim)
  if (!is.null(ylim))
    p1$plot <- p1$plot + ggplot2::coord_cartesian(ylim=ylim, xlim = xlim)
  
  # if (length(chosen.btas) == 1) {
    # p1$plot <- p1$plot + ggplot2::ggtitle(paste0(gsub('_', ' ', plot.title),
                                                 # '\np_value = ',p_value))
  # } else {
    # p1$plot <- p1$plot + ggplot2::ggtitle(paste0(gsub('_', ' ', plot.title)))
  # }
  
  p1$plot <- p1$plot +
    ggplot2::labs(colour = paste0("p-value = ", format(p_value)))
  
  p1$plot <- p1$plot + ggplot2::theme(
    legend.key = ggplot2::element_blank(),
    # legend.title = ggplot2::element_text(colour = "grey10",
                                         # size = 10),
    legend.background = ggplot2::element_rect(colour = "gray")
  )
  
  if (legend.outside == TRUE)
    p1$plot <- p1$plot +
    ggplot2::theme(legend.key.size = ggplot2::unit(20,"points"))
  else
    p1$plot <- p1$plot + ggplot2::theme(legend.position = c(1,1),
                                        legend.justification = c(1, 1),
                                        legend.key.size = ggplot2::unit(20,"points"))
  
  # return p-value, plot and km object
  # return(list(pvalue = p_value, plot = p1, km = km, index = prognostic.index))
  return(list(plot = p1, km = km, index = prognostic.index))
}