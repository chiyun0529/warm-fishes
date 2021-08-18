library(mizerExperimental) # for projectToSteady() 
library(mizer)
library(tidyverse)
library(plotly)
library(tictoc) 
library(shiny) 
library(shinyWidgets) 
library(parallel)
library(optimParallel)
library(lme4)
library(gmailr)


require(cowplot)
require(gridExtra)

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

plotSummary <- function (x, y, power = 2, wlim = c(.1,NA), short = F, ...) {
  xlim = c(NA,10^log10(max(x@params@species_params$w_inf)))
  font_size = 7
  
  # need to display the legend at the bottom and only p1 has the background so using that one
  
  p1 <- plotSpectra(x, power = power, wlim = wlim, ...)
  p1 <- p1  + scale_x_continuous(limits = xlim, trans = "log10", name = "Individual size [g]") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          text = element_text(size=font_size),
          legend.position = "bottom", legend.key = element_rect(fill = "white"))+
    guides(color = guide_legend(nrow=2))
  
  mylegend<-g_legend(p1) # save the legend
  p1 <- p1 + theme(legend.position = "none") # now remove it from the plot itself
  
  p7 <- plotBiomass(x)
  p7 <- p7 + theme(legend.position = "none",
                   text = element_text(size=font_size),)
  
  if(short)
  {
    p1 <- p1 +theme(axis.title.x=element_text(),
                    axis.text.x=element_text(),
                    axis.ticks.x=element_line())
    p10 <- plot_grid(p1, p7, mylegend,  
                     rel_heights = c(3,3,1),
                     # rel_widths = c(2,1),
                     ncol = 1, align = "v", axis = "l")
    
  } else {
    
    
    p2 <- plotFeedingLevel(x, include_critical = F, ...)
    p2 <- p2 + scale_x_continuous(limits = xlim, trans = "log10") +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            text = element_text(size=font_size),
            legend.position = "none")
    p3 <- plotPredMort(x, ...)
    p3 <- p3 + scale_x_continuous(limits = xlim, trans = "log10") +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            text = element_text(size=font_size),
            legend.position = "none")
    p4 <- plotFMort(x, ...)
    p4 <- p4 + scale_x_continuous(limits = xlim, trans = "log10", name = "Individual size [g]") +
      theme(legend.position = "none",
            text = element_text(size=font_size),)
    
    # yeild and ssb |
    
    bm <- getBiomassFrame(x)
    plot_dat <- filter(bm, Year == max(unique(bm$Year)))
    # plot_dat$w_inf <- x@params@species_params$w_inf
    yieldDat <- getYield(x)
    plot_dat$yield <- yieldDat[dim(yieldDat)[1],]
    plot_dat$Year <- NULL
    plot_dat <- melt(plot_dat,"Species")
    p5 <- ggplot(plot_dat) +
      geom_bar(aes(x = Species,y = value, fill = Species, alpha = variable), stat = "identity", position = position_dodge()) +
      coord_cartesian(ylim = c(0.5*min(plot_dat$value),NA)) +
      # geom_text(aes(x = Species, y = value, label = Species), check_overlap = T)+
      
      # geom_point(aes(x = w_inf, y = Biomass, color = species)) +
      # geom_point(aes(x = w_inf, y = yield, color = species), shape = "+", size = 5) +
      # scale_x_continuous(name = "SSB and Yield") +
      scale_y_continuous(trans = "log10", name = "SSB and Yield") + #, limits = c(0.5*min(plot_dat$value),NA)) +
      scale_fill_manual(name = "Species", values = x@params@linecolour) +
      scale_alpha_manual(name = "Stat", values = c(1, 0.5), labels = c("SSB","Yield")) +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            text = element_text(size=font_size),
            legend.position = "bottom", legend.key = element_rect(fill = "white"))
    
    mylegend<-g_legend(p5) # save the legend
    p5 <- p5 + theme(legend.position = "none")
    
    # try yield divided by total biomass, we want to see the difference
    
    # r0
    
    plot_dat <- as.data.frame(getRDI(x@params)/getRDD(x@params))
    plot_dat$species <- factor(rownames(plot_dat),x@params@species_params$species)
    colnames(plot_dat)[1] <- "ratio"
    plot_dat$w_inf <- as.numeric(x@params@species_params$w_inf)
    
    # trying to have bars at their w_inf but on a continuous scale
    plot_dat$label <- plot_dat$species
    plot_dat2 <- plot_dat
    plot_dat2$ratio <- 0
    plot_dat2$label <- NA
    plot_dat <- rbind(plot_dat,plot_dat2)
    
    
    p6 <- ggplot(plot_dat) +
      geom_line(aes(x = w_inf, y = ratio, color = species), size = 15, alpha = .8) +
      geom_text(aes(x = w_inf, y = ratio, label = label),position = position_stack(vjust = 0.5), angle = 30, size = 3)+
      scale_color_manual(name = "Species", values = x@params@linecolour) +
      scale_y_continuous(name = "RDI/RDD") +
      scale_x_continuous(name = "Asymptotic size (g)", trans = "log10") +
      theme(legend.position = "none",
            text = element_text(size=font_size),) 
    
    p10 <- plot_grid(p1,p2,p3,p4, p5, p6, p7, mylegend, byrow = F,  
                     # rel_heights = c(1,1,1,2),
                     rel_widths = c(2,1),
                     nrow = 4, align = "v", axis = "l")
  }
  # p <- grid.arrange(p10,mylegend, nrow=2,heights=c(9.5,0.5))
  return(p10)
}

plotPredObsYield <-function(sim, dat, returnData = FALSE){
  ## check obs vs. predicted yield
  plot_dat <-melt(getYield(sim)[100,]/1e6)
  plot_dat$obs <- log10(dat)
  plot_dat$value <- log10(plot_dat$value)
  plot_dat$Species <-row.names(plot_dat)
  w_inf <- log10(sim@params@species_params$w_inf)
  names(w_inf) <- sim@params@species_params$species
  # window size
  winLim <- c(min(plot_dat$obs,plot_dat$value), max(plot_dat$obs,plot_dat$value))
  p <- ggplot(plot_dat) + # plot predicted and observed yields
    geom_point(aes(x = value, y = obs, color = Species, size = Species)) +
    scale_size_manual(values = w_inf) +
    scale_color_manual(values = sim@params@linecolour) +
    geom_abline(color = "black", slope = 1, intercept = 0, linetype = "dashed", alpha = .5) + 
    scale_x_continuous(name = "log10 Predicted Yield", limits = winLim) + 
    scale_y_continuous(name = "log10 Observed Yield", limits = winLim) +
    theme(legend.position = "bottom", legend.key = element_rect(fill = "white"))
  if(returnData) return(plot_dat) else return(p)
  } 

getError <- function(vary,params,dat,env=state,data_type="catch", tol = 0.1,timetorun=10) {
  
  no.species<-length(params@species_params$species)
  #env$params@species_params$R_max[]<-10^vary[1:12]
  params@species_params$R_max[]<-10^vary[1:no.species]
  
  params <- setParams(params)
  # run to steady state and update params
  # env$params<- projectToSteady(env$params, distance_func = distanceSSLogN,
  #                 tol = tol, t_max = 200,return_sim = F)
  params<- projectToSteady(params, distance_func = distanceSSLogN,
                           tol = tol, t_max = 200,return_sim = F)
  
  # create sim object 
  
  sim <- project(params, effort = 1, t_max = timetorun) #Change t_max to determine how many years the model runs for
  
  # 
  # sim <- project(env$params, effort = 1, t_max = timetorun) #Change t_max to determine how many years the model runs for
  # 
  # env$params <-sim@params
  # 
  
  ## what kind of data and output do we have?
  if (data_type=="SSB") {
    output <-getSSB(sim)[timetorun,]   #could change to getBiomass if using survey, also check units.
  }
  
  if (data_type=="catch") {
    output <-getYield(sim)[timetorun,]/1e6 
    #' using n . w . dw so g per year per volume (i.e. North Sea since kappa is set up this way). 
    #'The data are in tonnes per year so converting to tonnes.
  }
  
  pred <- log(output)
  dat  <- log(dat)
  # sum of squared errors, here on log-scale of predictions and data (could change this or use other error or likelihood options)
  discrep <- pred - dat
  discrep <- (sum(discrep^2))
  
  # can use a strong penalty on the error to ensure we reach a minimum of 10% of the data (biomass or catch) for each species
  # if(any(pred < 0.1*dat)) discrep <- discrep + 1e10
  
  return(discrep)
}

plotFmsy <- function(params, effortRes = 20, returnData = F, speciesData = NULL)
{
  # make one gear per species so we can bary the effort per species
  gear <- gear_params(params)
  gear$gear <- params@species_params$species
  gear_params(params) <- gear
  # catchability <- params@species_params$catchability
  catchability <- gear$catchability
  xlim <- 1.5 # maximum effort* catchability / xaxis limit
  # we want to vary effort value so we get a scale from 0 to 1 of effort * catchability per species
  # the "species" arg allows to run the function for only one species, which should be faster but it means "species" must also contain the result of every other species (so it's a two object list)
  if(!is.null(speciesData))
  {
    speciesName <- speciesData[[1]] # which species are we changing?
    plot_dat <- speciesData[[2]] # plot_dat of all species
    plot_dat <- filter(plot_dat, species!= speciesName) # remove previous result of the concerned species
    iSpecies <- which(params@species_params$species == speciesName)
    counter = 0 # sim counter
    # determine effort range
    effortMax <- round(xlim/catchability[iSpecies],1)+.1
    SpDat <- NULL
    
    effortSeq <- exp(seq(0,log(effortMax+1), length.out =  effortRes)) -1 
    effortSeq <- effortSeq[effortSeq<effortMax] # creating an exponentially increasing effort sequence
    for(iEffort in effortSeq)
    {
      effort_vec <- rep(1,dim(params@species_params)[1]) # all effort set to one
      effort_vec[iSpecies] <- iEffort # except that one which varies
      if(!counter )
      {
        tempSim <- project(params, effort = effort_vec, t_max = 20)
        counter <- 1
      } else  tempSim <- project(params, effort = effort_vec, t_max = 10, initial_n = tempSim@n[dim(tempSim@n)[1],,], 
                                 initial_npp = tempSim@n_pp[dim(tempSim@n_pp)[1],])
      #catch
      yieldDat <- getYield(tempSim)
      SpDat <- rbind(SpDat,c(yieldDat[dim(yieldDat)[1],iSpecies],iEffort))
    }
    SpDat <- as.data.frame(SpDat)
    SpDat$species <- params@species_params$species[iSpecies]
    SpDat$V2 <- SpDat$V2*catchability[iSpecies] # so V2 is effort * catchability
    colnames(SpDat) <- c("yield","effort","species")
    plot_dat <- rbind(plot_dat,SpDat)
    
    
    
  } else {
    plot_dat <- NULL
    for(iSpecies in 1:dim(params@species_params)[1])
    {
      counter = 0 # sim counter
      # determine effort range
      effortMax <- round(xlim/catchability[iSpecies],1)+.1
      SpDat <- NULL
      
      effortSeq <- exp(seq(0,log(effortMax+1), length.out =  effortRes)) -1 # every .1 takes 2 min to run, evry .2 takes 1 min but lesser resolution
      effortSeq <- effortSeq[effortSeq<effortMax] # creating an exponentially increasing effort sequence
      for(iEffort in effortSeq)
      {
        effort_vec <- rep(1,dim(params@species_params)[1]) # all effort set to one
        effort_vec[iSpecies] <- iEffort # except that one which varies
        if(!counter )
        {
          tempSim <- project(params, effort = effort_vec, t_max = 20)
          counter <- 1
        } else  tempSim <- project(params, effort = effort_vec, t_max = 10, initial_n = tempSim@n[dim(tempSim@n)[1],,], 
                                   initial_npp = tempSim@n_pp[dim(tempSim@n_pp)[1],])
        #catch
        yieldDat <- getYield(tempSim)
        SpDat <- rbind(SpDat,c(yieldDat[dim(yieldDat)[1],iSpecies],iEffort))
      }
      SpDat <- as.data.frame(SpDat)
      SpDat$species <- params@species_params$species[iSpecies]
      SpDat$V2 <- SpDat$V2*catchability[iSpecies] # so V2 is effort * catchability
      colnames(SpDat) <- c("yield","effort","species")
      plot_dat <- rbind(plot_dat,SpDat)
      
    }
  }
  
  plot_dat$species <- factor(plot_dat$species, levels = params@species_params$species)
  # colnames(plot_dat) <- c("yield","effort","species")
  if(!is.null(speciesData)) p <- ggplot(filter(plot_dat, species == speciesName)) else p <- ggplot(plot_dat)
  p <- p + geom_line(aes(x = effort , y = yield, color = species))+
    facet_wrap(species~., scales = "free") +
    scale_x_continuous(limits= c(0,xlim),name = "fishing mortality rate")+#, limits = c(1e10,NA))+
    scale_y_continuous(trans = "log10") +
    scale_color_manual(name = "Species", values = params@linecolour) +
    theme(legend.position = "none", legend.key = element_rect(fill = "white"),
          panel.background = element_blank(), panel.grid.minor = element_line(color = "gray"),
          strip.background = element_blank())
  if(returnData) return(plot_dat) else return(p)
}

plotGrowthCurves2 <- function (object, 
                               species = NULL, 
                               max_age = 20, 
                               percentage = FALSE, 
                               species_panel = FALSE, 
                               highlight = NULL,
                               returnData = F) 
{
  if (is(object, "MizerSim")) {
    params <- object@params
    t <- dim(object@n)[1]
    params@initial_n[] <- object@n[t, , ]
    params@initial_n_pp <- object@n_pp[t, ]
  }
  else if (is(object, "MizerParams")) {
    params <- validParams(object)
  }
  species <- valid_species_arg(params, species)
  ws <- getGrowthCurves(params, species, max_age, percentage)
  plot_dat <- reshape2::melt(ws)
  plot_dat$Species <- factor(plot_dat$Species, params@species_params$species)
  plot_dat$legend <- "model"
  if (all(c("a", "b", "k_vb") %in% names(params@species_params))) {
    if ("t0" %in% names(params@species_params)) {
      t0 <- params@species_params$t0
    }
    else {
      t0 <- 0
    }
    VBdf <- data.frame(species = params@species_params$species, 
                       w_inf = params@species_params$w_inf, a = params@species_params$a, 
                       b = params@species_params$b, k_vb = params@species_params$k_vb, 
                       t0 = t0)
    VBdf$L_inf <- (VBdf$w_inf/VBdf$a)^(1/VBdf$b)
    plot_dat2 <- plot_dat
    plot_dat2$value <- apply(plot_dat, 1, function(x) {
      sel <- VBdf$species == x[1]
      bodymass <- VBdf$w_inf[sel] * (1 - exp(-VBdf$k_vb[sel] * 
                                             (as.numeric(x[2]) - VBdf$t0[sel])))
    })
    plot_dat2$legend <- "von Bertalanffy"
    plot_dat <- rbind(plot_dat, plot_dat2)
  }
  p <- ggplot(filter(plot_dat, legend == "model")) + geom_line(aes(x = Age, 
                                                                   y = value, colour = Species, linetype = Species, size = Species))
  y_label <- if (percentage) "Percent of maximum size" else "Size [g]"
  linesize <- rep(0.8, length(params@linetype))
  names(linesize) <- names(params@linetype)
  linesize[highlight] <- 1.6
  p <- p + scale_x_continuous(name = "Age [Years]") + scale_y_continuous(name = y_label) + 
    scale_colour_manual(values = params@linecolour) + scale_linetype_manual(values = params@linetype) + 
    scale_size_manual(values = linesize)
  if (!percentage) {
    if (length(species) == 1) {
      idx <- which(params@species_params$species == species)
      w_inf <- params@species_params$w_inf[idx]
      p <- p + geom_hline(yintercept = w_inf, colour = "grey") + 
        annotate("text", 0, w_inf, vjust = -1, label = "Maximum")
      w_mat <- params@species_params$w_mat[idx]
      p <- p + geom_hline(yintercept = w_mat, linetype = "dashed", 
                          colour = "grey") + annotate("text", 0, w_mat, 
                                                      vjust = -1, label = "Maturity")
      if ("von Bertalanffy" %in% plot_dat$legend) 
        p <- p + geom_line(data = filter(plot_dat, legend == 
                                           "von Bertalanffy"), aes(x = Age, y = value))
    }
    else if (species_panel) {
      p <- ggplot(plot_dat) + 
        geom_line(aes(x = Age, y = value, colour = legend)) + 
        scale_x_continuous(name = "Age [years]") +
        scale_y_continuous(name = "Size [g]") +
        facet_wrap(.~Species, scales = "free") +
        geom_hline(aes(yintercept = w_mat),
                   data = tibble(Species = as.factor(object@params@species_params$species[]),
                                 w_mat = object@params@species_params$w_mat[]),
                   linetype = "dashed", colour = "grey") +
        geom_hline(aes(yintercept = w_inf),
                   data = tibble(Species = as.factor(object@params@species_params$species[]),
                                 w_inf = object@params@species_params$w_inf[]),
                   linetype = "solid", colour = "grey") +
        theme(panel.background = element_blank(), panel.grid.minor = element_line(color = "gray"),
              strip.background = element_blank(), legend.key = element_blank())+
        scale_color_discrete(name = "Growth", labels = c("Modelled","von Bertalanffy"))
    }
  }
  if(returnData) return(plot_dat) else return(p)
}

plotDiet2 <- function (sim, species = NULL, xlim = c(1,NA), returnData = F) 
{
  params <- sim@params
  # if (is.integer(species)) {
  #     species <- params@species_params$species[species]
  # }
  
  # diet <- getDiet(params)[params@species_params$species == 
  #     species, , ]
  # prey <- dimnames(diet)$prey
  # prey <- factor(prey, levels = rev(prey))
  # plot_dat <- data.frame(Proportion = c(diet), w = params@w, 
  #     Prey = rep(prey, each = length(params@w)))
  # plot_dat <- plot_dat[plot_dat$Proportion > 0, ]
  # 
  # ggplot(plot_dat) + geom_area(aes(x = w, y = Proportion, fill = Prey)) + 
  #     scale_x_log10(limits = xlim) + labs(x = "Size [g]") + 
  #   scale_fill_manual(values = sim@params@linecolour) +
  #   ggtitle(species)
  
  
  diet <- getDiet(params)
  plot_dat <- melt(diet)
  plot_dat <- plot_dat[plot_dat$value > 0, ]
  colnames(plot_dat) <- c("Predator", "size", "Prey", "Proportion")
  
  if(is.null(species)) p <- ggplot(plot_dat) + facet_wrap(.~Predator, scales = "free") else p <- ggplot(filter(plot_dat, Predator == species))
  
  p <- p +
    geom_area(aes(x = size, y = Proportion, fill = Prey))+
    scale_x_continuous(limits = c(1,NA), name = "Size [g]", trans = "log10") + 
    scale_fill_manual(values = sim@params@linecolour)+
    theme(legend.position = "right", legend.key = element_rect(fill = "white"),
          panel.background = element_blank(), panel.grid.minor = element_line(color = "gray"),
          strip.background = element_blank())
  
  if(returnData) return(plot_dat) else return(p)
  
}


shiny_erepro <- function(input, dat = NULL) {
  params_shiny <- input
  no_sp<-length(params_shiny@species_params$species)
  ui=fluidPage(
    
    # Application title
    titlePanel("erepro calibration"),
    
    fluidRow(
      column(4, wellPanel(
        #sliderInput("kappa", "log10 Resource Carrying Capacity:", min = 8, max = 12, value = log10(params_shiny@resource_params$kappa),
        #            step = 0.1),
        #   sliderInput("Rmax", "log10 Maximum Recruitment:", min = 1, max = 12, value = 12,
        #              step = 0.1),
        sliderInput("erepro", "log10 Reproductive Efficiency:", min = -8, max = 1, value = log10(params_shiny@species_params$erepro[1]),
                    step = 0.1)
      )),
      column(6,
             plotOutput("plot1", width = 600, height = 600),
             if(!is.null(dat)) plotOutput("plot2", width = 600, height = 600)
      ))
    
    
    
  )
  server = function(input, output) {
    
    #    sim <- observeEvent(input$erepro, {
    #       print("yo")
    #       params_shiny@species_params$erepro <- rep(10^input$erepro,12)
    #       params_shiny <- setParams(params_shiny)
    #       sim_shiny <- project(params_shiny, effort = 1, t_max = 50)
    # })
    #     output$plot1 <- renderPlot({
    #       plot(sim())
    #     })
    
    output$plot1 <- renderPlot({
      # set up params using values given, need check and change parameter values so units work in days units
      params_shiny@species_params$erepro <- rep(10^input$erepro,no_sp)
      # params@species_params$Rmax <- rep(10^input$Rmax,12)
      params_shiny <- setParams(params_shiny)#,kappa=10^input$kappa)
      # run without fishing
      sim_shiny <- project(params_shiny, effort = 1, t_max = 50)
      plot(sim_shiny)
    })
    
    if(!is.null(dat))
    {
      output$plot2 <- renderPlot({
        # set up params using values given, need check and change parameter values so units work in days units
        params_shiny@species_params$erepro <- rep(10^input$erepro,no_sp)
        # params@species_params$Rmax <- rep(10^input$Rmax,12)
        params_shiny <- setParams(params_shiny)#,kappa=10^input$kappa)
        # run without fishing
        sim_shiny <- project(params_shiny, effort = 1, t_max = 100)
        plotPredObsYield(sim_shiny,dat)
        # plotBiomass(sim_shiny)
      })
    }
  }
  shinyApp(ui, server)
}



#' shiny app to tweak erepro per species and see the effects on Fmsy
#'
#' @export
shiny_fmsy <- function(params,dat) {
  
  
  params_shiny <- params
  FmsyDat <- dat
  
  ui=fluidPage(
    
    # Application title
    titlePanel("North Sea Fmsy"),
    
    fluidRow(
      column(4, wellPanel(
        # sliderInput("kappa", "log10 Resource Carrying Capacity:", min = 8, max = 12, value = log10(params_optim@resource_params$kappa),
        #             step = 0.1),
        #   sliderInput("Rmax", "log10 Maximum Recruitment:", min = 1, max = 12, value = 12,
        #              step = 0.1),
        sliderInput("erepro", "log10 Reproductive Efficiency:", min = -8, max = 1, value = log10(params_shiny@species_params$erepro[1]),
                    step = 0.1),
        sliderTextInput(
          inputId = "species",
          label = "Species name",
          choices = params_shiny@species_params$species,
          selected = params_shiny@species_params$species[1],
          grid = T
        ),
      )),
      
      column(6,
             plotOutput("distPlot", width = 600, height = 600)
      ))
    
    
    
  )
  
  server = function(input, output) {
    
    output$distPlot <- renderPlot({
      # set up params using values given, need check and change parameter values so units work in days units
      params_shiny@species_params$erepro[which(params_shiny@species_params$species == input$species)] <-10^(input$erepro)
      
      plotFmsy(params_shiny, speciesData = list(input$species,FmsyDat), effortRes = 10)
    })
  }
  # add species slider to display jsut one species at a time
  
  shinyApp(ui = ui, server = server)
  
}




#' shiny app to tweak gamma per species and see the effects on growth
#'
#' @export

shiny_gamma <- function(params, dat = NULL)
{
  
  params_shiny <- params
  
  ui=fluidPage(
    titlePanel("Gamma calibration"),
    fluidRow(
      column(4,
             wellPanel(
               actionButton("run", "Run the simulation")
             ),
             wellPanel(
               "Ap intitial value:",
               textOutput("ApGamma"),
               numericInput("gamma1", "Search volume:", min = .1*params_shiny@species_params$gamma[1], max = 10*params_shiny@species_params$gamma[1],
                            value = params_shiny@species_params$gamma[1]),
               "As intitial value:",
               textOutput("AsGamma"),
               numericInput("gamma2", "Search volume:", min = .1*params_shiny@species_params$gamma[2], max = 10*params_shiny@species_params$gamma[2],
                            value = params_shiny@species_params$gamma[2], step = 0.01*params_shiny@species_params$gamma[2]
               ),
               "Ad intitial value:",
               textOutput("AdGamma"),
               numericInput("gamma3", "Search volume:", min = .1*params_shiny@species_params$gamma[3], max = 10*params_shiny@species_params$gamma[3],
                            value = params_shiny@species_params$gamma[3], step = 0.01*params_shiny@species_params$gamma[3]
               ),
               "Am intitial value:",
               textOutput("AmGamma"),
               numericInput("gamma4", "Search volume:", min = .1*params_shiny@species_params$gamma[4], max = 10*params_shiny@species_params$gamma[4],
                            value = params_shiny@species_params$gamma[4], step = 0.01*params_shiny@species_params$gamma[4]
               ),
               "Ch intitial value:",
               textOutput("ChGamma"),
               numericInput("gamma5", "Search volume:", min = .1*params_shiny@species_params$gamma[5], max = 10*params_shiny@species_params$gamma[5],
                            value = params_shiny@species_params$gamma[5], step = 0.01*params_shiny@species_params$gamma[5]
               ),
               "Gm intitial value:",
               textOutput("GmGamma"),
               numericInput("gamma6", "Search volume:", min = .1*params_shiny@species_params$gamma[6], max = 10*params_shiny@species_params$gamma[6],
                            value = params_shiny@species_params$gamma[6], step = 0.01*params_shiny@species_params$gamma[6]
               ),
               "Gc intitial value:",
               textOutput("GcGamma"),
               numericInput("gamma7", "Search volume:", min = .1*params_shiny@species_params$gamma[7], max = 10*params_shiny@species_params$gamma[7],
                            value = params_shiny@species_params$gamma[7], step = 0.01*params_shiny@species_params$gamma[7]
               ),
               "Lf intitial value:",
               textOutput("LfGamma"),
               numericInput("gamma8", "Search volume:", min = .1*params_shiny@species_params$gamma[8], max = 10*params_shiny@species_params$gamma[8],
                            value = params_shiny@species_params$gamma[8], step = 0.01*params_shiny@species_params$gamma[8]
               ),
               "La intitial value:",
               textOutput("LaGamma"),
               numericInput("gamma9", "Search volume:", min = .1*params_shiny@species_params$gamma[9], max = 10*params_shiny@species_params$gamma[9],
                            value = params_shiny@species_params$gamma[9], step = 0.01*params_shiny@species_params$gamma[9]
               ),
               "Ma intitial value:",
               textOutput("MaGamma"),
               numericInput("gamma10", "Search volume:", min = .1*params_shiny@species_params$gamma[10], max = 10*params_shiny@species_params$gamma[10],
                            value = params_shiny@species_params$gamma[10], step = 0.01*params_shiny@species_params$gamma[10]
               ),
               "Mb intitial value:",
               textOutput("MbGamma"),
               numericInput("gamma11", "Search volume:", min = .1*params_shiny@species_params$gamma[11], max = 10*params_shiny@species_params$gamma[11],
                            value = params_shiny@species_params$gamma[11], step = 0.01*params_shiny@species_params$gamma[11]
               ),
               "Mo intitial value:",
               textOutput("MoGamma"),
               numericInput("gamma12", "Search volume:", min = .1*params_shiny@species_params$gamma[12], max = 10*params_shiny@species_params$gamma[12],
                            value = params_shiny@species_params$gamma[12], step = 0.01*params_shiny@species_params$gamma[12]
               ),
               "Pd intitial value:",
               textOutput("PdGamma"),
               numericInput("gamma13", "Search volume:", min = .1*params_shiny@species_params$gamma[13], max = 10*params_shiny@species_params$gamma[13],
                            value = params_shiny@species_params$gamma[13], step = 0.01*params_shiny@species_params$gamma[13]
               ),
               "Pa intitial value:",
               textOutput("PaGamma"),
               numericInput("gamma14", "Search volume:", min = .1*params_shiny@species_params$gamma[14], max = 10*params_shiny@species_params$gamma[14],
                            value = params_shiny@species_params$gamma[14], step = 0.01*params_shiny@species_params$gamma[14]
               ),
               "Re intitial value:",
               textOutput("ReGamma"),
               numericInput("gamma15", "Search volume:", min = .1*params_shiny@species_params$gamma[15], max = 10*params_shiny@species_params$gamma[15],
                            value = params_shiny@species_params$gamma[15], step = 0.01*params_shiny@species_params$gamma[15]
               ),
               "Sa intitial value:",
               textOutput("SaGamma"),
               numericInput("gamma16", "Search volume:", min = .1*params_shiny@species_params$gamma[16], max = 10*params_shiny@species_params$gamma[16],
                            value = params_shiny@species_params$gamma[16], step = 0.01*params_shiny@species_params$gamma[16]
               ),
               "Uc intitial value:",
               textOutput("UcGamma"),
               numericInput("gamma17", "Search volume:", min = .1*params_shiny@species_params$gamma[17], max = 10*params_shiny@species_params$gamma[17],
                            value = params_shiny@species_params$gamma[17], step = 0.01*params_shiny@species_params$gamma[17]
               ),
               "Ut intitial value:",
               textOutput("UtGamma"),
               numericInput("gamma18", "Search volume:", min = .1*params_shiny@species_params$gamma[18], max = 10*params_shiny@species_params$gamma[18],
                            value = params_shiny@species_params$gamma[18], step = 0.01*params_shiny@species_params$gamma[18]
               ),
               "Za intitial value:",
               textOutput("ZaGamma"),
               numericInput("gamma19", "Search volume:", min = .1*params_shiny@species_params$gamma[19], max = 10*params_shiny@species_params$gamma[19],
                            value = params_shiny@species_params$gamma[19], step = 0.01*params_shiny@species_params$gamma[19]
               )
             )
      ),
      column(6,
             plotOutput("plot1", width = 600, height = 600),
             plotOutput("plot2", width = 600, height = 600),
             plotOutput("plot3", width = 600, height = 600),
             if(!is.null(dat)) plotOutput("plot4", width = 600, height = 600)
      )
    )
    
  )
  
  server = function(input, output) {
    output$ApGamma <- renderText({
      params_shiny@species_params$gamma[1]
    })
    output$AsGamma <- renderText({
      params_shiny@species_params$gamma[2]
    })
    output$AdGamma <- renderText({
      params_shiny@species_params$gamma[3]
    })
    output$AmGamma <- renderText({
      params_shiny@species_params$gamma[4]
    })
    output$ChGamma <- renderText({
      params_shiny@species_params$gamma[5]
    })
    output$GmGamma <- renderText({
      params_shiny@species_params$gamma[6]
    })
    output$GcGamma <- renderText({
      params_shiny@species_params$gamma[7]
    })
    output$LfGamma <- renderText({
      params_shiny@species_params$gamma[8]
    })
    output$LaGamma <- renderText({
      params_shiny@species_params$gamma[9]
    })
    output$MaGamma <- renderText({
      params_shiny@species_params$gamma[10]
    })
    output$MbGamma <- renderText({
      params_shiny@species_params$gamma[11]
    })
    output$MoGamma <- renderText({
      params_shiny@species_params$gamma[12]
    })
    output$PdGamma <- renderText({
      params_shiny@species_params$gamma[13]
    })
    output$PaGamma <- renderText({
      params_shiny@species_params$gamma[14]
    })
    output$ReGamma <- renderText({
      params_shiny@species_params$gamma[15]
    })
    output$SaGamma <- renderText({
      params_shiny@species_params$gamma[16]
    })
    output$UcGamma <- renderText({
      params_shiny@species_params$gamma[17]
    })
    output$UtGamma <- renderText({
      params_shiny@species_params$gamma[18]
    })
    output$ZaGamma <- renderText({
      params_shiny@species_params$gamma[19]
    })
    # reactive expression
    sim <- eventReactive(input$run, {
      params_shiny@species_params$gamma <- c(input$gamma1,input$gamma2,input$gamma3,input$gamma4,input$gamma5,input$gamma6,
                                             input$gamma7,input$gamma8,input$gamma9,input$gamma10,input$gamma11,input$gamma12,
                                             input$gamma13,input$gamma14,input$gamma15,input$gamma16,input$gamma17,input$gamma18,
                                             input$gamma19)
      params_shiny <- setParams(params_shiny)
      # print(params_shiny@species_params$gamma) # check if everything is going well
      sim_shiny <- project(params_shiny, effort = 1, t_max = 100)
    })
    
    output$plot1 <- renderPlot({
      plotGrowthCurves2(sim(), species_panel = T)
    })
    
    output$plot2 <- renderPlot({
      plotFeedingLevel2(sim(), include_critical = T)
    })
    
    output$plot3 <- renderPlot({
      plotCalibration(sim())
    })
    if(!is.null(dat))
    {output$plot4 <- renderPlot({
      plotPredObsYield(sim(),dat)
    })
    }
  }
  
  shinyApp(ui, server)
}


#' @export
shiny_kappa <- function(params, dat = NULL)
{
  params_shiny <- params
  
  ui=fluidPage(
    
    # Application title
    titlePanel("kappa calibration"),
    
    fluidRow(
      column(4, wellPanel(
        sliderInput("kappa", "log10 Resource Carrying Capacity:", min = 8, max = 12, value = log10(params_shiny@resource_params$kappa),
                    step = 0.1),
      )),
      column(6,
             plotOutput("plot1", width = 600, height = 600),
             if(!is.null(dat)) plotOutput("plot2", width = 600, height = 600)
      ))
    
    
    
  )
  
  server = function(input, output) {
    
    output$plot1 <- renderPlot({
      
      params_shiny <- setParams(params_shiny,kappa=10^input$kappa)
      sim_shiny <- project(params_shiny, effort = 1, t_max = 100)
      plotGrowthCurves2(sim_shiny, species_panel = T)
    })
    
    
    if(!is.null(dat))
    {
      output$plot2 <- renderPlot({
        params_shiny <- setParams(params_shiny,kappa=10^input$kappa)
        sim_shiny <- project(params_shiny, effort = 1, t_max = 100)
        plotPredObsYield(sim_shiny,dat)
      })
    }
  }
  
  
  shinyApp(ui, server)
}

checkGrowth<-function(sim){
  
  sim_loop<-sim
  
  ws <- getGrowthCurves(sim_loop)
  plot_dat <- reshape2::melt(ws)
  
  VBdf <- data.frame(species = sim_loop@params@species_params$species, 
                     w_inf = sim_loop@params@species_params$w_inf, 
                     k_vb = sim_loop@params@species_params$k_vb, 
                     t0 = 0)
  
  plot_dat2<-plot_dat
  
  plot_dat2$value<-apply(plot_dat, 1, function(x) {
    sel <- VBdf$species == x[1]
    bodymass <- VBdf$w_inf[sel] * (1 - exp(-VBdf$k_vb[sel] * (as.numeric(x[2]) - 0.001)))
  })
  
  plot_dat<-rbind(plot_dat,plot_dat2)
  
  plot_dat$legend<-c(rep("model",950),rep("von Bertanlaffy",950))
  
  p <- ggplot(plot_dat) + 
    geom_line(aes(x = Age, y = value, colour = legend)) + 
    scale_x_continuous(name = "Age [years]") +
    scale_y_continuous(name = "Size [g]") +
    facet_wrap(.~Species, scales = "free") +
    geom_hline(aes(yintercept = w_mat),
               data = tibble(Species = as.factor(sim_loop@params@species_params$species[]),
                             w_mat = sim_loop@params@species_params$w_mat[]),
               linetype = "dashed", colour = "grey") +
    geom_hline(aes(yintercept = w_inf),
               data = tibble(Species = as.factor(sim_loop@params@species_params$species[]),
                             w_inf = sim_loop@params@species_params$w_inf[]),
               linetype = "solid", colour = "grey") +
    theme(panel.background = element_blank(), panel.grid.minor = element_line(color = "gray"),
          strip.background = element_blank(), legend.key = element_blank())+
    scale_color_discrete(name = "Growth", labels = c("Modelled","von Bertalanffy"))
  
  p
  
}
