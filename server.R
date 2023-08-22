library(shiny) 
library(survival) #provides tools for survival analysis
library(tidyverse) #data manipulation and coding
library(ggsurvfit) #has an useful plot function
library(ggplot2) #best graphing library, tied with Plotly
library(plotly) #great for interactive charts
library(broom) #for exporting to well organized dfs
library(shinyjs) #dynamic of some components
library(arsenal) #among many things, creates good tables
library(xtable)


#This is a list with all the available distributions in survfit
available_distributions = names(survreg.distributions)

#Color palette
colorsPalette <- c(CompPFS= "#000000", CompOS = "#E69F00", PFS = "#56B4E9", OS = "#009E73")#, "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#Loading all the data sets
#stack them in a single data set
#filter if necessary
#code will be cleaner in later stages because of this
#eventually, for larger applications, data will probably come from SQL tables in this format
dataCompPFS = read.csv("./data/COMP_PFS.csv", sep = ";", header = TRUE, col.names = c("time", "survival"))
dataCompPFS$intervention = "Non-intervention"
dataCompPFS$state = "PFS"
dataCompOS = read.csv("./data/COMP_OS.csv", sep = ";", header = TRUE, col.names = c("time", "survival"))
dataCompOS$intervention = "Non-intervention"
dataCompOS$state = "OS"
dataPFS = read.csv("./data/INT_PFS.csv", sep = ";", header = TRUE, col.names = c("time", "survival"))
dataPFS$intervention = "Intervention"
dataPFS$state = "PFS"
dataOS = read.csv("./data/INT_OS.csv", sep = ";", header = TRUE, col.names = c("time", "survival"))
dataOS$intervention = "Intervention"
dataOS$state = "OS"
data = rbind(dataCompPFS,dataCompOS,dataPFS,dataOS) %>%
  arrange(time)

#models
#All models will be fitted to every data set as to prevent
#a lot of estimation from an user trying to visualize and compare them
# I tried to make this very general, to avoid code repetition, using do.call
# I will come back to it later and do the basic for now.
#Done is better than perfect

#Weibull models
weibull_estimate <- function(df){
  params = lm(log(-log(df$survival))~log(df$time))$coefficients
  k = params[2]
  lambda = exp(-params[1]/k)
  return(c(k = as.numeric(k), lambda = as.numeric(lambda)))
}

weibull_cost = function(params){
  s_x = 1 - pweibull(df$time, params[1],params[2])
  return(sum((s_x-df$survival)^2))
}

#Comp PFS
df = dataCompPFS
model_weibull_compPFS = weibull_estimate(df)

#PFS
df = dataPFS
model_weibull_PFS =  weibull_estimate(df)

#Comp OS
df = dataCompOS
model_weibull_compOS = weibull_estimate(df)

#OS
df = dataOS
model_weibull_OS = weibull_estimate(df)

fig_weibull_PFS = ggplot(data = dataCompPFS,aes(time, survival, color = "CompPFS"))+
  geom_step()+
  geom_step(data = dataPFS, aes(time,survival, color = "PFS"))+
  geom_line(data = dataCompPFS, aes(time,1-pweibull(time,model_weibull_compPFS[1], model_weibull_compPFS[2]), color = "CompPFS"))+
  geom_line(data = dataPFS, aes(time,1-pweibull(time, model_weibull_PFS[1],model_weibull_PFS[2]), color = "PFS"))+
  scale_color_manual(values = colorsPalette)+
  ggtitle("Weibull")+
  labs(color ="Intervention")

fig_weibull_OS = ggplot(data = dataCompOS,aes(time, survival, color = "CompOS"))+
  geom_step()+
  geom_step(data = dataOS, aes(time,survival, color = "OS"))+
  geom_line(data = dataCompOS, aes(time,1-pweibull(time,model_weibull_compOS[1], model_weibull_compOS[2]), color = "CompOS"))+
  geom_line(data = dataOS, aes(time,1-pweibull(time, model_weibull_OS[1],model_weibull_OS[2]), color = "OS"))+
  scale_color_manual(values = colorsPalette)+
  ggtitle("Weibull")+
  labs(color ="Intervention")

#Log-logistic models
loglogistic_estimate <- function(df){
  params = lm(log(1/df$survival -1 )~log(df$time))$coefficients
  beta = params[2]
  alpha = exp(-params[1]/beta)
  return(c(alpha = as.numeric(alpha), lambda = as.numeric(beta)))
}

pllogis = function(x,alpha,beta){
  return (1/(1+(x/alpha)^(-beta)))
}

#Comp PFS
df = dataCompPFS
model_loglogistic_compPFS = loglogistic_estimate(df)

#PFS
df = dataPFS
model_loglogistic_PFS =  loglogistic_estimate(df)

#Comp OS
df = dataCompOS
model_loglogistic_compOS = loglogistic_estimate(df)

#OS
df = dataOS
model_loglogistic_OS = loglogistic_estimate(df)

fig_loglogistic_PFS = ggplot(data = dataCompPFS,aes(time, survival, color = "CompPFS"))+
  geom_step()+
  geom_step(data = dataPFS, aes(time,survival, color = "PFS"))+
  geom_line(data = dataCompPFS, aes(time,1-pllogis(time,model_loglogistic_compPFS[1], model_loglogistic_compPFS[2]), color = "CompPFS"))+
  geom_line(data = dataPFS, aes(time,1-pllogis(time, model_loglogistic_PFS[1],model_loglogistic_PFS[2]), color = "PFS"))+
  scale_color_manual(values = colorsPalette)+
  ggtitle("loglogistic")+
  labs(color ="Intervention")

fig_loglogistic_OS = ggplot(data = dataCompOS,aes(time, survival, color = "CompOS"))+
  geom_step()+
  geom_step(data = dataOS, aes(time,survival, color = "OS"))+
  geom_line(data = dataCompOS, aes(time,1-pllogis(time,model_loglogistic_compOS[1], model_loglogistic_compOS[2]), color = "CompOS"))+
  geom_line(data = dataOS, aes(time,1-pllogis(time, model_loglogistic_OS[1],model_loglogistic_OS[2]), color = "OS"))+
  scale_color_manual(values = colorsPalette)+
  ggtitle("loglogistic")+
  labs(color ="Intervention")

#Exponential models
exponential_estimate <- function(df){
  params = lm(log(df$survival)~ -1 + df$time)$coefficients
  lambda = -params[1]
  return(c(lambda = as.numeric(lambda)))
}

exponential_cost = function(params){
  s_x = 1 - pexp(df$time, params[1],params[2])
  return(sum((s_x-df$survival)^2))
}

#Comp PFS
df = dataCompPFS
model_exponential_compPFS = exponential_estimate(df)

#PFS
df = dataPFS
model_exponential_PFS =  exponential_estimate(df)

#Comp OS
df = dataCompOS
model_exponential_compOS = exponential_estimate(df)

#OS
df = dataOS
model_exponential_OS = exponential_estimate(df)

fig_exponential_PFS = ggplot(data = dataCompPFS,aes(time, survival, color = "CompPFS"))+
  geom_step()+
  geom_step(data = dataPFS, aes(time,survival, color = "PFS"))+
  geom_line(data = dataCompPFS, aes(time,1-pexp(time,model_exponential_compPFS[1]), color = "CompPFS"))+
  geom_line(data = dataPFS, aes(time,1-pexp(time, model_exponential_PFS[1]), color = "PFS"))+
  scale_color_manual(values = colorsPalette)+
  ggtitle("Exponential")+
  labs(color ="Intervention")

fig_exponential_OS = ggplot(data = dataCompOS,aes(time, survival, color = "CompOS"))+
  geom_step()+
  geom_step(data = dataOS, aes(time,survival, color = "OS"))+
  geom_line(data = dataCompOS, aes(time,1-pexp(time,model_exponential_compOS[1]), color = "CompOS"))+
  geom_line(data = dataOS, aes(time,1-pexp(time, model_exponential_OS[1]), color = "OS"))+
  scale_color_manual(values = colorsPalette)+
  ggtitle("Exponential")+
  labs(color ="Intervention")

#rayleigh models
prayleigh = function(x,sigma){
  return(1-exp(-(x^2)/(2*sigma^2)))
}

rayleigh_estimate <- function(df){
  params = mean(exp(log(-log(df$survival))-2*log(df$time))/2-log(sqrt(2))) 
  return(c(sigma = params))
}

rayleigh_cost = function(params){
  s_x = 1 - pexp(df$time, params[1],params[2])
  return(sum((s_x-df$survival)^2))
}

#Comp PFS
df = dataCompPFS
model_rayleigh_compPFS = rayleigh_estimate(df)

#PFS
df = dataPFS
model_rayleigh_PFS =  rayleigh_estimate(df)

#Comp OS
df = dataCompOS
model_rayleigh_compOS = rayleigh_estimate(df)

#OS
df = dataOS
model_rayleigh_OS = rayleigh_estimate(df)

fig_rayleigh_PFS = ggplot(data = dataCompPFS,aes(time, survival, color = "CompPFS"))+
  geom_step()+
  geom_step(data = dataPFS, aes(time,survival, color = "PFS"))+
  geom_line(data = dataCompPFS, aes(time,1-prayleigh(time,model_rayleigh_compPFS[1]), color = "CompPFS"))+
  geom_line(data = dataPFS, aes(time,1-prayleigh(time, model_rayleigh_PFS[1]), color = "PFS"))+
  scale_color_manual(values = colorsPalette)+
  ggtitle("Rayleigh")+
  labs(color ="Intervention")

fig_rayleigh_OS = ggplot(data = dataCompOS,aes(time, survival, color = "CompOS"))+
  geom_step()+
  geom_step(data = dataOS, aes(time,survival, color = "OS"))+
  geom_line(data = dataCompOS, aes(time,1-prayleigh(time,model_rayleigh_compOS[1]), color = "CompOS"))+
  geom_line(data = dataOS, aes(time,1-prayleigh(time, model_rayleigh_OS[1]), color = "OS"))+
  scale_color_manual(values = colorsPalette)+
  ggtitle("Rayleigh")+
  labs(color ="Intervention")

#main server code
shinyServer(function(input, output) {
  
    #get data only when the button is pressed
    output$histPlot = renderPlot({
      if (input$dataset == "All data"){
      #Produce a survival plot with all the data sets
        fig <-  ggplot(data = dataCompPFS, aes(time,survival, color = "CompPFS")) +
          geom_step() +
          geom_step(data = dataCompOS, aes(time, survival, color = "CompOS")) +
          geom_step(data = dataOS, aes(time, survival, color = "OS")) +
          geom_step(data = dataPFS, aes(time, survival, color = "PFS")) +
          scale_color_manual(values = colorsPalette)+
          labs(color = "Data sets")+
          ggtitle("All data")
      } else if (input$dataset == "Progression-free survival"){
        fig <- ggplot(data = dataCompPFS, aes(time, survival, color = "CompPFS"))+
          geom_step()+
          geom_step(data = dataPFS, aes(time, survival, color = "PFS"))+
          scale_color_manual(values = colorsPalette)+
          labs(color = "Intervention")+
          ggtitle("Progression-free survival")
      } else if (input$dataset == "Overall survival"){
        fig <- ggplot(data = dataCompOS, aes(time, survival, color = "CompOS"))+
          geom_step()+
          geom_step(data = dataOS, aes(time, survival, color = "OS"))+
          scale_color_manual(values = colorsPalette)+
          labs(color = "Intervention")+
          ggtitle("Overall survival")
      } else {
        print("Unknown option for data set")
      }
      return (fig)
    })
    
    output$modelPlotPFS <- renderPlot({
      if(input$dist == "Weibull"){
        return(fig_weibull_PFS)
      } else if(input$dist == "Log-logistic"){
        return(fig_loglogistic_PFS)
      } else if(input$dist == "Exponential"){
        return(fig_exponential_PFS)
      }
    })
    
    output$modelPlotOS <- renderPlot({
      if(input$dist == "Weibull"){
        return(fig_weibull_OS)
      } else if(input$dist == "Log-logistic"){
        return (fig_loglogistic_OS)
      } else if(input$dist == "Exponential"){
        return(fig_exponential_OS)
      }
    })
    
    output$modelTablePFS <- renderUI({
      if(input$dist == "Weibull"){
        tags$table(
          tags$tr(
            tags$th("Parameter"),
            tags$th("Estimate"),
            tags$th("Standard Error"),
          ),
          tags$tr(
            tags$td("Non-intervention")
          ),
          tags$tr(
            tags$td("k"),
            tags$td(model_weibull_compPFS["k"]),
            tags$td("TODO")
          ),
          tags$tr(
            tags$td("lambda"),
            tags$td(model_weibull_compPFS["lambda"]),
            tags$td("TODO")
          ),
          tags$tr(
            tags$td("Intervention")
          ),
          tags$tr(
            tags$td("k"),
            tags$td(model_weibull_PFS["k"]),
            tags$td("TODO")
          ),
          tags$tr(
            tags$td("lambda"),
            tags$td(model_weibull_PFS["lambda"]),
            tags$td("TODO")
          )
        )  
      } else if(input$dist == "Log-logistic"){
        
      } else if(input$dist == "Exponential"){
        
      }
    })
    
    output$modelTableOS <- renderTable({
      
    })
    
    #exploratory data analysis
    #Histograms
    #output$histPlot <- renderPlotly({
    #  df = data()
    #  event = df %>% filter(status == 1)
    #  censored = df %>% filter(status == 0)
    #  fig <- df %>%
    #    plot_ly(x=event$time, type = "histogram", name = "Event", alpha = 0.5) %>%
    #    add_histogram(x = censored$time, name = "Censored") %>%
    #    layout(barmode = "overlay")
    #  return (fig)
    #})
    
    tableControls = tableby.control(test = TRUE, numeric.stats = c("meansd","medianq1q3","range"))
    tableLabels = c(time = "Time")
    
    #descriptive statistics table
    output$descriptiveTablePFS <- renderTable({
      df = data %>% filter(state == "PFS")
      return(
        as.data.frame(
          summary(
            tableby(intervention ~ time, data = df, control = tableControls),
            text = 'html', Translations = tableLabels, pfootnote = TRUE)
            )
          )
    }, sanitize.text.function = identity)
    
    #Descriptive table for the OS variable
    output$descriptiveTableOS <- renderTable({
      df = data %>% filter(state == "OS")
      return(
        as.data.frame(
          summary(
            tableby(intervention ~ time, data = df, control = tableControls),
            text = 'html', labelTranslations = tableLabels, pfootnote = TRUE,
            title = "Overall survival")
            )
          )
    }, sanitize.text.function = identity)
    
    #Changing modes between Descriptive and Models
    observeEvent(input$radioBtnMode, {
      if(input$radioBtnMode == "Descriptive analysis"){
        hide("dist")
      } else if(input$radioBtnMode == "Model analysis"){
        show("dist")  
      }
    })
    
    output$mainPanelContent <- renderUI({
      if(input$radioBtnMode == "Descriptive analysis"){
          ui = tags$div(id = "EDApanel",
            plotOutput("histPlot"),
            if(input$dataset %in% c("All data", "Progression-free survival")){
              "Progression-free survival"
            },
            if(input$dataset %in% c("All data", "Progression-free survival")){
              tableOutput("descriptiveTablePFS")
            },
            if(input$dataset %in% c("All data", "Overall survival")){
              "Overall survival"
            },
            if(input$dataset %in% c("All data", "Overall survival")){
              tableOutput("descriptiveTableOS")
            }
          )
      } else if(input$radioBtnMode == "Model analysis"){
        if(input$dist != ""){
          ui = tags$div(id = "modelPanel",
            if(input$dataset %in% c("All data", "Progression-free survival")){
              plotOutput("modelPlotPFS")
            },
            if(input$dataset %in% c("All data", "Overall survival")){
              plotOutput("modelPlotOS")
            },
            if(input$dataset %in% c("All data", "Progression-free survival")){
              uiOutput("modelTablePFS")
            },
            if(input$dataset %in% c("All data", "Overall survival")){
              tableOutput("modelTableOS")
            }
          )
        }
      }
    })
      #  df = data()
    #  tableControls = tableby.control(test = FALSE,
    #                                  numeric.stats = c("meansd", "medianq1q3", "range"))
    #  return(as.data.frame(summary(tableby(status_label ~time, data = df, 
    #                                       control = tableControls), 
    #                               text = 'html')))
    #}, sanitize.text.function = identity)
    #
    ##at risk table
    #output$atRiskTable <- renderTable({
    #  xtable(atRiskTable(data()))
    #}, rownames = TRUE, colnames = FALSE)
    
    
    #model analysis
    #observeEvent(input$submitBtn,{
    #  hide("exploratoryElements", anim = TRUE, animType =  'fade')
    #  show("modelElements")
    #})
    
    
    
    
    
    
    
    
    #output$modelsPlot <- renderPlot({
    #  df = data()
    #  option = input$dist
    #  n = nrow(df)
    #  quantiles = seq(0.98,0.01,length.out = n)
    #  m_weibull = model_weibull()
    #  pred_weibull = predict(m_weibull, type = 'quantile', p = quantiles)[1,]
    #  df_weibull = data.frame(x = pred_weibull, y = 1-quantiles)
    #  
    #  m_exponential = model_exponential()
    #  pred_exponential = predict(m_exponential, type = 'quantile', p = quantiles)[1,]
    #  df_exponential = data.frame(x = pred_exponential, y = 1-quantiles)
    #  
    #  m_logistic = model_logistic()
    #  pred_logistic = predict(m_logistic, type = 'quantile', p = quantiles)[1,]
    #  df_logistic = data.frame(x = pred_logistic, y = 1-quantiles)
    #  
    #  m_rayleigh = model_rayleigh()
    #  pred_rayleigh = predict(m_rayleigh, type = 'quantile', p = quantiles)[1,]
    #  df_rayleigh = data.frame(x = pred_rayleigh, y = 1-quantiles)
    #  
    #  km_model = survfit(Surv(time,status)~1, data = df) %>% tidy()
    #  if (option == "All distributions"){
    #    fig <- km_model %>%
    #     ggplot(aes(x = time, y = estimate)) +
    #     geom_point()+
    #     geom_line()+
    #     geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2)+
    #     geom_line(data = df_weibull, aes(x,y), color = "red", lwd = 1)+
    #     geom_line(data = df_exponential, aes(x,y), color = "blue", lwd = 1)+
    #     geom_line(data = df_logistic, aes(x,y), color = "green", lwd = 1)+
    #     geom_line(data = df_rayleigh, aes(x,y), color = "brown", lwd = 1)  
    #  } else{
    #    fig <- km_model %>%
    #     ggplot(aes(x = time, y = estimate)) +
    #     geom_point()+
    #     geom_line()+
    #     geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2)
    #  }
    #  
    #  return(fig)
    #})
    
    
})

