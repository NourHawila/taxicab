library(shiny)
library(MCMCpack)
library(ggplot2)
library(gridExtra)
library(partitions)
library(shinyjs)
library(tidyverse)
library(ggpubr)
library(cowplot)
library(stringr)

textInputRow<-function (inputId, label, value = "")
{
  div(style="display:inline-block",
      tags$label(label, `for` = inputId),
      tags$input(id = inputId, type = "text", value = value,class="input-small"))
}

# Define UI for dataset viewer app ----
ui <- fluidPage(# App title ----
                titlePanel("Bayesian Inference with the Taxicab Problem"),

                # Sidebar layout with input and output definitions ----
                sidebarLayout(
                  # Sidebar panel for inputs ----
                  sidebarPanel(width=5,
                               h3("Data"),
                               textInput("x", h5("Enter observed numbers (separated by commas): "), 1729),

                               h3("Bayesian Parameters"),
                               textInput("Nmin",
                                         h5("Enter smallest feasible value for N (Nmin):"),
                                         "50"),
                               textInput("Nmax",
                                         h5("Enter largest feasible value for N (Nmax):"),
                                         "10000"),

                              #selecting prior
                               radioButtons(
                                 "prior",
                                 label = "Select the prior for N",
                                 choiceNames = list(
                                   HTML("<b>Uniform prior</b> (equal probabilities)"),
                                   HTML("<b> Decaying prior</b> (decreasing probabilities)")
                                 ),
                                 choiceValues = c("uniform","decaying"),
                                 selected = "uniform"
                               ),
                               numericInput("alpha", h5("alpha (α):"), 0.2),
                  ),



                  # Main panel for displaying outputs ----
                  mainPanel(width=7,
                            htmlOutput("num.taxicabs"),
                            htmlOutput("max"),
                            htmlOutput("summaries"),
                            h4("Graphing the frequentist and Bayesian estimates"),
                            plotOutput('Plots'),
                  )
                ))



# Define server logic to summarize and view selected dataset ----
server <- function(input, output) {

  prior.ind.uniform=reactive({
    out<-FALSE
    if(is.element("uniform",input$prior)){out<-TRUE}
    return(out)
  })

  prior.ind.decaying=reactive({
    out<-FALSE
    if(is.element("decaying",input$prior)){out<-TRUE}
    return(out)
  })

  #splitting input data
  k <- reactive({
    x = input$x
    outp= strsplit(x,",")
    return(length(outp[[1]]))
  })

  #max observed
  m <- reactive({
    x = input$x
    outp= strsplit(x,",")[[1]]
    return(max(as.numeric(outp),na.rm=T))
  })



  #posterior probability (N|m,k)
  pt3 <- reactive({
    m=m()
    Nmin = as.numeric(input$Nmin)
    Nmax = as.numeric(input$Nmax)
    k=as.numeric(k())



    myprob=rep(0,Nmax)
    for(i in 1:Nmax){
      myprob[i]=((m/i)^k-((m-1)/i)^k)*(m<i)
    }
    likelihood=myprob
    prior=rep(0,Nmax)
    posterior=rep(0,Nmax)

    if(prior.ind.uniform()) {
      for(i in Nmin:Nmax)
      {prior[i]= 1/(Nmax-Nmin+1)}
    }
    if(prior.ind.decaying()) {
      for(i in Nmin:Nmax) {
        prior[i]=1/i
      }
      prior=prior/sum(prior)
    }

    posterior=prior*likelihood
    if(sum(posterior)!=0) {posterior=posterior/sum(posterior)}


    df = data.frame(prior=prior,likelihood=likelihood,posterior=posterior)
    #View(df)

    prior.mean= 0
    for(i in Nmin:Nmax){
      prior.mean = prior.mean + i*prior[i]
    }

    post.mean=(Nmin:Nmax)[which.min(abs(posterior-rep(mean(posterior[posterior>0]),Nmax)))][1]
    post.median=sum((cumsum(posterior)<0.5)==T)+1

    summaries = data.frame(Legend = factor(c("MLE","Unbiased","Prior mean" ,"Posterior mean", "Posterior median"),
                                           levels=c("MLE","Unbiased", "Prior mean","Posterior mean", "Posterior median")),
                                    vals = c(m(),2*m(),prior.mean,post.mean,post.median))

    samp.ub=sum((cumsum(posterior)<1-input$alpha)==T)+1
    samp.ub.plot=sum((cumsum(posterior)<0.95)==T)+1

    gg=ggplot(df)+
      geom_line(aes(x=1:Nmax,y=posterior),size=1)+
      geom_vline(data=summaries,mapping=aes(xintercept=vals, colour=Legend),size=1)+
      theme(text = element_text(size=15)) + xlab("") +
               ylab("Posterior probability") + scale_x_continuous(limits=c(m(), samp.ub.plot))+
      geom_segment(aes(x = m(), xend = samp.ub, y = 0, yend = 0),
                  inherit.aes = FALSE,
                  size = 1.5,
                  color = "blue")+ annotate(geom="text", x=samp.ub, y=0,label=paste0(" ",round((1-input$alpha)*100),"% CI"), color="blue",size=4.5,hjust=0)
    gg
  })



  #Output
  output$num.taxicabs <- renderPrint({
    h4(HTML(paste("Number of taxicab(s): ",(k()),collapse=", ")))
  })

  output$max <- renderPrint({
    h4(HTML(paste("Maximum observed taxicab number: ",(m()),collapse=", ")))
  })

  output$Plots = renderPlot({
    ptlist <- list(pt3())
    grid.arrange(grobs=ptlist,nrow=length(ptlist))
  })

  output$summaries <- renderPrint({
    m=m()
    Nmin = as.numeric(input$Nmin)
    Nmax = as.numeric(input$Nmax)
    k=as.numeric(k())

    myprob=rep(0,Nmax)
    for(i in 1:Nmax){
      myprob[i]=((m/i)^k-((m-1)/i)^k)*(m<i)
    }

    likelihood=myprob
    prior=rep(0,Nmax)
    posterior=rep(0,Nmax)
    if(prior.ind.uniform()) {
      for(i in Nmin:Nmax)
      {prior[i]= 1/(Nmax-Nmin+1)}
    }
    if(prior.ind.decaying()) {
      for(i in Nmin:Nmax) {
        prior[i]=1/i
      }
      prior=prior/sum(prior)
    }

    posterior=prior*likelihood
    if(sum(posterior)!=0) {posterior=posterior/sum(posterior)}
    prior.mean= 0
    for(i in Nmin:Nmax){
      prior.mean = prior.mean + i*prior[i]
    }
    post.mean=(Nmin:Nmax)[which.min(abs(posterior-rep(mean(posterior[posterior>0]),Nmax)))][1]
    post.median=sum((cumsum(posterior)<0.5)==T)+1
    samp.ub=sum((cumsum(posterior)<1-input$alpha)==T)+1

    HTML(paste("Maximum likelihood estimate: ", m(),"<br>","Approx. unbiased estimate: ",round(m()+m()/k(),3),"<br>",
               "Prior mean: ",round(prior.mean),"<br>","Posterior mean: ",round(post.mean),"<br>","Posterior median: ",round(post.median),
               "<br>","Credible interval lower bound: ",m(),"<br>","Credible interval upper bound: ",samp.ub))
  })

}

# Create Shiny app ----
shinyApp(ui, server)
