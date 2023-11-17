## Test yourself SOLUTIONS 
1.  We can modify the different bits of code in the ui and server functions to change the display and functionality of the code. Copy the above functions into this app. Run and play around!
   _Increase the number of repeats,
   change the number of bins,
   and test different inputs._





2.Try different functions. Currently we are using the random "normal" generator. Make one for the Poisson distribution, and two others (see which here). Make sure you change the inputs needed accordingly (ie add more input boxes, change text).

```
ui <- fluidPage( 
  
  titlePanel("Random numbers: chisquare distribution"),   # Title panel
  
  # Inputs
  sidebarLayout(
    sidebarPanel(
      textInput("one", "Mean", value=0),
      textInput("two", "Standard deviation", value=1),
      textInput("three", "Number of repeats", value=10000),
      textInput("caption", "X-axis label"),      
      actionButton("add", "Run"),
      checkboxInput(inputId = "obs", label = strong("Show individual observations"), value = FALSE),
      sliderInput("bins", "Number of bins:",  min = 1, max = 100, value = 5) 
    ),
    
    # Result display
    mainPanel(
      textOutput("check"),
      plotOutput("distPlot"),
      tableOutput("table")
    )
  )
)


server <- function(input,output,session) {
  # Reactive input 
  observeEvent( input$add,{
    x <- as.numeric(input$one)
    y <- as.numeric(input$two)
    title <- as.character(input$caption)
    
    # Reactive expressions
    n <- x+y
    d <- round(abs(rchisq(n)  ) * n ) ### Changed function here!!! 
    m <- mean(d)
    sdd <- sd(d)
    bins <- seq(min(d), max(d), length.out = input$bins + 1)
      
    # Reactive output
    output$check <- renderPrint("Done!")
    output$distPlot <- renderPlot({
      
      # Draw the histogram with the specified number of bins
      hist(d, breaks = bins, freq=F, col = 'darkgray', border = 'white', main=title)
      ld <- density(d)
      lines(ld, col=4) 
      abline(v=m, lty=2, lwd=3, col=4)      
      if (input$obs) { rug(d) }
      
    })
    output$table <- renderTable({
      data = cbind(n,m, sdd)
      return( data ) 
    })
    
  })
)

```  
 
