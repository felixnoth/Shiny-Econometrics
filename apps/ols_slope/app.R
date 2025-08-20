#| echo: false
library(shiny)
library(ggplot2)
library(MASS)
ui <- fluidPage(
  titlePanel("Effect of Covariance and Variance on Linear Regression"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("var_x", "Variance of X:", min = 0.1, max = 10, value = 1, step = 0.1),
      sliderInput("cov_xy", "Covariance between X and Y:", min = -10, max = 10, value = 1, step = 0.1),
      numericInput("var_y", "Variance of Y (fixed for now):", value = 1, min = 0.1),
      numericInput("n", "Number of Data Points:", value = 100, min = 10)
    ),
    mainPanel(
      plotOutput("regPlot"),
      verbatimTextOutput("regSummary")
    )
  )
)
server <- function(input, output, session) {
  data <- reactive({
    Sigma <- matrix(c(input$var_x, input$cov_xy,
                      input$cov_xy, input$var_y), nrow = 2)
    set.seed(123)
    dat <- as.data.frame(mvrnorm(n = input$n, mu = c(0, 0), Sigma = Sigma))
    colnames(dat) <- c("x", "y")
    return(dat)
  })
  
  output$regPlot <- renderPlot({
    dat <- data()
    model <- lm(y ~ x, data = dat)
    ggplot(dat, aes(x = x, y = y)) +
      geom_point(color = "#4682B4", alpha = 0.6) +
      geom_smooth(method = "lm", se = FALSE, color = "#B22222", size = 1.2) +
      labs(title = "Linear Regression Line",
           subtitle = paste("Slope:", round(coef(model)[2], 3)),
           x = "X", y = "Y") +
      theme_minimal()
  })
}

shinyApp(ui = ui, server = server)