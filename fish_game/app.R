#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
# age-structured population dynamics ----


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("魚とりゲーム"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          HTML("<p>ルール</p>"),
          HTML("<ol> <li> １～１４年目まで、まいとし何びきの魚をとるかきめてください
                     <li> 15年目に50匹以上の魚をのこして、14年で合計200匹以上の魚をとったら勝ちです </ol>"),
          HTML("＊魚は、とってしまっても、卵から毎年新しい魚(白丸の魚)が生まれてきます<br><hr>"),  
          conditionalPanel(
                condition = "input.start == false",
                checkboxInput("start", "ゲームスタート")
            ),
          conditionalPanel(
                condition = "input.start == true  & input.year1 == false",
                #condition = "input.start == true  & input.catch1 < 0",
                sliderInput("catch1",
                        "1年目にとる魚の数",
                        min =01,
                        max = 120,
                        value = 0),
                #submitButton("けってい", icon("refresh")),
                checkboxInput("year1", "けってい")
                
            ),
            conditionalPanel(
                condition = "input.year1 == true & input.year2 == false",
                sliderInput("catch2",
                            "2年目にとる魚の数",
                            min = 0,
                            max = 120,
                            value = 0),
                checkboxInput("year2", "けってい")
            ),
            conditionalPanel(
                condition = "input.year2 == true & input.year3 == false",
                sliderInput("catch3",
                            "3年目にとる魚の数",
                            min = 0,
                            max = 120,
                            value = 0),
                checkboxInput("year3", "けってい")
            ),
            conditionalPanel(
                condition = "input.year3 == true & input.year4 == false",
                sliderInput("catch4",
                            "4年目にとる魚の数",
                            min = 0,
                            max = 120,
                            value = 0),
                checkboxInput("year4", "けってい")
            ),
            conditionalPanel(
                condition = "input.year4 == true & input.year5 == false",
                sliderInput("catch5",
                            "5年目にとる魚の数",
                            min = 0,
                            max = 120,
                            value = 0),
                checkboxInput("year5", "けってい")
            ),
            conditionalPanel(
                condition = "input.year5 == true & input.year6 == false",
                sliderInput("catch6",
                            "6年目にとる魚の数",
                            min = 0,
                            max = 120,
                            value = 0),
                checkboxInput("year6", "けってい")
            ),     
            conditionalPanel(
                condition = "input.year6 == true & input.year7 == false",
                sliderInput("catch7",
                            "7年目にとる魚の数",
                            min = 0,
                            max = 120,
                            value = 0),
                checkboxInput("year7", "けってい")
            ),   
            conditionalPanel(
                condition = "input.year7 == true & input.year8 == false",
                sliderInput("catch8",
                            "8年目にとる魚の数",
                            min = 0,
                            max = 120,
                            value = 0),
                checkboxInput("year8", "けってい")
            ),   
            conditionalPanel(
                condition = "input.year8 == true & input.year9 == false",
                sliderInput("catch9",
                            "9年目にとる魚の数",
                            min = 0,
                            max = 120,
                            value = 0),
                checkboxInput("year9", "けってい")
            ),   
            conditionalPanel(
                condition = "input.year9 == true & input.year10 == false",
                sliderInput("catch10",
                            "10年目にとる魚の数",
                            min = 0,
                            max = 120,
                            value = 0),
                checkboxInput("year10", "けってい")
            ),   
            conditionalPanel(
                condition = "input.year10 == true & input.year11 == false",
                sliderInput("catch11",
                            "11年目にとる魚の数",
                            min = 0,
                            max = 120,
                            value = 0),
                checkboxInput("year11", "けってい")
            ),   
            conditionalPanel(
                condition = "input.year11 == true & input.year12 == false",
                sliderInput("catch12",
                            "12年目にとる魚の数",
                            min = 0,
                            max = 120,
                            value = 0),
                checkboxInput("year12", "けってい")
            ),   
            conditionalPanel(
                condition = "input.year12 == true & input.year13 == false",
                sliderInput("catch13",
                            "13年目にとる魚の数",
                            min = 0,
                            max = 120,
                            value = 0),
                checkboxInput("year13", "けってい")
            ),   
            conditionalPanel(
                condition = "input.year13 == true & input.year14 == false",
                sliderInput("catch14",
                            "14年目にとる魚の数",
                            min = 0,
                            max = 120,
                            value = 0),
                checkboxInput("year14", "けってい")
            ),               
            #actionButton("year1", "けってい", icon = NULL, width = NULL)
          # Button
          conditionalPanel(
            condition = "input.year14 == true",
            textInput("name", "おなまえ（できれば入れてください）", "name"),
            downloadButton("data", "結果を保存")
          ),  
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    output$distPlot <- renderPlot({

        catch_number_given <- c(input$catch1,input$catch2,input$catch3,input$catch4,
                                input$catch5,input$catch6,input$catch7,input$catch8,
                                input$catch9,input$catch10,input$catch11,input$catch12,
                                input$catch13,input$catch14)
        year_decision <- c(input$year1,input$year2,input$year3,input$year4,
                           input$year5,input$year6,input$year7,input$year8,
                           input$year9,input$year10,input$year11,input$year12,
                           input$year13,input$year14)

        catch_number_given[year_decision==FALSE] <- NA
        
        res <- pop_dyn(Fc=rep(0,100),R0=R0,M=M,nage=nage,nyear=nyear,init_number=init_number,
                       interactive=TRUE,debug_mode=FALSE,shiny_mode=TRUE,
                       fish_image="fish_isaki.png",xngrid=15,yngrid=8,
                       catch_number_given=catch_number_given)      
        
        if(input$year14==TRUE){
          final_color <- ifelse(sum(catch_number_given)>200 & sum(res$nfish[,15])>50, "yellow", "gray")        
          final_image <- ifelse(sum(catch_number_given)>200 & sum(res$nfish[,15])>50, "win.png", "lose.png")        
          res <- pop_dyn(Fc=rep(0,100),R0=R0,M=M,nage=nage,nyear=nyear,init_number=init_number,
                       interactive=TRUE,debug_mode=FALSE,shiny_mode=TRUE,
                       fish_image="fish_isaki.png",xngrid=15,yngrid=8,
                       catch_number_given=catch_number_given, background_color = final_color)
          res$g1 <- res$g1 +
            geom_image(data=tibble(x=15/2*10,y=8/2*10,img=final_image),
                       aes(x=x, y=y, image=img),asp=3, size=1) 
        }
        gridExtra::grid.arrange(res$g1,res$g2)
    })
  
    # Downloadable csv of selected dataset ----
    output$data <- downloadHandler(
      filename = function() {
        tmp <- str_sub(as.character(Sys.time()),1,16) %>%
          stringr::str_replace(" ", "-") %>%
          stringr::str_replace(":", "-") 
        print(tmp)
        stringr::str_c("fish_game", tmp, ".csv", sep = "")
      },
      content = function(filename) {
        catch_number_given <- c(input$catch1,input$catch2,input$catch3,input$catch4,
                                input$catch5,input$catch6,input$catch7,input$catch8,
                                input$catch9,input$catch10,input$catch11,input$catch12,
                                input$catch13,input$catch14)
        print(catch_number_given)
        print(nyear)
        res <- pop_dyn(Fc=rep(0,100),R0=R0,M=M,nage=nage,nyear=nyear,init_number=init_number,
                       interactive=TRUE,debug_mode=FALSE,shiny_mode=TRUE,
                       fish_image="fish_isaki.png",xngrid=15,yngrid=8,
                       catch_number_given=catch_number_given)
        char_name <- as.character(input$name)
        res <- as_tibble(t(res$nfish)) %>%
          mutate(catch=res$catch_number, name=char_name)
        write.csv(x=res, file=filename)        
      }
    )        

}

library(shiny)
library(png)
library(grid)
library(ggimage)
library(ggthemes)
library(tidyverse)

source("functions.R", encoding="UTF-8", local=FALSE)
shinyApp(ui = ui, server = server)
