# Week 10: So shiny! (Python version)
## Objectives 
Let's make something interactive! 
Shiny is an R/python package that makes it easy to build interactive web applications (apps) straight from R/Python. This lesson will get you started building Shiny apps right away. You will learn how to: 
- Build a user interface
- Add control widgets
- Display reactive output
- Use python scripts and data
- Use reactive expressions
 
## Setting up
You will be creating a shiny app, so there is no need for a notebook this week.  

## What's Shiny, precious, eh?
- Shiny is a package that makes it easy to build interactive web apps.
- Shiny apps are interactive ways to display/showcase or work with your data. See here: 
  
Start off by installing the shiny package: 
```
pip install shiny
```
Then create a package from the terminal. 
```
shiny create  
```



### Reactive == Interactive 
- Reactivity is the ability of a program to compute outputs from user inputs
- The program 'responds' to user defined changes
- Three components: 
1. Reactive Inputs
- Triggers that cause values to be set in the environment 
- actions from buttons, selections, uploads, forms, etc.  
2. Reactive Outputs
- Outputs displayed on the browser
- graphs, plots, tables, etc. 
3. Reactive Expressions
- expressions that take in the reactive input and convert it to a reactive output 

 
### Run!
```
shiny run --reload --launch-browser app_dir/app.py
```

#### Widgets (input)  
- go in the ui function 

function | widget
--- | --- 
actionButton | Action Button
checkboxGroupInput | A group of check boxes
checkboxInput | A single check box
dateInput | A calendar to aid date selection
dateRangeInput | A pair of calendars for selecting a date range
fileInput | A file upload control wizard
helpText | Help text that can be added to an input form
numericInput | A field to enter numbers
radioButtons | A set of radio buttons
selectInput | A box with choices to select from
sliderInput | A slider bar
submitButton | A submit button
textInput | A field to enter text



#### Output functions 
- go in the ui function 


 Output function   |   Creates  
  --- | ---  
 dataTableOutput   |   DataTable 
 htmlOutput   |   raw HTML  
 imageOutput   |   image  
 plotOutput   |   plot 
 ableOutput   |   table  
 textOutput   |   text 
 uiOutput   |   raw HTML  
 verbatimTextOutput   |   text  

#### Render functions 
- go in the server function 


render function   |   Creates
---   |   ---  
renderDataTable   |   DataTable
renderImage   |   images (saved as a link to a source file)
renderPlot   |   plots
renderPrint   |   any printed output
renderTable   |   data frame, matrix, other table like structures
renderText   |   character strings
renderUI   |   a Shiny tag object or HTML


 

## Resources and other 
https://shiny.posit.co/py/docs/overview.html

Back to the [homepage](../README.md)

