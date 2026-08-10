library(shiny)
library(shinyjs)
library(DT)
library(dplyr)
library(shinydashboard)
library(rxode2)
library(mlr3misc)
library(shinybusy)
library(spsComps)
library(CGNM)


js <- '
$(document).on("keyup", function(e) {
  if(e.keyCode == 13){
    Shiny.onInputChange("keyPressed", Math.random());
  }
});
'

callback <- c(
  "var tbl = $(table.table().node());",
  "var id = tbl.closest('.datatables').attr('id');",
  "table.on('autoFill', function(e, datatable, cells){",
  "  var out = Array(cells.length);",
  "  for(var i=0; i<cells.length; ++i){",
  "    var c = cells[i][0];",
  "    var value = c.set === null ? '' : c.set;",
  # null causes problem in R
  "    out[i] = {row: c.index.row+1, col: c.index.column, value: value};",
  # if you want a color for the autofilled cells:
  #"    $(table.cell(c.index.row, c.index.column).node())",
  #"      .css('background-color', 'yellow');",
  "  }",
  "  Shiny.setInputValue(id + '_cells_filled:DT.cellInfo', out);",
  "  table.rows().invalidate();",
  # this updates the column type
  "});"
)


dt_output = function(title, id) {
  fluidRow(column(
    12, h1(paste0('Table ', sub('.*?([0-9]+)$', '\\1', id), ': ', title)),
    hr(), DTOutput(id)
  ))
}
render_dt = function(data, editable = 'cell', server = TRUE, ...) {
  renderDT(data, selection = 'none', server = server, editable = editable, ...)
}


ui <- function(req) {
  dashboardPage(
    dashboardHeader(title = "Shiny CGNM"),
    dashboardSidebar(
      textInput("runNameText", "run name", value = "", width = NULL, placeholder = NULL),
      sidebarMenu(
        menuItem(
          "Step1: specify ODE",
          tabName = "ODE_specify",
          icon = icon("glass-whiskey")
        ),
        menuItem(
          "Step2: give dosing information",
          tabName = "DoseInfo",
          icon = icon("glass-whiskey")
        ),
        menuItem(
          "Step3: specify parameter info",
          tabName = "ParameterInfo",
          icon = icon("glass-whiskey")
        ),
        menuItem(
          "Step4: give Observed data",
          tabName = "ObservedData",
          icon = icon("glass-whiskey")
        ),
        menuItem(
          "Step5: generate code",
          tabName = "generateCODE",
          icon = icon("glass-whiskey")
        ),
        menuItem(
          "Step6: generate simulation code",
          tabName = "generate_simulation_CODE",
          icon = icon("glass-whiskey")
        )
      ),
      wellPanel(
        style = "background: darkgrey",
        h4("Model Variables"),
        textOutput("ODE_lhs"),
        h4("Parameters"),
        textOutput("ODE_para"),
        textOutput("DOSE_para")
      ),
      wellPanel(
        style = "background: grey",
        h4("Load and Save"),
        h6(
          "All GUI inputs can be saved and use later by saving and loading .rds file"
        ),
        fileInput(
          "openSavedStatus",
          label = "Load saved .rds file",
          placeholder = ".rds file"
        ),
        downloadButton(
          "saveCurrentStatus",
          label = "Save model as .rds file",
          icon = shiny::icon("download"),
          style = "color: black; background-color: white; border-color: grey"
        )
      ),
      h6(
        "To see an example model implemented using this GUI you can press the button below. (PBPK model for the DDI between pitavastatin and rifampicin, as published in Yoshikado et al.)"
      ),
      actionButton("loadExample_button", "Load example model")



    ),
    dashboardBody(tabItems(
      tabItem(
        tabName = "ODE_specify",
        add_busy_spinner(spin = "fading-circle"),
        tags$script(js),
        textAreaInput(
          "ODE_text",
          "Specify the system of ODEs (follow rxODE syntax)",
          width = "100%",
          height = "500",
          resize = "vertical"
        ),
        actionButton("compileODE", "Compile ODEs"),
        htmlOutput("ODE_compileMessage"),
        h5("example ode:"),
        h6("d/dt(A_admin)=-ka*A_admin"),
        h6("d/dt(C_central)=1/V_central*(ka*A_admin-CL*C_central-Q*C_central+Q*C_pref)"),
        h6("d/dt(C_pref)=1/V_pref*(Q*C_central-Q*C_pref)"),
      ),

      tabItem(
        tabName = "DoseInfo",
        add_busy_spinner(spin = "fading-circle"),
        numericInput(
          "numDoses",
          "Number of dose events",
          1,
          min = 1,
          step = 1
        ),
        #  dt_output('client-side processing (editable = "cell")', 'dose_dt')
        DTOutput("dose_dt"),
        verbatimTextOutput("table"),
        fileInput("dose_csv", "You can import dose data from a csv file by uploarding it below:"),
        h6("Dose csv file should contain the following columns: ID, dose, dosing.to, start.time and rate."),
        downloadButton('download_dose_csv',"Download dosing information as csv file"),


      ),

      tabItem(
        tabName = "ParameterInfo",
        add_busy_spinner(spin = "fading-circle"),
        DTOutput("parameterInfo_dt"),
        h4(
          "Lower_bound and Upper_bound are optional so can be left empty.  All table entries needs to be numerical values."
        ),
        h4(
          "VaryByID=0: all individials share the same parameter, VaryByID=1: each individiaul has its own parameter."
        ),
        fileInput("para_csv", "You can import parameter information from a csv file by uploarding it below:"),
        h6("Parameter information csv file should contain the following columns: ParameterName, Initial_lower_range, Initial_upper_range, Lower_bound, Upper_bound, and Unit. Only the rows with relevant ParameterName to the current model will be read in."),
        downloadButton('download_para_csv',"Download parameter information as csv file"),
        #   verbatimTextOutput("parameterInfo_table"),

      ),

      tabItem(
        tabName = "ObservedData",
        add_busy_spinner(spin = "fading-circle"),
        numericInput(
          "numObservations",
          "Number of observations",
          1,
          min = 1,
          step = 1
        ),
        DTOutput("ObservedData_dt"),
        h4("ResidualError_model: 0 additive model (default), 1 proportional model."),
        fileInput("obs_csv", "You can import observation data from a csv file by uploarding it below:"),
        h6("Observation csv file should contain the following columns: ID, time, Observation_expression, Observed_value, and Memo."),
        downloadButton('download_obs_csv',"Download observation as csv file"),


        #   verbatimTextOutput("ObservedData_table"),

      ),

      tabItem(
        tabName = "generateCODE",
        add_busy_spinner(spin = "fading-circle"),
        shinyjs::useShinyjs(),
        h3("Conduct test according to the order below."),
        h5("If preceeding test fails most likely the following test will crash the GUI, so if the error is detected fix it before proceeding to the next test."),
        h4("Test 1: check the code can be put together, also check if the code make sense"),
        h6("tip: focusing to see if any unexpected NA values used in the code or not"),
        actionButton("makeCode_button", "Make model code"),
        verbatimTextOutput("test1_output"),
        h4("Test 2: check if your model can be evaluated, if error occurs, fix model code."),
        h6("tip: check online resources on rxode2. If model cannot be evaluated it is a problem independent of CGNM."),
        actionButton("testModelCode_button", "Test model code"),
        verbatimTextOutput("testModelCode_output"),
        h4("Test 3: check if CGNM runs, if error occurs, check all the observation values and initial lower and upper bounds of the parameters."),
        actionButton("testCGNMCode_button", "Test CGNM code"),
        h6("tip: if model evaluation test (Test 2) runs without an error and cannot spot error for CGNM inputs yourself, you may contact the author of CGNM with rds file and generated R-script which can be obtained by pressing the button below."),
        htmlOutput("testCode_output"),
        verbatimTextOutput("testCode_output2"),
        downloadButton("download_CGNMcode_button", label = "Download R script to run CGNM"),
        downloadButton("download_CGNMcode_parallel_mac_button", label = "Download R script to run CGNM in parallel in MAC"),
        downloadButton("download_CGNMcode_parallel_win_button", label = "Download R script to run CGNM in parallel in Windows"),
        downloadButton("download_posthocMiddleout_button", label = "Download R script for posthoc middle out"),

      ),

      tabItem(
        tabName = "generate_simulation_CODE",
        add_busy_spinner(spin = "fading-circle"),
        shinyjs::useShinyjs(),
        h3("Provide dosing information for simulation"),
        numericInput(
          "sim_numDoses",
          "Number of dose events to simulate",
          1,
          min = 1,
          step = 1
        ),
        DTOutput("sim_dose_dt"),
        actionButton("copySimDoseFromDose_button", "Copy original dose events to simulation dose events"),
        downloadButton("download_sim_dose_csv", label = "Download simulation Dose as .csv"),
        fileInput("sim_dose_csv", "You can import simulation dose from a csv file by uploarding it below:"),
        h6("Dose csv file should contain the following columns: ID, dose, dosing.to, start.time and rate."),
        h3("Provide data skelton for simulation"),
        numericInput(
          "sim_numObservations",
          "Number of observations",
          1,
          min = 1,
          step = 1
        ),
        DTOutput("sim_ObservedData_dt"),
        actionButton("copySimObsFromObs_button", "Copy original observation time points to simulation time points"),
        downloadButton('download_sim_obs_csv',"Download simulation timepoints as csv file"),
        fileInput("sim_obs_csv", "You can import simulation timepoints from a csv file by uploarding it below:"),
        h6("Simulation time point csv file should contain the following columns: ID, time, Observation_expression, and Memo."),
        h3("Simulation code"),
        actionButton("sim_makeCode_button", "Generate simulation code"),
        verbatimTextOutput("simulationCode_output"),
        downloadButton("download_simulationCode_button", label = "Download R script to run Simulation")

      )

    ))
  )
}


#shinyApp(ui = ui,
#         server = server,
#         enableBookmarking = "server")
