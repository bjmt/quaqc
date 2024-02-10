#!/usr/bin/env Rscript
#
#   quaqc: QUick Atac-seq Quality Control
#   Copyright (C) 2024  Benjamin Jean-Marie Tremblay
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <https://www.gnu.org/licenses/>.

#-------------------------------------------------------------------------------
# Choose whether to launch app in interactive sessions.

# launchApp <- !interactive()
launchApp <- TRUE

#-------------------------------------------------------------------------------
# Load packages

if (launchApp) message("Loading dependencies ...")

pkgs <- c("shiny"
  , "shinydashboard"
  , "ggplot2"
  , "jsonlite"
  , "DT"
)

checkPkg <- function(x) {
  if (!requireNamespace(x, quietly = TRUE)) x else c()
}

missing.packages <- c()
for (pkg in pkgs) missing.packages <- c(missing.packages, checkPkg(pkg))

if (length(missing.packages)) {
  message("Missing ", length(missing.packages), " packages.")
  answer <- menu(c("Yes", "No"), title = "Do you wish to install missing packages?")
  if (answer == 1) {
    message("Installing packages.")
    install.packages(missing.packages)
  } else {
    message("Please install the following packages:")
    message("    c(", paste0(paste0("\"", missing.packages, "\""), collapse = ", "), ")")
    quit(save = "no")
  }
}

for (pkg in pkgs) suppressPackageStartupMessages(library(pkg, character.only = TRUE))

#-------------------------------------------------------------------------------
# Functions

readQuaqcReport <- function(x) {
  f <- gzfile(x, "rt")
  on.exit(close(f))
  fromJSON(readLines(f), simplifyDataFrame = FALSE)
}

multiReadQuaqcReports <- function(x) {
  lapply(x, function(y) {
    tryCatch(readQuaqcReport(y), error = function(e) NULL)
  })
}

checkReport <- function(x) {
  if (any(!c("quaqc_version", "quaqc_reports", "quaqc_time_end") %in% names(x))) {
    FALSE
  } else {
    TRUE
  }
}

reportToDf <- function(f, x) {
  if (is.null(x) || !checkReport(x)) {
    data.frame(row.names = NULL, check.names = FALSE,
      File = f,
      Status = "BAD INPUT",
      Version = NA,
      Samples = NA,
      Date = NA
    )
  } else {
    data.frame(row.names = NULL, check.names = FALSE,
      File = f,
      Status = "OK",
      Version = x$quaqc_version,
      Samples = length(x$quaqc_reports),
      Date = x$quaqc_time_end
    )
  }
}

#-------------------------------------------------------------------------------
# App UI

ui <- dashboardPage(
  dashboardHeader(title = "quaqc"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Run quaqc / load results", tabName = "load"),
      menuItem("Explore results", tabName = "results"),
      menuItem("Install / configure quaqc", tabName = "config"),
      menuItem("About", tabName = "about")
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "load",
        box(width = NULL, status = "primary", solidHeader = TRUE,
          title = "Load quaqc JSON reports (.json, .json.gz)",
          fileInput("report_files", NULL, multiple = TRUE,
            accept = c(".json", ".json.gz")
          ),
          dataTableOutput("report_files_table")
        ),
        box(width = NULL, status = "primary", solidHeader = TRUE,
          title = "Parameters of loaded JSON report",
          dataTableOutput("report_parameters")
        ),
        box(width = NULL, status = "primary", solidHeader = TRUE,
          title = "Run quaqc"
        )
      ),
      tabItem(tabName = "results"),
      tabItem(tabName = "config"),
      tabItem(tabName = "about")
    )
  )
)

#-------------------------------------------------------------------------------
# App server

server <- function(input, output, session) {

  output$report_files_table <- renderDataTable(datatable({
    req(input$report_files)
    reports <- multiReadQuaqcReports(input$report_files$datapath)
    reportRows <- mapply(function(x, y) reportToDf(x, y),
      input$report_files$name, reports, SIMPLIFY = FALSE)
    reportTable <- do.call(rbind, reportRows)
    reportTable[order(reportTable$Status, decreasing = TRUE), ]
  }, selection = list(mode = "single", selected = 1), rownames = FALSE,
    options = list(dom = "t", columnDefs = list(list(className = "dt-center", targets = "_all"))))
  )

}

#-------------------------------------------------------------------------------
# Launch app

if (launchApp) {
  message("All OK. Launching App.")
  runApp(list(ui = ui, server = server), launch.browser = TRUE, quiet = FALSE)
}

