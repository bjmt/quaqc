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
#

quaqc.R_version <- "0.1"
quaqc.R_year <- "2024"

warning("Note that quaqc.R is not yet complete! Use at your own risk.")

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
      Title = "",
      Status = "BAD INPUT",
      Version = NA,
      Samples = NA,
      Date = NA
    )
  } else {
    data.frame(row.names = NULL, check.names = FALSE,
      File = f,
      Title = x$quaqc_run_title,
      Status = "OK",
      Version = x$quaqc_version,
      Samples = length(x$quaqc_reports),
      Date = x$quaqc_time_end
    )
  }
}

reportToParams <- function(x, i) {
  if (is.null(i)) {
    data.frame(row.names = NULL, check.names = FALSE,
      Flag = c("--mitochondria", "--plastids", "--peaks", "--tss",
        "--target-names", "--target-list", "--blacklist", "--rg-names",
        "--rg-list", "--use-secondary", "--use-nomate", "--use-dups",
        "--use-chimeric", "--use-dovetails", "--no-se", "--mapq",
        "--min-qlen", "--min-flen", "--max-qlen", "--max-flen", "--use-all",
        "--max-depth", "--max-qhist", "--max-fhist", "--tss-size", 
        "--tss-qlen", "--tss-tn5", "--omit-gc", "--omit-depth",
        "--fast", "--lenient", "--nfr", "--nbr", "--footprint", "--chip",
        "--output-dir", "--output-ext", "--no-output", "--json",
        "--keep", "--keep-dir", "--keep-ext", "--threads", "--title",
        "--continue", "--verbose"),
      Value = NA
    )
  } else {
    r <- x[[i]]$quaqc_params
    data.frame(row.names = NULL, check.names = FALSE,
      Flag = c("--mitochondria", "--plastids", "--peaks", "--tss",
        "--target-names", "--target-list", "--blacklist", "--rg-names",
        "--rg-list", "--use-secondary", "--use-nomate", "--use-dups",
        "--use-chimeric", "--use-dovetails", "--no-se", "--mapq",
        "--min-qlen", "--min-flen", "--max-qlen", "--max-flen", "--use-all",
        "--max-depth", "--max-qhist", "--max-fhist", "--tss-size", 
        "--tss-qlen", "--tss-tn5", "--omit-gc", "--omit-depth",
        "--fast", "--lenient", "--nfr", "--nbr", "--footprint", "--chip",
        "--output-dir", "--output-ext", "--no-output", "--json",
        "--keep", "--keep-dir", "--keep-ext", "--threads", "--title",
        "--continue", "--verbose"),
      Value = NA
    )
  }
}

#-------------------------------------------------------------------------------
# App UI

ui <- dashboardPage(
  dashboardHeader(title = "quaqc"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Load results", tabName = "load"),
      menuItem("Explore results", tabName = "results"),
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
        )
      ),
      tabItem(tabName = "results"),
      tabItem(tabName = "about")
    )
  )
)

#-------------------------------------------------------------------------------
# App server

server <- function(input, output, session) {

  reports <- reactiveValues(
    files = NULL, data = NULL, selected = NULL
  )

  output$report_files_table <- renderDataTable(datatable({
    req(input$report_files)
    reports$files <- input$report_files
    reports$data <- multiReadQuaqcReports(input$report_files$datapath)
    reportRows <- mapply(function(x, y) reportToDf(x, y),
      reports$files$name, reports$data, SIMPLIFY = FALSE)
    reportTable <- do.call(rbind, reportRows)
    reportTable[order(reportTable$Status, decreasing = TRUE), ]
  }, selection = list(mode = "single", selected = 1), rownames = FALSE,
    options = list(dom = "t", columnDefs = list(list(className = "dt-center", targets = "_all"))))
  )

  output$report_parameters <- renderDataTable(datatable({
    # req(input$report_files)
    reportToParams(reports$data, input$report_files_table_rows_selected)
  }, rownames = FALSE, selection = "none",
    options = list(ordering = FALSE, dom = "t", pageLength = 999,
      columnDefs = list(list(className = "dt-left", targets = "_all"))))
  )

}

#-------------------------------------------------------------------------------
# Launch app

if (launchApp) {
  message("All OK. Launching App.")
  runApp(list(ui = ui, server = server), launch.browser = TRUE, quiet = FALSE)
}

