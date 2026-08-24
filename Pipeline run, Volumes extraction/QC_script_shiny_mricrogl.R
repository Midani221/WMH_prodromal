# install.packages("shiny")

library(shiny)

# =========================
# 1. EDIT THESE PATHS
# =========================

ROOT <- "Folder path"   # parent folder containing sub-xxx_ses-BL folders
# Parent folder containing sub-xxx_ses-BL folders

MRICROGL_EXE <- "Path/MRIcroGL.exe"
# Change this to the real MRIcroGL.exe path on your computer

QC_CSV <- file.path(ROOT, "wmh_qc_results.csv")
SCRIPT_DIR <- file.path(ROOT, "_mricrogl_scripts")

FLAIR_REGEX <- "^FLAIR_biascorr_brain_to_MNI_nonlin\\.nii(\\.gz)?$"
MASK_REGEX  <- "^results2mni_nonlin_combined\\.nii(\\.gz)?$"


# =========================
# 2. FILE DISCOVERY
# =========================

if (!dir.exists(ROOT)) {
  stop("ROOT folder does not exist. Please check ROOT.")
}

if (!dir.exists(SCRIPT_DIR)) {
  dir.create(SCRIPT_DIR, recursive = TRUE)
}

find_case_files <- function(folder) {
  nii_files <- list.files(
    folder,
    pattern = "\\.nii(\\.gz)?$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  
  flair_candidates <- nii_files[
    grepl(FLAIR_REGEX, basename(nii_files), ignore.case = TRUE)
  ]
  
  mask_candidates <- nii_files[
    grepl(MASK_REGEX, basename(nii_files), ignore.case = TRUE)
  ]
  
  data.frame(
    case = basename(folder),
    folder = folder,
    flair = if (length(flair_candidates) > 0) flair_candidates[1] else NA_character_,
    mask  = if (length(mask_candidates) > 0) mask_candidates[1] else NA_character_,
    stringsAsFactors = FALSE
  )
}

case_dirs <- list.dirs(ROOT, recursive = FALSE, full.names = TRUE)

if (length(case_dirs) == 0) {
  stop("No subfolders found inside ROOT.")
}

cases <- do.call(rbind, lapply(case_dirs, find_case_files))
cases <- cases[!is.na(cases$flair) & !is.na(cases$mask), ]
cases <- cases[order(cases$case), ]

if (nrow(cases) == 0) {
  stop("No valid cases found. Check ROOT, FLAIR_REGEX, and MASK_REGEX.")
}

cat("Found", nrow(cases), "valid cases.\n")
print(head(cases[, c("case", "flair", "mask")], 5))


# =========================
# 3. QC HELPERS
# =========================

load_existing_qc <- function() {
  if (file.exists(QC_CSV)) {
    read.csv(QC_CSV, stringsAsFactors = FALSE)
  } else {
    data.frame()
  }
}

save_qc_row <- function(row) {
  old <- load_existing_qc()
  
  if (nrow(old) == 0 && ncol(old) == 0) {
    write.csv(row, QC_CSV, row.names = FALSE)
    return(invisible(TRUE))
  }
  
  if (nrow(old) > 0 && "case" %in% names(old)) {
    old <- old[old$case != row$case, ]
  }
  
  all_cols <- union(names(old), names(row))
  
  for (nm in setdiff(all_cols, names(old))) {
    old[[nm]] <- NA
  }
  
  for (nm in setdiff(all_cols, names(row))) {
    row[[nm]] <- NA
  }
  
  old <- old[, all_cols, drop = FALSE]
  row <- row[, all_cols, drop = FALSE]
  
  out <- rbind(old, row)
  write.csv(out, QC_CSV, row.names = FALSE)
  
  invisible(TRUE)
}

safe_value <- function(x) {
  if (length(x) == 0 || is.na(x)) "" else x
}


# =========================
# 4. MRICROGL SCRIPT HELPERS
# =========================

py_quote <- function(x) {
  x <- gsub("\\\\", "/", x)
  x <- gsub('"', '\\"', x, fixed = TRUE)
  paste0('"', x, '"')
}

kill_existing_mricrogl <- function() {
  if (.Platform$OS.type == "windows") {
    for (proc in c("MRIcroGL.exe", "mricrogl.exe")) {
      suppressWarnings(
        try(
          system2(
            "taskkill",
            args = c("/IM", proc, "/F"),
            stdout = FALSE,
            stderr = FALSE
          ),
          silent = TRUE
        )
      )
    }
    Sys.sleep(0.4)
  }
}

make_mricrogl_script <- function(case_row,
                                 script_file,
                                 mask_min = 0.5,
                                 mask_max = 1,
                                 opacity = 50,
                                 overlay_color = "red") {
  
  flair <- normalizePath(case_row$flair, winslash = "/", mustWork = TRUE)
  mask  <- normalizePath(case_row$mask,  winslash = "/", mustWork = TRUE)
  
  lines <- c(
    "import gl",
    "",
    "gl.resetdefaults()",
    "gl.backcolor(0, 0, 0)",
    "",
    "# Load FLAIR background",
    sprintf("gl.loadimage(%s)", py_quote(flair)),
    "",
    "# Load binary WMH mask as jagged/non-smoothed overlay",
    "gl.overlayloadsmooth(False)",
    sprintf("gl.overlayload(%s)", py_quote(mask)),
    sprintf("gl.minmax(1, %.4f, %.4f)", mask_min, mask_max),
    sprintf("gl.colorname(1, %s)", py_quote(overlay_color)),
    sprintf("gl.opacity(1, %d)", as.integer(opacity)),
    "",
    "# Start near MNI origin; you can move crosshairs freely in MRIcroGL",
    "gl.orthoviewmm(0, 0, 0)"
  )
  
  writeLines(lines, script_file)
  script_file
}


# =========================
# 5. SHINY UI
# =========================

ui <- fluidPage(
  titlePanel("WMH QC Dashboard + MRIcroGL Viewer"),
  
  sidebarLayout(
    sidebarPanel(
      width = 4,
      
      selectInput(
        "case_idx",
        "Subject/session",
        choices = cases$case,
        selected = cases$case[1]
      ),
      
      fluidRow(
        column(6, actionButton("prev_case", "Previous", width = "100%")),
        column(6, actionButton("next_case", "Next", width = "100%"))
      ),
      
      br(),
      
      actionButton(
        "open_mricrogl",
        "Open current case in MRIcroGL",
        width = "100%",
        class = "btn-primary"
      ),
      
      br(), br(),
      
      actionButton(
        "save_next_open",
        "Save, next, and open MRIcroGL",
        width = "100%",
        class = "btn-success"
      ),
      
      hr(),
      
      checkboxInput(
        "close_existing",
        "Close previous MRIcroGL window before opening next",
        value = TRUE
      ),
      
      numericInput(
        "mask_min",
        "Mask lower threshold",
        value = 0.5,
        min = 0,
        max = 10,
        step = 0.1
      ),
      
      numericInput(
        "mask_max",
        "Mask upper threshold",
        value = 1,
        min = 0,
        max = 10,
        step = 0.1
      ),
      
      sliderInput(
        "overlay_opacity",
        "Mask opacity",
        min = 0,
        max = 100,
        value = 50,
        step = 5
      ),
      
      selectInput(
        "overlay_color",
        "Mask color",
        choices = c("red", "hot", "blue", "green", "actc"),
        selected = "red"
      ),
      
      hr(),
      
      radioButtons(
        "quality",
        "Quality",
        choices = c("Good", "Maybe", "Poor"),
        selected = character(0),
        inline = TRUE
      ),
      
      checkboxGroupInput(
        "reasons",
        "Reasons",
        choices = c(
          "Missed dWMHs",
          "Missed pWMHs",
          "Poor deskulling",
          "Poor brain selection",
          "Motion artefact",
          "Extra neuroimaging abnormality",
          "Other"
        )
      ),
      
      radioButtons(
        "estimation",
        "Estimation",
        choices = c("Overestimation", "Decent", "Underestimation"),
        selected = character(0)
      ),
      
      textAreaInput(
        "comments",
        "Comments",
        value = "",
        rows = 4
      ),
      
      fluidRow(
        column(6, actionButton("save", "Save QC", width = "100%")),
        column(6, actionButton("save_next", "Save and next", width = "100%"))
      )
    ),
    
    mainPanel(
      width = 8,
      
      h4("Current case"),
      verbatimTextOutput("case_info"),
      
      hr(),
      
      h4("Workflow"),
      p("1. Select a subject/session."),
      p("2. Click 'Open current case in MRIcroGL'."),
      p("3. Review FLAIR + mask overlay smoothly in MRIcroGL."),
      p("4. Return here, fill QC, then use 'Save, next, and open MRIcroGL'."),
      
      hr(),
      
      h4("MRIcroGL launch status"),
      verbatimTextOutput("launch_status")
    )
  )
)


# =========================
# 6. SHINY SERVER
# =========================

server <- function(input, output, session) {
  
  launch_message <- reactiveVal("No case launched yet.")
  
  selected_index <- reactive({
    req(input$case_idx)
    i <- match(input$case_idx, cases$case)
    validate(need(!is.na(i), "Selected case not found."))
    i
  })
  
  selected_case <- reactive({
    cases[selected_index(), ]
  })
  
  load_qc_for_case <- function(case_name) {
    qc <- load_existing_qc()
    
    if (nrow(qc) > 0 && "case" %in% names(qc)) {
      existing <- qc[qc$case == case_name, ]
    } else {
      existing <- data.frame()
    }
    
    if (nrow(existing) > 0) {
      q <- safe_value(existing$quality[1])
      e <- safe_value(existing$estimation[1])
      cmt <- safe_value(existing$comments[1])
      r <- safe_value(existing$reasons[1])
      
      updateRadioButtons(
        session,
        "quality",
        selected = if (q == "") character(0) else q
      )
      
      if (r != "") {
        existing_reasons <- unlist(strsplit(r, ";\\s*"))
        existing_reasons <- existing_reasons[existing_reasons != ""]
      } else {
        existing_reasons <- character(0)
      }
      
      updateCheckboxGroupInput(
        session,
        "reasons",
        selected = existing_reasons
      )
      
      updateRadioButtons(
        session,
        "estimation",
        selected = if (e == "") character(0) else e
      )
      
      updateTextAreaInput(session, "comments", value = cmt)
      
    } else {
      updateRadioButtons(session, "quality", selected = character(0))
      updateCheckboxGroupInput(session, "reasons", selected = character(0))
      updateRadioButtons(session, "estimation", selected = character(0))
      updateTextAreaInput(session, "comments", value = "")
    }
  }
  
  observeEvent(input$case_idx, {
    load_qc_for_case(input$case_idx)
  }, ignoreInit = FALSE)
  
  output$case_info <- renderText({
    row <- selected_case()
    
    paste0(
      "Case: ", row$case, "\n",
      "Index: ", selected_index(), " / ", nrow(cases), "\n\n",
      "FLAIR:\n", row$flair, "\n\n",
      "Mask:\n", row$mask, "\n\n",
      "QC CSV:\n", QC_CSV
    )
  })
  
  output$launch_status <- renderText({
    launch_message()
  })
  
  open_case_in_mricrogl <- function(i) {
    if (!file.exists(MRICROGL_EXE)) {
      msg <- paste0(
        "MRIcroGL executable not found:\n",
        MRICROGL_EXE,
        "\n\nFix MRICROGL_EXE at the top of the script."
      )
      launch_message(msg)
      showNotification("MRIcroGL executable not found. Check MRICROGL_EXE.", type = "error")
      return(FALSE)
    }
    
    row <- cases[i, ]
    
    script_file <- file.path(
      SCRIPT_DIR,
      paste0(gsub("[^A-Za-z0-9_\\-]", "_", row$case), "_mricrogl.py")
    )
    
    make_mricrogl_script(
      case_row = row,
      script_file = script_file,
      mask_min = input$mask_min,
      mask_max = input$mask_max,
      opacity = input$overlay_opacity,
      overlay_color = input$overlay_color
    )
    
    if (isTRUE(input$close_existing)) {
      kill_existing_mricrogl()
    }
    
    mricrogl <- normalizePath(MRICROGL_EXE, winslash = "/", mustWork = TRUE)
    script_norm <- normalizePath(script_file, winslash = "/", mustWork = TRUE)
    
    system2(
      mricrogl,
      args = shQuote(script_norm),
      wait = FALSE
    )
    
    msg <- paste0(
      "Opened in MRIcroGL:\n",
      row$case,
      "\n\nScript:\n",
      script_norm
    )
    
    launch_message(msg)
    showNotification(paste("Opened", row$case, "in MRIcroGL"), type = "message")
    
    TRUE
  }
  
  save_current_qc <- function() {
    row0 <- selected_case()
    
    row <- data.frame(
      case = row0$case,
      folder = row0$folder,
      flair = row0$flair,
      mask = row0$mask,
      quality = if (length(input$quality) == 0) "" else input$quality,
      reasons = paste(input$reasons, collapse = "; "),
      estimation = if (length(input$estimation) == 0) "" else input$estimation,
      comments = input$comments,
      mask_min = input$mask_min,
      mask_max = input$mask_max,
      overlay_opacity = input$overlay_opacity,
      overlay_color = input$overlay_color,
      qc_time = as.character(Sys.time()),
      stringsAsFactors = FALSE
    )
    
    save_qc_row(row)
    showNotification(paste("Saved QC for", row0$case), type = "message")
    
    TRUE
  }
  
  observeEvent(input$open_mricrogl, {
    open_case_in_mricrogl(selected_index())
  })
  
  observeEvent(input$prev_case, {
    i <- selected_index()
    if (i > 1) {
      updateSelectInput(session, "case_idx", selected = cases$case[i - 1])
    }
  })
  
  observeEvent(input$next_case, {
    i <- selected_index()
    if (i < nrow(cases)) {
      updateSelectInput(session, "case_idx", selected = cases$case[i + 1])
    }
  })
  
  observeEvent(input$save, {
    save_current_qc()
  })
  
  observeEvent(input$save_next, {
    save_current_qc()
    
    i <- selected_index()
    if (i < nrow(cases)) {
      updateSelectInput(session, "case_idx", selected = cases$case[i + 1])
    }
  })
  
  observeEvent(input$save_next_open, {
    save_current_qc()
    
    i <- selected_index()
    
    if (i < nrow(cases)) {
      next_i <- i + 1
      updateSelectInput(session, "case_idx", selected = cases$case[next_i])
      open_case_in_mricrogl(next_i)
    } else {
      showNotification("Already at last case.", type = "warning")
    }
  })
}

shinyApp(ui, server)