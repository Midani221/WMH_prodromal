# install.packages("shiny")

library(shiny)

# =========================
# 1. EDIT THESE PATHS
# =========================

ROOT <- "Path for folder"   # parent folder containing sub-xxx_ses-BL folders
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
MAX_OVERLAYS <- 3


# =========================
# 2. FILE DISCOVERY
# =========================

if (!dir.exists(ROOT)) {
  stop("ROOT folder does not exist. Please check ROOT.")
}

if (!dir.exists(SCRIPT_DIR)) {
  dir.create(SCRIPT_DIR, recursive = TRUE)
}

parse_subject_session <- function(folder_name) {
  # Expected format: sub-xxxxxx_ses-Vxx
  m <- regexec("^(sub-[^_]+)_ses-([^_]+)$", folder_name, ignore.case = TRUE)
  hit <- regmatches(folder_name, m)[[1]]
  
  if (length(hit) == 3) {
    return(list(subject = hit[2], session = hit[3]))
  } else {
    return(list(subject = NA_character_, session = NA_character_))
  }
}

session_rank <- function(session) {
  s <- toupper(session)
  
  if (s %in% c("BL", "BASELINE")) {
    return(0)
  }
  
  # V06 -> 6, V10 -> 10
  num <- suppressWarnings(as.numeric(gsub("[^0-9]", "", s)))
  
  if (!is.na(num)) {
    return(num)
  }
  
  return(9999)
}

find_visit_files <- function(folder) {
  folder_name <- basename(folder)
  parsed <- parse_subject_session(folder_name)
  
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
    folder_name = folder_name,
    subject = parsed$subject,
    session = parsed$session,
    session_rank = session_rank(parsed$session),
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

visits <- do.call(rbind, lapply(case_dirs, find_visit_files))

# Keep only folders that match sub-xxx_ses-yyy and have both FLAIR and mask
visits <- visits[
  !is.na(visits$subject) &
    !is.na(visits$session) &
    !is.na(visits$flair) &
    !is.na(visits$mask),
]

if (nrow(visits) == 0) {
  stop("No valid visit folders found. Check ROOT, FLAIR_REGEX, and MASK_REGEX.")
}

make_subject_case <- function(df) {
  df <- df[order(df$session_rank, df$session), ]
  
  # Need at least 2 sessions for this workflow
  if (nrow(df) < 2) {
    return(NULL)
  }
  
  n_use <- if (is.infinite(MAX_OVERLAYS)) {
    nrow(df)
  } else {
    min(nrow(df), MAX_OVERLAYS)
  }
  
  use <- df[seq_len(n_use), ]
  baseline <- use[1, ]
  
  data.frame(
    case = baseline$subject,
    subject = baseline$subject,
    baseline_session = baseline$session,
    baseline_flair = baseline$flair,
    overlay_sessions = paste(use$session, collapse = "; "),
    overlay_masks = paste(use$mask, collapse = "|"),
    overlay_folders = paste(use$folder, collapse = "|"),
    n_sessions_available = nrow(df),
    n_overlays = nrow(use),
    all_sessions = paste(df$session, collapse = "; "),
    stringsAsFactors = FALSE
  )
}

cases_list <- lapply(split(visits, visits$subject), make_subject_case)
cases_list <- Filter(Negate(is.null), cases_list)

if (length(cases_list) == 0) {
  stop("No subjects with at least 2 valid sessions found.")
}

cases <- do.call(rbind, cases_list)
cases <- cases[order(cases$subject), ]

cat("Found", nrow(cases), "subjects with at least 2 sessions.\n")
print(head(cases[, c("subject", "baseline_session", "overlay_sessions")], 10))


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
  
  if (nrow(old) > 0 && "subject" %in% names(old)) {
    old <- old[old$subject != row$subject, ]
  } else if (nrow(old) > 0 && "case" %in% names(old)) {
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
  if (length(x) == 0 || is.na(x[1])) "" else as.character(x[1])
}

split_masks <- function(x) {
  if (length(x) == 0 || is.na(x) || x == "") {
    character(0)
  } else {
    unlist(strsplit(x, "\\|"))
  }
}

split_sessions <- function(x) {
  if (length(x) == 0 || is.na(x) || x == "") {
    character(0)
  } else {
    unlist(strsplit(x, ";\\s*"))
  }
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
                                 overlay_color_1 = "red",
                                 overlay_color_2 = "blue") {
  
  flair <- normalizePath(case_row$baseline_flair, winslash = "/", mustWork = TRUE)
  masks <- split_masks(case_row$overlay_masks)
  sessions <- split_sessions(case_row$overlay_sessions)
  
  if (length(masks) == 0) {
    stop("No overlay masks found for this subject.")
  }
  
  color_pool <- c(
    overlay_color_1,
    overlay_color_2,
    "green",
    "hot",
    "actc",
    "blue",
    "red"
  )
  
  colors <- rep(color_pool, length.out = length(masks))
  
  lines <- c(
    "import gl",
    "",
    "gl.resetdefaults()",
    "gl.backcolor(0, 0, 0)",
    "",
    "# Load baseline FLAIR background",
    sprintf("gl.loadimage(%s)", py_quote(flair)),
    "",
    "# Load WMH masks as jagged/non-smoothed overlays",
    "gl.overlayloadsmooth(False)",
    ""
  )
  
  for (k in seq_along(masks)) {
    mask_path <- normalizePath(masks[k], winslash = "/", mustWork = TRUE)
    sess_label <- if (k <= length(sessions)) sessions[k] else paste0("overlay_", k)
    
    lines <- c(
      lines,
      sprintf("# Overlay %d: %s", k, sess_label),
      sprintf("gl.overlayload(%s)", py_quote(mask_path)),
      sprintf("gl.minmax(%d, %.4f, %.4f)", k, mask_min, mask_max),
      sprintf("gl.colorname(%d, %s)", k, py_quote(colors[k])),
      sprintf("gl.opacity(%d, %d)", k, as.integer(opacity)),
      ""
    )
  }
  
  lines <- c(
    lines,
    "# Start near MNI origin; you can move crosshairs freely in MRIcroGL",
    "gl.orthoviewmm(0, 0, 0)"
  )
  
  writeLines(lines, script_file)
  script_file
}


# =========================
# 5. SHINY UI
# =========================

case_choices <- setNames(
  as.character(seq_len(nrow(cases))),
  paste0(
    cases$subject,
    " | baseline FLAIR: ",
    cases$baseline_session,
    " | overlays: ",
    cases$overlay_sessions
  )
)

ui <- fluidPage(
  titlePanel("WMH QC Dashboard: Baseline FLAIR + Two Visit Overlays"),
  
  sidebarLayout(
    sidebarPanel(
      width = 4,
      
      selectInput(
        "case_idx",
        "Subject",
        choices = case_choices,
        selected = "1"
      ),
      
      fluidRow(
        column(6, actionButton("prev_case", "Previous", width = "100%")),
        column(6, actionButton("next_case", "Next", width = "100%"))
      ),
      
      br(),
      
      actionButton(
        "open_mricrogl",
        "Open current subject in MRIcroGL",
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
        "overlay_color_1",
        "Visit 1 mask color",
        choices = c("red", "hot", "blue", "green", "actc"),
        selected = "red"
      ),
      
      selectInput(
        "overlay_color_2",
        "Visit 2 mask color",
        choices = c("blue", "green", "hot", "red", "actc"),
        selected = "blue"
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
          "Longitudinal mismatch",
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
      
      h4("Current subject"),
      verbatimTextOutput("case_info"),
      
      hr(),
      
      h4("Workflow"),
      p("1. Select a subject."),
      p("2. Click 'Open current subject in MRIcroGL'."),
      p("3. MRIcroGL will load the earliest session FLAIR as the background."),
      p("4. The first two session masks will be loaded as overlays with different colors."),
      p("5. Return here, fill QC, then use 'Save, next, and open MRIcroGL'."),
      
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
  
  launch_message <- reactiveVal("No subject launched yet.")
  
  selected_index <- reactive({
    req(input$case_idx)
    i <- as.integer(input$case_idx)
    validate(need(!is.na(i), "Selected subject not found."))
    i
  })
  
  selected_case <- reactive({
    cases[selected_index(), ]
  })
  
  load_qc_for_subject <- function(subject_name) {
    qc <- load_existing_qc()
    
    if (nrow(qc) > 0 && "subject" %in% names(qc)) {
      existing <- qc[qc$subject == subject_name, ]
    } else if (nrow(qc) > 0 && "case" %in% names(qc)) {
      existing <- qc[qc$case == subject_name, ]
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
    row <- selected_case()
    load_qc_for_subject(row$subject)
  }, ignoreInit = FALSE)
  
  output$case_info <- renderText({
    row <- selected_case()
    masks <- split_masks(row$overlay_masks)
    sessions <- split_sessions(row$overlay_sessions)
    
    overlay_text <- paste0(
      paste0(
        "Overlay ", seq_along(masks), " — session ", sessions, ":\n",
        masks
      ),
      collapse = "\n\n"
    )
    
    paste0(
      "Subject: ", row$subject, "\n",
      "Index: ", selected_index(), " / ", nrow(cases), "\n\n",
      "All available sessions:\n", row$all_sessions, "\n\n",
      "Baseline FLAIR session:\n", row$baseline_session, "\n\n",
      "Baseline FLAIR:\n", row$baseline_flair, "\n\n",
      "Overlay masks:\n", overlay_text, "\n\n",
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
      paste0(gsub("[^A-Za-z0-9_\\-]", "_", row$subject), "_mricrogl.py")
    )
    
    make_mricrogl_script(
      case_row = row,
      script_file = script_file,
      mask_min = input$mask_min,
      mask_max = input$mask_max,
      opacity = input$overlay_opacity,
      overlay_color_1 = input$overlay_color_1,
      overlay_color_2 = input$overlay_color_2
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
      row$subject,
      "\n\nBaseline FLAIR session:\n",
      row$baseline_session,
      "\n\nOverlay sessions:\n",
      row$overlay_sessions,
      "\n\nScript:\n",
      script_norm
    )
    
    launch_message(msg)
    showNotification(paste("Opened", row$subject, "in MRIcroGL"), type = "message")
    
    TRUE
  }
  
  save_current_qc <- function() {
    row0 <- selected_case()
    
    row <- data.frame(
      case = row0$case,
      subject = row0$subject,
      baseline_session = row0$baseline_session,
      baseline_flair = row0$baseline_flair,
      overlay_sessions = row0$overlay_sessions,
      overlay_masks = row0$overlay_masks,
      n_sessions_available = row0$n_sessions_available,
      n_overlays = row0$n_overlays,
      all_sessions = row0$all_sessions,
      quality = if (length(input$quality) == 0) "" else input$quality,
      reasons = paste(input$reasons, collapse = "; "),
      estimation = if (length(input$estimation) == 0) "" else input$estimation,
      comments = input$comments,
      mask_min = input$mask_min,
      mask_max = input$mask_max,
      overlay_opacity = input$overlay_opacity,
      overlay_color_1 = input$overlay_color_1,
      overlay_color_2 = input$overlay_color_2,
      qc_time = as.character(Sys.time()),
      stringsAsFactors = FALSE
    )
    
    save_qc_row(row)
    showNotification(paste("Saved QC for", row0$subject), type = "message")
    
    TRUE
  }
  
  observeEvent(input$open_mricrogl, {
    open_case_in_mricrogl(selected_index())
  })
  
  observeEvent(input$prev_case, {
    i <- selected_index()
    if (i > 1) {
      updateSelectInput(session, "case_idx", selected = as.character(i - 1))
    }
  })
  
  observeEvent(input$next_case, {
    i <- selected_index()
    if (i < nrow(cases)) {
      updateSelectInput(session, "case_idx", selected = as.character(i + 1))
    }
  })
  
  observeEvent(input$save, {
    save_current_qc()
  })
  
  observeEvent(input$save_next, {
    save_current_qc()
    
    i <- selected_index()
    if (i < nrow(cases)) {
      updateSelectInput(session, "case_idx", selected = as.character(i + 1))
    }
  })
  
  observeEvent(input$save_next_open, {
    save_current_qc()
    
    i <- selected_index()
    
    if (i < nrow(cases)) {
      next_i <- i + 1
      updateSelectInput(session, "case_idx", selected = as.character(next_i))
      open_case_in_mricrogl(next_i)
    } else {
      showNotification("Already at last subject.", type = "warning")
    }
  })
}

shinyApp(ui, server)
