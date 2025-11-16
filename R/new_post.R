# Create a new R Markdown draft for a Jekyll post
new_post <- function(title = NULL,
                     date = Sys.Date(),
                     categories = NULL,
                     tags = NULL,
                     drafts_dir = "_drafts") {
  
  # If the user didn't supply a title, ask interactively
  if (is.null(title)) {
    title <- rstudioapi::showPrompt("New Post", "Enter post title:")
    if (is.null(title) || title == "") return(invisible())
  }
  
  # Helper function:
  # Convert a title like "My New Post!" → "my-new-post"
  slugify <- function(x) {
    x <- tolower(x)
    x <- gsub("[^a-z0-9]+", "-", x)  # replace non-alphanumerics with '-'
    x <- gsub("(^-|-$)", "", x)      # trim leading/trailing dashes
    x
  }
  
  slug <- slugify(title)
  filename <- paste0(as.character(date), "-", slug, ".Rmd")
  path <- file.path(drafts_dir, filename)
  
  # Make sure the drafts folder exists
  dir.create(drafts_dir, showWarnings = FALSE)
  
  # Every post gets its own image directory
  img_dir <- file.path("assets", "img", paste0(as.character(date), "-", slug))
  dir.create(img_dir, , showWarnings = FALSE)

  # --- YAML front matter ---
  # preserve_yaml = TRUE ensures Jekyll sees this untouched
  yaml <- c(
    "---",
    paste0('title: "', title, '"'),
    paste0("date: ", as.character(date)),
    paste0("categories: [", paste(categories %||% "", collapse=", "), "]"),
    paste0("tags: [", paste(tags %||% "", collapse=", "), "]"),
    "output:",
    "  md_document:",
    "    variant: gfm",       # GitHub-style markdown (works with Jekyll)
    "    preserve_yaml: true",# Keep YAML at top when knitting
    "---",
    ""
  )
  
  template <- c(
    yaml,
    "```{r setup, include=FALSE}",
    "library(knitr)",
    sprintf("knitr::opts_chunk$set(fig.path = '%s', fig.width = 6, fig.height = 4)",
            paste0("../", img_dir, "/")),
    "```"
  )

  # A simple starter template
  body <- c(
    paste0("# ", title),
    "",
    "## Introduction",
    "",
    "```{r Introduction}",
    "library(tidyverse)",
    "```",
    "",
    "## Section 1",
    "",
    "```{r Section-1}",
    "library(tidyverse)",
    "```",
    "",
    "## Section 2",
    "",
    "```{r Section-2}",
    "library(tidyverse)",
    "```",
    "",
    "## Section 3",
    "",
    "```{r Section-3}",
    "library(tidyverse)",
    "```",
    ""
  )
  
  # Write new .Rmd file
  writeLines(c(template, body), con = path)
  
  # Open the file in RStudio for editing
  rstudioapi::navigateToFile(path)
  
  message("✔ New draft created at ", path)
  invisible(path)
}
