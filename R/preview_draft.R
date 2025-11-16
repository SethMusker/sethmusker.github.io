preview_draft <- function() {
  # Detect active file
  context <- rstudioapi::getActiveDocumentContext()
  file <- context$path

  if (!nzchar(file)) {
    stop("No active file detected.")
  }

  if (!grepl("_drafts/", file)) {
    stop("Active file is not inside the _drafts/ directory.")
  }

  # Confirm knitting
  answer <- utils::menu(
    c("Yes — knit before previewing", "No — preview without knitting"),
    title = "Knit draft before preview?"
  )

  knit_it <- (answer == 1)

  # Knit only if user agreed
  if (knit_it) {
    message("Knit requested by user…")

    # # where to put assets
    post_slug <- sub(".*/_drafts/(.*)\\.R?md$", "\\1", file)

    assets_dir <- file.path("assets/img", post_slug)
    if (!dir.exists(assets_dir)) dir.create(assets_dir, recursive = TRUE)

    # Knit with correct figure directory
    rmarkdown::render(
      file,
      output_format = "github_document",
      output_dir = dirname(file),
      envir = new.env(parent = globalenv()),
      output_options = list(
        fig_path = file.path(assets_dir, "figure-")
      )
    )
  }

  # Preview site with live rebuilding
  message("Launching preview server using rmarkdown::serve_site()…")
  system2("/home/seth/gems/bin/bundle", "exec jekyll serve --drafts")
  # Open browser
  browseURL("http://127.0.0.1:4000/")
}
