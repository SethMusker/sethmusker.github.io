# finalise a post after knitting:
# - move the .md file from _drafts/ → _posts/
# - move image folders into assets/img/
# - update figure paths inside the markdown
finalise_post <- function(title = NULL,
                          date = Sys.Date(),
                          drafts_dir = "_drafts",
                          posts_dir = "_posts",
                          assets_dir = "assets/img") {
  
  # If no title passed, ask interactively
  if (is.null(title)) {
    title <- rstudioapi::showPrompt("finalise Post", "Enter post title:")
    if (is.null(title) || title == "") return(invisible())
  }
  
  # Same slugify helper as above
  slugify <- function(x) {
    x <- tolower(x)
    x <- gsub("[^a-z0-9]+", "-", x)
    x <- gsub("(^-|-$)", "", x)
    x
  }
  
  slug <- slugify(title)
  date_str <- as.character(date)
  base <- paste0(date_str, "-", slug)
  
  # File paths
  md_src <- file.path(drafts_dir, paste0(base, ".md"))  # knitted file
  md_dst <- file.path(posts_dir,  paste0(base, ".md"))  # final location
  
  # Paths for figures created by knitting
  # fig_src_dir <- file.path(drafts_dir, paste0(base, "_files"))
  # fig_dst_dir <- file.path(assets_dir, base)
  
  # Must knit the file first, otherwise .md won't exist
  if (!file.exists(md_src)) {
    stop("Markdown not found: Knit the .Rmd first")
  }
  
  # Ensure required folders exist
  dir.create(posts_dir, showWarnings = FALSE)
  # dir.create(assets_dir, showWarnings = FALSE)
  
  # Move the markdown into _posts/
  file.rename(md_src, md_dst)
  
  # Move figures (if any)
  # if (dir.exists(fig_src_dir)) {
  #   dir.create(fig_dst_dir, recursive = TRUE)
  #   file.copy(fig_src_dir, fig_dst_dir, recursive = TRUE)
  #   unlink(fig_src_dir, recursive = TRUE)  # delete old folder
  # }
  
  # Rewrite figure paths inside the markdown file
  # knitr uses paths like: base_files/figure-gfm/plot.png
  # Jekyll needs them like: /assets/img/base/figure-gfm/plot.png
  # md <- readLines(md_dst)
  
  # md <- gsub(
  #   paste0(base, "_files/figure-gfm/"),
  #   paste0("/", fig_dst_dir, "/figure-gfm/"),
  #   md
  # )
  
  # writeLines(md, md_dst)
  
  # Open the final markdown post in RStudio
  rstudioapi::navigateToFile(md_dst)
  
  message("✔ Post finalised: ", md_dst)
  # message("✔ Figures moved to: ", fig_dst_dir)
  
  invisible(md_dst)
}

# Utility operator:
# a %||% b means: return a if non-null/non-empty, else b
`%||%` <- function(a, b) if (is.null(a) || length(a)==0) b else a
