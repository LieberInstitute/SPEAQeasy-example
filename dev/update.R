## Update the sile

styler::style_dir(usethis::proj_path("dev"), transformers = biocthis::bioc_style())
styler::style_dir(
    usethis::proj_path(""),
    transformers = biocthis::bioc_style(),
    filetype = "Rmd"
)
