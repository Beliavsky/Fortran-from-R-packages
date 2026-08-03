library(showtext)
library(tidyverse)
library(ggpath)
library(ggtext)
library(glue)
library(cropcircles)
library(magick)


image <- image_read("devel/logo.webp")

font_add("an", regular = "devel/ArchivoNarrow-Regular.ttf")
font_add("an", regular = "devel/ArchivoNarrow-Bold.ttf")

pkg_name <- "R4GoodPersonalFinances"
pkg_name <- "R4Good\nPersonalFinances"

txt <- "#2780E3"
img_cropped <- hex_crop(
  images = image,
  border_colour = txt,
  border_size = 24
)
img_cropped

offset <- 0.004
lineheight <- 0.75
moving_up <- 0.03
moving_left <- -0.02
size <- 10.5
angle <- 30

ggplot() +
  geom_from_path(aes(0.5, 0.5, path = img_cropped)) +
  annotate("text",
    x = 0.48 - offset + moving_left,
    y = 0.01 - offset + moving_up,
    size = size,
    label = pkg_name,
    family = "an",
    colour = txt,
    fontface = "bold",
    angle = angle,
    hjust = -offset, lineheight = lineheight
  ) +
  annotate("text",
    x = 0.48 + moving_left,
    y = 0.01 + moving_up,
    size = size,
    label = pkg_name,
    family = "an",
    color = "white",
    fontface = "bold",
    angle = angle,
    hjust = 0, lineheight = lineheight
  ) +
  annotate("text",
    x = 0.06,
    y = 0.26,
    size = 8,
    label = "R4Good.Academy",
    family = "an",
    colour = txt,
    fontface = "bold",
    angle = -angle,
    hjust = 0,
    lineheight = lineheight
  ) +
  annotate("text",
    x = 0.06 + offset,
    y = 0.26 + offset,
    size = 8,
    label = "R4Good.Academy",
    family = "an",
    colour = "white",
    fontface = "bold",
    angle = -angle,
    hjust = 0,
    lineheight = lineheight
  ) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_void() +
  coord_fixed()

# ggsave(
#   "devel/r4gpf-hex.png",
#   scale = 2.3,
#   bg = "transparent"
# )

# usethis::use_logo("devel/r4gpf-hex.png")
