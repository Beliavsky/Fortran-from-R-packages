### define font
library("showtext")
font_add_google("Martel", "my_font")
showtext_auto()

### build plot
library("ao")
himmelblau <- function(x) {
  ((x[1] - 2)^2 + (x[2] - 2) + 11)^2 + ((x[1] - 2) + (x[2] - 2)^2 + 7)^2 + 10
}
init <- c(1.6, 1.2)
out <- ao(
  f = himmelblau, p = init, partition = list(1, 2),
  base_optimizer = optimizer_optim(lower = -5, upper = 5, method = "L-BFGS-B"),
  tolerance = 0.1
)
library("ggplot2")
df <- expand.grid(
  x = seq(1.57, 1.8, length.out = 200),
  y = seq(1.17, 1.4, length.out = 200)
)
df$z <- apply(df, 1, himmelblau)
p <- ggplot() +
  geom_contour(
    data = df, aes(x = x, y = y, z = z),
    breaks = c(
      apply(out$sequence[c("p1", "p2")], 1, himmelblau), 169.65, 169.619
    ),
    colour = "#93032E"
  ) +
  geom_point(
    data = out$sequence, aes(x = p1, y = p2),
    colour = "#E94F37"
  ) +
  geom_line(
    data = out$sequence, aes(x = p1, y = p2),
    colour = "#E94F37"
  ) +
  theme_void() +
  scale_x_log10()
plot(p)

### build sticker
library("hexSticker")
sticker_file <- sticker(
  ### image
  subplot = p,
  s_x = 1,
  s_y = 0.8,
  s_width = 1.7,
  s_height = 0.8,
  ### package name
  package = "ao",
  p_x = 1,
  p_y = 1.5,
  p_color = "#000000",
  p_family = "my_font",
  p_fontface = "plain",
  p_size = 40,
  ### sticker
  h_size = 1.2,
  h_fill = "#F6F7EB",
  h_color = "#000000",
  spotlight = FALSE,
  l_x = 0.9,
  l_y = 1.4,
  l_width = 2,
  l_height = 1,
  l_alpha = 0.8,
  white_around_sticker = FALSE,
  ### URL
  url = "loelschlaeger.de/ao",
  u_x = 1,
  u_y = 0.1,
  u_color = "#000000",
  u_family = "my_font",
  u_size = 7,
  u_angle = 30,
  ### save file
  filename = "sticker/ao_sticker.png",
  asp = 1,
  dpi = 300
)
plot(sticker_file)
