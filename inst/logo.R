library(ggplot2)
library(hexSticker)
library(sysfonts)
library(showtext)
library(ggtext)
library(systemfonts)
library(prospectr)

library(magick)

img <- image_read("inst/zebra_handdrawn-wb.png")
print(img)          # displays it
image_info(img)     # check it has alpha channel (colorspace/format info)

#   Fonts  
font_add(
  family = "helvetica_neue",
  regular = "/home/leo/.local/share/fonts/HelveticaNeue/HelveticaNeueLTStd-Lt.otf",
  bold = "/home/leo/.local/share/fonts/HelveticaNeue/HelveticaNeueLTStd-Bd.otf",
  italic = "/home/leo/.local/share/fonts/HelveticaNeue/HelveticaNeueLTStd-LtIt.otf",
  bolditalic = "/home/leo/.local/share/fonts/HelveticaNeue/HelveticaNeueLTStd-BdIt.otf"
)
font_add(
  family = "helvetica_neue_bdex",
  regular = "/home/leo/.local/share/fonts/HelveticaNeue/HelveticaNeueLTStd-BdEx.otf"
)
font_add(
  family = "helvetica_neue_ltex",
  regular = "/home/leo/.local/share/fonts/HelveticaNeue/HelveticaNeueLTStd-LtEx.otf"
)
showtext_auto()

#   Palette  
palette <- c(
  primary = "#64B445",
  secondary1 = "#289A93",
  secondary2 = "#4DB9D2",
  secondary3 = "#4F719A",
  accent1 = "#CF554E",
  accent2 = "#E08B55",
  accent3 = "#EAC473"
)


#   Sticker — spotlight OFF, will be added manually as centered subview  
s <- sticker(
  img,
  package = "",
  s_x = 1,
  s_y = 1.05,
  s_width = 1.4,
  s_height = 1.4,
  h_fill = "#4f719a",
  h_size = 1.4,
  h_color = "#3B82F6", ##64B445
  filename = "proximetricsr_hex.png",
  dpi = 300,
  url = "https://buchi-nir.io/",
  # u_x = 0.94,
  # u_y = 0.085,
  u_family = "helvetica_neue_ltex",
  u_color = "white",
  u_size = 3,
  spotlight = FALSE
)

s <- s +
  geom_richtext(
    aes(x = 1, y = 0.50),
    label = "<span style='color:#FFFFFF;'>proximetrics</span><span style='color:#FFFFFF00;'>R</span>",
    family = "helvetica_neue_bdex",
    fontface = "plain",
    size = 12,
    fill = NA,
    label.color = NA
  ) +
  geom_richtext(
    aes(x = 1.02, y = 0.42),
    label = "<span style='color:#FFFFFF00;'>proximetrics</span><span style='color:#64B445;'>R</span>",
    family = "helvetica_neue_bdex",
    fontface = "plain",
    size = 12,
    fill = NA,
    # hjust = 0,
    angle = 8,
    label.color = NA
  ) 

s

ggsave("man/figures/logo.png", s, width = 43.94, height = 50.8, units = "mm", dpi = 300, bg = "transparent")


