# A new way of making maps

library(data.table)
library(ggplot2)
library(here)
library(sp)
library(raster)
library(viridis)

## Plot a map

# Load in Carlos's island outline
bi.outline <- raster::shapefile(here("data/clean/BI_maps/bioko/bioko_admin1.shp"))
plot(bi.outline)

# Load in Carlos's centroids
bi.centroids <- raster::shapefile(here("data/raw/mapArea_centroids/mapareagrid_centroids.shp"))

#geom_polygon(data = bioko, aes(x = long, y = lat, group = group), color = "black", fill="grey", size = 0.25)

# Load in the population data
pop.data <- fread(here("data/clean/aggregated_2015_2018_travel_data.csv"))



areas.inh <- merge(bi.centroids, pop.data[year == 2018, .(areaId, ad2, pop)], by.x = "OBJECTID", by.y = "areaId", all = FALSE)

# Look! It's a beautiful population plot!
example.plot <- ggplot() + 
  geom_polygon(data = bi.outline, aes(x = long, y = lat, group = group), color = "black", fill="grey", size = 0.25) + 
  geom_point(data = data.frame(areas.inh), mapping = aes(x = coords.x1, y = coords.x2, color = log10(pop)), 
             shape = 15, size = 2.1) + 
  scale_color_gradient(name = "Log(Population)", low = "#FFFFFF", high = "#5519C4", 
                       breaks = c(0,1,2,3,4), 
                       labels = c(1, 10, bquote(10^2), bquote(10^3), bquote(10^4))) +
  #scale_color_viridis(name="Log(Population)", limits=c(0, 5), option="plasma") + # set scale bar parameters
  # move the scale bar
  coord_fixed() + # aspect ratio
  theme(aspect.ratio = 1, 
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.background = element_blank(),
        legend.title = element_text(size =18),
        legend.text = element_text(size = 14),
        legend.position = c(.25, .8))
example.plot

ggsave("/Users/dtcitron/Documents/Presentation Materials/ASTMH_2019/Figures/example_plot.pdf", 
       plot = example.plot,
       width = 6, height = 6, units = "in"
      )

## Plot a time series
mat_mean <- fread(here("data/ASTMH/data/baseline_pr_means.csv"))
mat_sd <- fread(here("data/ASTMH/data/baseline_pr_sds.csv"))

h <- melt(mat_mean[areaId %in% c(220,335,502,1175,2199,2457)],
          id.vars = c("time", "areaId"), 
          measure.vars = c("s","i","p"),
          value.name = "fraction")
h.sd <- melt(mat_sd[areaId %in% c(220,335,502,1175,2199,2457)],
             id.vars = c("time", "areaId"), 
             measure.vars = c("s","i","p"),
             value.name = "fraction.sd")
h <- merge(h, h.sd, by = c("time", "areaId", "variable"))

h[, lower.bound := max(fraction - fraction.sd, 0), by = c("areaId", "variable", "time")]

ggplot(data = h[areaId == 335 & variable == "i"]) + 
  geom_errorbar(mapping = aes(x = time, 
                              ymin = lower.bound,
                              ymax = fraction + fraction.sd, 
                              color = variable), 
                alpha = .2) + 
  scale_x_continuous(name = "Time (Years)", 
                     breaks = c(365, 365*2, 365*3, 365*4, 365*5), labels = c(1:5),
                     expand = c(0,0)) + 
  scale_y_continuous(name = "Prevalence", 
                     expand = c(0,0),
                     limits = c(0,.25)) + #ylab("Prevalence") +
  geom_point(mapping = aes(x = time, y = fraction), color = "black", shape = 16, size = .1) +
  theme(panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(color = "black", size = 1),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        panel.grid = element_line(color = NULL),
        legend.position = c(.8, .8))

