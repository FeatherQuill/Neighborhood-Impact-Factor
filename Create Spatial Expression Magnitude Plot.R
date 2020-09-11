#After settings are set, run code.

#Generate the graphs type "graph_gen(all_data)"


#[OPTIONS]:
#Saving Directory name (This saves in your working directory)
save_dir <- "images"

#Cell types (Enter Uppercase) (CASE-SENSITIVE)
ct_1 <- "CD34"
ct_2 <- "CEBPa"
ct_rad <- "HA"

#File/Dataset Name
f <- "objectNormalized_Nuclei.csv"



#DATA SETUP

#Initial Comments for Data:
#The script has been used on CSV formatted data to extract the specific column names of 
# ("ImageNumber", "ObjectNumber", 
# "noquote(sprintf("Children_%s_stain_Count", ct_1))", 
# "noquote(sprintf("Children_%s_stain_Count", ct_2))", 
# "Children_HA_stain_Count", 
# "Location_Center_X", 
# "Location_Center_Y", 
# "Math_normalized_cd34", 
# "Math_normalized_cebpa", 
# "Math_normalized_ha")
# It does not matter how many columns you have, as long as it has those specific names.
# Furthermore, if there are any NaN or NA spaces, those will be replaced with a zero.


#These libraries might have to be installed as it might not come standard
library(ggplot2)
library(data.table)
library(dplyr)
library(stringr)

# Cebpa+/cd34+    purple
# Cebpa+/cd34-     blue
# Cebpa-/cd34+     red
# Ceboa-/cd34-      gray


# Directory Creation
dr <- sprintf("%s/%s", noquote(getwd()), noquote(save_dir))
dir.create(dr)

# noquote(sprintf("Children_%s_stain_Count", ct_1))
# noquote(sprintf("Children_%s_stain_Count", ct_2))
# sprintf("Math_normalized_%s", tolower(ct_1))
# sprintf("Math_normalized_%s", tolower(ct_2))
# sprintf("%s+/%s+", ct_1, ct_2)

chld_ct1 <- sprintf("Children_%s_stain_Count", ct_1)
chld_ct2 <- sprintf("Children_%s_stain_Count", ct_2)
mth_ct1 <- sprintf("Math_normalized_%s", tolower(ct_1))
mth_ct2 <- sprintf("Math_normalized_%s", tolower(ct_2))
mth_rad <- sprintf("Math_normalized_%s", tolower(ct_rad))
population <- "population" #Used in Part 5

data_initial <- fread(sprintf("%s", f), select= c(
  "ImageNumber", "ObjectNumber", 
  chld_ct1, chld_ct2, 
  "Children_HA_stain_Count", 
  "Location_Center_X", 
  "Location_Center_Y", 
  mth_ct1, mth_ct2, mth_rad), stringsAsFactors = FALSE)

data_initial[is.na(data_initial)] <- 0

#This is the CD34+/CEBPa+ population. Includes both high and low HA from the dataframe.
#Setup so each is equivalent HighCD34_HighHA is paired with HighCEBPa_HighHA and LowHA with LowHA.

a <- subset(data_initial, data_initial[[noquote(sprintf("Children_%s_stain_Count", ct_1))]] > 0 
            & data_initial[[noquote(sprintf("Children_%s_stain_Count", ct_2))]] > 0)

a$population <- c(rep(sprintf("%s+/%s+", ct_1, ct_2), nrow(a)))

#This population is setup so the low CD34 is treated as a minus in the CD34-/CEBPa+ setup.

b <- subset(data_initial, data_initial[[noquote(sprintf("Children_%s_stain_Count", ct_1))]] == 0 
            & data_initial[[noquote(sprintf("Children_%s_stain_Count", ct_2))]] > 0)

b$population <- c(rep(sprintf("%s-/%s+", ct_1, ct_2), nrow(b)))


#Same as above but for CEBPa-

c <- subset(data_initial, data_initial[[noquote(sprintf("Children_%s_stain_Count", ct_1))]] > 0 
            & data_initial[[noquote(sprintf("Children_%s_stain_Count", ct_2))]] == 0)

c$population <- c(rep(sprintf("%s+/%s-", ct_1, ct_2), nrow(c)))

#CD34-/CEBPa-

d <- subset(data_initial, data_initial[[noquote(sprintf("Children_%s_stain_Count", ct_1))]] == 0 
            & data_initial[[noquote(sprintf("Children_%s_stain_Count", ct_2))]] == 0)

d$population <- c(rep(sprintf("%s-/%s-", ct_1, ct_2), nrow(d)))

#Combining all of the different populations for splitting into images.

all_data <- rbind(a, b, c, d)



#col_scale
pop_order <- unique(all_data$population)
col_scale <- scale_fill_manual(limits = c(pop_order), values = c
                               ("#cc00cc", # ct_1+/ct_2+ PURPLE/MAGENTA
                               "#0000ff", # ct_1-/ct_2+ BLUE
                               "#ff0000", # ct_1+/ct_2- RED
                               "#d3d3d3" # ct_1-/ct_2- GREY
                               ))


                               



#Include Debugging Plot

tes <- ggplot(a) + aes(a$ImageNumber, a$ObjectNumber) + geom_point()
j <- tes + ggtitle("DEBUGGING")
print(j)


#Graphing Function
graph_gen <- function(x, na.rm = TRUE, ...)
{
  imgs <- unique(all_data$ImageNumber)
  
  for (i in seq_along(imgs))
  {
    # #XY Location Creation
    p <- ggplot(subset(all_data, all_data$ImageNumber == imgs[i]),
                aes(x=Location_Center_X, y=Location_Center_Y, fill = population)) + col_scale
    l <- p + geom_point(aes(size=sqrt((subset(all_data, all_data$ImageNumber == imgs[i])[[mth_rad]])/pi)), pch = 21) +
      guides(col = guide_legend(override.aes = list(shape = 20, size = 10))) +
      scale_size_continuous(range = c(2.5,7.5),
                            #guide = guide_legend(title = mth_rad, title.theme = element_text(size=25, angle = 0)),
                            guide = FALSE) +
      guides(fill = guide_legend(override.aes = list(size=15),
                                 keywidth = 3, keyheight = 3, title = "Population", title.theme = element_text(size=25, angle = 0)
      )) +
      theme(plot.title = element_text(size = 25, face = "bold"),
            axis.text = element_text(size=20),
            axis.title = element_text(size=20, face= "bold"),
            legend.text = element_text(size=20) ) +
      ggtitle(sprintf("Image %s %s and %s Location Plot", 1, ct_1, ct_2))
    ggsave(l, file=sprintf("img_%s_XY.png", i), path = dr, scale = 3, dpi = 800)

    #ct_rad vs ct1
    p <- ggplot(subset(all_data, all_data$ImageNumber == imgs[i]),
                aes_string(x=mth_rad, y=mth_ct1, fill = population))
    s <- p + geom_point(pch = 21, size = 3) +
      col_scale + ggtitle(sprintf("Image %s %s vs %s Plot", i, ct_rad, ct_1)) +
      guides(fill = guide_legend(override.aes = list(size=15),
                                 keywidth = 3, keyheight = 3, title = "Population", title.theme = element_text(size=25, angle = 0))) +
      theme(plot.title = element_text(size = 25, face = "bold"),
            axis.text = element_text(size=20),
            axis.title = element_text(size=20, face= "bold"),
            legend.text = element_text(size=20) )
    ggsave(s, file=sprintf("img_%s_%s_%s_scatter.png", i, ct_rad, ct_1), path = dr, scale = 2.5, dpi = 800)


    #ct_rad vs ct2
    p <- ggplot(subset(all_data, all_data$ImageNumber == imgs[i]),
                aes_string(x=mth_rad, y=mth_ct2, fill = population))
    s <- p + geom_point(pch = 21, size = 3) +
      col_scale + ggtitle(sprintf("Image %s %s vs %s Plot", i, ct_rad, ct_2)) +

    guides(fill = guide_legend(override.aes = list(size=15),
                               keywidth = 3, keyheight = 3, title = "Population", title.theme = element_text(size=25, angle = 0))) +
      theme(plot.title = element_text(size = 25, face = "bold"),
            axis.text = element_text(size=20),
            axis.title = element_text(size=20, face= "bold"),
            legend.text = element_text(size=20) )

    ggsave(s, file=sprintf("img_%s_%s_%s_scatter.png", i, ct_rad, ct_2), path = dr, scale = 2.5, dpi = 800)

    #Part 5 Cell_type 1 vs Cell_type 2
    p <- ggplot(subset(all_data, all_data$ImageNumber == imgs[i]),
                aes_string(x=mth_ct1, y=mth_ct2, fill = population))
    s <- p + geom_point(pch = 21, size = 3) +
      col_scale + ggtitle(sprintf("Image %s %s vs %s Plot", i, ct_1, ct_2)) +

    guides(fill = guide_legend(override.aes = list(size=15),
                               keywidth = 3, keyheight = 3, title = "Population", title.theme = element_text(size=25, angle = 0))) +
      theme(plot.title = element_text(size = 25, face = "bold"),
            axis.text = element_text(size=20),
            axis.title = element_text(size=20, face= "bold"),
            legend.text = element_text(size=20) )

    ggsave(s, file=sprintf("img_%s_%s_%s_scatter.png", i, ct_1, ct_2), path = dr, scale = 2.5, dpi = 800)

    #Part 6 Histograms
    p <- ggplot(subset(all_data, all_data$ImageNumber == imgs[i] & all_data[[mth_rad]] != 0),
                aes_string(mth_rad, fill = population))
    h <- p + geom_histogram(binwidth = 0.1, colour="black", lwd=1) + col_scale +
      scale_x_continuous(breaks=seq(0.05,max(all_data[[mth_rad]]), 0.1)) +
      ggtitle(sprintf("Image %s HA Histogram", i)) +

    guides(fill = guide_legend(override.aes = list(size=15),
                               keywidth = 3, keyheight = 3, title = "Population", title.theme = element_text(size=25, angle = 0))) +
      theme(plot.title = element_text(size = 25, face = "bold"),
            axis.text = element_text(size=20),
            axis.title = element_text(size=20, face= "bold"),
            legend.text = element_text(size=20) )

    ggsave(h, file=sprintf("img_%s_%s_%s_hist.png", i, ct_1, ct_2), path = dr, scale = 2.5, dpi = 800)

  }
  
}

