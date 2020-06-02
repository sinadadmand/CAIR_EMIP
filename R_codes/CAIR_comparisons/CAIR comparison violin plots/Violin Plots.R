## 0. Loading packages ##
#############################################################################################
if (!require("plotly")) install.packages("plotly")
if (!require("brunnermunzel")) install.packages("brunnermunzel")
if (!require("data.table")) install.packages("data.table")
if (!require("processx")) install.packages("processx")
library (plotly)
library (brunnermunzel)
library (data.table)
library (processx)
###########################
## 1. Preparing the data ##
#############################################################################################
directory <- getwd()
df <- fread(paste(directory, 'Supplementary/runCAIR/Complete proteome CAIRs.csv', sep = '/'),
            select = c(1,2,3,5))
colnames(df) <- c("org", "cair", "tax", "code")
each_plot_width = 70
assign('num_of_plots', 0)
#######################################################################################
## 2. Creating a function to draw violin plots and attach Brunner Munzel test results##
#############################################################################################
violinplot <- function(test_num, grp1, grp2){
    num_of_plots <<- num_of_plots + 1
    BMtest <- brunnermunzel.test(x = df$cair[grp1],
                                 y = df$cair[grp2])
    if (BMtest$p.value < 0.0001){
      sig <- "<b>***</b>"
    } else if (BMtest$p.value < 0.001){
      sig <- "<b>**</b>"
    } else if (BMtest$p.value < 0.05){
      sig <- "<b>*</b>"
    } else {
      sig <- ""
    }
    
    test <- paste(test_num, sig, ' ES: ',
                  formatC(BMtest$estimate, format = 'g', digits = 2),
                  '\nCI: ',
                  formatC(BMtest$conf.int[1], format = 'g', digits = 2),
                  ' - ',
                  formatC(BMtest$conf.int[2], format = 'g', digits = 2),
                  sep = '')
    
    p <- plot_ly(type = 'violin') %>%
      add_trace(y = df$cair[grp1],
                x0 = test,
                side = 'negative',
                width = each_plot_width,
                opacity = 0.85,
                pointpos = -0.4,
                points = 'suspectedoutliers',
                fillcolor = "#4DBBD5",
                box = list(visible = T, width = 0.2),
                meanline = list(visible = T, color="#91D1C2"),
                line = list (width = 0.3, color = "#29636D"),
                marker = list (size = 1.5, color="#29636D",
                               line = list (outlierwidth = 0.15),
                               symbol="diamond")) %>% 
      
      add_trace(y = df$cair[grp2],
                x0 = test,
                side = 'positive',
                width = each_plot_width,
                opacity = 0.85,
                pointpos = 0.4,
                points = 'suspectedoutliers',
                fillcolor = "#00A087",
                box = list(visible = T, width = 0.2),
                meanline = list(visible = T, color="#91D1C2"),
                line = list (width = 0.3, color = "#005141"),
                marker = list (size = 1.5, color="#005141",
                               line = list (outlierwidth = 0.15),
                               symbol="diamond")) %>% 
      
      layout(font = list(family = "Arial, sans-serif", size = 6.5),
            yaxis = list(zeroline = F, color = "#444444",
            tickfont = list(family = "Arial, sans-serif", size = 6.5)),
            showlegend = F)
}

################################
## 3. Runing comparison tests ##
#############################################################################################
# WARNING: The tests contating a group with less than 10
# observations are #NOT to be executed and thus #HASHED
#.......................................................
# The below codes are based on the tree of life
#.......................................................
#library(plyr)
#count(df, vars = "code")
fig_rows = 10
fig <- subplot(
  violinplot('T1', df$code<=56, df$code>=57),
  violinplot('T2', df$code<=55, df$code==56),
  violinplot('T3', df$code<=51, df$code<=55 & df$code>=52),
  violinplot('T4', df$code<=42, df$code<=51 & df$code>=43),
  violinplot('T5', df$code<=34, df$code<=42 & df$code>=35),
  violinplot('T6', df$code<=18, df$code<=34 & df$code>=19),
  violinplot('T7', df$code<=10, df$code<=18 & df$code>=11),
  violinplot('T8', df$code<=6, df$code<=10 & df$code>=7),
  violinplot('T9', df$code<=4, df$code<=6 & df$code>=5),
  violinplot('T10', df$code==1, df$code<=4 & df$code>=2),
# violinplot('T11', df$code<=3 & df$code>=2, df$code==4),
# violinplot('T12', df$code==2, df$code==3),
# violinplot('T13', df$code==5, df$code==6),
  violinplot('T14', df$code==7, df$code<=10 & df$code>=8),
  violinplot('T15', df$code==8, df$code<=10 & df$code>=9),
  violinplot('T16', df$code==9, df$code==10),
  violinplot('T17', df$code==11, df$code<=18 & df$code>=12),
  violinplot('T18', df$code==12, df$code<=18 & df$code>=13),
  violinplot('T19', df$code==13, df$code<=18 & df$code>=14),
  violinplot('T20', df$code==14, df$code<=18 & df$code>=15),
  violinplot('T21', df$code<=16 & df$code>=15, df$code<=18 & df$code>=17),
# violinplot('T22', df$code==15, df$code==16),
# violinplot('T23', df$code==17, df$code==18),
# violinplot('T24', df$code==19, df$code<=34 & df$code>=20),
# violinplot('T25', df$code==20, df$code<=34 & df$code>=21),
  violinplot('T26', df$code<=22 & df$code>=21, df$code<=34 & df$code>=23),
# violinplot('T27', df$code==21, df$code==22),
  violinplot('T28', df$code<=32 & df$code>=23, df$code<=34 & df$code>=33),
  violinplot('T29', df$code==23, df$code<=32 & df$code>=24),
# violinplot('T30', df$code==24, df$code<=32 & df$code>=25),
# violinplot('T31', df$code==25, df$code<=32 & df$code>=26),
  violinplot('T32', df$code<=29 & df$code>=26, df$code<=32 & df$code>=30),
  violinplot('T33', df$code<=27 & df$code>=26, df$code<=29 & df$code>=28),
# violinplot('T34', df$code==26, df$code==27),
# violinplot('T35', df$code==28, df$code==29),
  violinplot('T36', df$code==30, df$code<=32 & df$code>=31),
# violinplot('T37', df$code==31, df$code==32),
  violinplot('T38', df$code==33, df$code==34),
  violinplot('T39', df$code<=36 & df$code>=35, df$code<=42 & df$code>=37),
  violinplot('T40', df$code==35, df$code==36),
  violinplot('T41', df$code==37, df$code<=42 & df$code>=38),
  violinplot('T42', df$code<=39 & df$code>=38, df$code<=42 & df$code>=40),
# violinplot('T43', df$code==38, df$code==39),
  violinplot('T44', df$code==40, df$code<=42 & df$code>=41),
# violinplot('T45', df$code==41, df$code==42),
  violinplot('T46', df$code<=49 & df$code>=43, df$code<=51 & df$code>=50),
  violinplot('T47', df$code<=48 & df$code>=43, df$code==49),
# violinplot('T48', df$code<=47 & df$code>=43, df$code==48),
  violinplot('T49', df$code<=45 & df$code>=43, df$code<=47 & df$code>=46),
# violinplot('T50', df$code<=44 & df$code>=43, df$code==45),
  violinplot('T51', df$code==43, df$code==44),
# violinplot('T52', df$code==46, df$code==47),
# violinplot('T53', df$code==50, df$code==51),
  violinplot('T54', df$code<=53 & df$code>=42, df$code<=55 & df$code>=54),
  violinplot('T55', df$code==52, df$code==53),
  violinplot('T56', df$code==54, df$code==55),
  violinplot('T57', df$code<=66 & df$code>=57, df$code>=67),
# violinplot('T58', df$code==57, df$code<=66 & df$code>=58),
  violinplot('T59', df$code<=60 & df$code>=58, df$code<=66 & df$code>=61),
  violinplot('T60', df$code==58, df$code<=60 & df$code>=59),
# violinplot('T61', df$code==59, df$code==60),
  violinplot('T62', df$code==61, df$code<=66 & df$code>=62),
# violinplot('T63', df$code==62, df$code<=66 & df$code>=63),
  violinplot('T64', df$code<=64 & df$code>=63, df$code<=66 & df$code>=65),
# violinplot('T65', df$code==63, df$code==64),
  violinplot('T66', df$code==65, df$code==66),
  violinplot('T67', df$code==67, df$code>=68),
  violinplot('T68', df$code<=75 & df$code>=68, df$code<=92 & df$code>=76),
  violinplot('T69', df$code<=70 & df$code>=68, df$code<=75 & df$code>=71),
  violinplot('T70', df$code==68, df$code<=70 & df$code>=69),
# violinplot('T71', df$code==69, df$code==70),
# violinplot('T72', df$code==71, df$code<=75 & df$code>=72),
  violinplot('T73', df$code<=73 & df$code>=62, df$code<=75 & df$code>=74),
# violinplot('T74', df$code==72, df$code==73),
  violinplot('T75', df$code==74, df$code>=75),
  violinplot('T76', df$code==76, df$code>=77),
  violinplot('T77', df$code==77, df$code>=78),
# violinplot('T78', df$code==78, df$code>=79),
  violinplot('T79', df$code<=87 & df$code>=79, df$code<=92 & df$code>=88),
  violinplot('T80', df$code==79, df$code<=87 & df$code>=80),
  violinplot('T81', df$code<=86 & df$code>=80, df$code==87),
# violinplot('T82', df$code==80, df$code<=86 & df$code>=81),
  violinplot('T83', df$code<=82 & df$code>=81, df$code<=86 & df$code>=83),
# violinplot('T84', df$code==81, df$code==82),
# violinplot('T85', df$code==83, df$code<=86 & df$code>=84),
# violinplot('T86', df$code==84, df$code<=86 & df$code>=85),
# violinplot('T87', df$code==85, df$code==86),
  violinplot('T88', df$code<=90 & df$code>=88, df$code<=92 & df$code>=91),
# violinplot('T89', df$code<=89 & df$code>=88, df$code==90),
# violinplot('T90', df$code==88, df$code==89),
# violinplot('T91', df$code==91, df$code==92),
  shareY = T, nrows = fig_rows)

##########################
## 4. Displaying output ##
#############################################################################################
orca(p = fig, file = "Fig1_Violinplots5.svg", format = 'svg',
     width = round(num_of_plots/fig_rows, 0)*each_plot_width,
     height = fig_rows*each_plot_width)
