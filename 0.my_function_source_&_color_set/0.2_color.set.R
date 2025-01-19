#--------------------------------
#>>>> Set color
#--------------------------------
colors.time.2=c('#FF0000','#FFA500','#ffff00','limegreen','#30d5c8','dodgerblue','#4B0082',"#8DD3C7","#8DD3C8")
names(colors.time.2)=c('9ss','12ss','15ss','18ss','21ss','24ss','27ss','E10.5',"E10.5ss")

colors.time=c('#FF0000','#FFA500','#ffff00','limegreen','#30d5c8','dodgerblue','#4B0082')
names(colors.time)=c('ss9','ss12','ss15','ss18','ss21','ss24','ss27')

colors.time = c(colors.time, colors.time.2)

color.mg1 = c("#90C9E7","#219EBC","#4ad2d4","#DAADF7","#DF3AC9")
names(color.mg1) = c("MG.1","MG.1-Stomach","Stomach",
                     "MG.1-Small.intestine.1","Small.intestine.1")
color.mg2 = c("#737ADE","#5252F5","#2A15B9","#F9AEEC","#DF3AC9")
names(color.mg2) = c("MG.2","MG.2-Small.intestine.2","Small.intestine.2",
                     "MG.2-Small.intestine.1","Small.intestine.1")
color.mg3 = c("#4169e1",
              "#5000b8",
              '#66cdaa','#66cdaa','#66cdaa',
              '#008b8b','#33FFAD','#318fac','#6FEEFF','#A01BA7',
              '#A10BA7',
              '#0047ab',"#2A15B9",
              '#9370db',"#DF3AC9")
names(color.mg3) = c("MG.3",
                     "MG.3.P",
                     "MG.3.A/M","MG.3.A","MG.3.M",
                     "DP","EP.1","EP.2","EP","Pancreas",
                     "Small.intestine",
                     "MG.3.P-Small.intestine.2","Small.intestine.2",
                     "MG.3.P-Small.intestine.1","Small.intestine.1")
mg.3.color = c(color.mg1,color.mg2,color.mg3)


color.al3 = c("#F7BF63","#C7946E","#9e7a00","#D1BA98","#8b4513","#8b4514",
              "#FAAFDC","#DF3AC9","#821500")
names(color.al3) = c("AL.3","AL.3-VP","VP","AL.3-EHBD/VP","EHBD","Extrahepatic biliary tract",
                     "AL.3-Small.intestine.1","Small.intestine.1","AL.3-Liver")

color.al12 = c("#DC7E95","#E9D62C","#F7B543","#F0983F","#DEABC0","#f2855e")
names(color.al12) = c("AL.1","AL.2","AL.1-Liver","AL.2-Liver","AL.1/2-Liver","Liver")
al.3.color = c(color.al12,color.al3)


color.fg5 = c("#FA84AC","#D23A7C","#a12559","#6694AD","#3D71A2")
names(color.fg5) = c("FG.5","FG.5-Thyroid","Thyroid",
                     "FG.5-Pharynx.organ.3","Pharynx.organ.3")
color.fg4 = c("#7f4a88","#B274BC","#C97D66","#BF4F87","#8EA0DC",#"#4ad2d4",
              "#C0732E")
names(color.fg4) = c("FG.4","FG.4-Lung/Stomach","FG.4-Liver",
                     "FG.4-Lung", "FG.4-Stomach",#"Stomach",
                     "FG.4-Pharynx.organ.4")

color.fg3 = c("#ff938b","#FF3099","#C61549","#8EA0Dc",
              "#cc7320","#BA7A74")
names(color.fg3) = c("FG.3","FG.3-Lung","Lung","Pharynx.organ.5",
                     "FG.3-Pharynx.organ.4","Pharynx.organ.4")

color.fg1 = c("#f0bec1","#ffbbff",# "#A01e8E",
              "#FFAEDC","#DF7A6B","#FF6262","#A42613",
              "#91b48e","#5C0F36",
              "#E68686","#cca4e3")
names(color.fg1) = c("FG.1","Esophagus",# "FG.6",
                     "FG.2","FG.1-Pharynx.organ.1","FG.2-Pharynx.organ.1","Pharynx.organ.1",
                     "FG.1-Pharynx.organ.2","Pharynx.organ.2",
                     "FG.1-Esophagus","FG.6" # "Esophagus"
                     )

fg.3.color = c(color.fg1,#color.fg2,
               color.fg3,color.fg4,color.fg5)


color.hg2 = c("#057748","#BCCCB8","#8CDB83","#B9C78D","#9acd32")
names(color.hg2) = c("HG.2","HG.2-Large.intestine.1","Large.intestine.1",
                     "HG.2-Large.intestine.3","Large.intestine.3")

color.hg1 = c("#1abb9c","#5DB7AE","#2A15B9",
              "#37A254","#8CDB83",
              "#69b076","#36686B",
              "#274635")
names(color.hg1) = c("HG.1","HG.1-Small.intestine.2","Small.intestine.2",
                     "HG.1-Large.intestine.1","Large.intestine.1",
                     "HG.1-Large.intestine.2","Large.intestine.2",
                     "HG.1-Sm.2/Lar.1")
hg.3.color = c(color.hg1,color.hg2)


color_plus = c( # "#DF3AC9",
  "#A497FF","#964E5C","#D26910","#A8A865","#F9f788",
  "#3DBFF0",'#CC00ff',"#ddcf44",
  "#BC69BB","#143F9A","#240cd8","#274635",
  "#9AD2A2","#C5fF9B")
names(color_plus) = c(
  "FG.4-MG.1","FG.3/4",'FG.4-AL.1/2/3',"AL.3-MG.2","AL.1/2",
  "FG.4-MG.1/3","AL.3-MG.1/2/3","FG.4-AL.1/2",
  "MG.1/3","MG.1/2/3","MG.2/3","HG.1/2",
  "MG.2/3-HG.1","HG.1/2")

color_plus.add = c("#d9333f","#ffcf00",'#0e4bef','LimeGreen')
names(color_plus.add) = c("FG","AL","MG","HG")

cluster.endoderm.color.v5 = c(fg.3.color,al.3.color,
                              mg.3.color,hg.3.color,
                              color.na,color_plus,color_plus.add)


color.lineage.re = c(
  "#BF8B50","#EE4C97","#BC3C29","#274635",
           "#0072B5","#0072B6","#8358CA","#E77C2F","#419379",'#ffc110','#ffc111',
           "#50d2f0","#FCCDE5")
           
names(color.lineage.re) = c(
  "Pax9","Sox2","Nepn",'Wnt5b',
  "Nkx2-3","Nkx2_3","Hhex",'Pdx1',"Mnx1","Mnx1GFP","Mnx1\nGFP",
  "Ngn3Cre","Ngn3GFP")

colors.type.lineage = c("#ffdf12","#8F05dc","#ffdf12","#8F05dc","#ffdf12","#8F05dc") #,"#8DD3C7")
names(colors.type.lineage) = c("Ngn3GFPHM","Ngn3GFPHZ","Neurog3GFPHM","Neurog3GFPHZ",
                               "Ngn3GFP_HM","Ngn3GFP_HZ") #,"Mnx1")

color.lineage = c(color.lineage.re, colors.type.lineage)

colors.geneset = c(
  "#d76364", "#9dc3e7", "#f1d77e", "#b1ce46",
           "#63e398", "#9394e7", "#5f97d2", "#14517c",
           "#ef7a6d", "#f7e1ed", "#c497b2", "#f8f3f9")
names(colors.geneset) = c(1:12) # 4 3 9 7
           
           
colors.num = c(
  "#F1D2BC", "#D8B694", "#D49261", "#CC5A35", "#B03015", "#5C2519",
           "#A78693", "#9F86CF", "#9970D7", "#4A346F",
           "#C3F1E6", "#6AD4C2", "#8CCFA3", "#40908B", "#2A595B", "#92745C",
           "#AFD5F1", "#2F9BE1", "#3971CB", "#2348A2", "#273F7B", "#141F45",
           "#ECDEE2", "#D190D8", "#C66A7A", "#A54050", "#612F4E", "#3F2AA0",
           "#F4E8BC", "#F9F388", "#D1A955", "#C7AA89", "#715738",
           "#CAE3F1", "#A9D7EC", "#A8C2D5", "#436B8D", "#807D80", "#92533E",
           "#FFFCCB", "#E1E993", "#A3B256", "#82A843", "#C9AE85", "#B69277",
           "#EDF4DC", "#DBEBCA", "#B6DABE", "#88BB9C", "#7DAF75", "#577D62")
names(colors.num) = c(0,1:(length(colors.num)-1))
           

library("ggplot2")

p.add = 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks =  element_blank(),
        axis.text = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        plot.title = element_blank(),
        aspect.ratio=1)


p_add = 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks =  element_blank(),
        axis.text = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        #plot.title = element_blank(),
        aspect.ratio=1)

p_add_leg = 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks =  element_blank(),
        axis.text = element_blank(),
        # legend.title = element_blank(),
        # legend.text = element_text(size = 10),
        #plot.title = element_blank(),
        aspect.ratio=1)






