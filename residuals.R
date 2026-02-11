
# load("output") it contains fitmodel, the output of the file "example_JEBS.R"
# the following codes help us to create box-plots of the standardized residuals.

fitmod<-fitmodel

# fitmod$vecres provides the standardized residuals.
# These can be grouped by Political orientation (SR, NR, IND, ND, SD)
# across the four strata defined by the combination of Race and Education
# (Black–Low, Black–High, White–Low, White–High).
# Creating boxplots by strata (see Figure 8 in the paper).

stdres11 <- fitmod$vecres[1:25,3] # SR B L
stdres21 <- fitmod$vecres[26:50,3] # NR 
stdres31 <- fitmod$vecres[51:75,3] # IND 
stdres41 <- fitmod$vecres[76:100,3] # ND 
stdres51 <- fitmod$vecres[101:125,3] # SD 

stdres12 <- fitmod$vecres[126:150,3] # SR B H
stdres22 <- fitmod$vecres[151:175,3] # NR  
stdres32 <- fitmod$vecres[176:200,3] # IND 
stdres42 <- fitmod$vecres[201:225,3] # ND  
stdres52 <- fitmod$vecres[226:250,3] # SD  

stdres13 <- fitmod$vecres[251:275,3] # SR W L
stdres23 <- fitmod$vecres[276:300,3] # NR  
stdres33 <- fitmod$vecres[301:325,3] # IND 
stdres43 <- fitmod$vecres[326:350,3] # ND  
stdres53 <- fitmod$vecres[351:375,3] # SD  

stdres14 <- fitmod$vecres[376:400,3] # SR W H
stdres24 <- fitmod$vecres[401:425,3] # NR  
stdres34 <- fitmod$vecres[426:450,3] # IND 
stdres44 <- fitmod$vecres[451:475,3] # ND  
stdres54 <- fitmod$vecres[476:500,3] # SD  

# Cleaning a few exceptional cases

stdres11 <- stdres11[!stdres11 > 20] # SR B L
stdres21 <- stdres21[!stdres21 > 20] # NR 
stdres31 <- stdres31[!stdres31 > 20] # IND 
stdres41 <- stdres41[!stdres41 > 20] # ND 
stdres51 <- stdres51[!stdres51 > 20] # SD 

stdres12 <- stdres12[!stdres12 > 20] # SR B H
stdres22 <- stdres22[!stdres22 > 20] # NR  
stdres32 <- stdres32[!stdres32 > 20] # IND 
stdres42 <- stdres42[!stdres42 > 20] # ND  
stdres52 <- stdres52[!stdres52 > 20] # SD  

stdres13 <- stdres13[!stdres13 > 20] # SR W L
stdres23 <- stdres23[!stdres23 > 20] # NR  
stdres33 <- stdres33[!stdres33 > 20] # IND 
stdres43 <- stdres43[!stdres43 > 20] # ND  
stdres53 <- stdres53[!stdres53 > 20] # SD  

stdres14 <- stdres14[!stdres14 > 20] # SR W H
stdres24 <- stdres24[!stdres24 > 20] # NR  
stdres34 <- stdres34[!stdres34 > 20] # IND 
stdres44 <- stdres44[!stdres44 > 20] # ND  
stdres54 <- stdres54[!stdres54 > 20] # SD 

################## plots
par(mfrow=c(2,2))

# stratum B L
mat<-cbind(stdres11,stdres21,stdres31,stdres41,stdres51)
boxplot.matrix(mat,col = c(scales::alpha("darkblue",1), scales::alpha("steelblue3",1),
                           scales::alpha("grey",1), scales::alpha("darkgreen",0.3),
                           scales::alpha("darkgreen",1)), 
               names=c("S Rep", "N Rep", "Ind", "N Demo", "S Demo"), main = "Black-Low",cex.names = 1.1, 
               cex.axis = 1.2, cex.lab = 2, cex.main = 1.5)
# stratum B H
mat<-cbind(stdres12,stdres22,stdres32,stdres42,stdres52)
boxplot.matrix(mat,col = c(scales::alpha("darkblue",1), scales::alpha("steelblue3",1),
                           scales::alpha("grey",1), scales::alpha("darkgreen",0.3),
                           scales::alpha("darkgreen",1)), 
               names=c("S Rep", "N Rep", "Ind", "N Demo", "S Demo"), main = "Black-High",cex.names = 1.1, 
               cex.axis = 1.2, cex.lab = 2, cex.main = 1.5)
# stratum W L
mat<-cbind(stdres13,stdres23,stdres33,stdres43,stdres53)
boxplot.matrix(mat,col = c(scales::alpha("darkblue",1), scales::alpha("steelblue3",1),
                           scales::alpha("grey",1), scales::alpha("darkgreen",0.3),
                           scales::alpha("darkgreen",1)), 
               names=c("S Rep", "N Rep", "Ind", "N Demo", "S Demo"), main = "White-Low",cex.names = 1.1, 
               cex.axis = 1.2, cex.lab = 2, cex.main = 1.5)
# stratum W H
mat<-cbind(stdres14,stdres24,stdres34,stdres44,stdres54)
boxplot.matrix(mat,col = c(scales::alpha("darkblue",1), scales::alpha("steelblue3",1),
                           scales::alpha("grey",1), scales::alpha("darkgreen",0.3),
                           scales::alpha("darkgreen",1)), 
               names=c("S Rep", "N Rep", "Ind", "N Demo", "S Demo"), main = "White-High",cex.names = 1.1, 
               cex.axis = 1.2, cex.lab = 2, cex.main = 1.5)

