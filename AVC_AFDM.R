#Vider la mémoire 
rm(list=ls()) 


#Packages utilis?s
library('xlsx')                                   
library('ade4')  
library('glue')   
library('FactoMineR')  
library('factoextra')  
library('PCAmixdata')
library('tidyverse')
library('ggsci')


#Importation des donn?es "AVC_data.csv" 
AVC_donnees <- read.table(file.choose(),fill = TRUE, header = TRUE,sep=",",dec=".", row.names =1,stringsAsFactors = TRUE) 
print(AVC_donnees)
AVC_donnee=na.omit(AVC_donnees) 
#Le jeu de donn?es est bien rang?. Chaque colonne est une variable et chaque ligne est une observation.



#Affichage des donn?es 
View(AVC_donnees) 



#Structure des donn?es 
str(AVC_donnees) 
glue("{ncol(AVC_donnees)} variables et {nrow(AVC_donnees)} observations")



#Affichage des 20 premières lignes 
print(AVC_donnees[1:20,]) 



#Combien de patients ont un AVC ?
data_column <- AVC_donnees %>%
  mutate(
    stroke = ifelse(stroke == 0, 'no stroke', 'stroke'),
    hypertension = ifelse(hypertension == 0, 'no hypertension', 'hypertension'),
    heart_disease  = ifelse(heart_disease == 0, 'no heart disease', 'heart disease')
  )

total_observation <- nrow(data_column)
disease_tables <- data_column %>%
  group_by(stroke) %>%
  dplyr::count() %>% 
  mutate(proportion = n/total_observation * 100)
print(disease_tables)



#La proportion des malades par age et par sexe
total_observation <- nrow(AVC_donnees)
t1 <- AVC_donnees %>%
  group_by(gender, stroke) %>%
  dplyr::count() %>%
  mutate(
    percentage = n/total_observation * 100)

t2 <- AVC_donnees %>%
  group_by(gender, stroke) %>%
  summarise(
    age_mean = mean(age),
    age_sd = sd(age),
    n = n(),
    .groups = 'drop'
  )

t3 <- inner_join(t1, t2)
t3



#Fonction pour centrage-r?duction 
centrage_reduction <- function(x){ 
  n <- length(x) 
  m <- mean(x) 
  v <- (n-1)/n*var(x)  
  return((x-m)/sqrt(v)) 
} 



#Appliquer la fonction sur les variables continues 
AVC_cr_varCont <- data.frame(lapply(subset(AVC_donnees,select=c(2,3,4,8,9,11)),centrage_reduction))  
#print(AVC_cr_varCont) 
summary(AVC_cr_varCont) 



#Codage disjonctif complet 
AVC_disjonctif <- acm.disjonctif(subset(AVC_donnees,select=c(1,5,6,7,10))) 
View(AVC_disjonctif) 



#Fonction pour pond?ration des indicatrices 
fonction_ponderation <- function(x){  
  m <- mean(x)  
  return(x/sqrt(m)) 
} 



#Appliquer la pond?ration sur les indicatrices 
AVC_disjonctif_pond <- data.frame(lapply(AVC_disjonctif,fonction_ponderation)) 
print(AVC_disjonctif_pond) 



#Donn?es transform?es envoy?es ? l?'ACP 
AVC_pour_acp <- cbind(AVC_cr_varCont,AVC_disjonctif_pond)  
rownames(AVC_pour_acp) <- rownames(AVC_donnees)  
print(round(AVC_pour_acp,3)) 



#ACP avec le package ade4 
acp_AVC <- dudi.pca(AVC_pour_acp,center=T, scale=F, scannf=F) 
print(acp_AVC)



#Valeurs propres 
print(round(acp_AVC$eig,5)) 



#Coordonn?es ACP des variables : Gkh.  pour les quanti -> corr?lations avec les facteurs
print(acp_AVC$co[,1:2]) 



#pour les quali -> calculs suppl?mentaires n?cessaires R?cup?rer coord. acp des modalit?s 
moda <- acp_AVC$co[7:21,1:2] 
print(moda) 



#Fr?quence des modalit?s 
freq_moda <- colMeans(AVC_disjonctif) 
print(freq_moda) 
print(AVC_disjonctif) 



#Calcul des moyennes conditionnelles sur les 2 premiers facteurs 
coord_moda <- moda[,1]*sqrt(acp_AVC$eig[1]/freq_moda) 
coord_moda <- cbind(coord_moda,moda[,2]*sqrt(acp_AVC$eig[2]/freq_moda))  
print(coord_moda) 



#Coordonn?es des individus 
print(round(acp_AVC$li[,1:2],5)) 



#Carr? des corr?lations 1er facteur 
r2 <- acp_AVC$co[1:6,1]^2 
print(r2) 



#Carr? du rapport de corr?lation, var. qualitatives 
eta2 <- NULL 
eta2[1] <- sum(acp_AVC$co[7:9,1]^2)  
eta2[2] <- sum(acp_AVC$co[10:11,1]^2)  
eta2[3] <- sum(acp_AVC$co[12:16,1]^2) 
eta2[4] <- sum(acp_AVC$co[17:18,1]^2) 
eta2[5] <- sum(acp_AVC$co[19:21,1]^2) 
print(eta2) 



#Valeurs ? sommer 
criteres <- c(eta2,r2)  
names(criteres) <- colnames(AVC_donnees)  
print(criteres) 



#Crit?re de l'AFDM au 1er facteur 
lambda1 <- sum(criteres) 
print(lambda1) 



#confrontation avec r?sultat (v.p.) de l'ACP sur variables transform?es  au 1er facteur 
print(acp_AVC$eig[1]) 



#----------------AFDM avec FactoMineR-------------------------- 
#Lancement de la proc?dure sur 10000 individus 
afdm_AVC_donnees <- FAMD(AVC_donnees[1:10000,], ncp = 2)



#Lancement de la proc?dure sur tous les individus (analyse g?n?ral) 
afdm_AVC_donnees <- FAMD(AVC_donnees, graph = FALSE)
#Affichage des r?sultats 
print(summary(afdm_AVC_donnees)) 



#visualisation des donn?es 
# Tracage des variables
viz1 <- fviz_famd_var(afdm_AVC_donnees, repel = TRUE)
# Contribution de la première dimension
viz2 <- fviz_contrib(afdm_AVC_donnees, "var", axes = 1)
# Contribution de la deuxième dimension
viz3 <- fviz_contrib(afdm_AVC_donnees, "var", axes = 2)
viz1
viz2
viz3



#############################
quali.var <- get_famd_var(afdm_AVC_donnees, "quali.var") 
quanti.var <- get_famd_var(afdm_AVC_donnees, "quanti.var") 
fviz_famd_var(afdm_AVC_donnee, "quanti.var", col.var = "contrib",  
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
              repel = TRUE) 



##########################################
fviz_famd_var(afdm_AVC_donnees, "quali.var", col.var = "contrib", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
              repel = TRUE 
) 



##########################################"
#Affichage du nuage de donn?es qui represente les differentes variables
famd.stroke <- fviz_mfa_ind(afdm_AVC_donnees, 
                            habillage = "stroke", # color by groups
                            geom = c('point'),
                            palette = pal_simpsons("springfield", alpha = 0.6)(16),
                            title = "Stroke Status"
) 
famd.stroke



famd.smoking <- fviz_mfa_ind(afdm_AVC_donnees, 
                             habillage = "smoking_status", # color by groups
                             geom = c('point'),
                             palette = pal_simpsons("springfield", alpha = 0.6)(16),
                             title = "Smoking Status"
)
famd.smoking



famd.married <- fviz_mfa_ind(afdm_AVC_donnees, 
                             habillage = "ever_married", # color by groups
                             geom = c('point'),
                             palette = pal_simpsons("springfield", alpha = 0.6)(16),
                             title = "Martial Status"
)
famd.married



famd.work_type <- fviz_mfa_ind(afdm_AVC_donnees, 
                               habillage = "work_type", # color by groups
                               geom = c('point'), 
                               palette = pal_simpsons("springfield", alpha = 0.6)(16),
                               title = "Work Status"
)
famd.work_type



famd.residence <- fviz_mfa_ind(afdm_AVC_donnees, 
                               habillage = "Residence_type", # color by groups
                               geom = c('point'), # rempove labels
                               palette = pal_simpsons("springfield", alpha = 0.6)(16), # use a color blind friendly palette 
                               title = "Residence type"
)
famd.residence



famd.heart_disease <- fviz_mfa_ind(afdm_AVC_donnees, 
                                   habillage = "heart_disease", # color by groups
                                   geom = c('point'), 
                                   palette = pal_simpsons("springfield", alpha = 0.6)(16),
                                   title = "Heart Disease Status"
)
famd.heart_disease



famd.hypertension <- fviz_mfa_ind(afdm_AVC_donnees, 
                                  habillage = "hypertension", # color by groups
                                  geom = c('point'),
                                  palette = pal_simpsons("springfield", alpha = 0.6)(16),
                                  title = "Hypertension Status"
)
famd.hypertension

all_FAMD <- (famd.smoking + famd.married) / (famd.work_type + famd.residence) / (famd.heart_disease + famd.hypertension)
plot(all_FAMD)


