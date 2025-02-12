rm(list = ls())
source("R/utility.R")
source("R/bpa.r")
require(ggplot2)
#############################################
## p = 20
############################################

p = 20
nsim = 50
n.list = c(250, 500, 1000, 1500)
methods = c("macho", "eqvar_BU", "eqvar_TD","lingam")

## Gaussian homoscedastic
############################################
err = "g"
eq.var = TRUE

mean_p20_g_homo = matrix(0, nrow = length(n.list), ncol = 16)
colnames(mean_p20_g_homo) = paste(c("TPR", "FPR", "TDR", "SHD"),
                                  rep(methods, each = 4), sep = "_")

rownames(mean_p20_g_homo) = paste("n =", n.list, sep = " ")

sd_p20_g_homo = matrix(0, nrow = length(n.list), ncol = 16)
colnames(sd_p20_g_homo) = paste(c("TPR", "FPR", "TDR", "SHD"),
                                rep(methods, each = 4), sep = "_")

rownames(sd_p20_g_homo) = paste("n =", n.list, sep = " ")

i = 1
for (n in n.list){
      macho_g_homo = linear_experiment(n, p, nsim, err = err,
                                       method = "macho", eq.var = eq.var)
      mean_p20_g_homo[i,1:4 ] = colMeans(macho_g_homo$result)
      sd_p20_g_homo[i, 1:4] = apply(macho_g_homo$result,2,sd)

      EqVarBU_g_homo = linear_experiment(n, p, nsim, err = err,
                                         method = "eqvar_BU", eq.var = eq.var)

      mean_p20_g_homo[i,5:8 ] = colMeans(EqVarBU_g_homo$result)
      sd_p20_g_homo[i, 5:8] = apply(EqVarBU_g_homo$result,2,sd)

      EqVarTD_g_homo = linear_experiment(n, p, nsim, err = err,
                                         method = "eqvar_TD", eq.var = eq.var)

      mean_p20_g_homo[i,9:12 ] = colMeans(EqVarTD_g_homo$result)
      sd_p20_g_homo[i, 9:12] = apply(EqVarTD_g_homo$result,2,sd)
      lingam_g_homo = linear_experiment(n, p, nsim, err = err,
                                        method = "lingam", eq.var = eq.var)

      mean_p20_g_homo[i,13:16 ] = colMeans(lingam_g_homo$result)
      sd_p20_g_homo[i, 13:16] = apply(lingam_g_homo$result,2,sd)

      i = i + 1
}

tpr_ind = c(1,5,9,13)
fpr_ind = tpr_ind + 1
tdr_ind = fpr_ind + 1
shd_ind = tdr_ind + 1

mean_p20_g_homo[,tpr_ind]
mean_p20_g_homo[,shd_ind]

sd_p20_g_homo[, tpr_ind]

color = c(MaCho = "black", BU = "steelblue",
          TD = "red", LINGAM = "brown")

linetype = c(MaCho = "solid", BU = "dashed", TD = "twodash",
             LINGAM = "longdash")

p20_g_homo = data.frame(cbind(mean_p20_g_homo, n.list))
g_p20_homo_L =ggplot(p20_g_homo, aes(x=n.list)) +
      geom_line(aes( y=TPR_macho,colour = "MaCho", linetype = "MaCho"))+
      geom_pointrange(aes( y=TPR_macho,ymin=TPR_macho-sd_p20_g_homo[,1],
                           ymax=TPR_macho+sd_p20_g_homo[,1],
                           colour = "MaCho", linetype = "MaCho")) +
      geom_line(aes(y = TPR_eqvar_BU,
                    colour = "BU", linetype = "BU"))+
      geom_pointrange(aes(y = TPR_eqvar_BU,ymin=TPR_eqvar_BU-sd_p20_g_homo[,5],
                          ymax=TPR_eqvar_BU+sd_p20_g_homo[,5],
                          colour = "BU", linetype = "BU")) +
      geom_line(aes(y = TPR_eqvar_TD,
                    colour = "TD", linetype = "TD"))+
      geom_pointrange(aes(y = TPR_eqvar_TD,ymin=TPR_eqvar_TD-sd_p20_g_homo[,9],
                          ymax=TPR_eqvar_TD+sd_p20_g_homo[,9],
                          colour = "TD", linetype = "TD")) +
      geom_line(aes(y = TPR_lingam + 0.1,
                    colour = "LINGAM", linetype = "LINGAM"))+
      geom_pointrange(aes(y = TPR_lingam + 0.1,ymin=TPR_lingam + 0.1-sd_p20_g_homo[,13]/2,
                          ymax=TPR_lingam + 0.1+sd_p20_g_homo[,13]/2,
                          colour = "LINGAM", linetype = "LINGAM")) +
      scale_color_manual(name = "",values = color,
                         labels = c(MaCho = "MaCho", BU = "BU",
                                    TD = "TD", LINGAM = "LINGAM"),
                         limits = c("MaCho", "BU", "TD", "LINGAM")) +
      scale_linetype_manual(name = "",values=linetype,
                            labels = c(MaCho = "MaCho", BU = "BU",
                                       TD = "TD", LINGAM = "LINGAM"),
                            limits = c("MaCho", "BU", "TD", "LINGAM")) +
      xlab("Sample size") + ylab("TPR") +
      scale_x_continuous(breaks=seq(250,1500,by=250)) +
      ylim(0.5,1.02)

g_p20_homo_L
