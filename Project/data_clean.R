rm(list=ls())

library(raster)
library(foreign)
library(dplyr)
library(INSP)
library(ggplot2)
library(INLA)
library(GGally)
library(tidyr)
library(surveillance)
library(ape)
library(spdep)
library(knitr)
library(INLA)

setwd("~/Documents/Classes/MPSpatialDemog/Project")

read_vital <- function(x){
    f_ <- paste0("~/Documents/MXU5MR/nacimientos/data/inegi/NACIM", x,".dbf")
    read.dbf(f_, as.is=TRUE)
}

tidymap <- fortify(mx.sp.df) %>% mutate(id=as.numeric(id)) %>%
    left_join(mx.sp.df@data %>% mutate(id=0:(nrow(mx.sp.df)-1)))

base_plot <- function(datur, var_, custom=NA, continuous=F, inverse=F, 
                      interval=T, bins=8){
    if(interval){
        datur$QVAR <- cut_interval(datur[,var_], bins)
    }
    else{
        datur$QVAR <- cut_number(datur[,var_], bins)
    }
    if(!is.na(custom[1])){
        datur$QVAR <- cut(datur[,var_], custom)
    }
    datur$VAR <- datur[,var_]
    inv <- ifelse(inverse, 1, -1)
    bmap <- tidymap %>% left_join(datur) %>%
        ggplot(aes(x=long, y=lat)) +
        theme_classic()
    if(!continuous){
        bmap <- bmap + geom_polygon(aes(group=group, fill = QVAR)) +
            scale_fill_brewer(palette = "Spectral", direction=inv) +
            scale_color_brewer(palette = "Spectral", direction=inv) + 
            theme(axis.line = element_blank(),
                  legend.justification=c(1,2.9),legend.position=c(.3, .95),
                  legend.title=element_blank(), axis.text=element_blank(),
                  axis.ticks=element_blank(), axis.title=element_blank(),
                  legend.key.size=unit(.6, "in"),
                  legend.text=element_text(size=20),
                  title=element_text(size=26))
    }
    else{
        bmap <- bmap + 
            geom_polygon(aes(group=group, fill = VAR)) +
            scale_fill_distiller(palette = "Spectral", direction=inv) + 
            theme(axis.line = element_blank(),
                  legend.justification=c(.8,3.7),legend.position=c(.3, .95),
                  legend.title=element_blank(), axis.text=element_blank(),
                  axis.ticks=element_blank(), axis.title=element_blank(),
                  legend.key.size=unit(.8, "in"),
                  legend.text=element_text(size=20),
                  title=element_text(size=26))
    }
    return(bmap + geom_path(aes(group=group), size=.1))
}

f2n <- function(x){
    as.numeric(as.character(x))
}

DFpop <- "./data/ITER_NALDBF10.dbf" %>%
    read.dbf(as.is=TRUE) %>% 
    filter(LOC=="0000" & MUN != "000") %>%
    select(POBTOT, POBFEM, P_12A14_F, P_15A17_F, P_18A24_F, P_15A49_F, PHOG_IND,
           P_15YMAS_F, ENTIDAD, MUN, P15YM_AN_F, P_12YMAS_F, P_12YMAS_M,
           POCUPADA_F, POCUPADA_M, P_0A2, P_3A5) %>%
    mutate(P_15A49_F=f2n(P_15A49_F), P_15A17_F=f2n(P_15A17_F), 
           P_18A24_F=f2n(P_18A24_F), POBTOT=f2n(POBTOT), POBFEM=f2n(POBFEM),
           P_15YMAS_F=f2n(P_15YMAS_F), PHOG_IND=f2n(PHOG_IND), 
           P15YM_AN_F=f2n(P15YM_AN_F), P_0A2=f2n(P_0A2), P_3A5=f2n(P_3A5),
           P_12YMAS_F=f2n(P_12YMAS_F), P_12YMAS_M=f2n(P_12YMAS_M),
           POCUPADA_F=f2n(POCUPADA_F), POCUPADA_M=f2n(POCUPADA_M),
           P_12A14_F=f2n(P_12A14_F)) %>%
    mutate(P_25A49_F=P_15A49_F - P_15A17_F - P_18A24_F) %>% # New age group
    mutate(Pr_IND=PHOG_IND/POBTOT) %>% # percentage of pop that is indeginous
    mutate(Pr_LIT=1 - (P15YM_AN_F/P_15YMAS_F)) %>% # percentage literate
    mutate(Pr_WORK_F=POCUPADA_F/P_12YMAS_F) %>% # working Female Pr
    mutate(Pr_WORK_M=POCUPADA_M/P_12YMAS_M) %>% # working male Pr
    mutate(P_CHILD=P_0A2 + P_3A5) %>% # Number of Children
    mutate(CWR=P_CHILD/P_15A49_F) %>% # Child Woman Ratio (traditional)
    mutate(GEOID=paste0(ENTIDAD, MUN))

png("./plots/CWR_discrete.png", width = 560*3, height = 480*3)
base_plot(DFpop, "CWR") + 
    labs(title="Child Woman Ratio 2010")
dev.off()

png("./plots/CWR_continuous.png", width = 560*3, height = 480*3)
base_plot(DFpop, "CWR", continuous=T) + 
    labs(title="Child Woman Ratio 2010")
dev.off()

png("./plots/indigenous.png", width = 560*3, height = 480*3)
base_plot(DFpop, "Pr_IND") + 
    labs(title="Indigenous Population Percentage 2010")
dev.off()

png("./plots/literacy.png", width = 560*3, height = 480*3)
base_plot(DFpop, "Pr_LIT", inverse=T) + 
    labs(title="Literate Population Percentage 2010")
dev.off()

png("./plots/working_female.png", width = 560*3, height = 480*3)
base_plot(DFpop, "Pr_WORK_F", inverse=T) + 
    labs(title="Percentage Working Female 2010")
dev.off()

png("./plots/working_male.png", width = 560*3, height = 480*3)
base_plot(DFpop, "Pr_WORK_M", inverse=T) + 
    labs(title="Percentage Working Male 2010")
dev.off()

simpreg <- glm(
    P_CHILD ~ Pr_LIT + Pr_WORK_F + Pr_WORK_M + Pr_IND + offset(log(P_15A49_F)),
    family=poisson, data=DFpop)

summary(simpreg)

png("./plots/simple_model_residuals_map.png", width = 560*3, height = 480*3)
base_plot(DFpop %>% mutate(resid=simpreg$residuals), "resid") + 
    labs(title="Simple Model Residual Plot")
dev.off()

png("./plots/simple_model_residuals_hist.png", width = 560*3, height = 480*3)
DFpop %>% mutate(Res=simpreg$residuals, ResQ=cut_interval(Res, 8)) %>%
    ggplot(aes(x=Res, fill=ResQ)) + 
    geom_histogram(bins=50) + 
    scale_fill_brewer(palette = "Spectral", direction=-1, guide=F) + 
    theme_minimal() + 
    labs(title="Distribution of 2010 CWR Residuals", y="Count",
         x="") +
    theme(title=element_text(size=26), axis.text=element_text(size=20))
dev.off()

DFbirths <- bind_rows(read_vital(10), read_vital(11)) %>% 
    filter(ANO_NAC==2010 & EDAD_MADN !=99) %>% 
    select(ENT_RESID, MUN_RESID, EDAD_MADN) %>%
    mutate(GEOID=paste0(sprintf("%02d", ENT_RESID), sprintf("%03d", MUN_RESID)))

png("./plots/mother_age_hist.png", width = 560*3, height = 480*3)
DFbirths %>% ggplot(aes(x=EDAD_MADN)) + 
    geom_bar() + 
    labs(title="Age of Mother Giving Birth 2010", x="Age", y="Count") +
    theme(title=element_text(size=26), axis.text=element_text(size=20))
dev.off()

png("./plots/age_map.png", width = 560*3, height = 480*3)
DFbirths %>% group_by(GEOID) %>% summarize(AVG_MADN=mean(EDAD_MADN)) %>%
    as.data.frame %>%
    base_plot("AVG_MADN") + 
    labs(title="Average Age of Mother")
dev.off()

png("./plots/age_map_var.png", width = 560*3, height = 480*3)
DFbirths %>% group_by(GEOID) %>% summarize(AVG_MADN_VAR=sd(EDAD_MADN)) %>%
    as.data.frame %>%
    base_plot("AVG_MADN_VAR", interval=F) + 
    labs(title="SD Age of Mother")
dev.off()

my_fn <- function(data, mapping, ...){
    p <- ggplot(data = data, mapping = mapping) + 
        geom_point() + 
        geom_smooth(method=loess, fill="red", color="red", ...) +
        geom_smooth(method=lm, fill="blue", color="blue", ...)
    p
}

png("./plots/cov_corr_plot.png", width = 560*2.1, height = 480*2.1)
DFbirths %>% group_by(GEOID) %>% 
    summarize(AVG_MADN_VAR=sd(EDAD_MADN), AVG_MADN=mean(EDAD_MADN)) %>%
    as.data.frame %>% left_join(DFpop) %>% 
    select(Pr_LIT, Pr_WORK_F, Pr_WORK_M, Pr_IND, 
           CWR, AVG_MADN, AVG_MADN_VAR) %>%
    ggpairs(lower = list(continuous = my_fn))
dev.off()

temp_ <- DFpop %>% 
    select(GEOID, P_12A14_F, P_15A17_F, P_18A24_F, P_25A49_F) %>%
    gather(Age, Population, -GEOID) %>% arrange(GEOID, Age) %>% 
    left_join(DFpop %>% select(Pr_LIT, Pr_WORK_F, Pr_WORK_M, Pr_IND, GEOID))

DF <- DFbirths %>% 
    mutate(Age=cut(EDAD_MADN, c(11,14,17,24,49), labels=unique(temp_$Age))) %>%
    filter(!is.na(Age)) %>%
    group_by(GEOID, Age) %>% summarize(Births=n()) %>% as.data.frame %>%
    right_join(temp_) %>% mutate(Births=ifelse(!is.na(Births), Births, 0)) %>%
    mutate(RawFR=Births/Population) %>%
    mutate(CHPOAX=strtrim(GEOID,2) == "20" | strtrim(GEOID,2) == "07")

DF %>% group_by(Age) %>% summarize(P=sum(Population), B=sum(Births)) %>%
    mutate(R=B/P) %>% pander

for(a in unique(DF$Age)[2:4]){
    p <- DF %>% filter(Age==a) %>%
        base_plot("RawFR", interval=FALSE) + 
        labs(title=paste0("Raw Fertility Rate Age: ", a))
    ggsave(paste0("./plots/RawFR_", a, ".png"), width=444.5/.9, height=381/.9, 
           plot=p, units="mm")
}

nullModel <- glm(Births ~ Age, family=poisson, offset=log(Population), data=DF)

glmModel <- glm(Births ~ Age + Pr_LIT + Pr_WORK_F + Pr_WORK_M + Pr_IND,
                family=poisson, offset=log(Population), data=DF)

# Estimation of deviance explained by age

1 - (nullModel$deviance / nullModel$null.deviance)

# Deviance not explained by age
nullModel$deviance / nullModel$null.deviance

# total Deviance explained by covariates and age
1 - (glmModel$deviance / glmModel$null.deviance)

# remainder residual explained
1 -  glmModel$deviance/nullModel$deviance

summary(glmModel)$coefficients %>% round(4) %>% kable("markdown")

cutvalsBasic <- c(-10, -.75, -.5, -.25, 0, .25, .5, .75, 40)

for(a in unique(DF$Age)[1:4]){
    p <- DF %>% mutate(Residuals=glmModel$residuals) %>% filter(Age==a) %>%
        base_plot("Residuals", custom=cutvalsBasic) + 
        labs(title=paste0("Residuals Age: ", a))
    ggsave(paste0("./plots/Residuals_", a, ".png"), width=444.5/.9, 
           height=381/.9, plot=p, units="mm")
}

graph <- poly2adjmat(mx.sp.df)
for(a in unique(DF$Age)[1:4]){
    Moran.I((mx.sp.df@data %>% left_join(
        DF %>% mutate(Rez=glmModel$residuals) %>% filter(Age==a)))$Rez, 
        graph) %>% print
}

MEXlistw <- poly2nb(mx.sp.df, queen=TRUE) %>% nb2listw

DFrez <- bind_rows(lapply(unique(DF$Age), function(a){
    DF %>% mutate(Residuals=glmModel$residuals) %>% 
        filter(Age==a) %>%
        right_join(mx.sp.df@data) %>% 
        mutate(pLocal=localmoran(Residuals, MEXlistw)[,"Pr(z > 0)"]) %>% 
        mutate(pSig=pLocal < .05) %>%
        mutate(laggedResidual=lag.listw(MEXlistw, Residuals)) %>% 
        mutate(Neighbors=sapply(MEXlistw$weights, length)) %>%
        mutate(Cluster=ifelse(
            Residuals < 0 & laggedResidual < 0, "Low", "High")) %>%
        mutate(Cluster=ifelse(pSig, Cluster, NA)) %>%
        mutate(pAlpha=ifelse(pSig, .4, .35))
}))


DFrez %>%
    ggplot(aes(x=Residuals, y=laggedResidual, alpha=pAlpha,
               size=log(Population), color=Cluster)) + 
    geom_point() + 
    xlim(c(-1, 1)) + ylim(c(-1, 1)) + 
    labs(title="Lisa Examination of Clusters", 
         x="Residuals", y="Lagged Residuals") + 
    scale_alpha_continuous(guide=F) + 
    geom_hline(yintercept=0, linetype=2, alpha=.6) + 
    geom_vline(xintercept=0, linetype=2, alpha=.6) + 
    facet_wrap(~Age)
ggsave("./plots/LisaResiduals.png")


DFrez %>% left_join(tidymap, by="GEOID") %>% 
    ggplot(aes(x=long, y=lat)) +
    theme_classic() + 
    geom_polygon(aes(group=group, fill = Cluster)) + 
    facet_wrap(~Age) + 
    theme(axis.line = element_blank(),
          legend.title=element_blank(), axis.text=element_blank(),
          axis.ticks=element_blank(), axis.title=element_blank(),
          title=element_text(size=26)) + 
    labs(title="Significant Cluster By Location")
ggsave("./plots/LisaResidualsMap.png")

DFinla <- DFrez %>% 
    mutate(Residuals=NULL, laggedResidual=NULL, 
           Cluster=NULL, pLocal=NULL, pSig=NULL) %>%
    mutate(ID1=rep(1:nrow(mx.sp.df@data), 4), ID2=ID1, ID3=ID1) %>%
    mutate(age15=Age=="P_15A17_F", age18=Age=="P_18A24_F",
           age25=Age=="P_25A49_F")

funcform <- Births ~ Age + Pr_LIT + Pr_WORK_F + Pr_WORK_M + Pr_IND + 
    f(ID1, age15, model="besag", graph=graph, 
      hyper=list(prec=list(prior="loggamma", param=c(1, 10)))) + 
    f(ID2, age18, model="besag", graph=graph, 
      hyper=list(prec=list(prior="loggamma", param=c(1, 10)))) +
    f(ID3, age25, model="besag", graph=graph, 
      hyper=list(prec=list(prior="loggamma", param=c(1, 10))))

funcform2 <- Births ~ Age + Pr_LIT + Pr_WORK_F + Pr_WORK_M + Pr_IND + 
    f(ID1, age15, model="bym", graph=graph, 
      hyper=list(
          prec.spatial=list(prior="loggamma", param=c(1, 10)),
          prec.unstruct=list(prior="loggamma", param=c(1, 20)))) + 
    f(ID2, age18, model="bym", graph=graph, 
      hyper=list(
          prec.spatial=list(prior="loggamma", param=c(1, 10)),
          prec.unstruct=list(prior="loggamma", param=c(1, 20)))) +
    f(ID3, age25, model="bym", graph=graph, 
      hyper=list(
          prec.spatial=list(prior="loggamma", param=c(1, 10)),
          prec.unstruct=list(prior="loggamma", param=c(1, 20))))

inlaModel <- inla(funcform, family="poisson", E=Population, data=DFinla,
                  control.predictor=list(link = 1))
#inlaModel2 <- inla(funcform2, family="poisson", E=Population, data=DFinla)

cVars <- c("Pr_LIT", "Pr_WORK_F", "Pr_WORK_M", "Pr_IND")
betas <- matrix(inlaModel$summary.fixed[cVars, "mean"], nrow=4)

cutvals <- c(-1, -.75, -.5, -.25, 0, .25, .5, .75, 1.25)

DFeffects <- bind_rows(
    DFinla %>% filter(age15==1) %>% 
        select(Pr_LIT, Pr_WORK_F, Pr_WORK_M, Pr_IND) %>% as.matrix %>%
        `%*%`(betas) %>% c %>% `-`(., mean(.)) %>%
        mutate(mx.sp.df@data, Param="Fixed", Values=.),
    mx.sp.df@data %>% 
        mutate(Param="RE1", Values=inlaModel$summary.random$ID1$mean,
               ParamSD=inlaModel$summary.random$ID1$sd),
    mx.sp.df@data %>% 
        mutate(Param="RE2", Values=inlaModel$summary.random$ID2$mean,
               ParamSD=inlaModel$summary.random$ID2$sd),
    mx.sp.df@data %>% 
        mutate(Param="RE3", Values=inlaModel$summary.random$ID3$mean,
               ParamSD=inlaModel$summary.random$ID3$sd)) %>%
    mutate(ValueQ=cut(Values , cutvals))

summary(glmModel)$coefficients %>% kable("markdown")

summary(inlaModel)$fixed %>% as.data.frame %>% select(-kld, -mode) %>%
    kable("markdown")

hist(DFinla$Births)

DFeffects %>% filter(Param=="RE2") %>%
    left_join(filter(DFinla, age18) %>% select(GEOID, Population)) %>%
    ggplot(aes(x=log(Population), y=ParamSD)) + geom_point() + 
    labs(x="Log Population", y="RE Std. Dev.", 
         title="Random Effects By Population Size") + 
    theme(title=element_text(size=26))
ggsave("./plots/REandPopSize.png")

DFeffects %>% ggplot(aes(x=Values, fill=ValueQ)) +
    geom_histogram() + 
    scale_fill_brewer(palette="Spectral", direction=-1) +
    facet_wrap(~Param)
ggsave("./plots/ParamEffectDistribution.png")

DFeffects %>% select(GEOID, Param, Values) %>%
    reshape(idvar='GEOID', timevar='Param', direction='wide') %>% 
    select(-GEOID) %>% ggpairs + 
    labs(title="Correlation of Parameter Effects") + 
    theme(title=element_text(size=26))
ggsave("./plots/ParamEffectCorrelation.png")

tidymap %>% left_join(DFeffects, by="GEOID") %>%
    ggplot(aes(x=long, y=lat)) +
    theme_classic() + 
    geom_polygon(aes(group=group, fill = ValueQ)) + 
    facet_wrap(~Param) + 
    scale_fill_brewer(palette="Spectral", direction=-1) +
    theme(axis.line = element_blank(),
          legend.title=element_blank(), axis.text=element_blank(),
          axis.ticks=element_blank(), axis.title=element_blank(),
          title=element_text(size=26)) + 
    labs(title="Parameter Effects (Zero Centered)")
ggsave("./plots/ParamEffectMap.png")

png("./plots/inla_estimate15.png", width = 560*3, height = 480*3)
DFinla %>% mutate(pred=inlaModel$summary.fitted.values$mean) %>%
    filter(age15) %>% base_plot("pred") + 
    labs(title="Estimated Fertility Rates: 15-17")
dev.off()

png("./plots/inla_estimate18.png", width = 560*3, height = 480*3)
DFinla %>% mutate(pred=inlaModel$summary.fitted.values$mean) %>%
    filter(age18) %>% base_plot("pred") + 
    labs(title="Estimated Fertility Rates: 18-24")
dev.off()

png("./plots/inla_estimate25.png", width = 560*3, height = 480*3)
DFinla %>% mutate(pred=inlaModel$summary.fitted.values$mean) %>%
    filter(age25) %>% base_plot("pred") + 
    labs(title="Estimated Fertility Rates: 25-49")
dev.off()

### Interactive map of final mean estimates for middle age group
mx.sp.df@data$data <- DFinla %>% 
    mutate(pred=inlaModel$summary.fitted.values$mean) %>%
    filter(age18) %>% .$pred

# find that crazy outlier
spdf2leaf(mx.sp.df, label="Fertility")

# Whats its GEOID
mx.sp.df@data %>% filter(NOM_MUN=="Coyame del Sotol")

# whats the pop size for dat guy
DFinla %>% filter(GEOID=="08015")

# what about RE uncertainty?
DFeffects %>% filter(GEOID=="08015")
