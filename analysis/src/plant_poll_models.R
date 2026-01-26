
remove_subset_formula <- function(form){
    char.form <- as.character(form)
    no.sub <-
        gsub("\\| subset\\(Weights[:alpha:⁠]*\\)",
             "", char.form[2])
    form.out <- formula(paste(no.sub, "~", char.form[3]))
    return(form.out)
}


if(site.or.lat ==  "lat"){

    ## flower diversity
    formula.flower.div1 <- formula(scale(MeanFloralDiversity) |
                                  subset(Weights) ~
                                    scale(Lat) + 
                                    (1|Site)
                                  )

    formula.flower.div2 <- formula(scale(MeanFloralDiversity) |
                                     subset(Weights) ~
                                     scale(Cumulative_Precip) + 
                                     scale(Lat) +
                                     (1|Site)
    )
    formula.flower.div3 <- formula(scale(MeanFloralDiversity) |
                                     subset(Weights) ~
                                     scale(Area) +                              
                                     (1|Site)
    )

    ## bee diversity
    formula.bee.div1 <- formula(scale(Net_BeeDiversity) |
                               subset(Weights)~
                                   scale(Lat) +
                                   (1|Site)
                               )
    formula.bee.div2 <- formula(scale(Net_BeeDiversity) |
                                 subset(Weights)~
                                 scale(Cumulative_Precip) +
                                  scale(Lat) +
                                 (1|Site)
    )
    formula.bee.div3 <- formula(scale(Net_BeeDiversity) |
                                 subset(Weights)~
                                 scale(Area) +
                                 (1|Site)
    )
    formula.bee.div4 <- formula(scale(Net_BeeDiversity) |
                                 subset(Weights)~
                                 scale(MeanFloralDiversity) +
                                 scale(Cumulative_Precip) +
                                 scale(Lat) +
                                 scale(Area) +
                                 (1|Site)
    )
    ## bombus abund
    formula.bombus.abund1 <- formula(scale(Net_BombusAbundance) |
                                    subset(Weights)~
                                         scale(Lat) +
                                        (1|Site)
                                    )
    formula.bombus.abund2 <- formula(scale(Net_BombusAbundance) |
                                      subset(Weights)~
                                      scale(Cumulative_Precip) +
                                       scale(Lat) +
                                      (1|Site)
    )
    formula.bombus.abund3 <- formula(scale(Net_BombusAbundance) |
                                      subset(Weights)~
                                      scale(Area) +
                                      (1|Site)
    )
    formula.bombus.abund4 <- formula(scale(Net_BombusAbundance) |
                                      subset(Weights)~
                                      scale(MeanFloralDiversity) +
                                      scale(Cumulative_Precip) +
                                      scale(Lat) +
                                      scale(Area) +
                                      (1|Site)
    )
    ## HB abund
    formula.HB.abund1 <- formula(scale(Net_HBAbundance) |
                                subset(Weights)~
                                   scale(Lat) +
                                    (1|Site)
                                )
    formula.HB.abund2 <- formula(scale(Net_HBAbundance) |
                                  subset(Weights)~
                                  scale(Cumulative_Precip) +
                                  scale(Lat) +
                                  (1|Site)
    )
    formula.HB.abund3 <- formula(scale(Net_HBAbundance) |
                                  subset(Weights)~
                                  scale(Area) +
                                  (1|Site)
    )
    formula.HB.abund4 <- formula(scale(Net_HBAbundance) |
                                  subset(Weights)~
                                  scale(MeanFloralDiversity) +
                                  scale(Cumulative_Precip) +
                                  scale(Lat) +
                                  scale(Area) +
                                  (1|Site)
    )

} else{

    ## flower diversity
    formula.flower.div <- formula(MeanFloralDiversity |
                                  subset(Weights) ~
                                      Cumulative_Precip +
                                      Site
                                  )

    ## bee diversity
    formula.bee.div <- formula(Net_BeeDiversity |
                               subset(Weights) ~
                                   MeanFloralDiversity +
                                   Cumulative_Precip +
                                   Site
                               )
    ## bombus abund
    formula.bombus.abund <- formula(Net_BombusAbundance |
                                    subset(Weights)~
                                        MeanFloralDiversity +
                                        Cumulative_Precip +
                                        Site
                                    )
    ## HB abund
    formula.HB.abund <- formula(Net_HBAbundance |
                                subset(Weights) ~
                                    Site
                                )


}

## **********************************************************
## convert formulas to brms forma
## **********************************************************


bf.fdiv.lat <- bf(formula.flower.div1, family="student")
bf.fdiv.cp <- bf(formula.flower.div2, family="student")
bf.fdiv.a <- bf(formula.flower.div3, family="student")
bf.bombusabund.lat <- bf(formula.bombus.abund1, family="student")
bf.bombusabund.cp <- bf(formula.bombus.abund2, family="student")
bf.bombusabund.a <- bf(formula.bombus.abund3, family="student")
bf.bombusabund.fd <- bf(formula.bombus.abund4, family="student")
bf.HBabund.lat <- bf(formula.HB.abund1, family="student")
bf.HBabund.cp <- bf(formula.HB.abund2, family="student")
bf.HBabund.a <- bf(formula.HB.abund3, family="student")
bf.HBabund.fd <- bf(formula.HB.abund4, family="student")
bf.bdiv.lat <- bf(formula.bee.div1)
bf.bdiv.cp <- bf(formula.bee.div2)
bf.bdiv.a <- bf(formula.bee.div3)
bf.bdiv.fd <- bf(formula.bee.div4)
