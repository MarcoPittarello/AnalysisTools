#' Validation plots for models: lm, glm, glmm, gam,gamm
#'
#' @description 
#' @param model a model object, e.g. lm, glm, lme, glmmTMB, gam, lmer, etc
#' @param residuals.model e.g. residuals(model,type ="pearson") . Any kind of "residuals" 
#' function able to extract residuals from a model
#' @param database.model a data.frame used as data input in the model
#' @param response.var the name of response variable between quotes, e.g. "Y"
#' @param explanatory.var a vector with listed the explanatory variables used in the model, e.g.
#' c("X1","X2","X3")
#' @param coordinate.name a vector with listed the names of the column reporting the coordinates.
#' Default is NULL.
#' @return plots for graphical check of model
#' @import tidyverse, ape, vegan, geoR
#' @examples xx
#' @details  xx
#' @references xx
#' @export

modelCheck_plots<-function(model,residuals.model,database.model,response.var,
                            explanatory.var,coordinate.name=NULL){
  db.check1<-data.frame(
    fitted= fitted(model),
    residuals=residuals.model
  )
  db.check<-cbind(db.check1,
                  database.model[,c(explanatory.var)],
                  database.model[,c(coordinate.name)],
                  database.model[,paste0(response.var,collapse = "'")])
  colnames(db.check)[ncol(db.check)]<-"y"
  
  #normality of residuals
  res<-qplot(sample = residuals, data = db.check)+
    stat_qq_line()+
    stat_qq() +
    theme(axis.line = element_line(colour = "black"))+
    theme(panel.background = element_blank())+
    ylab("Sample Quantiles")+
    xlab("Theoretical Quantiles")+
    ggtitle("Normal Q-Q Plot")
  
  #heterogeneity
  het<-ggplot(db.check,aes(x = fitted,y=residuals))+
    geom_point(shape = 1, colour = "grey40", fill = "white", size = 1, stroke = 0.5)+ 
    geom_smooth(method = "loess")+
    theme(axis.line = element_line(colour = "black"))+
    theme(panel.background = element_blank())+
    ylab("Residuals")+
    xlab("Fitted values")+
    ggtitle("Heteroscedasticity")+
    geom_hline(yintercept=0)
  
  #response vs fitted
  
  resp.vs.fitt<-ggplot(db.check,aes(x = fitted,y=y))+
    geom_point(shape = 1, colour = "grey40", fill = "white", size = 1, stroke = 0.5)+
    geom_smooth(method = "loess")+
    theme(axis.line = element_line(colour = "black"))+
    theme(panel.background = element_blank())+
    ylab("Response")+
    xlab("Fitted values")+
    ggtitle("Response vs. Fitted Value")
  
  # independence 
  i.list<-list()
  
  levels<-unique(explanatory.var)
  
  for (i in levels){
    if (is.numeric(database.model[,i])){
      
      i.list[[i]]<-ggplot(db.check,aes_string(x = "residuals",y=i))+
        geom_point(shape=1)+geom_smooth(method = "loess")+
        theme(axis.line = element_line(colour = "black"))+
        theme(panel.background = element_blank())+
        ylab(noquote(i))+
        xlab("Residuals")+
        ggtitle("Independence")
      
    }
    else {
      i.list[[i]]<-ggplot(db.check,aes_string(x = "residuals",y=i))+
        geom_boxplot()+ 
        theme(axis.line = element_line(colour = "black"))+
        theme(panel.background = element_blank())+
        ylab(noquote(i))+
        xlab("Residuals")+
        ggtitle("Independence")
      
    }
  }
  
  
  if (length(coordinate.name)>0){
    ## spatial autocorrelation
    newdat<-cbind(residuals.model,
                  database.model[,c(coordinate.name)])
    colnames(newdat)[2]<-"x"
    colnames(newdat)[3]<-"y"
    
    dist1 <- dist(cbind(newdat$x, newdat$y))
    
    Moran<-ape::Moran.I(x = residuals.model, w = as.matrix(dist1))
    
    dist3 <- dist(residuals.model)
    Mantel<-vegan::mantel(dist3, dist1)
    
    geodat <- geoR::as.geodata(data.frame(newdat$x, newdat$y
                                          ,residuals.model))
    geodat.v1 <- variog(geodat, uvec=1:1000, max.dist=10000, option='bin')
    
    vario.plot<-data.frame(distance=geodat.v1$u,semivariance=geodat.v1$v)
    variog<-ggplot(vario.plot,aes(x=distance,y=semivariance))+geom_point()+
      ylim(0,max(vario.plot$semivariance))+
      theme(axis.line = element_line(colour = "black"))+
      theme(panel.background = element_blank())+
      ylab("semivariance")+
      xlab("Distance")+
      ggtitle("Spatial autocorrelation")+
      labs(caption =paste0("P-value Moran test: ",
                           round(Moran$p.value,3),
                           " | P-value Mantel test: ",
                           round(Mantel$signif,3)))
    
    i.list1<-list(variog,res,het,resp.vs.fitt)
    i.list2<-append(i.list1,i.list)
    
    ggpubr::ggarrange(plotlist = i.list2)
    
  } else if (is.null(coordinate.name)){
    i.list1<-list(res,het,resp.vs.fitt)
    i.list2<-append(i.list1,i.list)
    
    ggpubr::ggarrange(plotlist = i.list2)
  }
  
  
}


