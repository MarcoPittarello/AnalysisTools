#' Find the best spatial correlation structure in GLMM run with 'lme' function in 'nlme' package.
#'
#' @description When through the function (\link[AnalysisTool]{modelCheck_plots}) spatial autocorrelation is found in the residuals of a model performed
#' with (\link[nlme]{gls}) or (\link[nlme]{lme}), it may be necessary to specify the structure of this spatial correlation through the 'correlation' argument.\cr
#' Among the different correlation structures (see (\link[nlme]{CorStruct})), those supported by the present function are:
#'  - corGaus	 --> Gaussian spatial correlation
#'  - corExp	 --> Exponential spatial correlation
#'  - corSpher --> Spherical spatial correlation
#'  - corLin   --> Linear spatial correlation.\cr
#'  
#'  See @examples  section for operative details
#'  
#' @param data a data.frame used as data input in the model
#' @param X column within the dataset containing the coordinates X (e.g.: dataset$X)
#' @param Y column within the dataset containing the coordinates Y (e.g.: dataset$Y)
#' @param var.res column within the dataset containing the response variable or model residuals obtained with residuals(model,type="pearson") see @examples  
#' @param range value of the 'range' identifiable from the semivariogram obtained with (\link[AnalysisTool]{modelCheck_plots}). 
#' @param sill value of the 'sill' identifiable from the semivariogram obtained with (\link[AnalysisTool]{modelCheck_plots}). 
#' @param nugget value of the 'nugget' identifiable from the semivariogram obtained with (\link[AnalysisTool]{modelCheck_plots}). 
#' @param spatial.formula a one sided formula of the form ~ S1+...+Sp, or ~ S1+...+Sp | g, specifying spatial covariates S1 through Sp and, optionally, a grouping factor g. For more details see the 'form' argument of e.g. (\link[nlme]{corGaus}).
#' @return Linear, exponential, gaussian and spherical models plotted on the semivariogram.\cr 
#' It is also provided in the console the string to copy and paste in 'correlation' argument in the functions: (\link[nlme]{gls}), (\link[nlme]{lme}) of 'nlme' package"
#' @import vegan geoR nlme
#' @examples 
#' 1- define from (\link[AnalysisTool]{modelCheck_plots}) which are the values of range, sill, and nugget. 
#' For a definition of range, sill, and nugget terminology see \href{https://sakai.unc.edu/access/content/group/2842013b-58f5-4453-aa8d-3e01bacbfc3d/public/Ecol562_Spring2012/images/lectures/lecture31/fig1d.png}{here} \cr
#' 
#' 2- run 'best_Spatial_CorrStruct' function:
#' best_Spatial_CorrStruct(data=Dataset, 
#' X = Dataset$X,
#' Y = Dataset$Y,
#' var.res = residuals(model,type="pearson"),
#' range = 5,
#' sill = 5,
#' nugget = 4,
#' spatial.formula = "~ X + Y")
#' @references 
#' 
#' * Gałecki, A., & Burzykowski, T. (2013). Linear mixed-effects model. In Linear mixed-effects models using R (pp. 245-273). Springer, New York, NY.
#' * \href{https://sakai.unc.edu/access/content/group/2842013b-58f5-4453-aa8d-3e01bacbfc3d/public/Ecol562_Spring2012/docs/lectures/lecture32.htm}{Jack Weiss tutorial (2012). Fitting theoretical models to the empirical semivariogram}
#' 
#' @export

best_Spatial_CorrStruct<-function(data,X,Y,var.res,range, sill,nugget,spatial.formula){
  newdat <- data.frame(X,Y,var.res)
  geodat <- as.geodata(newdat)
  dist1 <- dist(cbind(X, Y))
  geodat.v1 <- variog(geodat)
  plot(geodat.v1,main="Semivariogram")
  
  Dataset<-data
  
  partial.sill<-sill-nugget#sill-nugget
  form.formula<-as.formula(spatial.formula)
  ### correlazione spaziale lineare -----------
  correlation_structure_lin <- corLin(form =form.formula, nugget=T, metric="euclidean", fixed=F)
  correlation_structure_lin_Iniz <- Initialize(correlation_structure_lin, data = Dataset)
  ols1 <- variofit(geodat.v1, ini=c(partial.sill,range), nugget=nugget, cov.model='linear')
  
  ### correlazione spaziale esponenziale  -----------
  correlation_structure_Exp <- corExp(form = form.formula, nugget=T, metric="euclidean", fixed=F)
  correlation_structure_Exp_Iniz <- Initialize(correlation_structure_Exp, data = Dataset)
  ols2 <- variofit(geodat.v1, ini=c(partial.sill,range), nugget=nugget, cov.model='exponential')
  
  
  ### correlazione spaziale gaussiana  -----------
  correlation_structure_Gaus <- corGaus(form = form.formula, nugget=T, metric="euclidean", fixed=F)
  correlation_structure_Gaus_Iniz <- Initialize(correlation_structure_Gaus, data = Dataset)
  ols3 <- variofit(geodat.v1, ini=c(partial.sill,range), nugget=nugget, cov.model='gaussian')
  
  ### correlazione spaziale sferica  -----------
  correlation_structure_Spher <- corSpher(form =form.formula, nugget=T, metric="euclidean", fixed=F)
  correlation_structure_Spher_Iniz <- Initialize(correlation_structure_Spher, data = Dataset)
  ols4 <- variofit(geodat.v1, ini=c(partial.sill,range), nugget=nugget, cov.model='spherical')
  
  lines(ols1, col='gray20',lwd = 1.2)
  lines(ols2, col="red",lwd = 1.2)
  lines(ols3, col="green",lwd = 1.2)
  lines(ols4, col="blue",lwd = 1.2)
  legend('bottomright', c('Linear','Exponential', 'Gaussian', 'Spherical'), 
         col=c('gray20',"red","green","blue"), lty=rep(1,4), cex=1, bty='n')
  
  print("                                                                              ")
  print("                                                                              ")
  
  print("Based on the fits of the models (see figure), copy the strings of the best spatial correlation model, execute them, and copy the object 'corr.structure' into the 'correlation' argument in the functions 'gls','lme', etc. of the 'nlme' package")
  print("                                                                              ")
  print("### LINEAR -----------------------------------------------------------------")
  print(paste("correlation_structure_lin <- corLin(form =",spatial.formula, ",nugget=T, metric='euclidean', fixed=F)"))
  print("corr.structure <- Initialize(correlation_structure_lin, data = )")
  
  print("                                       ")
  print("### EXPONENTIAL -----------------------------------------------------------------")
  print(paste("correlation_structure_Exp <- corLin(form =",spatial.formula, ",nugget=T, metric='euclidean', fixed=F)"))
  print("corr.structure <- Initialize(correlation_structure_Exp, data = )")
  
  print("                                                                              ")
  print("### GAUSSIAN -----------------------------------------------------------------")
  print(paste("correlation_structure_Gaus <- corLin(form =",spatial.formula, ",nugget=T, metric='euclidean', fixed=F)"))
  print("corr.structure <- Initialize(correlation_structure_Gaus, data = )")
  
  print("                                                                              ")
  print("### SPHERICAL -----------------------------------------------------------------")
  print(paste("correlation_structure_Spher <- corLin(form =",spatial.formula, ",nugget=T, metric='euclidean', fixed=F)"))
  print("corr.structure <- Initialize(correlation_structure_Spher, data = )")
  
}
