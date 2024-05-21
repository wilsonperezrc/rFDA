library(fda.usc)

#####################################################################################
# Function that obtains the FPCs given a % variability
#   Parameters: 
#     - obj_fd: Functional object
#     - porcentaje_explica: % of variability to be explained with n FPCs
#   Returns:
#     - obj_fpc: Resulting FPC object
#####################################################################################
get_harmonics_fpc<-function(obj_fd, porcentaje_explica){
  obj_fpc <- pca.fd(obj_fd, nharm=dim(obj_fd$coefs)[1])
  num_harmonic_calculation<-0
  obj_fpc$values
  for (i in obj_fpc$varprop) {
    if(sum(obj_fpc$varprop[1:(num_harmonic_calculation+1)])<porcentaje_explica)
      num_harmonic_calculation<-num_harmonic_calculation+1
  }
  num_harmonic_calculation<-num_harmonic_calculation+1
  print(paste('Num components:',num_harmonic_calculation))
  print(paste('Sum FPCA:',sum(obj_fpc$varprop[1:num_harmonic_calculation])))
  
  obj_fpc <- pca.fd(obj_fd, nharm=num_harmonic_calculation)
  obj_fpc$sum_harmonics<-sum(obj_fpc$varprop[1:(num_harmonic_calculation)])
  return(obj_fpc)
}

#####################################################################################
# Function that obtains the FPCs under the restriction LAMBDA > 1
#   Parameters: 
#     - obj_fd: Functional object
#   Returns:
#     - obj_fpc: Resulting FPC object
#####################################################################################
get_harmonics_fpc_lambda<-function(obj_fd){
  obj_fpc <- pca.fd(obj_fd, nharm=25)
  num_harmonic_calculation<-length(which (obj_fpc$values[c(1:25)]>1))
  obj_fpc <- pca.fd(obj_fd, nharm=num_harmonic_calculation)
  obj_fpc$sum_harmonics<-sum(obj_fpc$varprop[1:(num_harmonic_calculation)])
  return(obj_fpc)
}

#####################################################################################
# Function that obtains the 2 FPCs
#   Parameters: 
#     - obj_fd: Functional object
#   Returns:
#     - obj_fpc: Resulting FPC object
#####################################################################################
get_harmonics_2C<-function(obj_fd){
  pc_obj <- pca.fd(obj_fd, nharm=2)
  pc_obj$sum_harmonics<-sum(pc_obj$varprop[1:2])
  
  print(paste('Num components:',2))
  print(paste('Sum FPCA-2:',sum(pc_obj$varprop[1:2])))
  
  return(pc_obj)
}

#####################################################################################
# Function that obtains the outlier observations using the Mahalanobis method
#   Parameters: 
#     - pc_obj : PC object
#     - p_value: Confidence level to consider the observations as outliers
#     - obj_fd : Functional object
#     - metodo : Mahalanobis method to be executed
#               -> b:  Classic Mahalanobis
#               -> r:  Robust Mahalanobis
#               -> 50: Robust Mahalanobis with the 50th quartile
#               -> 75: Robust Mahalanobis with the 75th quartile
#   Returns:
#     - pos_out: Vector of positions of the outlier observations
#####################################################################################
get_atypical_mahalnobis<-function(pc_obj,p_value,obj_fd, metodo, grafico='F', tipo=0){
  if (tipo==0)
    data<-as.data.frame(pc_obj$scores)
  if (tipo==1)
    data<-as.data.frame(t(obj_fd$coefs))
  
  # Set the confidence level and calculate the critical chi-square value
  chi_square_value <- qchisq(p_value, df = ncol(data))
  print (chi_square_value)

  if (metodo=='b'){
    data$mahalnobis<- mahalanobis(data, colMeans(data), cov(data))
  } else if(metodo=='r'){
    cov_robusta <- cov.rob(data, method  = 'mcd')
    data$mahalnobis <- mahalanobis(data, center = cov_robusta$center, cov = cov_robusta$cov)
  } else if(metodo=='50'){
    cov_robusta <- cov.mcd(data,quantile.used =nrow(data)*0.5)
    data$mahalnobis <- mahalanobis(data, center = cov_robusta$center, cov = cov_robusta$cov)
  } else {
    cov_robusta <- cov.mcd(data,quantile.used =nrow(data)*0.75)
    data$mahalnobis <- mahalanobis(data, center = cov_robusta$center, cov = cov_robusta$cov)
  }
  
  data$outlier <- (data$mahalnobis>chi_square_value)
  
  data$color <- ifelse(data$outlier, "red", "gray80")
  pos_out<-which(data$outlier == TRUE)
  print(paste('Total atypical',length(pos_out)))
  print(pos_out)
  
  # Plot the functional curves
  if(grafico=="T"){
    par(mfrow=c(1,3))
    plot(obj_fd, xlab="", ylab="", col=data$color,lwd = 1,  axes =F,
         main ='Curvas SJR',
         lty=1, ylim=(c((min(obj_fd$coefs)),(max(obj_fd$coefs)))))
    par(new = TRUE)
    plot(obj_fd[pos_out,], col='red',xlab='Año', ylab="SJR", axes =F,  lty=1,lwd = 1, 
         ylim=(c((min(obj_fd$coefs)),(max(obj_fd$coefs)))))
    axis(1, at=seq(0,24,5), labels=c('2000','2005','2010','2015','2020'),las=1)
    axis(2)
    box()
    
    # Plot outliers in components
    colores<-c(1:length(pos_out))
    forma<-c((1):(length(pos_out)))
    
    if(dim(pc_obj$scores)[2]>1){
      V1<-pc_obj$scores[,1]
      V2<-pc_obj$scores[,2]
      ylim2_max<-round(max(V2),2)
      ylim2_min<-round(min(V2),2)
      xlim2_max<-round(max(V1),2)
      xlim2_min<-round(min(V1),2)
      
      plot(V1, V2, col="gray55",pch = 21, xlab="PC 1", ylab="PC 2", main ="PC-1 vs PC-2", 
           ylim = c(ylim2_min,ylim2_max), xlim = c(xlim2_min,xlim2_max))
      par(new = TRUE)
      plot(V1[pos_out], V2[pos_out], col=data$color[pos_out], bg = data$color[pos_out],pch = 20, xlab="", 
           ylab="",cex = 1,  ylim = c(ylim2_min,ylim2_max), xlim = c(xlim2_min,xlim2_max),lwd = 1)
    }  
    print(sort(data$mahalnobis))
    
    plot(sort(data$mahalnobis), col=sort(replace(data$color, data$color=='gray80', 'gray55')), pch = 20
         , xlab="Observation", ylab="Distance", main ="Mahalanobis Distance",log = "y")
    abline(h=chi_square_value)
    text(x = length(data$color)/5, y = chi_square_value, labels = paste("\u03C7\u00B2",'=', round(chi_square_value , 3)), pos = 3)
  }
  par(mfrow=c(1,1))
  return(pos_out)
}

#####################################################################################
# Function that obtains the outlier observations using the fbplot method
#   Parameters: 
#     - pc_obj : PC object
#     - obj_fd : Functional object
#     - disc_array : Discrete data matrix
#   Returns:
#     - pos_out: Vector of positions of the outlier observations
#####################################################################################
get_atipicos_fplot<-function(pc_obj,obj_fd,disc_array, grafico=FALSE){
  out<-fbplot(disc_array, col = "blue", main = 'boxplot', method='MBD', plot = FALSE)
  pos_out<-out$outpoint
  print(pos_out)
  
  if(grafico==TRUE){
    # Plot the functional curves
    par(mfrow=c(1,3))
    plot(obj_fd, xlab="", ylab="", col='gray80',lwd = 2,  axes =F,
         main ='Curvas SJR',lty=1, ylim=(c((min(obj_fd$coefs)),(max(obj_fd$coefs)))))
    par(new = TRUE)
    plot(obj_fd[pos_out,], col='red',xlab='Año', ylab="SJR", axes =F,  lty=1, 
         ylim=(c((min(obj_fd$coefs)),(max(obj_fd$coefs)))))
    axis(1, at=seq(0,24,5), labels=c('2000','2005','2010','2015','2020'),las=1)
    axis(2)
    box()

    # Plot outliers in components
    data<-as.data.frame(pc_obj$scores)
    
    if(length(data$V2)!=0){
      ylim2_max<-round(max(data$V2),2)
      ylim2_min<-round(min(data$V2),2)
      xlim2_max<-round(max(data$V1),2)
      xlim2_min<-round(min(data$V1),2)
      
      plot(data$V1, data$V2, col="gray75",pch = 21, xlab="PC 1", ylab="PC 2", main ="PC-1 vs PC-2", 
           ylim = c(ylim2_min,ylim2_max), xlim = c(xlim2_min,xlim2_max))
      par(new = TRUE)
      plot(data$V1[pos_out], data$V2[pos_out], col='red', bg = 'red', xlab="", pch = 20,
           ylab="",cex = 1,  ylim = c(ylim2_min,ylim2_max), xlim = c(xlim2_min,xlim2_max))
    }
  }
  par(mfrow=c(1,1))
  if(length(out$outpoint)==0)
    return(-1)
  else
    return(pos_out)
}


#####################################################################################
# Function that obtains the outlier observations using the Rainbow library
#   Parameters: 
#     - m1     : Model matrix
#     - metodo : Method to be executed
#               -> robMah    :  Integrated Square Forecast Errors
#               -> HUoutliers:  Robust Mahalanobis on the discrete matrix
#     - num_pc : Number of PCs to use, only defined for HUoutliers, default is 2
#   Returns:
#     - outliers: Vector of positions of the outlier observations
#####################################################################################
foutliers_out<-function(disc_array, metodo, num_pc){
  outliers<-c()
  result <- tryCatch({
    obj_sfts=sfts(ts(as.numeric(t(disc_array)),start=c(1,1),frequency=dim(disc_array)[2]), xname="Anio", yname="SJR")
    if (num_pc==0)
      out_v<-foutliers(data=obj_sfts, method = metodo)
    else{
      out_v<-foutliers2(data=obj_sfts, method = metodo, order =num_pc)
    }
    outliers<-out_v$outliers
  }, error = function(err) {  
    outliers<-c(-1)
  }
  )
  print(outliers)
  return(outliers)
}

