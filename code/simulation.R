library(fda.usc)
library(fdaoutlier)

#####################################################################################
# Function that obtains the metrics of the simulation
#   Parameters: 
#     - outliers_simulation: Vector of the original outlier observations of the model
#     - outliers_method:     Vector of the outlier observations of the applied method
#     - dim_df:              Number of observations of the model
#     - num_fpc:             Number of FPCs used in the applied method
#     - sum_arm:             Percentage of variability accumulated by the num_fpc
#   Returns:
#     - df_i: The list with the resulting data
#             [1]-> (Number of correct outliers identified by the method) / (Number of outliers in the model)
#             [2]-> (Number of incorrect outliers identified by the method) / (Number of non-outlier observations in the model)
#             [3] -> Vector of outlier observations
#             [4] -> Number of FPCs used in the applied method
#             [5] -> Percentage of variability accumulated by the num_fpc
#####################################################################################
get_metrics<-function(outliers_simulation,outliers_method,dim_df,num_fpc,sum_arm){
  outliers_common<-intersect(outliers_simulation, outliers_method)
  outliers_wrong<-setdiff(outliers_method,outliers_simulation)
  
  if(outliers_method[1]==-1 || length(outliers_method)==0){
    pc_i<-0
    pf_i<-0
  } else{
    pc_i<-round(length(outliers_common)/length(outliers_simulation), 3)
    pf_i<-round(length(outliers_wrong)/(dim_df-length(outliers_simulation)), 3)
  }
  df_i <- data.frame(pc = pc_i, pf = pf_i, atypical=toString(outliers_method),num_cp=num_fpc,sum=sum_arm)
  return(df_i)
}

#####################################################################################
# Function that obtains the metrics of the simulation using the HDR method
#   Parameters: 
#     - model: Matrix of the model to be executed
#####################################################################################
HDR_out<-function(modelo){
  disc_array2<-modelo$data
  SJRTs=sfts(ts(as.numeric(t(disc_array2)),start=c(1,1),frequency=dim(modelo$data)[2]), xname="x",yname="y")
  par(mfrow=c(1,2))
  fboxplot(data= SJRTs, plot.type = "bivariate", type = "hdr", projmethod = "PCAproj", alpha=c(0.2,0.5), ncol=4)
  fboxplot(data= SJRTs, plot.type = "functional", type = "hdr", projmethod = "PCAproj", alpha=c(0.2,0.5), ncol=4)
  
  fboxplot(data= SJRTs, plot.type = "bivariate", type = "bag", projmethod = "PCAproj" , alpha=c(0.2,0.5), ncol=4)
  fboxplot(data= SJRTs, plot.type = "functional", type = "bag", projmethod = "PCAproj", alpha=c(0.2,0.5), ncol=4)
  par(mfrow=c(1,1))
}

HDR_out_v2<-function(modelo){
  disc_array2<-modelo
  outiers<-''
  SJRTs=sfts(ts(as.numeric(t(disc_array2)),start=c(1,1),frequency=dim(modelo)[2]), xname="x",yname="y")
  par(mfrow=c(1,2))
  outiers$hdr<-fboxplot2(data= SJRTs, plot.type = "bivariate", type = "hdr", alpha = c(0.07, 0.5))
  outiers$bag<-fboxplot2(data= SJRTs, plot.type = "bivariate", type = "bag", factor = 2.5)
  
  par(mfrow=c(1,1))
  return (outiers)
}


#####################################################################################
# Function that executes each method 
#   m1: matrix of the model to be executed
#     m1$data: matrix of observations
#     m1$true_outliers: vector of positions of the outlier observations
#####################################################################################
df_result <- data.frame()
run_modelos<-function(m1){
  df_result_j <- data.frame()
  disc_array<-m1$data
  outliers_simulation<-m1$true_outliers
  
  # create bspline funcional 
  dim_x<-dim(m1$data)[2]
  basisobj_tasa <- create.bspline.basis(c(1,dim_x), nbasis=25, norder=4)
  disc_array<-matrix(as.numeric(disc_array), ncol = ncol(disc_array))
  obj_fd <- Data2fd(y=t(disc_array), argvals=c(1:dim_x), basisobj=basisobj_tasa)
  
  ##################################### 
  # Mahalanobis analysis
  #####################################
  
  # Call the function that obtains the n FPCs that sum up to 95% of variability
  pc_obj_95<-get_harmonics_fpc(obj_fd,0.95)
  num_arm_1<-length(pc_obj_95$varprop)
  
  # Process to obtain the outlier observations using Mahalanobis, FPCs that sum up to 95% of variability, and a confidence level x
  # Mahalanobis classic
  outliers_method_1<-get_atypical_mahalnobis(pc_obj_95,0.975 ,obj_fd,'b')
  outliers_method_2<-get_atypical_mahalnobis(pc_obj_95,0.99 ,obj_fd,'b')
  outliers_method_3<-get_atypical_mahalnobis(pc_obj_95,0.999,obj_fd,'b')
  # Mahalanobis robust
  outliers_method_4<-get_atypical_mahalnobis(pc_obj_95,0.975 ,obj_fd,'r')
  outliers_method_5<-get_atypical_mahalnobis(pc_obj_95,0.99 ,obj_fd,'r')
  outliers_method_6<-get_atypical_mahalnobis(pc_obj_95,0.999,obj_fd,'r')
  # Mahalanobis robust - MCD 50
  outliers_method_7<-get_atypical_mahalnobis(pc_obj_95,0.975 ,obj_fd,'50')
  outliers_method_8<-get_atypical_mahalnobis(pc_obj_95,0.99 ,obj_fd,'50')
  outliers_method_9<-get_atypical_mahalnobis(pc_obj_95,0.999,obj_fd,'50')
  # Mahalanobis robust - MCD 75
  outliers_method_10<-get_atypical_mahalnobis(pc_obj_95,0.975 ,obj_fd,'75')
  outliers_method_11<-get_atypical_mahalnobis(pc_obj_95,0.99 ,obj_fd,'75')
  outliers_method_12<-get_atypical_mahalnobis(pc_obj_95,0.999,obj_fd,'75')
  
  #############
  # Call the function that obtains the n FPCs that sum up to 99% of variability
  pc_obj_99<-get_harmonics_fpc(obj_fd,0.99)
  num_arm_2<-length(pc_obj_99$varprop)
  
  # Process to obtain the outlier observations using Mahalanobis, FPCs that sum up to 99% of variability, and a confidence level x
  # Mahalanobis classic
  outliers_method_13<-get_atypical_mahalnobis(pc_obj_99,0.975 ,obj_fd,'b')
  outliers_method_14<-get_atypical_mahalnobis(pc_obj_99,0.99 ,obj_fd,'b')
  outliers_method_15<-get_atypical_mahalnobis(pc_obj_99,0.999,obj_fd,'b')
  # Mahalanobis robust
  outliers_method_16<-get_atypical_mahalnobis(pc_obj_99,0.975 ,obj_fd,'r')
  outliers_method_17<-get_atypical_mahalnobis(pc_obj_99,0.99 ,obj_fd,'r')
  outliers_method_18<-get_atypical_mahalnobis(pc_obj_99,0.999,obj_fd,'r')
  # Mahalanobis robust - MCD 50
  outliers_method_19<-get_atypical_mahalnobis(pc_obj_99,0.975 ,obj_fd,'50')
  outliers_method_20<-get_atypical_mahalnobis(pc_obj_99,0.99 ,obj_fd,'50')
  outliers_method_21<-get_atypical_mahalnobis(pc_obj_99,0.999,obj_fd,'50')
  # Mahalanobis robust - MCD 75
  outliers_method_22<-get_atypical_mahalnobis(pc_obj_99,0.975 ,obj_fd,'75')
  outliers_method_23<-get_atypical_mahalnobis(pc_obj_99,0.99 ,obj_fd,'75')
  outliers_method_24<-get_atypical_mahalnobis(pc_obj_99,0.999,obj_fd,'75')
  
  
  #############
  # Analysis with 2C
  pc_obj_fpca_2C<-get_harmonics_2C(obj_fd)
  outliers_method_25<-get_atypical_mahalnobis(pc_obj_fpca_2C,0.95 ,obj_fd,'b')
  outliers_method_26<-get_atypical_mahalnobis(pc_obj_fpca_2C,0.95 ,obj_fd,'r')
  outliers_method_27<-get_atypical_mahalnobis(pc_obj_fpca_2C,0.95 ,obj_fd,'50')
  outliers_method_28<-get_atypical_mahalnobis(pc_obj_fpca_2C,0.95 ,obj_fd,'75')
  
  
  #####################################
  # Analysis to be performed using the fbplot method
  #####################################
  outliers_method_29<-get_atypical_fplot(pc_obj_95,obj_fd,t(disc_array))
  
  
  #####################################
  # Analysis to be performed using the Mahalanobis method from the Rainbow library
  #####################################
  outliers_method_30<-foutliers_out(m1$data,'robMah',0)
  
  
  #####################################
  # Analysis to be performed using the Integrated Square Forecast Errors method from the Rainbow library
  #####################################
  outliers_method_31<-foutliers_out(m1$data,'HUoutliers',2)
  outliers_method_32<-c(0)
  outliers_method_33<-c(0)
  
  outliers_method_34<-c(0)
  outliers_method_35<-c(0)
  
  
  #####################################
  # Analysis to be performed using the Mahalanobis method with LAMBDA > 1
  #####################################
  pc_obj_lambda<-get_harmonics_fpc_lambda(obj_fd)
  num_arm_lambda_1<-length(pc_obj_lambda$varprop)
  
  #####################################
  # Analysis to be performed using the Mahalanobis method with obj_fd coefficients
  #####################################
  # Mahalanobis classic
  outliers_method_48<-get_atypical_mahalnobis(pc_obj_95,0.975 ,obj_fd,'b',T,1)
  outliers_method_49<-get_atypical_mahalnobis(pc_obj_95,0.99 ,obj_fd,'b',T,1)
  outliers_method_50<-get_atypical_mahalnobis(pc_obj_95,0.999,obj_fd,'b',T,1)
  # Mahalanobis robust
  outliers_method_51<-get_atypical_mahalnobis(pc_obj_95,0.975 ,obj_fd,'r',T,1)
  outliers_method_52<-get_atypical_mahalnobis(pc_obj_95,0.99 ,obj_fd,'r',T,1)
  outliers_method_53<-get_atypical_mahalnobis(pc_obj_95,0.999,obj_fd,'r',T,1)
  # Mahalanobis robust - MCD 50
  outliers_method_54<-get_atypical_mahalnobis(pc_obj_95,0.975 ,obj_fd,'50',T,1)
  outliers_method_55<-get_atypical_mahalnobis(pc_obj_95,0.99 ,obj_fd,'50',T,1)
  outliers_method_56<-get_atypical_mahalnobis(pc_obj_95,0.999,obj_fd,'50',T,1)
  # Mahalanobis robust - MCD 75
  outliers_method_57<-get_atypical_mahalnobis(pc_obj_95,0.975 ,obj_fd,'75',T,1)
  outliers_method_58<-get_atypical_mahalnobis(pc_obj_95,0.99 ,obj_fd,'75',T,1)
  outliers_method_59<-get_atypical_mahalnobis(pc_obj_95,0.999,obj_fd,'75',T,1)
  
  #####################################
  # Analysis to be performed using the f-HDR boxplot and F-bagplot methods 
  #####################################
  outliers_met<-HDR_out_v2(m1$data)
  
  
  print('==========FIN MODELADO=======') 
  sum_arm<-round(pc_obj_95$sum_harmonics,4)
  print(sum_arm)
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_1,dim(m1$data)[1],num_arm_1,sum_arm))
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_2,dim(m1$data)[1],num_arm_1,sum_arm))
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_3,dim(m1$data)[1],num_arm_1,sum_arm))
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_4,dim(m1$data)[1],num_arm_1,sum_arm))
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_5,dim(m1$data)[1],num_arm_1,sum_arm))
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_6,dim(m1$data)[1],num_arm_1,sum_arm))
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_7,dim(m1$data)[1],num_arm_1,sum_arm))
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_8,dim(m1$data)[1],num_arm_1,sum_arm))
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_9,dim(m1$data)[1],num_arm_1,sum_arm))
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_10,dim(m1$data)[1],num_arm_1,sum_arm))
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_11,dim(m1$data)[1],num_arm_1,sum_arm))
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_12,dim(m1$data)[1],num_arm_1,sum_arm))
  sum_arm2<-round(pc_obj_99$sum_harmonics,4)
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_13,dim(m1$data)[1],num_arm_2,sum_arm2))
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_14,dim(m1$data)[1],num_arm_2,sum_arm2))
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_15,dim(m1$data)[1],num_arm_2,sum_arm2))
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_16,dim(m1$data)[1],num_arm_2,sum_arm2))
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_17,dim(m1$data)[1],num_arm_2,sum_arm2))
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_18,dim(m1$data)[1],num_arm_2,sum_arm2))
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_19,dim(m1$data)[1],num_arm_2,sum_arm2))
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_20,dim(m1$data)[1],num_arm_2,sum_arm2))
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_21,dim(m1$data)[1],num_arm_2,sum_arm2))
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_22,dim(m1$data)[1],num_arm_2,sum_arm2))
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_23,dim(m1$data)[1],num_arm_2,sum_arm2))
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_24,dim(m1$data)[1],num_arm_2,sum_arm2))
  sum_arm3<-round(pc_obj_fpca_2C$sum_harmonics,4)
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_25,dim(m1$data)[1],2,sum_arm3))
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_26,dim(m1$data)[1],2,sum_arm3))
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_27,dim(m1$data)[1],2,sum_arm3))
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_28,dim(m1$data)[1],2,sum_arm3))
  
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_29,dim(m1$data)[1],num_arm_1,sum_arm))
  
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_30,dim(m1$data)[1],0,0))
  
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_31,dim(m1$data)[1],2,sum_arm3))
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_32,dim(m1$data)[1],num_arm_1,sum_arm))
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_33,dim(m1$data)[1],num_arm_2,sum_arm2))
  
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_34,dim(m1$data)[1],num_arm_1,sum_arm))
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_35,dim(m1$data)[1],num_arm_2,sum_arm2))
  
  sum_arm4<-round(pc_obj_lambda$sum_harmonics,4)
  
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_48,dim(m1$data)[1],-1,-1))
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_49,dim(m1$data)[1],-1,-1))
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_50,dim(m1$data)[1],-1,-1))
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_51,dim(m1$data)[1],-1,-1))
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_52,dim(m1$data)[1],-1,-1))
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_53,dim(m1$data)[1],-1,-1))
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_54,dim(m1$data)[1],-1,-1))
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_55,dim(m1$data)[1],-1,-1))
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_56,dim(m1$data)[1],-1,-1))
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_57,dim(m1$data)[1],-1,-1))
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_58,dim(m1$data)[1],-1,-1))
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_method_59,dim(m1$data)[1],-1,-1))
  
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_met$hdr,dim(m1$data)[1],-1,-1))
  df_result_j<-rbind(df_result_j, get_metrics(outliers_simulation,outliers_met$bag,dim(m1$data)[1],-1,-1))
  
  df_result_j<-rbind(df_result_j, data.frame(pc = '-', pf = '-', atypical=toString(outliers_simulation), num_cp=0,sum=0))
  return(df_result_j)
}

df_pc <- data.frame()
df_pf <- data.frame()

num_observaciones<-200
out_rate<-0.2

for(j in c(1,2,4,6,8,9)){
  df_result <- data.frame(matrix(nrow = 50, ncol = 0))
  df_pc <- data.frame()
  df_pf <- data.frame()
  
  for(i in c(1:100)){
    print("======================")
    print(paste("Executing: model",j, ", iteration", i))
    if(j==1) m <-simulation_model1(n = num_observaciones, plot = FALSE, outlier_rate=out_rate,seed=i)
    else if(j==2) m <-simulation_model2(n = num_observaciones, plot = FALSE, outlier_rate=out_rate,seed=i)
    else if(j==4) m <-simulation_model4(n = num_observaciones, plot = FALSE, outlier_rate=out_rate,seed=i)
    else if(j==6) m <-simulation_model6(n = num_observaciones, plot = FALSE, outlier_rate=out_rate,seed=i)
    else if(j==8) m <-simulation_model8(n = num_observaciones, plot = FALSE, outlier_rate=out_rate,seed=i)
    else if(j==9) m <-simulation_model9(n = num_observaciones, plot = FALSE, outlier_rate=out_rate,seed=i)
    
    dir_esc<-paste0("sim_",num_observaciones,"_",gsub("0\\.", "",out_rate),"\\")
    name_arc_1<-paste0(dir_esc,"m",j,"_",i,".rds")
    saveRDS(m, file = name_arc_1)
    
    name_arc_2<-paste0(dir_esc,"m",j,".csv")
    resultado_modelo<-run_modelos(m)
    df_result<-cbind(df_result, resultado_modelo)
    write.csv(df_result, name_arc_2, row.names = FALSE)
    
    df_pc<-rbind(df_pc, t(resultado_modelo$pc))
    df_pf<-rbind(df_pf, t(resultado_modelo$pf))
    
    write.csv(df_pc, paste0(dir_esc,"pc_m",j,".csv"), row.names = FALSE)
    write.csv(df_pf, paste0(dir_esc,"pf_m",j,".csv"), row.names = FALSE)
  }
}
