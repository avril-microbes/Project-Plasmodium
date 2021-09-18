# script to calculate total infected RBC and total merozoite population for lagged
# dynamics df. Accepts output of lagged model.
total_parasite <- function(dyn){
  dyn2 <- data.frame()
  I <- dyn %>% 
    dplyr::filter(variable == "I") %>% 
    dplyr::select(value)
  
  Ig <- dyn %>% 
    dplyr::filter(variable == "Ig") %>% 
    dplyr::select(value)
  
  M <- dyn %>% 
    dplyr::filter(variable == "M") %>% 
    dplyr::select(value)
  
  Mg <- dyn %>% 
    dplyr::filter(variable == "Mg") %>% 
    dplyr::select(value)
  
  total_I <- data.frame(time = dyn$time, variable = rep("total_I", length(dyn)), value = I+Ig)
  total_M <- data.frame(time = dyn$time, variable = rep("total_M", length(dyn)), value = M+Mg)
  
  dyn2 <- rbind(dyn, total_I, total_M)
}