C_Index = function(df, lifetime = "lifetime", metric = "h_t", forward = TRUE){
  df = df[order(df[lifetime]),]
  edges = 0
  
  concordant = 0
  
  dfe = df[df$d_i != 0,]
  
  for(i in 1:nrow(dfe)){
    
    t_i = dfe[i, lifetime]
    h_i = dfe[i, metric]
    d_i = dfe[i, "d_i"]
    
    # message(paste("i ==", i, "  ", "d_i == ", d_i))
    
    for(j in 1:nrow(df)){
      if(i == j){ next }
      t_j = df[j, lifetime]
      h_j = df[j, metric]
      d_j = df[j, "d_i"]
      
      if(t_j > t_i){
        edges = edges + 1
        
        if(forward == TRUE) {
          if(h_i > h_j){
            concordant = concordant + 1
          }
        } else {
          if(h_i < h_j){
            concordant = concordant + 1
          }
        }

        
      }
      
      
      
    }
    
  }
  message(paste("Concordant Pairs", concordant, 
                ", Total Edges", edges, 
                ", C-Index", (concordant / edges) ))
  return((concordant / edges))
}

