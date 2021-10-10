#### SSVS function used in the project report

SSVS = function(X, y, nsave = 1000, nburn = 1000, tau0 = 0.01, 
                tau1 = 10, S0 = 0.01, s0 = 0.01, scale = FALSE){
  
  # Modell preliminaries
  X = as.matrix(X)
  N = nrow(X)
  
  
  ## Gibbs preliminaries
  nsave = nsave                # S_1
  nburn = nburn                # S_0
  ntot = nsave+nburn           # S
  
  
  ## Prior preliminaries für SSVS
  tau0 = tau0    # Prior Varianz für den Fall Beta_j == 0
  # -> Tau0 muss klein sein, damit die Präzision hoch ist
  #    (siehe Kapitel 3.2)
  tau1 = tau1      # Prior Varianz für den Fall Beta_j != 0
  
  
  ## Prior auf Sigma 
  s0 = s0 # uninformativ
  S0 = S0 # uninformativ
  
  
  ## Konstruieren einer neuen Designmatrix durch Hinzufügen
  ## nicht relevanter Kovariaten
  Y = matrix(y)
  X = cbind(rnorm(N,0,10),rnorm(N,0,1), X)  # Designmatrix
  
  if (scale == TRUE){
    Y = scale(Y)
    X = scale(X)
  }
  
  colnames(X)[1] = "x1"
  colnames(X)[2] = "x2"
  N = nrow(Y)
  K = ncol(X)
  
  
  ## Berechnung der OLS Größen
  A.OLS = solve(crossprod(X))%*%crossprod(X,Y)
  SSE = crossprod(Y-X%*%A.OLS)
  SIG.OLS = SSE/(N-K)
  
  
  ## Gamma
  ## Startwert für Gamma wird festgelegt
  ## gamma == gamma_1, gamma_2, ..., gamma_K (siehe Kapitel 3.2)
  ## Gamma == 0 -> Variable nicht inkludiert
  ## Gamma == 1 -> Variable inkludiert
  
  gamma = matrix(1,K,1) 
  # -> Startwert von Gamma ist ein Einheitsvektor. D.h. es wird mit dem
  #    vollen Modell gestartet. Zu Beginn existiert also die Annahme, dass alle
  #    Variablen für das Modell wichtig sind.
  
  
  ## Startwert für Sigma
  sigma2.draw = as.numeric(SIG.OLS)    # Start mit OLS Schätzer
  
  
  ## Damit Gibbs verwendet werden kann, muss die Prior Varianz angepasst werden
  V.prior = diag(as.numeric(gamma*tau1+(1-gamma)*tau0))
  # -> Abhängig von gamma können wir also jetzt flexibel bestimmen, ob die Varianz 
  #    im Prior groß oder klein ist. Ist gamma == 0, wäre auch V.prior nahe 0. Ist 
  #    gamma == 1, wäre V.prior abhängig von tau_1 z.B. 100 oder 1000, 
  #    also sehr hoch.
  
  
  ## Storage matrizzen
  ALPHA.store = matrix(NA,nsave,K)     # Regressionskoeffizienten
  SIGMA.store = matrix(NA,nsave,1)     # Fehlervarianz
  Gamma.store = matrix(NA,nsave,K)     # Gamma
  
  
  ## Gibbs Sampler
  for (irep in 1:ntot){
    # Draw ALPHA given rest from multivariate normal
    # Ziehen von den Regressionskoeffizienten
    V.post = solve(crossprod(X)*1/sigma2.draw+diag(1/diag(V.prior)))
    A.post = V.post%*%(crossprod(X,Y)*1/sigma2.draw)
    A.draw = A.post+t(chol(V.post))%*%rnorm(K)
    
    # Draw indicators conditional on ALPHA
    # Hierbei wird jeder Koeffizient getestet
    # Wir schauen für jeden Koeffizienten, ob wir eine große oder kleine
    # Varianz im Prior brauchen!
    for (jj in 1:K){
      p0 = dnorm(A.draw[[jj]],0,sqrt(tau0))
      p1 = dnorm(A.draw[[jj]],0,sqrt(tau1))
      p11 = p1/(p0+p1)     
      # -> Posterior Wahrscheinlichkeit, dass gamma_j == 1 gilt
      
      if (p11>runif(1)) gamma[[jj]] = 1 else gamma[[jj]] = 0
      # -> Vergleich p11 mit Einheitsverteilung
      # -> If TRUE: gamma_j == 1, ELSE: gamma_j == 0
    }
    
    # Construct prior VC matrix conditional on gamma
    # Neuer V.prior auf Basis der aktualisierten gamma_j's
    V.prior = diag(as.numeric(gamma*tau1+(1-gamma)*tau0))
    
    # Simulate sigma2 from inverse Gamma
    S.post = crossprod(Y-X%*%A.draw)/2+S0
    s.post = S0+N/2
    sigma2.draw = 1/rgamma(1,s.post,S.post)  
    
    # Speichern der Größen nach burn-in 
    if (irep>nburn){
      ALPHA.store[irep-nburn,] = A.draw
      SIGMA.store[irep-nburn,] = sigma2.draw
      Gamma.store[irep-nburn,] = gamma
    }
    
  }
  
  
  ## Posterior Inclusion Probabilities
  PIP.mean = apply(Gamma.store,2,mean)
  
  A.mean = apply(ALPHA.store,2,mean)
  SIG.mean = apply(SIGMA.store,2,mean)
  
  
  ##### PLOTS ######
  ### PIP Plot
  colnames(X)[1] = "x1"
  colnames(X)[2] = "x2"
  
  PIPs = as.data.frame(PIP.mean)
  PIPs$predictor = colnames(X)
  PIPs$predictor = factor(PIPs$predictor, levels = PIPs$predictor)
  
  PIP_Plot = ggplot(PIPs, aes(x=predictor, y=PIP.mean)) +
    geom_segment(
      aes(x=predictor, xend=predictor, y=0, yend=PIP.mean), 
      color=ifelse(PIPs$PIP.mean > 0.5, "orange", "grey"), 
      size=ifelse(PIPs$PIP.mean > 0.5, 1.3, 0.7)
    ) +
    geom_point(
      color=ifelse(PIPs$PIP.mean > 0.5, "orange", "grey"), 
      size=ifelse(PIPs$PIP.mean > 0.5, 5, 2)
    ) +
    theme_light() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    ylim(c(0,1)) +
    labs(title = "Posterior Inclusion Probability") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(plot.title = element_text(size = 20, face = "bold"),
          axis.title.x = element_text(size = 17),
          axis.title.y = element_text(size = 17),
          legend.title=element_text(size=17), 
          legend.text=element_text(size=15),
          text = element_text(size = 15))
  PIP_Plot
  
  
  ### Alpha_store density für alle Koeffizienten bei denen PIP.mean > 0.5 ist
  ## Base R
  plot(density(ALPHA.store[,1]), xlim = c(-1, 3), col = "grey",
       main = "Density Plot", bty = "n", xlab = "")
  grid(NA, NULL, lty=3, lwd=1, col="lightgrey")
  abline(v=c(-1,0,1,2,3), lty=3, lwd=1, col="lightgrey")
  
  for (i in 1:10){
    if (PIPs$PIP.mean[i] > 0.5){
      lines(density(ALPHA.store[,i]), col = adjustcolor('orange',  blue.f = 0.8), lwd = 2)
      text(mean(ALPHA.store[,i]), 4, colnames(X)[i], srt = -30, cex = 0.7)
    } else {
      lines(density(ALPHA.store[,i]), col = adjustcolor('grey', alpha.f = 0.5))
    }
  }
  abline(v = c(0), lty = 2, lwd = 1, col = "firebrick")
  
  PIP_Density_Plot = recordPlot()
  
  
  ### All Densities
  ## GGPLOT Data
  Alphas = as.data.frame(ALPHA.store)
  colnames(Alphas) = colnames(X)
  
  ## Get GGPlot Format
  Alphas2 = stack(Alphas)
  
  # Density Plot
  All_Densities = ggplot(Alphas2, aes(x = values, y = ind, fill = ind)) +
    geom_density_ridges() +
    theme_ridges()+
    theme(legend.position = "none") +
    labs(title = "Density plot", y = "Covariate", x = "") 
  
  
  ## Plot2
  # Density excl ESCAUSx Plot
  Alphas3 = Alphas2[!Alphas2$ind == "EXCAUSx", ]
  
  Ex_Densities = ggplot(Alphas3, aes(x = values, y = ind, fill = ind)) +
    geom_density_ridges() +
    theme_ridges() +
    theme(legend.position = "none") +
    labs(title = "Density plot", y = "Covariate", x = "") 
  
  
  
  ##### Ergebnisse plotten ######
  beta.mean = apply(ALPHA.store, 2, mean)
  
  # Ts.plot
  Estimation = ggplot(FRED2, aes(x = sasdate)) +
    geom_line(aes(y = Y, color = "Ölpreis Data")) +
    geom_line(aes(y = as.numeric(X %*% beta.mean), color = "Posterior Model Estimation")) +
    labs(
      title = "Ölpreisentwicklung vs. Bayesianische Schätzung",
      y = "Preis in USD",
      x = "",
      colour = "Legende"
    ) +
    scale_color_manual(values = c("lightgrey", "orange")) +
    theme_minimal() +
    theme(
      legend.position = c(0.15, 0.9),
      legend.background = element_rect(linetype = "solid")
    )
  if (scale == TRUE){
    Estimation = Estimation + labs(
      title = "Ölpreisentwicklung vs. Bayesianische Schätzung (skalierte Werte)",
      y = "Preis in USD (skaliert)"
    )
  }
  
  ### FIT
  pred = as.numeric(X%*%beta.mean)
  results = as.data.frame(cbind(pred, Y))
  colnames(results) = c("Prediction", "Y")
  results$MSE = (results$Y - results$Prediction)^2 / N
  MSE = colSums(results)["MSE"]
  
  ### Return
  ssvs = list(PIP_Plot, All_Densities, Ex_Densities, Estimation, MSE, results)
  
}