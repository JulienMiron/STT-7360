## Section 8 -- exemple sur un seul jeu de données (normale bivariée, MAR)
suppressPackageStartupMessages({library(mvtnorm); library(mice)})
figdir <- "/home/claude/refonte/notes/figures"

mu <- c(1, 2); s <- c(0.5, 1); rho <- 0.5
Sigma <- diag(s) %*% matrix(c(1, rho, rho, 1), 2) %*% diag(s)
expit <- function(x) 1 / (1 + exp(-x))

## Mécanismes : on choisit au hasard la variable « pilote » j ;
## MCAR : l'autre variable manque avec prob. constante ;
## MAR  : l'autre variable manque avec prob. croissante en Y_j (observée) ;
## MNAR : Y_j elle-même manque avec prob. croissante en Y_j.
generer <- function(n, mecanisme = c("MCAR", "MAR", "MNAR"), a = -0.6, g = 2) {
  mecanisme <- match.arg(mecanisme)
  Y <- rmvnorm(n, mu, Sigma)
  Z <- sweep(sweep(Y, 2, mu), 2, s, "/")          # valeurs centrées réduites
  j <- sample(2, n, replace = TRUE)
  zj <- Z[cbind(1:n, j)]
  p <- switch(mecanisme, MCAR = rep(expit(a), n),
              MAR = expit(a + g * zj), MNAR = expit(a + g * zj))
  manque <- runif(n) < p
  col <- if (mecanisme == "MNAR") j else 3 - j
  Y[cbind(which(manque), col[manque])] <- NA
  colnames(Y) <- c("Y1", "Y2")
  Y
}

## (1) Données disponibles (H26 : « analyse de données complètes »)
cc <- function(Y) c(mu1 = mean(Y[, 1], na.rm = TRUE), mu2 = mean(Y[, 2], na.rm = TRUE),
                    s1 = sd(Y[, 1], na.rm = TRUE), s2 = sd(Y[, 2], na.rm = TRUE),
                    rho = cor(Y[, 1], Y[, 2], use = "complete.obs"))

## (2) Vraisemblance des données observées : l11 + l10 + l01
nll <- function(th, Y) {
  m <- th[1:2]; sd <- th[3:4]; r <- th[5]
  if (any(sd <= 0) || abs(r) >= 1) return(1e10)
  S <- diag(sd) %*% matrix(c(1, r, r, 1), 2) %*% diag(sd)
  o1 <- !is.na(Y[, 1]); o2 <- !is.na(Y[, 2])
  -(sum(dmvnorm(Y[o1 & o2, , drop = FALSE], m, S, log = TRUE)) +
      sum(dnorm(Y[o1 & !o2, 1], m[1], sd[1], log = TRUE)) +
      sum(dnorm(Y[!o1 & o2, 2], m[2], sd[2], log = TRUE)))
}

## (3) Algorithme EM pour la normale bivariée
em <- function(Y, tol = 1e-8, maxit = 500, trace = FALSE) {
  n <- nrow(Y); o1 <- !is.na(Y[, 1]); o2 <- !is.na(Y[, 2])
  m <- colMeans(Y, na.rm = TRUE); S <- diag(apply(Y, 2, var, na.rm = TRUE))
  for (it in 1:maxit) {
    ## Étape E : espérances conditionnelles de Y et de Y Y^T
    Yh <- Y; C <- matrix(0, 2, 2)              # C : correction des moments d'ordre 2
    b21 <- S[1, 2] / S[1, 1]; v21 <- S[2, 2] - S[1, 2]^2 / S[1, 1]
    b12 <- S[1, 2] / S[2, 2]; v12 <- S[1, 1] - S[1, 2]^2 / S[2, 2]
    Yh[!o2, 2] <- m[2] + b21 * (Y[!o2, 1] - m[1]); C[2, 2] <- sum(!o2) * v21
    Yh[!o1, 1] <- m[1] + b12 * (Y[!o1, 2] - m[2]); C[1, 1] <- sum(!o1) * v12
    ## Étape M : EMV en données complètes à partir des statistiques exhaustives attendues
    m.new <- colMeans(Yh)
    S.new <- (crossprod(sweep(Yh, 2, m.new)) + C) / n
    ll <- -nll(c(m.new, sqrt(diag(S.new)), cov2cor(S.new)[1, 2]), Y)
    if (trace) cat(sprintf("it %2d  mu2 = %.5f  log-vrais = %.6f\n", it, m.new[2], ll))
    ok <- max(abs(c(m.new - m, S.new - S))) < tol
    m <- m.new; S <- S.new
    if (ok) break
  }
  list(mu = m, Sigma = S, iter = it)
}

## (4) Imputation simple par régression (moyenne conditionnelle)
imput.reg <- function(Y) {
  o1 <- !is.na(Y[, 1]); o2 <- !is.na(Y[, 2]); cpl <- o1 & o2
  f21 <- lm(Y2 ~ Y1, data = as.data.frame(Y[cpl, ]))
  f12 <- lm(Y1 ~ Y2, data = as.data.frame(Y[cpl, ]))
  Yi <- Y
  Yi[!o2, 2] <- predict(f21, newdata = data.frame(Y1 = Y[!o2, 1]))
  Yi[!o1, 1] <- predict(f12, newdata = data.frame(Y2 = Y[!o1, 2]))
  Yi
}

## (5) Équations chaînées « à la main » (algorithme du document H26, nettoyé)
tirage.reg <- function(y, x, obs) {
  # Tire (tau2*, beta*) de leur loi a posteriori (a priori non informatif)
  X <- cbind(1, x[obs]); fit <- lm.fit(X, y[obs])
  nobs <- sum(obs); XtXi <- chol2inv(qr.R(fit$qr))
  tau2 <- sum(fit$residuals^2) / rchisq(1, df = nobs - 2)
  beta <- drop(fit$coefficients + t(chol(tau2 * XtXi)) %*% rnorm(2))
  list(beta = beta, tau2 = tau2)
}
mice.main <- function(Y, M = 5, cycles = 10) {
  o1 <- !is.na(Y[, 1]); o2 <- !is.na(Y[, 2])
  lapply(1:M, function(m) {
    Yc <- Y
    Yc[!o1, 1] <- sample(Y[o1, 1], sum(!o1), replace = TRUE)   # valeurs initiales
    Yc[!o2, 2] <- sample(Y[o2, 2], sum(!o2), replace = TRUE)
    for (k in 1:cycles) {
      p <- tirage.reg(Yc[, 2], Yc[, 1], o2)                       # Y2 | Y1
      Yc[!o2, 2] <- p$beta[1] + p$beta[2] * Yc[!o2, 1] + rnorm(sum(!o2), 0, sqrt(p$tau2))
      p <- tirage.reg(Yc[, 1], Yc[, 2], o1)                       # Y1 | Y2
      Yc[!o1, 1] <- p$beta[1] + p$beta[2] * Yc[!o1, 2] + rnorm(sum(!o1), 0, sqrt(p$tau2))
    }
    Yc
  })
}

## (6) Augmentation de données, a priori de Jeffreys (normale--Wishart inverse)
rinvwish <- function(df, Psi) solve(rWishart(1, df, solve(Psi))[, , 1])
aug.donnees <- function(Y, M = 5, burn = 50, espacement = 10) {
  n <- nrow(Y); o1 <- !is.na(Y[, 1]); o2 <- !is.na(Y[, 2])
  Yc <- imput.reg(Y); res <- list()
  for (it in 1:(burn + M * espacement)) {
    ## Étape P : tirer (Sigma, mu) de leur loi a posteriori sachant les données courantes
    xb <- colMeans(Yc); S <- cov(Yc)
    Sig <- rinvwish(n - 1, (n - 1) * S)
    m <- drop(rmvnorm(1, xb, Sig / n))
    ## Étape I : tirer les valeurs manquantes de leur loi conditionnelle
    Yc[!o1, 1] <- rnorm(sum(!o1), m[1] + Sig[1, 2] / Sig[2, 2] * (Yc[!o1, 2] - m[2]),
                        sqrt(Sig[1, 1] - Sig[1, 2]^2 / Sig[2, 2]))
    Yc[!o2, 2] <- rnorm(sum(!o2), m[2] + Sig[1, 2] / Sig[1, 1] * (Yc[!o2, 1] - m[1]),
                        sqrt(Sig[2, 2] - Sig[1, 2]^2 / Sig[1, 1]))
    if (it > burn && (it - burn) %% espacement == 0) res[[length(res) + 1]] <- Yc
  }
  res
}

## Règles de Rubin pour un paramètre scalaire
rubin <- function(est, v) {
  M <- length(est); qb <- mean(est); W <- mean(v); B <- var(est)
  Tt <- W + (1 + 1/M) * B; r <- (1 + 1/M) * B / W
  nu <- (M - 1) * (1 + 1/r)^2
  c(est = qb, W = W, B = B, T = Tt, se = sqrt(Tt), dl = nu, fmi = (1 + 1/M) * B / Tt)
}

if (sys.nframe() == 0) {
  set.seed(7360)
  n <- 200
  Y <- generer(n, "MAR")
  cat("Patrons (R1, R2) :\n"); print(table(R1 = 1 * !is.na(Y[, 1]), R2 = 1 * !is.na(Y[, 2])))
  cat("\nDonnées disponibles :\n"); print(round(cc(Y), 3))
  op <- optim(cc(Y), nll, Y = Y, method = "BFGS", hessian = TRUE)
  cat("\nEMV (optim) :\n"); print(round(op$par, 3))
  cat("Erreurs-types (inverse de la hessienne) :\n"); print(round(sqrt(diag(solve(op$hessian))), 3))
  cat("\nAlgorithme EM :\n")
  e <- em(Y, trace = TRUE)
  print(round(c(e$mu, sqrt(diag(e$Sigma)), cov2cor(e$Sigma)[1, 2]), 3))
  Yi <- imput.reg(Y)
  cat("\nImputation simple par régression :\n"); print(round(cc(Yi), 3))
  cat("ET naïve de mu2 :", round(sd(Yi[, 2]) / sqrt(n), 4), "\n")

  set.seed(1)
  imp1 <- mice.main(Y, M = 20)
  est <- sapply(imp1, function(D) mean(D[, 2])); v <- sapply(imp1, function(D) var(D[, 2]) / n)
  cat("\nMICE à la main, mu2 :\n"); print(round(rubin(est, v), 4))
  imp2 <- aug.donnees(Y, M = 20)
  est <- sapply(imp2, function(D) mean(D[, 2])); v <- sapply(imp2, function(D) var(D[, 2]) / n)
  cat("\nAugmentation de données, mu2 :\n"); print(round(rubin(est, v), 4))

  ## mice + pool
  D <- as.data.frame(Y)
  imp <- mice(D, m = 20, method = "norm", maxit = 10, seed = 7360, printFlag = FALSE)
  fit <- with(imp, lm(Y2 ~ 1))
  cat("\nmice + pool, mu2 :\n"); print(summary(pool(fit)))
  print(pool(fit))
  fit2 <- with(imp, lm(Y2 ~ Y1))
  print(summary(pool(fit2), conf.int = TRUE))
  cat("\nRéférence données complètes inconnues ; vraie valeur mu2 = 2\n")

  ## Figure : imputation simple vs multiple
  cairo_pdf(file.path(figdir, "08_imputation.pdf"), width = 8, height = 3.8)   # cairo : accents
  par(mfrow = c(1, 2), mar = c(4, 4, 2, 1), cex = 0.85)
  m2 <- is.na(Y[, 2]) & !is.na(Y[, 1])
  plot(Yi[!m2, 1], Yi[!m2, 2], pch = 16, col = "grey60", xlab = "Y1", ylab = "Y2",
       main = "Imputation simple (régression)", xlim = range(Yi[, 1]), ylim = range(c(Yi[, 2], imp1[[1]][, 2])))
  points(Yi[m2, 1], Yi[m2, 2], pch = 17, col = "firebrick")
  plot(imp1[[1]][!m2, 1], imp1[[1]][!m2, 2], pch = 16, col = "grey60", xlab = "Y1", ylab = "Y2",
       main = "Imputation stochastique (une des M)", xlim = range(Yi[, 1]), ylim = range(c(Yi[, 2], imp1[[1]][, 2])))
  points(imp1[[1]][m2, 1], imp1[[1]][m2, 2], pch = 17, col = "steelblue")
  dev.off()
}
