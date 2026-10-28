## Section 9 -- attrition MAR, GEE non pondérées vs IPCW
suppressPackageStartupMessages(library(geepack))
figdir <- "/home/claude/refonte/notes/figures"
expit <- function(x) 1 / (1 + exp(-x))

## Données : n sujets, visites t = 0,1,2,3, groupe X (0/1)
## Modèle marginal : E[Y_it | X_i] = 10 + 1 t - 0.5 X + 0.5 X t
## Attrition monotone : P(C_t = 1 | C_{t-1} = 0, Ybar_{t-1}, X) = expit(a0 + a1 (Y_{t-1} - 12) + a2 X)
beta.vrai <- c(10, 1, -0.5, 0.5)
simuler <- function(n = 400, a = c(-1.5, 0.6, -0.5)) {
  id <- rep(1:n, each = 4); t <- rep(0:3, n); X <- rep(rbinom(n, 1, 0.5), each = 4)
  b <- rep(rnorm(n, 0, 1.5), each = 4)
  e <- as.vector(apply(matrix(rnorm(4 * n), n), 1, cumsum)) * 0.7   # marche aléatoire par sujet
  Y <- beta.vrai[1] + beta.vrai[2] * t + beta.vrai[3] * X + beta.vrai[4] * X * t + b + e
  d <- data.frame(id, t, X, Y)
  d$Yprec <- ave(d$Y, d$id, FUN = function(y) c(NA, y[-4]))
  ## censure séquentielle
  d$present <- 1
  for (k in 1:3) {
    i <- d$t == k
    p.abandon <- expit(a[1] + a[2] * (d$Yprec[i] - 12) + a[3] * d$X[i])
    prec.present <- d$present[d$t == k - 1]
    d$present[i] <- prec.present * (runif(sum(i)) > p.abandon)
  }
  d
}

## Poids IPCW : modèles logistiques séquentiels ajustés sur l'ensemble à risque
poids.ipcw <- function(d) {
  d$risque <- d$t > 0 & ave(d$present, d$id, FUN = function(r) c(0, r[-4])) == 1
  ar <- d[d$risque, ]
  m.den <- glm(present ~ factor(t) + Yprec + X, family = binomial, data = ar)  # P(C_t=0 | C_{t-1}=0, Zbar_{t-1})
  m.num <- glm(present ~ factor(t) + X, family = binomial, data = ar)         # numérateur stabilisé
  d$p.den <- 1; d$p.num <- 1
  d$p.den[d$risque] <- predict(m.den, type = "response")
  d$p.num[d$risque] <- predict(m.num, type = "response")
  d$w <- ave(1 / d$p.den, d$id, FUN = cumprod)
  d$sw <- ave(d$p.num / d$p.den, d$id, FUN = cumprod)
  list(d = d, m.den = m.den)
}

ajuster <- function(d) {
  obs <- d[d$present == 1, ]
  f <- list(
    complet.ind = geeglm(Y ~ t * X, id = id, data = d, corstr = "independence"),
    nonpond.ind = geeglm(Y ~ t * X, id = id, data = obs, corstr = "independence"),
    nonpond.ech = geeglm(Y ~ t * X, id = id, data = obs, corstr = "exchangeable"),
    ipcw.ind    = geeglm(Y ~ t * X, id = id, data = obs, weights = w, corstr = "independence"),
    sipcw.ind   = geeglm(Y ~ t * X, id = id, data = obs, weights = sw, corstr = "independence"),
    ipcw.ech    = geeglm(Y ~ t * X, id = id, data = obs, weights = w, corstr = "exchangeable"))
  f
}

if (sys.nframe() == 0) {
  set.seed(7360)
  d <- simuler()
  cat("Sujets présents à chaque visite :\n"); print(tapply(d$present, d$t, sum))
  pw <- poids.ipcw(d); d <- pw$d
  cat("\nModèle de censure (dénominateur) :\n")
  print(round(summary(pw$m.den)$coefficients, 3))
  obs <- d[d$present == 1, ]
  cat("\nRésumé des poids (observations présentes) :\n")
  print(round(rbind(w = tapply(obs$w, obs$t, mean), sw = tapply(obs$sw, obs$t, mean),
                    w.max = tapply(obs$w, obs$t, max)), 2))
  cat("Somme des poids w par visite (cible n = 400) :\n"); print(round(tapply(obs$w, obs$t, sum), 1))
  f <- ajuster(d)
  tabc <- t(sapply(f, function(m) coef(m)))
  tabs <- t(sapply(f, function(m) summary(m)$coefficients[, "Std.err"]))
  cat("\nCoefficients :\n"); print(round(tabc, 3))
  cat("\nErreurs-types sandwich :\n"); print(round(tabs, 3))
  cat("\nSortie geeglm IPCW (indépendance) :\n")
  print(summary(f$ipcw.ind)$coefficients)

  ## Bootstrap par sujets (réestimation des poids à chaque réplique)
  set.seed(1)
  nb <- 300; ids <- unique(d$id)
  bt <- t(replicate(nb, {
    s <- sample(ids, replace = TRUE)
    db <- do.call(rbind, lapply(seq_along(s), function(k) { z <- d[d$id == s[k], 1:6]; z$id <- k; z }))
    db <- poids.ipcw(db)$d; ob <- db[db$present == 1, ]
    coef(geeglm(Y ~ t * X, id = id, data = ob, weights = w, corstr = "independence"))
  }))
  cat("\nET bootstrap (poids réestimés) vs sandwich (poids fixes) :\n")
  print(round(rbind(bootstrap = apply(bt, 2, sd), sandwich = tabs["ipcw.ind", ]), 4))

  ## Figure : moyennes par visite et poids
  cairo_pdf(file.path(figdir, "09_ipcw.pdf"), width = 8, height = 3.6)   # cairo : accents
  par(mfrow = c(1, 2), mar = c(4, 4, 2, 1), cex = 0.85)
  plot(NA, xlim = c(0, 3), ylim = c(9, 14), xlab = "Visite t", ylab = "Moyenne de Y",
       main = "Moyennes par visite et groupe")
  for (x in 0:1) {
    lines(0:3, beta.vrai[1] + beta.vrai[2] * 0:3 + beta.vrai[3] * x + beta.vrai[4] * x * 0:3,
          lwd = 2, col = "grey40", lty = x + 1)
    o <- obs[obs$X == x, ]
    lines(0:3, tapply(o$Y, o$t, mean), type = "b", pch = 1, col = "firebrick", lty = x + 1)
    lines(0:3, sapply(split(o, o$t), function(z) weighted.mean(z$Y, z$w)), type = "b", pch = 16,
          col = "steelblue", lty = x + 1)
  }
  legend("topleft", c("vraie moyenne", "observés (non pondérée)", "pondérée IPCW"),
         col = c("grey40", "firebrick", "steelblue"), pch = c(NA, 1, 16), lwd = c(2, 1, 1), bty = "n", cex = 0.8)
  boxplot(w ~ t, data = obs[obs$t > 0, ], xlab = "Visite t", ylab = "Poids IPCW", main = "Distribution des poids",
          col = "grey90", log = "y")
  dev.off()
}
