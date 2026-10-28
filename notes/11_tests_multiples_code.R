## STT-7360 H27 -- Section 11 : Tests multiples. Code R complet des exemples des notes.
## Exécuter depuis le dossier notes/ (les figures sont écrites dans figures/).
library(multtest)

## ---- Exemple papier (10 valeurs p) ----
## Exemple papier : 10 valeurs p
p <- c(0.001, 0.004, 0.006, 0.0072, 0.008, 0.028, 0.041, 0.120, 0.380, 0.750)
res <- mt.rawp2adjp(p, proc = c("Bonferroni", "SidakSS", "Holm", "Hochberg", "BH", "BY"))
print(round(res$adjp, 4))
print(colSums(res$adjp <= 0.05))

## ---- Exemple simulé avec dépendance (figure 11_nrejets) ----
set.seed(7360)
m <- 1000; n1 <- 10; n2 <- 10; taille.bloc <- 20; rho <- 0.5
m1 <- 100                                   # gènes différentiellement exprimés
delta <- c(rep(2, m1/2), rep(1.25, m1/2), rep(0, m - m1))
## dépendance : facteur commun par bloc de 20 gènes (corrélation 0,5)
gen <- function(n) {
  f <- matrix(rnorm(n * m / taille.bloc), n)[, rep(1:(m / taille.bloc), each = taille.bloc)]
  t(sqrt(rho) * f + sqrt(1 - rho) * matrix(rnorm(n * m), n))
}
X <- cbind(gen(n1), gen(n2) + delta)        # m x (n1 + n2)
cl <- rep(0:1, c(n1, n2))
h0 <- delta == 0
tw <- mt.teststat(X, cl, test = "t")        # Welch (variances inégales)
rawp <- 2 * pt(-abs(tw), df = n1 + n2 - 2)
procs <- c("Bonferroni", "Holm", "Hochberg", "SidakSS", "SidakSD", "BH", "BY", "ABH", "TSBH")
res <- mt.rawp2adjp(rawp, proc = procs, alpha = 0.05)
adjp <- res$adjp[order(res$index), ]
rej <- adjp <= 0.05
tab <- rbind(R = colSums(rej), V = colSums(rej & h0), S = colSums(rej & !h0))
print(tab)
cat("m0 estimé (ABH) :", res$h0.ABH, " ; (TSBH) :", res$h0.TSBH, "\n")
lambda <- 0.5
cat("Storey pi0(0,5) :", round(sum(rawp > lambda) / (m * (1 - lambda)), 3), "\n")
pdf("figures/11_nrejets.pdf", width = 6.5, height = 4.6, pointsize = 10)
par(mar = c(4, 4, 1, 1))
choix <- c("rawp", "Bonferroni", "Holm", "Hochberg", "BH", "BY", "TSBH_0.05")
mt.plot(adjp[, choix], tw, plottype = "rvsa", proc = choix,
        alpha = seq(0, 0.2, length = 401), leg = c(0.003, 262),
        col = c(1, 2, 2, 2, 4, 4, 3), lty = c(1, 1, 2, 3, 1, 3, 2), lwd = 2,
        ylim = c(0, 260))
abline(h = m1, col = "grey60", lty = 3)
dev.off()

## ---- Figure des seuils (11_seuils) ----
m <- length(rawp); alpha <- 0.05
ps <- sort(rawp); j <- 1:m; k <- 120
cairo_pdf("figures/11_seuils.pdf", width = 6.5, height = 4.3)
par(mar = c(4, 4.2, 1, 1))
plot(j[1:k], ps[1:k], log = "y", pch = ifelse(h0[order(rawp)][1:k], 1, 16), cex = 0.7,
     xlab = expression("Rang " * j), ylab = expression("Valeur p ordonnée " * p[(j)]),
     ylim = c(1e-7, 0.1))
lines(j[1:k], rep(alpha / m, k), col = "black", lwd = 2, lty = 2)
lines(j[1:k], alpha / (m - j[1:k] + 1), col = "red", lwd = 2)
lines(j[1:k], j[1:k] * alpha / m, col = "blue", lwd = 2)
lines(j[1:k], j[1:k] * alpha / (m * sum(1 / (1:m))), col = "blue", lwd = 2, lty = 3)
abline(v = c(30.5, 76.5), col = "grey70", lty = 3)
legend("bottomright", bty = "n", cex = 0.8,
       legend = c(expression("Bonferroni " * alpha / m), expression("Holm/Hochberg " * alpha/(m-j+1)),
                  expression("BH " * j * alpha / m), expression("BY " * j * alpha / (m * c(m))),
                  "H0 vraie", "H0 fausse"),
       col = c(1, 2, 4, 4, 1, 1), lty = c(2, 1, 1, 3, NA, NA), lwd = c(2, 2, 2, 2, NA, NA), pch = c(NA, NA, NA, NA, 1, 16))
dev.off()

## ---- Procédures conjointes (permutation) et MTP (bootstrap) ----
maxT <- mt.maxT(X, cl, test = "t", side = "abs", B = 2000)
minP <- mt.minP(X, cl, test = "t", side = "abs", B = 2000)
head(maxT, 3); head(minP, 3)
minP50 <- mt.minP(X, cl, test = "t", side = "abs", B = 50000)
c(maxT = sum(maxT$adjp <= 0.05), minP = sum(minP$adjp <= 0.05),
  holm.perm = sum(p.adjust(maxT$rawp, "holm") <= 0.05), minP.50000 = sum(minP50$adjp <= 0.05))
c(faux.maxT = sum(maxT$adjp <= 0.05 & h0[maxT$index]))
boot <- MTP(X = X, Y = cl, test = "t.twosamp.unequalvar", alternative = "two.sided",
            typeone = "fwer", alpha = 0.05, nulldist = "boot.cs", B = 1000,
            method = "sd.maxT", seed = 2027)
summary(boot)

## ---- Permutation quand l échangeabilité échoue ----
set.seed(2027); B <- 999; nsim <- 2000
taux <- function(n1, n2, s1, s2) {        # taux de rejet empirique au seuil 5 %
  G <- t(replicate(B, sample(n1 + n2) <= n1)) * 1   # B permutations (indicatrices)
  stats <- function(z, g) {               # différence de moyennes et statistique de Welch
    m1 <- (g %*% z) / n1; m2 <- ((1 - g) %*% z) / n2
    v1 <- (g %*% z^2 - n1 * m1^2) / (n1 - 1); v2 <- ((1 - g) %*% z^2 - n2 * m2^2) / (n2 - 1)
    cbind(dif = m1 - m2, welch = (m1 - m2) / sqrt(v1 / n1 + v2 / n2)) }
  g0 <- matrix(rep(c(1, 0), c(n1, n2)), 1)
  rej <- replicate(nsim, { z <- c(rnorm(n1, 0, s1), rnorm(n2, 0, s2))
    t0 <- stats(z, g0); tb <- stats(z, G)
    (1 + colSums(abs(tb) >= rep(abs(t0), each = B))) / (B + 1) <= 0.05 })
  setNames(rowMeans(rej), c("dif", "welch"))
}
rbind("n1=n2=15, s1=s2=1" = taux(15, 15, 1, 1), "n1=n2=15, s1=4, s2=1" = taux(15, 15, 4, 1),
      "n1=5, n2=25, s1=s2=1" = taux(5, 25, 1, 1), "n1=5, n2=25, s1=4, s2=1" = taux(5, 25, 4, 1),
      "n1=5, n2=25, s1=1, s2=4" = taux(5, 25, 1, 4))

## ---- Taux d erreur empiriques (simulation) ----
## Taux d'erreur empiriques : indépendance vs dépendance par blocs
set.seed(1)
m <- 200; m1 <- 20; n <- 10; nsim <- 2000; alpha <- 0.05
delta <- c(rep(1.5, m1), rep(0, m - m1)); h0 <- delta == 0
une.sim <- function(rho, taille.bloc = 20) {
  gen <- function(n) {
    f <- matrix(rnorm(n * m / taille.bloc), n)[, rep(1:(m / taille.bloc), each = taille.bloc)]
    sqrt(rho) * f + sqrt(1 - rho) * matrix(rnorm(n * m), n)
  }
  x1 <- gen(n); x2 <- gen(n) + rep(delta, each = n)
  t <- (colMeans(x2) - colMeans(x1)) / sqrt((apply(x1, 2, var) + apply(x2, 2, var)) / n)
  p <- 2 * pt(-abs(t), 2 * n - 2)
  sapply(c("bonferroni", "holm", "hochberg", "BH", "BY"), function(meth) {
    r <- p.adjust(p, meth) <= alpha; V <- sum(r & h0); R <- sum(r)
    c(FWER = V > 0, FDP = if (R > 0) V / R else 0, puiss = mean(r[!h0]))
  })
}
for (rho in c(0, 0.8)) {
  out <- replicate(nsim, une.sim(rho))
  cat("rho =", rho, "\n"); print(round(apply(out, 1:2, mean), 3))
}
