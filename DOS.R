E_fermi <- -0.180000 * 27.2114

load_pdos <- function(file) {
  data <- read.table(file, header = FALSE, comment.char = "#")
  E <- data$V2 * 27.2114 - E_fermi
  occ <- data$V3
  s <- data$V4
  p <- data$V5
  d <- data$V6
  return(list(E=E, occ=occ, s=s, p=p, d=d))
}
cs <- load_pdos("SSCP_single_point_test_1p2nm-k1-1.pdos")
br <- load_pdos("SSCP_single_point_test_1p2nm-k2-1.pdos")
pb <- load_pdos("SSCP_single_point_test_1p2nm-k3-1.pdos")
E_grid <- seq(-10, 10, length.out = 2000)
sigma <- 0.2   # un po’ più largo (meglio per Γ-point)

gaussian <- function(E, Ei) {
  exp(-(E - Ei)^2 / (2 * sigma^2))
}

build_dos <- function(E_levels, weights) {
  DOS <- rep(0, length(E_grid))
  for (i in 1:length(E_levels)) {
    DOS <- DOS + weights[i] * gaussian(E_grid, E_levels[i])
  }
  return(DOS)
}

DOS_cs <- build_dos(cs$E, cs$s + cs$p + cs$d)
DOS_br <- build_dos(br$E, br$s + br$p + br$d)
DOS_pb <- build_dos(pb$E, pb$s + pb$p + pb$d)
DOS_tot <- DOS_cs + DOS_br + DOS_pb

DOS_cs_s <- build_dos(cs$E, cs$s)
DOS_br_s <- build_dos(br$E, br$s)
DOS_pb_s <- build_dos(pb$E, pb$s)
DOS_cs_p <- build_dos(cs$E, cs$p)
DOS_br_p <- build_dos(br$E, br$p)
DOS_pb_p <- build_dos(pb$E, pb$p)
DOS_cs_d <- build_dos(cs$E, cs$d)
DOS_br_d <- build_dos(br$E, br$d)
DOS_pb_d <- build_dos(pb$E, pb$d)

#-------------------------------------------------------------------------------
windows()
plot(E_grid, DOS_tot, type="l", col="black", lwd=2,
     ylab="Energy [eV]", xlab="DOS [eV]")
title(expression("Density of states of CsPbBr"[3]))
lines(E_grid, DOS_cs_s, col="darkgreen", lwd=2)
lines(E_grid, DOS_cs_p, col="forestgreen", lwd=2)
lines(E_grid, DOS_cs_d, col="limegreen", lwd=2)
lines(E_grid, DOS_pb_s, col="darkred", lwd=2)
lines(E_grid, DOS_pb_p, col="firebrick", lwd=2)
lines(E_grid, DOS_pb_d, col="tomato", lwd=2)
lines(E_grid, DOS_br_s, col="navy", lwd=2)
lines(E_grid, DOS_br_p, col="royalblue", lwd=2)
lines(E_grid, DOS_br_d, col="skyblue", lwd=2)
legend("topright",
       legend=c("Total","Cs-s","Cs-p","Cs-d","Pb-s","Pb-p","Pb-d", "Br-s","Br-p","Br-d"),
       col=c("black","darkgreen","forestgreen","limegreen", "darkred","firebrick","tomato", "navy","royalblue","skyblue"),
       lty=1, lwd=2)
#-------------------------------------------------------------------------------
windows()
plot(E_grid, DOS_tot, type="l", col="black", lwd=2,
     ylab="Energy [eV]", xlab="DOS [eV]")
title(expression("Contribution of Cs on DOS of CsPbBr"[3]))
lines(E_grid, DOS_cs_s, col="darkgreen", lwd=2)
lines(E_grid, DOS_cs_p, col="forestgreen", lwd=2)
lines(E_grid, DOS_cs_d, col="limegreen", lwd=2)
abline(v = 0, col="orange", lwd=2)
legend("topright",
       legend=c("Total","Cs-s","Cs-p","Cs-d","Fermi energy"),
       col=c("black","darkgreen","forestgreen","limegreen", "orange"),
       lty=1, lwd=2)
#-------------------------------------------------------------------------------
windows()
plot(E_grid, DOS_tot, type="l", col="black", lwd=2,
     ylab="Energy [eV]", xlab="DOS [eV]")
title(expression("Contribution of Pb on DOS of CsPbBr"[3]))
lines(E_grid, DOS_pb_s, col="darkred", lwd=2)
lines(E_grid, DOS_pb_p, col="firebrick", lwd=2)
lines(E_grid, DOS_pb_d, col="tomato", lwd=2)
abline(v = 0, col="orange", lwd=2)
legend("topright",
       legend=c("Total","Pb-s","Pb-p","Pb-d","Fermi energy"),
       col=c("black","darkred","firebrick","tomato", "orange"),
       lty=1, lwd=2)
#-------------------------------------------------------------------------------
windows()
plot(E_grid, DOS_tot, type="l", col="black", lwd=2,
     ylab="Energy [eV]", xlab="DOS [eV]")
title(expression("Contribution of Br on DOS of CsPbBr"[3]))
lines(E_grid, DOS_br_s, col="navy", lwd=2)
lines(E_grid, DOS_br_p, col="royalblue", lwd=2)
lines(E_grid, DOS_br_d, col="skyblue", lwd=2)
abline(v = 0, col="orange", lwd=2)
legend("topright",
       legend=c("Total", "Br-s","Br-p","Br-d", "Fermi energy"),
       col=c("black","navy","royalblue","skyblue", "orange"),
       lty=1, lwd=2)
#-------------------------------------------------------------------------------
windows()
plot(E_grid, DOS_tot, type="l", col="black", lwd=2,
     ylab="Energy [eV]", xlab="DOS [eV]")
title(expression("Density of states of CsPbBr"[3]))
lines(E_grid, DOS_pb, col="red", lwd=2)
lines(E_grid, DOS_br, col="blue", lwd=2)
lines(E_grid, DOS_cs, col="limegreen", lwd=2)
abline(v = 0, col="orange", lwd=2)

pts <- locator(2)
VBM <- pts$x[1]
CBM <- pts$x[2]
band_gap <- CBM - VBM

cat("VBM =", VBM, "eV\n")
cat("CBM =", CBM, "eV\n")
cat("Band gap =", band_gap, "eV\n")

#points(pts$x, pts$y, col="black", pch=19, cex=1)
abline(v=VBM, col="brown", lty=2)
abline(v=CBM, col="violet", lty=2)
legend("topright",
       legend = c("Total", "Cs", "Pb", "Br", "Fermi energy", "VBM", "CBM",
                  bquote(E[gap] == .(sprintf("%.2f", band_gap)) ~ "eV")),
       col = c("black", "limegreen", "red", "blue", "orange", "brown", "violet", NA),
       lty = c(1, 1, 1, 1, 1, 2, 2, NA),
       lwd = c(2, 2, 2, 2, 2, 1, 1, NA),
       bty = "n")