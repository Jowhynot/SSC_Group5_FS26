energy_cubic_SVZ <-read.table("G5_bandstrucutre_cubic_CsPbBr3_SZV_1.bs")
energy_cubic_DVZ <-read.table("G5_bandstrucutre_cubic_CsPbBr3_DZV_1.bs")
energy_ortho_SVZ <-read.table("G5_bandstrucutre_ortho_CsPbBr3_SZV_1.bs")
energy_ortho_DVZ <-read.table("G5_bandstrucutre_ortho_CsPbBr3_DZV_1.bs")
energy_tetra_SVZ <-read.table("G5_bandstrucutre_tetragonal_CsPbBr3_SZV_1.bs")
energy_tetra_DVZ <-read.table("G5_bandstrucutre_tetragonal_CsPbBr3_DVZ_1.bs")
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
x_cubic_SVZ <- rep(1:141, each = 21)
y_cubic_SVZ <- energy_cubic_SVZ[,2]
occupation_band <- energy_cubic_SVZ[,3]

mat_cubic_SVZ <- matrix(y_cubic_SVZ, nrow = 141, ncol = 21, byrow = TRUE)
mat_VB_cubic_SVZ <- mat_cubic_SVZ[, 18:21]
mat_CB_cubic_SVZ <- mat_cubic_SVZ[, 1:17]
VB_cubic_SVZ <- min(mat_VB_cubic_SVZ)
CB_cubic_SVZ <- max(mat_CB_cubic_SVZ)
bandgap_cubic_SVZ <- VB_cubic_SVZ - CB_cubic_SVZ

x_cubic_DVZ <- rep(1:141, each = 37)
y_cubic_DVZ <- energy_cubic_DVZ[,2]
occupation_band <- energy_cubic_DVZ[,3]

mat_cubic_DVZ <- matrix(y_cubic_DVZ, nrow = 141, ncol = 37, byrow = TRUE)
mat_VB_cubic_DVZ <- mat_cubic_DVZ[, 18:37]
mat_CB_cubic_DVZ <- mat_cubic_DVZ[, 1:17]
VB_cubic_DVZ <- min(mat_VB_cubic_DVZ)
CB_cubic_DVZ <- max(mat_CB_cubic_DVZ)
bandgap_cubic_DVZ <- VB_cubic_DVZ - CB_cubic_DVZ
#-------------------------------------------------------------------------------
x_ortho_SVZ <- rep(1:301, each = 84)
y_ortho_SVZ <- energy_ortho_SVZ[,2]
occupation_band <- energy_ortho_SVZ[,3]
occupation_band

mat_ortho_SVZ <- matrix(y_ortho_SVZ, nrow = 301, ncol = 84, byrow = TRUE)
mat_VB_ortho_SVZ <- mat_ortho_SVZ[, 69:84]
mat_CB_ortho_SVZ <- mat_ortho_SVZ[, 1:68]
VB_ortho_SVZ <- min(mat_VB_ortho_SVZ)
CB_ortho_SVZ <- max(mat_CB_ortho_SVZ)
bandgap_ortho_SVZ <- VB_ortho_SVZ - CB_ortho_SVZ

x_ortho_DVZ <- rep(1:301, each = 88)
y_ortho_DVZ <- energy_ortho_DVZ[,2]
occupation_band <- energy_ortho_DVZ[,3]
occupation_band

mat_ortho_DVZ <- matrix(y_ortho_DVZ, nrow = 301, ncol = 88, byrow = TRUE)
mat_VB_ortho_DVZ <- mat_ortho_DVZ[, 69:88]
mat_CB_ortho_DVZ <- mat_ortho_DVZ[, 1:68]
VB_ortho_DVZ <- min(mat_VB_ortho_DVZ)
CB_ortho_DVZ <- max(mat_CB_ortho_DVZ)
bandgap_ortho_DVZ <- VB_ortho_DVZ - CB_ortho_DVZ
#-------------------------------------------------------------------------------
x_tetra_SVZ <- rep(1:221, each = 34)
y_tetra_SVZ <- energy_tetra_SVZ[,2]
occupation_band <- energy_tetra_SVZ[,3]
occupation_band

mat_tetra_SVZ <- matrix(y_tetra_SVZ, nrow = 221, ncol = 34, byrow = TRUE)
mat_VB_tetra_SVZ<- mat_tetra_SVZ[, 27:34]
mat_CB_tetra_SVZ <- mat_tetra_SVZ[, 1:26]
VB_tetra_SVZ <- min(mat_VB_tetra_SVZ)
CB_tetra_SVZ <- max(mat_CB_tetra_SVZ)
bandgap_tetra_SVZ <- VB_tetra_SVZ - CB_tetra_SVZ

x_tetra_DVZ <- rep(1:221, each = 47)
y_tetra_DVZ <- energy_tetra_DVZ[,2]
occupation_band <- energy_tetra_DVZ[,3]
occupation_band

mat_tetra_DVZ <- matrix(y_tetra_DVZ, nrow = 221, ncol = 47, byrow = TRUE)
mat_VB_tetra_DVZ<- mat_tetra_DVZ[, 27:37]
mat_CB_tetra_DVZ <- mat_tetra_DVZ[, 1:26]
VB_tetra_DVZ <- min(mat_VB_tetra_DVZ)
CB_tetra_DVZ <- max(mat_CB_tetra_DVZ)
bandgap_tetra_DVZ <- VB_tetra_DVZ - CB_tetra_DVZ
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
bandgap_cubic_SVZ
bandgap_cubic_DVZ
bandgap_ortho_SVZ
bandgap_ortho_DVZ
bandgap_tetra_SVZ
bandgap_tetra_DVZ
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
windows()
par(mfrow = c(1, 2))
plot(x_cubic_SVZ, y_cubic_SVZ, pch = 20, cex = 0.3, ylim = c(-25,30),
     xlab = "", ylab = "Energy (eV)", xaxt = "n",
     main = expression("Band structure of cubic CsPbBr"[3] ~ "(SZV)"))
abline(v = c(1, 21, 41, 61, 81, 101, 121, 141), col = "grey", lty = 2)
axis(1,
     at = c(1, 21, 41, 61, 81, 101, 121, 141),
     labels = c(expression(Gamma), "X", "M", expression(Gamma), "R", "X", "M", "R"))
points(81, VB_cubic_SVZ, col = "blue", pch = 19, cex = 1)
points(81, CB_cubic_SVZ, col = "red", pch = 19, cex = 1)
legend("topright",
       legend = c(
         "Valence Band Min",
         "Conduction Band Max",
         paste0("Band gap = ", round(abs(bandgap_cubic_SVZ), 3), " eV")
       ),
       col = c("blue", "red", "black"),
       pch = c(19, 19, NA),
       bty = "n"
)

plot(x_cubic_DVZ, y_cubic_DVZ, pch = 20, cex = 0.3, ylim = c(-25,30),
     xlab = "", ylab = "Energy (eV)", xaxt = "n",
     main = expression("Band structure of cubic CsPbBr"[3] ~ "(DZV)"))

abline(v = c(1, 21, 41, 61, 81, 101, 121, 141), col = "grey", lty = 2)
axis(1,
     at = c(1, 21, 41, 61, 81, 101, 121, 141),
     labels = c(expression(Gamma), "X", "M", expression(Gamma), "R", "X", "M", "R"))

points(81, VB_cubic_DVZ, col = "blue", pch = 19, cex = 1)
points(81, CB_cubic_DVZ, col = "red", pch = 19, cex = 1)

legend("topright",
       legend = c(
         "Valence Band Min",
         "Conduction Band Max",
         paste0("Band gap = ", round(abs(bandgap_cubic_DVZ), 3), " eV")
       ),
       col = c("blue", "red", "black"),
       pch = c(19, 19, NA),
       bty = "n"
)

par(mfrow = c(1, 1))
#-------------------------------------------------------------------------------
windows()
par(mfrow = c(1, 2))
plot(x_ortho_SVZ, y_ortho_SVZ, pch = 20, cex = 0.3, ylim = c(-25,15),
     xlab = "", ylab = "Energy (eV)", xaxt = "n", main = expression("Band structure of orthorhombic CsPbBr"[3] ~ "(SZV)"))
abline(v = c(1, 21, 41, 61, 81, 101, 121, 141, 161, 181, 201, 221, 241, 261, 281, 301), col = "grey", lty = 2)
axis(1,
     at = c(1, 21, 41, 61, 81, 101, 121, 141, 161, 181, 201, 221, 241, 261, 281, 301),
     labels = c(expression(Gamma), "X", "S", "Y", expression(Gamma), "Z", "U", "R", "T", "Z", "Y", "T", "U", "X", "S", "R"))
points(81, VB_ortho_SVZ, col = "blue", pch = 19, cex=1)
points(81, CB_ortho_SVZ, col = "red", pch = 19, cex=1)
legend("topright",
       legend = c(
         "Valence Band Min",
         "Conduction Band Max",
         paste0("Band gap = ", round(abs(bandgap_ortho_SVZ), 3), " eV")
       ),
       col = c("blue", "red", "black"),
       pch = c(19, 19, NA),
       bty = "n"
)

plot(x_ortho_DVZ, y_ortho_DVZ, pch = 20, cex = 0.3, ylim = c(-25,15),
     xlab = "", ylab = "Energy (eV)", xaxt = "n", main = expression("Band structure of orthorhombic CsPbBr"[3] ~ "(DZV)"))
abline(v = c(1, 21, 41, 61, 81, 101, 121, 141, 161, 181, 201, 221, 241, 261, 281, 301), col = "grey", lty = 2)
axis(1,
     at = c(1, 21, 41, 61, 81, 101, 121, 141, 161, 181, 201, 221, 241, 261, 281, 301),
     labels = c(expression(Gamma), "X", "S", "Y", expression(Gamma), "Z", "U", "R", "T", "Z", "Y", "T", "U", "X", "S", "R"))
points(81, VB_ortho_DVZ, col = "blue", pch = 19, cex=1)
points(81, CB_ortho_DVZ, col = "red", pch = 19, cex=1)
legend("topright",
       legend = c(
         "Valence Band Min",
         "Conduction Band Max",
         paste0("Band gap = ", round(abs(bandgap_ortho_DVZ), 3), " eV")
       ),
       col = c("blue", "red", "black"),
       pch = c(19, 19, NA),
       bty = "n"
)

par(mfrow = c(1, 1))
#-------------------------------------------------------------------------------
windows()
par(mfrow = c(1, 2))
plot(x_tetra_SVZ, y_tetra_SVZ, pch = 20, cex = 0.3, ylim = c(-25,15),
     xlab = "", ylab = "Energy (eV)", xaxt = "n", main = expression("Band structure of tetragonal CsPbBr"[3] ~ "(SZV)"))
abline(v = c(1, 21, 41, 61, 81, 101, 121, 141, 161, 181, 201, 221), col = "grey", lty = 2)
axis(1,
     at = c(1, 21, 41, 61, 81, 101, 121, 141, 161, 181, 201, 221),
     labels = c(expression(Gamma), "X", "M", expression(Gamma), "Z", "R", "A", "Z", "X", "R", "M", "A"))
points(61, VB_tetra_SVZ, col = "blue", pch = 19, cex=1)
points(61, CB_tetra_SVZ, col = "red", pch = 19, cex=1)
legend("topright",
       legend = c(
         "Valence Band Min",
         "Conduction Band Max",
         paste0("Band gap = ", round(abs(bandgap_tetra_SVZ), 3), " eV")
       ),
       col = c("blue", "red", "black"),
       pch = c(19, 19, NA),
       bty = "n"
)

plot(x_tetra_DVZ, y_tetra_DVZ, pch = 20, cex = 0.3, ylim = c(-25,15),
     xlab = "", ylab = "Energy (eV)", xaxt = "n", main = expression("Band structure of tetragonal CsPbBr"[3] ~ "(DZV)"))
abline(v = c(1, 21, 41, 61, 81, 101, 121, 141, 161, 181, 201, 221), col = "grey", lty = 2)
axis(1,
     at = c(1, 21, 41, 61, 81, 101, 121, 141, 161, 181, 201, 221),
     labels = c(expression(Gamma), "X", "M", expression(Gamma), "Z", "R", "A", "Z", "X", "R", "M", "A"))
points(61, VB_tetra_DVZ, col = "blue", pch = 19, cex=1)
points(61, CB_tetra_DVZ, col = "red", pch = 19, cex=1)
legend("topright",
       legend = c(
         "Valence Band Min",
         "Conduction Band Max",
         paste0("Band gap = ", round(abs(bandgap_tetra_DVZ), 3), " eV")
       ),
       col = c("blue", "red", "black"),
       pch = c(19, 19, NA),
       bty = "n"
)

par(mfrow = c(1, 1))