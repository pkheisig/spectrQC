library(devtools)
load_all()

# Load the reference matrix we built earlier
M <- as.matrix(data.table::fread("reference_matrix.csv", select = 2:52))
rownames(M) <- data.table::fread("reference_matrix.csv")$file

print("Dimensions of Reference Matrix (M):")
print(dim(M))

# Derive OLS unmixing matrix
W <- derive_unmixing_matrix(M, method = "OLS")

print("Dimensions of Unmixing Matrix (W):")
print(dim(W))

# Verification: Unmixing M itself should yield the Identity matrix (approx)
# Unmixed_Data = Raw_Data %*% t(W)
Identity_Check <- M %*% t(W)
print("Diagonal of (M %*% W^T) - should be close to 1:")
print(diag(Identity_Check))

# Save it
save_unmixing_matrix(W, "unmixing_matrix_test.csv")

# Test manual unmixing vs calc_residuals
# (Using OLS to compare static matrix vs dynamic calc)
library(flowCore)
files_raw <- list.files("raw", pattern = ".*.fcs$", full.names = TRUE)[1]
ff <- read.FCS(files_raw)
Y <- exprs(ff)[1:10, colnames(M)] # 10 cells

# Method 1: Matrix multiplication
unmixed_manual <- Y %*% t(W)

# Method 2: calc_residuals
unmixed_calc <- calc_residuals(ff[1:10, ], M, method = "OLS")
unmixed_calc_markers <- as.matrix(unmixed_calc[, rownames(M)])

print("Difference between manual and calc_residuals (should be ~0):")
print(max(abs(unmixed_manual - unmixed_calc_markers)))
