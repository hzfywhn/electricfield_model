source(file = "LatticeKriging.r")

# predict potentials based on line-of-sight electric fields at high latitudes

hemi <- "north"

# resolution in degree, defined as the spacing between basis functions at each level
res <- 5

# input and output files
input <- paste(hemi, "_in.nc", sep = "")
output_grid <- paste(hemi, "_out_grid.nc", sep = "")
output <- paste(hemi, "_out_", res, "d.nc", sep = "")

# whether mean will be calculated
pred_m <- TRUE

# predicting variance is computation expensive
# it should be turned on only for demonstration purposes
pred_se <- FALSE

# modeling domain, (-range, range), this is suggested to be larger than the output grid
range <- pi/3

# earth radius, used for re-scaling between electric field and potential
a <- 6.371e6

# read input data, definitions of variables are given in prepare_input.py
nc <- ncdf4::nc_open(filename = input)
time <- ncdf4::ncvar_get(nc = nc, varid = "time")
nrec <- ncdf4::ncvar_get(nc = nc, varid = "nrec")
xi <- ncdf4::ncvar_get(nc = nc, varid = "x")
yi <- ncdf4::ncvar_get(nc = nc, varid = "y")
azim <- ncdf4::ncvar_get(nc = nc, varid = "azim")
Elos <- ncdf4::ncvar_get(nc = nc, varid = "Elos")
Elos_sd <- ncdf4::ncvar_get(nc = nc, varid = "Elos_sd")
Elos_prior <- ncdf4::ncvar_get(nc = nc, varid = "Elos_prior")
ncdf4::nc_close(nc = nc)

# read the designated output grid
nc <- ncdf4::nc_open(filename = output_grid)
xo <- ncdf4::ncvar_get(nc = nc, varid = "x")
yo <- ncdf4::ncvar_get(nc = nc, varid = "y")
prior <- ncdf4::ncvar_get(nc = nc, varid = "prior") / a
ncdf4::nc_close(nc = nc)

# set up the basis functions
res <- res * pi/180
nres <- length(res)
x0 <- vector(mode = "list", length = nres)
connect <- vector(mode = "list", length = nres)
width <- array(dim = nres)
alpha <- array(dim = nres)

for (ires in 1: nres) {
# within each level, the basis functions are equally spaced
    xy <- seq(from = -range[ires], to = range[ires], by = res[ires])
    nxy <- length(xy)
    loc <- array(dim = c(nxy^2, 2))
    con <- array(dim = c(nxy^2, 4))
    for (j in 1: nxy) {
        for (i in 1: nxy) {
            idx <- (j-1)*nxy + i
            loc[idx, 2] <- xy[j]
            loc[idx, 1] <- xy[i]
            if (j >= 2) con[idx, 1] <- (j-2)*nxy + i
            if (i >= 2) con[idx, 2] <- (j-1)*nxy + i-1
            if (j <= nxy-1) con[idx, 3] <- j*nxy + i
            if (i <= nxy-1) con[idx, 4] <- (j-1)*nxy + i+1
        }
    }
    x0[[ires]] <- loc
    connect[[ires]] <- con
# increase the width of basis functions to allow overlapping
    width[ires] <- res[ires] * 2.5
# weighting of each level can be adjusted based on actual problems
    alpha[ires] <- 1 / nres
}

# center weight is usually 2^d+0.01 where d is dimension
centerweight <- 4.01

# initial guesss of the marginal variance, will be re-estimated later
rho <- 1

# create the grid to be predicted
coor <- expand.grid(xo, yo)
xnew <- cbind(coor[, 1], coor[, 2])

# save essential modeling parameters for diagnostics
ntime <- length(time)
nx <- length(xo)
ny <- length(yo)
lambda <- array(dim = ntime)
d <- array(dim = ntime)
rhoMLE <- array(dim = ntime)
like <- array(dim = ntime)
m <- array(dim = c(nx, ny, ntime))
se <- array(dim = c(nx, ny, ntime))
for (itime in 1: ntime) {
    n <- nrec[itime]
    x <- cbind(xi[1: n, itime], yi[1: n, itime])
    t <- azim[1: n, itime]
# direction cosine goes with cos(<v,x>), cos(<v,y>)
    dircos <- cbind(cos(t), sin(t))
    y <- Elos[1: n, itime]
    w <- Elos_sd[1: n, itime]
    Z <- Elos_prior[1: n, itime]
    Znew <- array(data = prior[, , itime])

    lk <- LatticeKriging(x, dircos, y, w, Z,
        x0, width, alpha, connect, centerweight, rho,
        xnew, Znew, pred_m, pred_se)

    lambda[itime] <- lk$lambda
    rhoMLE[itime] <- lk$rhoMLE
    like[itime] <- lk$like
    d[itime] <- lk$d
    if (pred_m) m[, , itime] <- array(data = lk$m, dim = c(nx, ny))
    if (pred_se) se[, , itime] <- array(data = lk$se, dim = c(nx, ny))
}

# save outputs
dim_time <- ncdf4::ncdim_def(name = "time", units = "", vals = time)
dim_x <- ncdf4::ncdim_def(name = "x", units = "", vals = xo)
dim_y <- ncdf4::ncdim_def(name = "y", units = "", vals = yo)

lambda_out <- ncdf4::ncvar_def(name = "lambda", units = "", dim = dim_time)
rho_out <- ncdf4::ncvar_def(name = "rho", units = "", dim = dim_time)
like_out <- ncdf4::ncvar_def(name = "like", units = "", dim = dim_time)
d_out <- ncdf4::ncvar_def(name = "d", units = "", dim = dim_time)

# these are mandatory outputs
vars_out <- list(lambda_out, rho_out, like_out, d_out)

# output mean and variance if flagged
if (pred_m) {
    m_out <- ncdf4::ncvar_def(name = "m", units = "", dim = list(dim_x, dim_y, dim_time))
    vars_out <- c(vars_out, list(m_out))
}
if (pred_se) {
    se_out <- ncdf4::ncvar_def(name = "se", units = "", dim = list(dim_x, dim_y, dim_time))
    vars_out <- c(vars_out, list(se_out))
}

nc <- ncdf4::nc_create(filename = output, vars = vars_out, force_v4 = TRUE)
ncdf4::ncvar_put(nc = nc, varid = lambda_out, vals = lambda)
ncdf4::ncvar_put(nc = nc, varid = rho_out, vals = rhoMLE)
ncdf4::ncvar_put(nc = nc, varid = like_out, vals = like)
ncdf4::ncvar_put(nc = nc, varid = d_out, vals = d)
if (pred_m) ncdf4::ncvar_put(nc = nc, varid = m_out, vals = m * a)
if (pred_se) ncdf4::ncvar_put(nc = nc, varid = se_out, vals = se * a)
ncdf4::nc_close(nc = nc)