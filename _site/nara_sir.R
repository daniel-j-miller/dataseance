library(nara); library(grid)

grid_n  <- 200
state   <- matrix(0L, grid_n, grid_n)        # 0=S, 1=I, 2=R
state[sample(length(state), 20)] <- 1L       # a few initial seeds

p_inf    <- 0.25     # per infected neighbour, per step
p_rec    <- 0.10     # I -> R  (mean infectious period ~10 steps)
p_wane   <- 0.008    # R -> S  (mean immunity ~125 steps)  <- the new bit
p_import <- 2e-5     # spontaneous infection of S cells (keeps it alive)
n_steps  <- 3000

col_int <- function(col) nr_new(1, 1, fill = col)[1]
pal <- c(col_int("grey90"), col_int("firebrick"), col_int("steelblue"))
lattice <- nr_new(grid_n, grid_n)

# scrolling time-series strip
TW <- 400; TH <- 80
trace <- nr_new(TW, TH, "white")
ty <- function(frac) round((1 - frac) * (TH - 2))

shift <- function(m, dr, dc) {
  out <- matrix(0L, nrow(m), ncol(m))
  r <- max(1, 1 + dr):min(nrow(m), nrow(m) + dr)
  c <- max(1, 1 + dc):min(ncol(m), ncol(m) + dc)
  out[r, c] <- m[r - dr, c - dc]; out
}

windows(width = 6, height = 7); dev.control("inhibit")

for (k in seq_len(n_steps)) {
  inf <- (state == 1L) * 1L
  nb  <- shift(inf, 1, 0) + shift(inf, -1, 0) + shift(inf, 0, 1) + shift(inf, 0, -1)
  u   <- runif(length(state))

  # all transitions decided from the *current* state, then applied together
  p_SI     <- 1 - (1 - p_inf)^nb * (1 - p_import)
  new_inf  <- state == 0L & u < p_SI
  new_rec  <- state == 1L & u < p_rec
  new_wane <- state == 2L & u < p_wane
  state[new_inf]  <- 1L
  state[new_rec]  <- 2L
  state[new_wane] <- 0L

  # --- render lattice
  lattice[] <- pal[state + 1L]

  # --- update trace: clear this column, then plot S, I, R fractions
  col  <- (k - 1) %% TW
  frac <- tabulate(state + 1L, nbins = 3) / length(state)
  nr_rect(trace, x = col, y = 0, w = 1, h = TH, fill = "white")
  nr_rect(trace, x = col, y = ty(frac[1]), w = 1, h = 2, fill = "grey50")
  nr_rect(trace, x = col, y = ty(frac[2]), w = 1, h = 2, fill = "firebrick")
  nr_rect(trace, x = col, y = ty(frac[3]), w = 1, h = 2, fill = "steelblue")
  nr_rect(trace, x = (col + 1) %% TW, y = 0, w = 1, h = TH, fill = "black")  # cursor

  dev.hold()
  grid.newpage()
  grid.raster(lattice, y = 0.6, height = 0.75, interpolate = FALSE)
  grid.raster(trace,   y = 0.1, width = 0.95, height = 0.15, interpolate = FALSE)
  dev.flush()
}