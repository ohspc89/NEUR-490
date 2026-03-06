# --- Packages ---
if (!requireNamespace("crqa", quietly=T)) install.packages("crqa")
if (!requireNamespace("dplyr", quietly=T)) install.packages("dplyr")
if (!requireNamespace("readr", quietly=T)) install.packages("readr")
if (!requireNamespace("Matrix", quietly = TRUE)) install.packages("Matrix")
if (!requireNamespace("ggplot2", quietly=T)) install.packages("ggplot2")
if (!requireNamespace("stringr", quietly=T)) install.packages("stringr")
if (!requireNamespace("furrr", quietly=T)) install.packages("furrr")
if (!requireNamespace("progressr", quietly=T)) install.packages("progressr")

library(crqa)
library(dplyr)
library(readr)
library(Matrix)
library(ggplot2)
library(stringr)
library(furrr)
library(purrr)
library(progressr)

data.path <- '../data/processed/txt_transformed'
data.output.path <- '../data/processed/crqa_results'
fig.output.path <- file.path(data.output.path, "figs")
if (!dir.exists(data.output.path)) dir.create(data.output.path, recursive=T)
if (!dir.exists(fig.output.path)) dir.create(fig.output.path, recursive=T)

pid     <- c("TD13", "TD17", "TD24", "TD27", "TD30", "TD31", "TD54", "TD55")
months  <- c("M3", "M4")

# Trial numbers to consider
trial_ids <- 2:7

# --- CRQA parameters for categorical data ---
# Following the crqa tutorial
#  - m (embed) = 1, τ (delay) = 1
#  - radius very small to enforce exact-match recurrence on symbols
#  - small Theiler window to reduce trivial adjacency on main diagonal
delay       <- 1
embed       <- 1
rescale     <- 0
radius      <- 0.001
normalize   <- 0
mindiagline <- 6
minvertline <- 2
tw          <- 2
whiteline   <- F
recpt       <- F
side        <- "both"
method      <- "crqa"
metric      <- "euclidean"
datatype    <- "categorical"

# --- Utility functions ---
# Recode symbols so ONLY desired matches can ever be equal ---
# Goal:
#  - Always count 1-1 => give both series the SAME label "ONE" for 1
#  - Optionally count 0-0 => give both series the SAME label "ZERO" for 0 (if include_zero)
#  - Make ALL OTHER CODES non-matching by giving DIFFERENT labels per stream
#    (e.g., presenter's 'other' -> "POTHER", infant's 'other' -> "IOTHER")
 
# NOTE:
#  - Codes like 2 (MO) and 4 (unclear) are now safely mapped to "POTHER" / "IOTHER"
#    so they NEVER match each other (or 1/0). This cleanly excludes 2-2, 4-4, etc.
map_presenter_strict <- function(x) ifelse(x == 1L, "ONE", "POTHER")
map_infant_strict    <- function(x) ifelse(x == 1L, "ONE", "IOTHER")

map_presenter_with_zero <- function(x) {
    case_when(
      x == 1L ~ "ONE",
      x == 0L ~ "ZERO",
      TRUE    ~ "POTHER"
    )
}

map_infant_with_zero <- function(x) {
    case_when(
      x == 1L ~ "ONE",
      x == 0L ~ "ZERO",
      TRUE    ~ "IOTHER"
    )
}

# Select which mapping to use
map_presenter <- map_presenter_strict
map_infant    <- map_infant_strict

# --- HELPERS ---

# Return a tibble with columns: id, month, trial, hand ("LA/RA"), role, file
# We rely on the naming convention, e.g.:
#   TD13-M3_A2R1_CC_RE_R_for_LA.csv (presenter-for-LA)
#   TD13-M3_A2R1_CC_RE_LA.csv       (infant LA)
#   TD13-M3_A2R1_CC_RE_R_for_RA.csv (presenter-for-RA)
#   TD13-M3_A2R1_CC_RE_RA.csv       (infant RA)
index_files <- function(all_files, id, mon, trial) {
    # Build identifiers
    id_mon_trial <- sprintf("^%s-%s_A%d.*\\.csv$", id, mon, trial)
    rel <- all_files[str_detect(all_files, id_mon_trial)]
    if (!length(rel)) return (tibble())

    # Classify role & hand
    tibble(file = rel) %>%
        mutate(
               id    = id,
               mon   = mon,
               trial = trial,
               hand = case_when(
                    str_detect(file, "RE_LA|_LA\\.csv$") ~ "LA",
                    str_detect(file, "RE_RA|_RA\\.csv$") ~ "RA",
                    TRUE ~ NA_character_
                    ),
               role = case_when(
                    str_detect(file, "R_for_LA|R_for_RA") ~ "presenter",
                    TRUE ~"infant"
                    )
               ) %>%
        filter(!is.na(hand))
}

# Shared, fixed level set used by BOTH series
# Include ZERO only if you intend to allow 0-0 matches
CAT_LEVELS <- c("ONE", "POTHER", "IOTHER")

# Categorical to numeric
encode_categorical <- function(x) {
    as.integer(factor(x, levels=CAT_LEVELS))
}

drp_from_rp_idx <- function(RP,
                            maxlag    = 50L,
                            theiler   = 0L,    # set RR=NA for |lag| <= theiler
                            as_percent = TRUE,
                            bin_ms    = 100L) {
  # Coerce to compressed sparse column for fast indexing
  RP <- as(RP, "dgCMatrix")
  n  <- nrow(RP)
  stopifnot(ncol(RP) == n)

  maxlag <- as.integer(min(maxlag, n - 1L))
  lags   <- (-maxlag):maxlag

  # For each lag k, index the k-th diagonal explicitly.
  # If k >= 0: positions (i, j) = (1:(n-k), (1+k):n)
  # If k <  0: let d = -k, positions (i, j) = ((1+d):n, 1:(n-d))
  band_rr <- vapply(lags, function(k) {
    if (abs(k) <= theiler) return(NA_real_)  # apply Theiler window at profile level

    if (k >= 0L) {
      if (n - k <= 0L) return(NA_real_)
      i <- seq.int(1L, n - k)
      j <- seq.int(1L + k, n)
    } else {
      d <- -k
      if (n - d <= 0L) return(NA_real_)
      j <- seq.int(1L, n - d)
      i <- j + d
    }

    # Extract those positions; RP[...] returns numeric 0/1 for dgCMatrix
    v <- RP[cbind(i, j)]
    if (!length(v)) return(NA_real_)
    mean(v != 0)   # proportion of 1s on that diagonal
  }, numeric(1))

  if (as_percent) band_rr <- 100 * band_rr
  df <- data.frame(lag = lags, RR = band_rr)

  df
}

plot_drp <- function(df) {
  p <- ggplot(df, aes(lag, RR)) +
    geom_hline(yintercept = 0, color = "grey80") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_line(na.rm = TRUE, linewidth = 0.9) +
    labs(
      title = "Diagonal Cross-Recurrence Profile (computed from RP)",
      x = sprintf("Lag (bins; 1 bin = %d ms)", bin_ms),
      y = if (as_percent) "Recurrence Rate (%)" else "Recurrence (proportion)"
    ) +
    theme_minimal(base_size = 12)

  # Annotate the peak if any non-NA values exist
  if (any(!is.na(df$RR))) {
    df2  <- df[!is.na(df$RR), ]
    peak <- df2[which.max(df2$RR), ]
    p <- p + annotate(
      "text", x = peak$lag * 0.9, y = peak$RR,
      label = sprintf("Peak: %d bins (~%d ms), RR = %.2f%s",
                      peak$lag, peak$lag * bin_ms, peak$RR, if (as_percent) "%" else ""),
      vjust = -0.8, size = 3.6
    )
  }
  p
}


# Diagonal cross-recurrence profile
compute_dcrp <- function(ts1_num, ts2_num, ws_bins=50,
                         delay=1, embed=1,
                         rescale=0, radius=0.001, normalize=0,
                         method="crqa", metric="euclidean",
                         datatype="categorical", bin_ms=100) {

    drp <- drpfromts(ts1=ts1_num, ts2=ts2_num, windowsize=ws_bins,
                     method=method, datatype=datatype,
                     delay=delay, embed=embed,
                     rescale=rescale, radius=radius, normalize=normalize,
                     metric=metric)

    # drpdfromts returns a list with 'profile', 'maxrec', and 'delay'
    prof <- as.numeric(drp$profile)
    lags <- seq(-ws_bins, ws_bins)
    if (length(prof) == length(lags)) {
        pos_idx <- which(lags > 0)
        neg_idx <- which(lags < 0)

        pos_auc <- if (length(pos_idx)) sum(prof[pos_idx]) else 0
        neg_auc <- if (length(neg_idx)) sum(prof[neg_idx]) else 0

        asymmetry <- if ((pos_auc + neg_auc) > 0) {
            (pos_auc - neg_auc) / (pos_auc + neg_auc)
        } else {
            NA_real_
        }

        # Peak lag (bins and ms)
        peak_idx <- which.max(prof)
        peak_idx_positive <- which.max(prof[51:length(prof)])
        peak_lag_bins <- lags[peak_idx]
        peak_lag_bins_positive <- lags[50+peak_idx_positive]
        peak_lag_ms <- peak_lag_bins * bin_ms
        peak_lag_ms_positive <- peak_lag_bins_positive * bin_ms
        peak_rr <- prof[peak_idx]
        peak_rr_positive <- prof[peak_idx_positive+50]

        list(
             df = tibble(lag_bins = lags,
                         lag_ms = lags * bin_ms,
                         rr = prof),
             peak_lag_bins = peak_lag_bins,
             peak_lag_bins_positive = peak_lag_bins_positive,
             peak_lag_ms = peak_lag_ms,
             peak_lag_ms_positive = peak_lag_ms_positive,
             peak_rr = peak_rr,
             peak_rr_positive = peak_rr_positive,
             pos_auc = pos_auc,
             neg_auc = neg_auc,
             asymmetry = asymmetry
             )
    } else {
        list(
             df = tibble(lag_bins = NA,
                         lag_ms = NA,
                         rr = NA),
             peak_lag_bins = NA,
             peak_lag_bins_positive = NA, 
             peak_lag_ms = NA,
             peak_lag_ms_positive = NA,
             peak_rr = NA,
             peak_rr_positive = NA,
             pos_auc = NA,
             neg_auc = NA,
             asymmetry = NA
             )
    }
}

# One CRQA run for a (presenter_file, infant_file) pair.
run_one_crqa <- function(presenter_file, infant_file, hand, id, mon, trial) {
    pf <- file.path(data.path, presenter_file)
    inf <- file.path(data.path, infant_file)

    if (!file.exists(pf) || !file.exists(inf)) {
        warning(sprintf("Missing file(s) for %s %s A%d %s: %s | %s",
                        id, mon, trial, hand, presenter_file, infant_file))
        return (NULL)
    }
    
    presenter <- read_csv(pf, show_col_types=F)
    infant    <- read_csv(inf, show_col_types=F)

    dat <- presenter %>%
        select(center_ms, presenter = binary_value) %>%
        inner_join(infant %>% select(center_ms, infant = binary_value), by="center_ms") %>%
        arrange(center_ms)

    # Coerce to integer then map to categorical symbols (character)
    ts_presenter <- map_presenter(as.integer(dat$presenter))
    ts_infant    <- map_infant(as.integer(dat$infant))

    ts1_num <- encode_categorical(ts_presenter)
    ts2_num <- encode_categorical(ts_infant)

    # Run CRQA (categorical exact match)
    ans <- crqa(ts1 = ts1_num,
                ts2 = ts2_num,
                delay = delay, embed = embed,
                rescale = rescale, radius = radius, normalize = normalize,
                mindiagline = mindiagline, minvertline = minvertline,
                tw = tw, whiteline = whiteline, recpt = recpt, side = side,
                method = method, metric = metric, datatype = datatype)

# Key metrics
#  - RR: Share of points in the cross-recurrence matrix that are matches
#        under your rule (exact symbol equality).
#        Nearly half of all 100ms bins are "jointly in the same task-relevant state".
#        Long stretches of 1-1 (and 0-0 if included). A high RR is expected
#        when both streams spend a lot of time in the same gross state
#  - DET: Fraction of recurrent points that form diagonal lines.
#  - L: On average, when the two series match, they stay matched for ~L/10 seconds
#       at a time. Longer L typically reflects more persistent coupling.
#  - ENTR (Entropy of diagonal lengths): Shannon entropy of the diagonal length distribution; moderate-high implies a diverse mix of run lengths, not just one stereotyped bout size.
#  - LAM (Laminarity): Fraction of recurrent points that form vertical lines. Very high LAM means a lot of "sticking" in one series while the other continues to yield matches.
    out <- as.data.frame(ans[c("RR", "DET", "L", "maxL", "ENTR", "LAM", "TT", "max_vertlength")])

    # DCRP: RR as a function of lag
    # Choose a lag window that makes sense for your behavior.
    # With 100 ms bins, ws_bins = 50 -> +/- 5s
    ws_bins <- 50
    # lag_RR <- drp_from_rp_idx(ans$RP, plot=F)
    # max_lag <- which.max(lag_RR[lag_RR$lag > 0, 'RR'])
    dcrp <- compute_dcrp(ts1_num, ts2_num, ws_bins=ws_bins,
                         delay=delay, embed=embed,
                         rescale=rescale, radius=radius, normalize=normalize,
                         method=method, metric=metric, datatype=datatype,
                         bin_ms=100)

    # Safer stem (remove .csv robustly even if name contains multiple dots)
    stem <- sub("\\.csv$", "", infant_file)
    dcrp_csv <- file.path(data.output.path, paste0(stem, "_dcrp.csv"))
    write_csv(dcrp$df, dcrp_csv)

    # peak lag and asymmetric
    out$peak_lag_bins <- dcrp$peak_lag_bins
    out$peak_lag_ms   <- dcrp$peak_lag_ms
    out$peak_rr       <- dcrp$peak_rr
    out$pos_auc       <- dcrp$pos_auc
    out$neg_auc       <- dcrp$neg_auc
    out$asymmetry     <- dcrp$asymmetry
    out$peak_lag_bins_positive <- dcrp$peak_lag_bins_positive
    out$peak_lag_ms_positive   <- dcrp$peak_lag_ms_positive
    out$peak_rr_positive       <- dcrp$peak_rr_positive

    out_name <- paste0(stem, "_crqa_results.csv")
    write.csv(out, file = file.path(data.output.path, out_name), row.names=F)
    message(sprintf("[OK] %s %s A%d %s -> %s", id, mon, trial, hand, out_name))
    rp.plot <- plot_rp(ans$RP,
                       xlabel = "Presenter (100 ms bins)",
                       ylabel = "Infant (100ms bins)",
                       title = "Categorical Cross-Recurrence Plot"
                       )
    rp_path <- file.path(fig.output.path, paste0(stem, "_rp.png"))
    ggplot2::ggsave(filename=rp_path, plot=rp.plot, dpi=300, width=5, height=5, units="in")
    # dev.off()

    # RR(lag) plots
    dcrp_plot_path <- file.path(fig.output.path, paste0(stem, "_dcrp.png"))
    ggplot(dcrp$df, aes(x=lag_ms, y=rr)) +
        geom_line(color = "#2c7fb8", linewidth=1) +
        geom_vline(xintercept=0, linetype="dashed", color="grey50") +
        geom_vline(xintercept=dcrp$peak_lag_ms, color="#d95f0e") +
        labs(x="Lag (ms; + infant follows)",
             y="Recurrence rate",
             title=sprintf("DCRP (peak: %d ms, RR=%.3f)", dcrp$peak_lag_ms,
                           dcrp$peak_rr)) +
        theme_minimal(base_size=12)
    ggsave(filename=dcrp_plot_path, dpi=300, width=6, height=4, units="in")
}

# --- MAIN ---
all_files <- list.files(data.path, recursive=F)

# Optional: parallel plan
library(furrr); plan(multisession, workers=parallel::detectCores()-1)

handlers(global = T); handlers("progress")

with_progress({
    p <- progressor(along = pid)

    walk(pid, function(id) {
         p(sprintf("ID %s", id))
         walk(months, function(mon) {
          walk(trial_ids, function(tr) {
           idx <- index_files(all_files, id, mon, tr)
           if (nrow(idx) == 0L) return(invisible())

           # For each hand, find the presenter and infant files
           for (h in c("LA", "RA")) {
               pair <- idx %>% filter(hand == h)
               if (nrow(pair) < 2L) next

               presenter_file <- pair %>% filter(role == "presenter") %>% pull(file) %>% first()
               infant_file    <- pair %>% filter(role == "infant") %>% pull(file) %>% first()

               if (is.na(presenter_file) || is.na(infant_file)) next

               # Run the CRQA for this pair
               try(
                   run_one_crqa(presenter_file, infant_file, h, id, mon, tr),
                   silent=F
                   )
                }
            })
        })
    })
})


message("Batch CRQA finished.")


# # --- 6) Diagonal Recurrence Profile (lag profile) ---
# # profile across +/- maxlag bins; with 10 Hz binning, 50 bins ~ 5s
# 
# # ---------------------------
# # drp_from_rp(): Build a DRP from ans$RP
# # ---------------------------
# # RP       : a binary (0/1) cross-recurrence matrix (dgCMatrix from crqa)
# # maxlag   : compute lags from -maxlag ... +maxlag
# # theiler  : set RR=NA for |lag| <= theiler to suppress trivial near-diagonal matches
# # as_percent: if TRUE, RR is in percent; if FALSE, RR is proportion [0..1]
# # plot     : if TRUE, returns a ggplot object; otherwise returns a data.frame
# # bin_ms   : bin size in milliseconds for labeling (e.g., 100 ms in your data)
# 
# # Compute DRP directly from a (sparse) RP by explicit indexing
# drp_from_rp_idx <- function(RP,
#                             maxlag    = 50L,
#                             theiler   = 0L,    # set RR=NA for |lag| <= theiler
#                             as_percent = TRUE,
#                             plot      = TRUE,
#                             bin_ms    = 100L) {
#   # Coerce to compressed sparse column for fast indexing
#   RP <- as(RP, "dgCMatrix")
#   n  <- nrow(RP)
#   stopifnot(ncol(RP) == n)
# 
#   maxlag <- as.integer(min(maxlag, n - 1L))
#   lags   <- (-maxlag):maxlag
# 
#   # For each lag k, index the k-th diagonal explicitly.
#   # If k >= 0: positions (i, j) = (1:(n-k), (1+k):n)
#   # If k <  0: let d = -k, positions (i, j) = ((1+d):n, 1:(n-d))
#   band_rr <- vapply(lags, function(k) {
#     if (abs(k) <= theiler) return(NA_real_)  # apply Theiler window at profile level
# 
#     if (k >= 0L) {
#       if (n - k <= 0L) return(NA_real_)
#       i <- seq.int(1L, n - k)
#       j <- seq.int(1L + k, n)
#     } else {
#       d <- -k
#       if (n - d <= 0L) return(NA_real_)
#       j <- seq.int(1L, n - d)
#       i <- j + d
#     }
# 
#     # Extract those positions; RP[...] returns numeric 0/1 for dgCMatrix
#     v <- RP[cbind(i, j)]
#     if (!length(v)) return(NA_real_)
#     mean(v != 0)   # proportion of 1s on that diagonal
#   }, numeric(1))
# 
#   if (as_percent) band_rr <- 100 * band_rr
#   df <- data.frame(lag = lags, RR = band_rr)
# 
#   if (!plot) return(df)
# 
#   p <- ggplot(df, aes(lag, RR)) +
#     geom_hline(yintercept = 0, color = "grey80") +
#     geom_vline(xintercept = 0, linetype = "dashed") +
#     geom_line(na.rm = TRUE, linewidth = 0.9) +
#     labs(
#       title = "Diagonal Cross-Recurrence Profile (computed from RP)",
#       x = sprintf("Lag (bins; 1 bin = %d ms)", bin_ms),
#       y = if (as_percent) "Recurrence Rate (%)" else "Recurrence (proportion)"
#     ) +
#     theme_minimal(base_size = 12)
# 
#   # Annotate the peak if any non-NA values exist
#   if (any(!is.na(df$RR))) {
#     df2  <- df[!is.na(df$RR), ]
#     peak <- df2[which.max(df2$RR), ]
#     p <- p + annotate(
#       "text", x = peak$lag * 0.9, y = peak$RR,
#       label = sprintf("Peak: %d bins (~%d ms), RR = %.2f%s",
#                       peak$lag, peak$lag * bin_ms, peak$RR, if (as_percent) "%" else ""),
#       vjust = -0.8, size = 3.6
#     )
#   }
#   p
# }
# 
# 
# lag_RR <- drp_from_rp_idx(ans$RP, plot=F)
# max_lag <- which.max(lag_RR[lag_RR$lag > 0, 'RR'])
# # estimated infant response delay (in seconds)
# expected_delay <- max_lag * 0.1
