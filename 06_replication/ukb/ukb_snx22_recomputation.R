## snx22_row.R (2026-10-20): append the root-checked SNX22 synonymous row (chunk 676) to a NEW copy of the S10 lookup table
e <- new.env()
nm <- load("${ADSP_ROOT}/revision/ukb/gw_coding_spa/rootcheck/UKB_AD_Coding_676.Rdata", envir = e)
x <- get(nm[1], envir = e); cat("objects:", nm, "class:", class(x), "length:", length(x), "\n")
want <- c("Gene name", "Chr", "Category", "#SNV", "cMAC", "STAAR-B", "Burden(1,25)", "Burden(1,1)"); hit <- NULL
grab <- function(r) {
  if (is.matrix(r) || is.data.frame(r)) {
    r <- as.data.frame(r, check.names = FALSE, stringsAsFactors = FALSE)
    if ("Gene name" %in% names(r)) { s <- r[r[["Gene name"]] == "SNX22" & r[["Category"]] == "synonymous", , drop = FALSE]; if (nrow(s)) return(s) }
  } else if (is.list(r)) { for (z in r) { s <- grab(z); if (!is.null(s)) return(s) } }
  NULL
}
hit <- grab(x)
stopifnot(!is.null(hit), nrow(hit) == 1, all(want %in% names(hit))); print(t(hit[, want]))
csv <- "${ADSP_ROOT}/revision/tables/supp_v3/S10_ukb_coding_lookup_rootcheck.csv"; lines <- readLines(csv)
stopifnot(lines[1] == "arm,Gene name,Chr,Category,#SNV,cMAC,STAAR-B,\"Burden(1,25)\",\"Burden(1,1)\"", !any(grepl("SNX22", lines)))
stopifnot(as.integer(hit[["#SNV"]]) == 68L, as.integer(hit[["cMAC"]]) == 944L, abs(as.numeric(hit[["STAAR-B"]]) - 0.0518593585025056) < 1e-12)
new <- paste(c("rootcheck", "SNX22", as.character(hit[["Chr"]]), "synonymous", "68", "944",
               format(as.numeric(hit[["STAAR-B"]]), digits = 15), format(as.numeric(hit[["Burden(1,25)"]]), digits = 15), format(as.numeric(hit[["Burden(1,1)"]]), digits = 15)), collapse = ",")
out <- sub("_rootcheck.csv$", "_rootcheck_snx22.csv", csv); writeLines(c(lines, new), out)
cat("wrote", out, "\n"); cat(tail(readLines(out), 2), sep = "\n")
