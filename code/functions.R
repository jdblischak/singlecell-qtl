# Utility function to force R to output a specific number of decimal places.
# Avoids Git thinking results have changed simply because of slight changes
# in insignificant digits.
# http://stackoverflow.com/a/12135122
specify_decimal <- function(x, k) format(round(x, k), nsmall = k)

# Show the upper left corner of a table
lcorner <- function(x, row = 5, col = 5) {
  stopifnot(class(x)[1] %in% c("data.frame", "data.table", "matrix", "tbl_df"),
            row > 0, col > 0)
  x[seq(row), seq(col)]
}

# Show the upper right corner of a table
rcorner <- function(x, row = 5, col = 5) {
  stopifnot(class(x)[1] %in% c("data.frame", "data.table", "matrix", "tbl_df"),
            row > 0, col > 0)
  n <- ncol(x)
  x[seq(row), seq(n - col, n)]
}
