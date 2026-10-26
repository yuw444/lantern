# ============================================================================
# Argument checks shared by the exported functions
# ============================================================================
#
# Each helper stops with a message naming the offending argument, e.g.
# "`vcf_path` does not exist: x.vcf". `arg` defaults to the expression the
# caller passed, which is the argument's name when called as
# .check_file(vcf_path).

.check_file <- function(x, arg = deparse(substitute(x))) {
  .check_string(x, arg)
  if (!file.exists(x))
    stop("`", arg, "` does not exist: ", x, call. = FALSE)
  invisible(x)
}

.check_string <- function(x, arg = deparse(substitute(x)), null_ok = FALSE) {
  if (null_ok && is.null(x)) return(invisible(x))
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x))
    stop("`", arg, "` must be a single non-empty character string",
         call. = FALSE)
  invisible(x)
}

.check_flag <- function(x, arg = deparse(substitute(x))) {
  if (!is.logical(x) || length(x) != 1L || is.na(x))
    stop("`", arg, "` must be TRUE or FALSE", call. = FALSE)
  invisible(x)
}

.check_count <- function(x, arg = deparse(substitute(x))) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 1 || x != round(x))
    stop("`", arg, "` must be a single positive whole number", call. = FALSE)
  invisible(x)
}

# A numeric matrix whose non-missing values are whole numbers in the range
# of `values` (e.g. 0:2), if given. range() is one pass with no copy, so this
# stays cheap on cohort-sized genotype matrices; the whole-number check only
# runs for double storage.
.check_matrix <- function(x, arg = deparse(substitute(x)), values = NULL) {
  if (!is.matrix(x) || !is.numeric(x))
    stop("`", arg, "` must be a numeric matrix", call. = FALSE)
  if (!is.null(values) && length(x) && !all(is.na(x))) {
    r <- range(x, na.rm = TRUE)
    whole <- is.integer(x) || all(x == trunc(x), na.rm = TRUE)
    if (r[1] < min(values) || r[2] > max(values) || !whole)
      stop("`", arg, "` must only contain the whole numbers ", min(values),
           " to ", max(values), " (or NA)", call. = FALSE)
  }
  invisible(x)
}

.check_same_dim <- function(x, y, arg_x = deparse(substitute(x)),
                            arg_y = deparse(substitute(y))) {
  if (!identical(dim(x), dim(y)))
    stop("`", arg_x, "` and `", arg_y, "` must have the same dimensions (",
         toString(dim(x)), " vs ", toString(dim(y)), ")", call. = FALSE)
  invisible(x)
}

# A named integer-like vector of distinct codes, e.g. c(AFR = 0L, EUR = 1L).
.check_codes <- function(x, arg = deparse(substitute(x)), min_length = 1L) {
  nm <- names(x)
  if (!is.numeric(x) || anyNA(x) || any(x != round(x)) || is.null(nm) ||
      any(!nzchar(nm)) || anyDuplicated(nm) || anyDuplicated(x))
    stop("`", arg, "` must be a named vector of distinct integer codes, ",
         "e.g. c(AFR = 0L, EUR = 1L)", call. = FALSE)
  if (length(x) < min_length)
    stop("`", arg, "` must contain at least ", min_length, " populations",
         call. = FALSE)
  invisible(x)
}
