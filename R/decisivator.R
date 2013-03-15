isdecisive <- function(taxa, s) {
  .Call("IsDecisive",taxa, s, length(taxa),length(s))
}