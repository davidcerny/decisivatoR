isdecisive <- function(taxa, s, unrooted=T) {
  if(unrooted==F) {
    print("Rooted tree choice...")
    ans <- .Call("IsDecisiveRooted",taxa, s, length(taxa),length(s))
  } else {
    print("Unrooted tree choice...")
    ans <- .Call("IsDecisiveUnrooted",taxa, s, length(taxa),length(s))
  }
  ans
}
