verGetter <- function(){
x <- utils::packageVersion("deldir")
sub("^([^.]*\\.[^.]*)\\.(.*)$", "\\1-\\2", x) 
}
