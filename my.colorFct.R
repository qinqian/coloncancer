# Color function to generate green-red heat maps
my.colorFct <- function(n = 50, low.col = 0.45, high.col=1, saturation = 1) { 
	if (n < 2) stop("n must be greater than 2")
	n1 <- n%/%2
	n2 <- n - n1
	c(hsv(low.col, saturation, seq(1,0,length=n1)), hsv(high.col, saturation, seq(0,1,length=n2))) 
}
