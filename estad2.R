alpha <- runif(1000000, 0, 1000)
beta <-  runif(1000000, 0, 1000000000)

alpha_1 <- seq(0.1,1000, by=0.1)

#Funcion que verifica si los valores propios de la matriz son todos positivos
#retorna TRUE si todos los valores propios son positivos, FALSE de lo contrario
verificarAux <- function(alpha, beta)
{
	matrix <- cbind(c(trigamma(alpha),1/beta), c(1/beta,alpha/beta^2));
	temp <- eigen(matrix,only.values = TRUE, symmetric = TRUE);
	ans <- sum(temp$values<=0);

	if(ans == 0)
	{
	return(TRUE);
	}
	else
	{
	return(FALSE);
	}
}


#Se vectoriza la funcion auxiliar
verificar <- Vectorize(verificarAux, SIMPLIFY = TRUE,
          USE.NAMES = TRUE)


main <- function()
{
	#cantidad de matrices con algun valor propio negativo
	temp <- sum(!verificar(alpha, beta))
	if(temp == 0)
	{
	cat("todos los valores corresponden a matrices definidas positivas (usando alpha y beta)\n");
	}
	else
	{
	cat("no todos los valores corresponden a matrices definidas positivas (usando alpha y beta)\n");
	}	

	temp <- trigamma(alpha_1)*alpha_1
	
	plot(log(alpha_1), log(temp), type="o",col="red", 
			xlab="log(alpha)",ylab="log(trigamma(alpha))", ylim=c(min(log(temp)),max(log(temp))));

	temp1 <- sum(temp<=0);


	temp <- sum(!verificar(alpha, beta))
	if(temp == 0)
	{
	cat("todos los valores corresponden a matrices definidas positivas (usando solo alpha) ");
	}
	else
	{
	cat("no todos los valores corresponden a matrices definidas positivas (usando solo alpha)");
	}	
}





avg <- 0
log_avg <- 0
alphas <- 0 
betas <- 0 

punto3 <- function(alpha, beta, n)
{
	x <- replicate(500, rgamma(n,alpha,scale=beta));
	y <- lapply(seq_len(ncol(x)), function(i) x[,i]);

	all_avg <- unlist(lapply(y, mean))
	alphas <<- unlist(lapply(y, punto2))
	betas <<- all_avg/alphas
}


punto2 <-function(muestra)
{
    avg <<-  mean(muestra)
    log_avg <<- mean(log(muestra))

    return(newton.raphson(f, 0.1,1e-5, 1000 ))
}


f <-function(x)
{
	return(log(avg)-log(x)+digamma(x)-log_avg)
}


newton.raphson <- function(f, a, tol = 1e-5, n = 1000) {

    require(numDeriv) # Package for computing f'(x)
    
    x0 <- a # Set start value to supplied lower bound
    k <- n # Initialize for iteration results
     
    for (i in 1:n) {
        dx <- genD(func = f, x = x0)$D[1] # First-order derivative f'(x0)
        x1 <- x0 - (f(x0) / dx) # Calculate next value x1
 		k[i] <- x1
        # Once the difference between x0 and x1 becomes sufficiently small, output the results.
        if (abs(x1 - x0) < tol) {
            root.approx <- tail(k, n=1)
            res <- root.approx
            return(res)
        }
        # If Newton-Raphson has not yet reached convergence set x1 as x0 and continue
        x0 <- x1
    }
    print('Too many iterations in method')
}





























