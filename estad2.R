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
	cat("todos los valores corresponden a matrices definidas positivas (mirando valores propios)\n");
	}
	else
	{
	cat("no todos los valores corresponden a matrices definidas positivas (mirando valores propios)\n");
	}	

	temp <- trigamma(alpha_1)*alpha_1
	
	plot(log(alpha_1), log(temp), type="o",col="red", 
			xlab="log(alpha)",ylab="log(trigamma(alpha)*alpha)", ylim=c(min(log(temp)),max(log(temp))));

	temp1 <- sum(temp<=0);


	temp <- sum(!verificar(alpha, beta))
	if(temp == 0)
	{
	cat("todos los valores corresponden a matrices definidas positivas (mirando trigamma*alpha) ");
	}
	else
	{
	cat("no todos los valores corresponden a matrices definidas positivas (mirando trigamma*alpha)");
	}	
}




avg <- 0
log_avg <- 0
alphas <- 0 
betas <- 0 
x <- 0
y <- 0
valores <- 0

#funcion que retorna el porcentaje de veces que el vector está en la elipse
main2 <- function(alpha, beta, n, epsilon)
{
	x <<- replicate(500, rgamma(n,alpha,scale=beta));
	y <<- lapply(seq_len(ncol(x)), function(i) x[,i]);

	all_avg <- unlist(lapply(y, mean))
	alphas <<- unlist(lapply(y, punto2))
	betas <<- all_avg/alphas

	temp <- rbind(alphas, betas)
	aux <- lapply(seq_len(ncol(temp)), function(i) temp[,i]);
	valores <<- n*unlist(lapply(aux, mult, alpha = alpha, beta=beta))
	qchisq <- qchisq(1-epsilon, 2)

	return(list( 'porcentaje' = sum(valores < qchisq)/500, 'error relativo alpha' = sd(alphas)/alpha, 'error relativo beta' = sd(betas)/beta))
}


gm_mean = function(a){prod(a)^(1/length(a))}

area <- function(n)
{

	return(gm_mean(areas(n)))
}

areas <- function(n)
{
	x <<- replicate(500, rgamma(n,1.5,scale=4));
	y <<- lapply(seq_len(ncol(x)), function(i) x[,i]);

	all_avg <- unlist(lapply(y, mean))
	alphas <<- unlist(lapply(y, punto2))
	betas <<- all_avg/alphas


	temp <- (mapply(valores, alphas, betas))
	impar <- seq(1,999, by=2)
	par <- seq(2,1000, by=2)

	eigen1 <- temp[impar]
	eigen2 <- temp[par]

	return(pi*sqrt(qchisq(1-0.1,2)/(n*eigen1))*sqrt(qchisq(1-0.1,2)/(n*eigen2)))
}


valores <- function(alpha,beta)
{
	return(eigen(matrix(cbind(trigamma(alpha),1/beta,1/beta,alpha/beta^2),ncol=2))$values)
}


#aplicamos newton raphson a f con una muestra dada por parametro
punto2 <-function(muestra)
{
    avg <<-  mean(muestra)
    log_avg <<- mean(log(muestra))

    return(newton.raphson(f, 0.1,1e-5, 1000 ))
}


#la funcion a la que le vamos a aplicar newton raphson
f <-function(x)
{
	return(log(avg)-log(x)+digamma(x)-log_avg)
}


#metodo de newton-raphson
newton.raphson <- function(f, a, tol = 1e-5, n = 1000) {

    require(numDeriv) 
    
    x0 <- a
    k <- n 
     
    for (i in 1:n) {
        dx <- genD(func = f, x = x0)$D[1] 
        x1 <- x0 - (f(x0) / dx) 
 		k[i] <- x1
        if (abs(x1 - x0) < tol) {
            root.approx <- tail(k, n=1)
            res <- root.approx
            return(res)
        }
        x0 <- x1
    }
    print('Muchas iteraciones')
}



#funcion auxiliar para ver si el vector x cumple la condicion de estar en la elipse
mult <- function(x, alpha, beta){

	temp <- matrix(cbind(trigamma(x[1]),1/x[2],1/x[2],x[1]/x[2]^2),ncol=2)

	return((x-c(alpha,beta))%*%temp%*%(x-c(alpha,beta)))
}













