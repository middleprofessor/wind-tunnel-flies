## ---------------------------------------------------------
library(geomorph) # ols fit


## ----array-functions, echo=FALSE--------------------------
matrix_2_array <- function(x){
	# x is a matrix of n rows and p*2 columns of shape data
  # in column order x1, y1, x2, y2, ..., xp, yp
	# the output, y, is in the format for plotshapes in the shapes package
	# dim = p, 2, n where p = no. landmarks and n = n_figures
	n <- dim(x)[1]
	p <- dim(x)[2]/2
	y <- array(0, c(p, 2, n))
	for(i in 1:n){
		l <- 0
		for(j in 1:p){
			l <- l + 1
			y[j, 1, i] <- x[i, l]
			l <- l + 1
			y[j, 2, i] <- x[i, l]
			}
		}
	return(y)
	}

array_2_matrix <- function(x){
	# converts array to matrix - see matrix2array
	n <- dim(x)[3]
	p <- dim(x)[1]
	y <- matrix(0, n, p*2)
	for(i in 1:n){
		l <- 0
		for(j in 1:p){
			l <- l + 1
			y[i, l] <- x[j, 1, i]
			l <- l + 1
			y[i, l] <- x[j, 2, i]
			}
		}
	return(y)
}



## ----tpr-functions, echo=FALSE----------------------------
tpr <- function(x, lm1, lm2){
	# input is an array
	n <- dim(x)[3]
	p <- dim(x)[1]
	tpr_mat <- x
	for(i in 1:n){
		xbxa <- (x[lm2, 1, i] - x[lm1, 1, i])
		xcxa <- (x[, 1, i] - x[lm1, 1, i])
		ybya <- (x[lm2, 2, i] - x[lm1, 2, i])
		ycya <- (x[, 2, i] - x[lm1, 2, i])
		tpr_mat[, 1, i] <- (xbxa*xcxa + ybya*ycya)/(xbxa^2 + ybya^2)
		tpr_mat[, 2, i] <- (xbxa*ycya - ybya*xcxa)/(xbxa^2 + ybya^2)
		}
	return(tpr_mat)
	}


## ----ols-functions, echo=FALSE----------------------------
mean_figure_old <- function(y_array){ # y is an array
  ybar <- matrix(0, dim(y_array)[1], dim(y_array)[2])
  for(i in 1:dim(y_array)[1]){
    for(j in 1:dim(y_array)[2]){
      ybar[i, j] <- mean(y_array[i, j,])
    }
  }
  return(ybar)
}

mean_figure <- function(y_array){ # y is an array
  ybar <- matrix(0, dim(y_array)[1], dim(y_array)[2])
  ybar[, 1] <- apply(y_array[,1,], 1, mean)
  ybar[, 2] <- apply(y_array[,2,], 1, mean)
  return(ybar)
}


centroid_size <- function(coords){
  # centroid size of a figure
  x <- coords[,1]
  y <- coords[,2]
  xbar <- mean(x)
  ybar <- mean(y)
  cs <- sqrt(sum((x-xbar)^2 + (y-ybar)^2))
  return(cs)
}

centroid_size_flat <- function(coords){
  # centroid size of a figure formatted x1, y1, x2, y2
  cols <- 1:length(coords)
  even_col <- cols[lapply(cols, "%%", 2) == 0]
  odd_col <- cols[lapply(cols, "%%", 2) != 0]
  x <- coords[even_col]
  y <- coords[odd_col]
  xbar <- mean(x)
  ybar <- mean(y)
  cs <- sqrt(sum((x-xbar)^2 + (y-ybar)^2))
  return(cs)
}

translate <- function(y, t){	#y is the row vector to translate by t
	for(j in 1:dim(y)[1]){
		y[j, ] <- y[j, ]-t
		}
	y
}

translate_ols <- function(y_i){	#y is the row vector to translate by t
  p <- nrow(y_i)
  centroid <- apply(y_i, 2, mean)
	y_i <- y_i - t(matrix(centroid, ncol = p, nrow = 2))
	return(y_i)
}

scale_ols <- function(y_i){	#y is the row vector to translate by t
  p <- nrow(y_i)
  cs <- centroid_size(y_i)
	y_i <- y_i/cs
	return(y_i)
}

rotate_ols <- function(ref, y_i){
  u <- t(ref) %*% y_i
  svd_u <- svd(u)
  v <- svd_u$v
  d <- svd_u$d
  v <- v * t(matrix(sign(d), nrow = 2, ncol = 2))
  h <- v %*% t(u)
  y_i <- y_i %*% h
  return(y_i)
}


ols <- function(y_array, by = "mean"){ #fits the figures in y_array to ref
  p <- nrow(y_array[,,1])
  if(by == "mean"){
    ref <- mean_figure(y_array) |>
      translate_ols() |>
      scale_ols()
  }
  y_ols <- array(0, dim = dim(y_array)) # create array for superimposed data
  n <- dim(y_array)[3]
  for(i in 1:n){ 		# number of figures in y
    y_i <- y_array[, , i]
    y_i <- translate_ols(y_i)
    y_i <- scale_ols(y_i)
    y_i <- rotate_ols(ref, y_i)
    y_ols[, , i] <- y_i
  }
  return(y_ols)
}

procrustes_d <- function(ref, y){
  d <- sqrt(sum((y - ref)^2))
  return(d)
}

gols <- function(y_array){
  # generalized ols
  # y_array is an array
  ybar <- mean_figure(y_array) |>
      translate_ols() |>
      scale_ols()
  tol <- 0.001
  done <- FALSE
  y_ols <- copy(y_array)
  d_old <- 1000
  while(done == FALSE){
    y_ols <- ols(y_ols)
    ystar <- mean_figure(y_ols)
    d <- procrustes_d(ybar, ystar)
    iter_d <- abs(d - d_old)
    iter_d
    if(iter_d < convergence){
      done <- TRUE
    } else{
      d_old <- d
      ybar <- ystar
    }
  }
  return(fit_data)
}



## ----rf-functions, echo=FALSE-----------------------------

median_figure<-function(y){
  ybar <- matrix(0, dim(y)[1], dim(y)[2])
  for(i in 1:dim(y)[1]){
    for(j in 1:dim(y)[2]){
      ybar[i, j] <- median(y[i, j,])
    }
  }
  return(ybar)
}

median_size <- function(y){
  #finds the median size of a figure
  n <- dim(y)[1] # number of landmarks
  p <- n*(n-1)/2 # number of interlandmark distances
  t <- numeric(p)
  l <- 0
  for(i in 1:(n-1)){
    for(j in (i + 1):n){
      l <- l + 1
      t[l] <- euclid(y[i, ], y[j, ])
    }
  }
  return(median(t))
}

twoDangle <- function(a, b){
	#a and b are vectors with start at 0, 0 and end at [1]=x,  [2]=y
	aa <- atan(a[2]/a[1])
	if(a[1]<0 & a[2]>0){
		aa <- pi+aa
	} else if (a[1]<0 & a[2]<0){
		aa <- pi+aa
	} else if (a[1]>0 & a[2]<0){
		aa <- 2*pi+aa
	}

	ab <- atan(b[2]/b[1])
	if(b[1]<0 & b[2]>0){
		ab  <- pi+ ab
	} else if(b[1]<0 & b[2]<0){
		ab  <- pi+ ab
	} else if(b[1]>0 & b[2]<0){
		ab  <- 2*pi+ ab
	}
	
	twoDangle <- ab-aa
	if(twoDangle>=pi){twoDangle <- twoDangle-2*pi}
	twoDangle
	}
	
rfscale <- function(x, y){			#find the resistant fit scale factor
	t <- numeric(dim(y)[1]-1)
	tau <- numeric(dim(y)[1])
	for(j in 1:dim(y)[1]){
		l <- 0
		for(k in 1:dim(y)[1]){
			if(j!=k){
				l <- l+1
				dx <- euclid(x[j, ],  x[k, ])
				dy <- euclid(y[j, ],  y[k, ])
				t[l] <- dx/dy
				}
			}
		tau[j] <- median(t)
		}
	rfscale <- 1/median(tau)
	}
	
rftranslate <- function(x, y){				#find the resistant fit translation parameter
	t <- numeric(dim(y)[1])
	out <- numeric(dim(y)[2])
	for(k in 1:dim(y)[2]){
		for(j in 1:dim(y)[1]){
			t[j] <- y[j, k]-x[j, k]
			}
		out[k] <- median(t)
		}
	rftranslate <- out
	}

rfrotate <- function(x, y){				#find the resistant fit rotation parameter
	t1 <- numeric(dim(y)[1]-1)
	t2 <- numeric(dim(y)[1])
	xvec <- numeric(dim(y)[2])
	yvec <- numeric(dim(y)[2])
	H <- matrix(0, 2, 2)
	for(j in 1:dim(y)[1]){
		l <- 0
		for(k in 1:dim(y)[1]){
			if(k!=j){
				l <- l+1
				for(m in 1:dim(y)[2]){
					xvec[m] <- x[k, m]-x[j, m]
					yvec[m] <- y[k, m]-y[j, m]
					}
				t1[l] <- twoDangle(xvec, yvec)
				}
			}
		t2[j] <- median(t1)
		}
	theta <- median(t2)
	H[1, 1] <-  cos(theta)
	H[1, 2] <-  -sin(theta)
	H[2, 1] <-  sin(theta)
	H[2, 2] <- cos(theta)
	rfrotate <- H
	}


rf <- function(x, y){ 			#fits the figures in y to x
	yt <- array(0, dim = dim(y))
	for(i in 1:dim(y)[3]){ 		# number of figures in y
		yy <- y[, , i]
		yy <- yy / (rfscale(x, yy))
		yy <- yy %*% rfrotate(x, yy)
		#t <- rftranslate(x, yy)
		# yt <- apply(yy, 1, function(x)(x-t)) # this seems to give the transpose of the correct result
		yy <- translate(yy,  rftranslate(x, yy))
		yt[, , i] <- yy
		}
	yt
	}

grf <- function(y){
  # generalized resistant fit
  # y is an array
  convergence <- 0.001
  #create arrays
  fit_data <- array(0, dim = dim(y))
  yorf <- array(0, dim = dim(y))
  x <- matrix(0, dim(y)[1], dim(y)[2])
  xstar <- matrix(0, dim(y)[1], dim(y)[2])
  #the code
  # I need to write an OLS code because the shapes function procGPA rotated everything by
  #45 degrees or so. Since my guppies started at TPR alignment, I don't need the 
  #data <- procGPA(y)$rotated
  fit_data <- y
  x <- median_figure(fit_data)
  x <- x/median_size(x)
  done <- FALSE
  while(done == FALSE){
    fit_data <- rf(x, fit_data)
    xstar <- median_figure(fit_data)
    xstar <- xstar/median_size(xstar)
    d <- abs(xstar - x)
    c <- median(c(d[, 1], d[, 2]))
    if(c <- convergence){
      done <- TRUE
    } else{
      x <- xstar
    }
  }
  return(fit_data)
}



## ----echo=FALSE-------------------------------------------
plot_shapes <- function(z, color_by = NULL){
	#x is an array as specified in the shapes package
	
	#vectorize x
	n <- dim(z)[3]
	p <- dim(z)[1]
	x <- numeric(n*p)
	y <- numeric(n*p)
	r <- 0
	for(i in 1:n){
		for(j in 1:p){
			r <- r+1
			x[r] <- z[j,1,i]
			y[r] <-	z[j,2,i]		
		}
	}
	# plot(x,y)
	plot_data <- data.table(x = x, y = y)
	gg <- ggplot(data = plot_data,
	             aes(x = x, y = y),
	             color = get(color_by)) +
	  geom_point() +
	  coord_fixed()
	return(gg)
}


## ----echo=FALSE-------------------------------------------
plot_2shapes <- function(z, color_col){
	#x is an array as specified in the shapes package
	
	#vectorize x
	n <- dim(z)[3]
	p <- dim(z)[1]
	x <- numeric(n*p)
	y <- numeric(n*p)
	r <- 0
	for(i in 1:n){
		for(j in 1:p){
			r <- r+1
			x[r] <- z[j,1,i]
			y[r] <-	z[j,2,i]		
		}
	}
	# plot(x,y)
	plot_data <- data.table(x = x,
	                        y = y,
	                        group = color_col)

	gg <- ggplot(data = plot_data,
	             aes(x = x,
	                 y = y,
	                 color = group)) +
	  geom_point() +
	  coord_fixed() +
	  theme_pubr() +
	  scale_color_manual(values = pal_okabe_ito_blue)

	return(gg)
}



## ----plot_difference_shapes, echo=FALSE-------------------
plot_difference_shapes <- function(
    ref, # ref figure
    diff, # difference from ref
    mag = 1,
    zoom = 1,
    start_at_ref = TRUE, # otherwise symmetric on either side of ref
    show_legend = TRUE,
    reference_label = "Reference",
    target_label = "Target"
){
  p <- length(ref)
  x_inc <- (1:p)[odd(1:p)]
  y_inc <- x_inc + 1
  x_mean <- ref[x_inc]
  y_mean <- ref[y_inc]
  x_diff <- diff[x_inc]
  y_diff <- diff[y_inc]
  

  
  fig_ref <- data.table(
    figure = reference_label,
    x = x_mean - 0.5 * x_diff * mag,
    y = y_mean - 0.5 * y_diff * mag
  )
  fig_target <- data.table(
    figure = target_label,
    x = x_mean + 0.5 * x_diff * mag,
    y = y_mean + 0.5 * y_diff * mag
  )
  
  # add last lm = first lm
  fig_ref <- rbind(fig_ref,
                   fig_ref[1,])
  fig_target <- rbind(fig_target,
                      fig_target[1,])
  
  # combine
  fig_dots <- rbind(fig_ref, fig_target)
  fig_dots[, figure := factor(figure,
                              levels = c(reference_label, target_label))]
  
  fig_lines <- fig_dots[-c(12,13,26,27)]

  gg <- ggplot(data = fig_dots,
         aes(x = x,
             y = y,
             color = figure)) +
    geom_point() +
    geom_path(data = fig_lines,
              aes(x = x,
                  y = y,
                  color = figure)) +
    coord_fixed() +
    xlab("X") +
    ylab("Y") +
    theme_void() +
    scale_color_manual(values = pal_okabe_ito) +
    NULL
  if(show_legend == TRUE){
      gg <- gg + theme(legend.position="top")
  }
  return(gg)
}


## ----plot_difference_vectors, echo=FALSE------------------
plot_difference_vectors <- function(
    ybar, # ref figure
    diff_vec, # difference from ref
    mag = 1,
    zoom = 1,
    start_at_ref = TRUE # otherwise symmetric on either side of ref
){
  p <- length(ybar)
  x_inc <- (1:p)[odd(1:p)]
  y_inc <- x_inc + 1
  
  if(start_at_ref == TRUE){
    vec_dt <- data.table(
      x = ybar[x_inc], # mean x
      y = ybar[y_inc], # mean y
      dx_1 = ybar[x_inc], # begin arrow x
      dy_1 = ybar[y_inc], # begin arrow y
      dx_2 = ybar[x_inc] + 1 * diff_vec[x_inc] * mag, # end arrow x
      dy_2 = ybar[y_inc] + 1 * diff_vec[y_inc] * mag # end arrow y
    )
  }else{
    vec_dt <- data.table(
      x = ybar[x_inc], # mean x
      y = ybar[y_inc], # mean y
      dx_1 = ybar[x_inc] - 0.5 * diff_vec[x_inc] * mag, # begin arrow x
      dy_1 = ybar[y_inc] - 0.5 * diff_vec[y_inc] * mag, # begin arrow y
      dx_2 = ybar[x_inc] + 0.5 * diff_vec[x_inc] * mag, # end arrow x
      dy_2 = ybar[y_inc] + 0.5 * diff_vec[y_inc] * mag # end arrow y
    )
  }
  gg <- ggplot(data = vec_dt,
         aes(x = x,
             y = y)) +
    geom_point() +
    geom_segment(aes(x = dx_1,
                     xend = dx_2,
                     y = dy_1,
                     yend = dy_2),
                 arrow = arrow(length = unit(0.06, "npc"))) +
    coord_fixed() +
    xlab("X") +
    ylab("Y") +
    theme_minimal() +
    NULL
  return(gg)
}


## ----output-as-R-file-------------------------------------
# highlight and run to put update into R folder
# knitr::purl("superimposition.Rmd")

