# general function plot
# dsim = object generated replicating the appropriate 
#         sim() function
# ipar = parameter of interest
#
f_plot1 = function(dsim, ipar, n, xR=c(-2.5,2.5), by=0.5, leg=T,
                   legend=c('true','biased','corrected'), loc='topleft'){
  
  # # test
  # dsim=dsim; ipar='E'; n=20; xR=c(-1.5,1.5); by=0.5
  
  # plot
  nam = unique( attr(dsim, 'dimnames')[[1]] )
  nam = nam[str_detect(nam , paste0('^', ipar))]
  
  # plot range
  yR = range( density(dsim[nam[1],], na.rm=T)$y, 
              density(dsim[nam[2],], na.rm=T)$y )
  
  dens( dsim[nam[1],], lwd=3, xaxt='n',
        xlab="posterior mean", xlim=xR, ylim=yR )
  axis(side=1, at=seq(xR[1], xR[2], by=by))
  mtext( paste0(ipar, ', n=', n), 3, adj=0, cex=2)
  abline( v=mean( dsim[nam[1],], na.rm=T ), lty=2, lwd=3)
  
  if(leg){
    legend(loc, legend=legend, col=1:3, 
           lty=rep(1,3), lwd=rep(3,3), bty='n')  
  }
  
  dens( dsim[nam[2],], lwd=3 , col=2 , add=TRUE )
  abline( v=mean( dsim[nam[2],], na.rm=T ), col=2, lty=2, lwd=3)
  
  if(length(nam)>2){
    dens( dsim[nam[3],], lwd=3 , col=3 , add=TRUE )
    abline( v=mean( dsim[nam[3],], na.rm=T ), col=3, lty=2, lwd=3)
  }
  
}





# general function plot
# dsim = object generated replicating the appropriate 
#         sim() function
# sX = measurement error
# a = alpha parameter for plot
#
f_plot2 = function(dsim, sX=0.1, a=0.5, xR=c(-6,6), yR=c(-6,6), 
                   colors=c('black','red')){
  
  # # test
  # dsim=d; sX=1; a=0.3; xR=c(-7,7); yR=c(-4,4);colors=c('black','red')
  
  # plot range
  if( is.null(xR) ){
    idx = str_detect( names(dsim), paste0('^A'))
    xR = range( dsim[,idx] )
    xR = c( floor( xR[1] ), ceiling( xR[2] ) ) 
  }
  
  if( is.null(yR) ){
    idx = str_detect( names(dsim), paste0('^D'))
    yR = range( c(dsim[,idx]) )
    yR = c( floor( yR[1] ), ceiling( yR[2] ) )
  }
  
  # plots
  ipar = names(dsim)[str_detect( names(dsim), 'true' )]
  
  if( ipar=='D_true' ){
    
    with(dsim, 
         {
           plot(A, D_true, col=col.alpha(colors[1], a), pch=19,
                xlim=xR, ylim=yR )
           points(A, D_obs, col=col.alpha(colors[2], a), pch=19)
           for(i in 1:nrow(dsim)){
             lines( rep(A[i], 2), c(D_true[i], D_obs[i]), 
                    col=col.alpha('red', a), lty=2)
           }
           b = coef( lm(D_obs ~ -1 + A) )
           abline( c(0, b) )
           mtext( paste0('sD: ', sX, ',  bAD: ', round(b, 3)), 
                  3, adj=0, cex=1.5, at=xR[1])
         }
    )
    
  } else if( ipar=='M_true' ){
    
    with(dsim, 
         {
           plot(M_true, D, col=col.alpha(colors[1], a), pch=19,
                xlim=xR, ylim=yR )
           points(M_obs, D, col=col.alpha(colors[2], a), pch=19)
           for(i in 1:nrow(dsim)){
             lines( c(M_true[i], M_obs[i]), rep(D[i], 2), 
                    col=col.alpha('red', a), lty=2)
           }
           b = coef( lm(D ~ -1 + A + M_obs) )['M_obs']
           abline( c(0, b) )
           mtext( paste0('sM: ', sX, ',  bMD: ', round(b, 3)), 
                  3, adj=0, cex=1.5, at=xR[1])
         }
    )
    
  }
  
}



# general function plot
# dsim = object generated replicating the appropriate 
#         sim() function
# sX = measurement error
# a = alpha parameter for plot
#
f_plot3 = function(dsim, sX=0.1, a=0.5, xR=c(-6,6), yR=c(-6,6),
                   colors=c('black','red')){
  
  # # test
  # dsim=d; sX=1; a=0.1; xR=c(-7,7); yR=c(-6,6); colors=c('black','red')
  
  # plot range
  if( is.null(xR) ){
    idx = str_detect( names(dsim), paste0('^A'))
    xR = range( dsim[,idx] )
    xR = c( floor( xR[1] ), ceiling( xR[2] ) ) 
  }
  
  if( is.null(yR) ){
    idx = str_detect( names(dsim), paste0('^D'))
    yR = range( c(dsim[,idx]) )
    yR = c( floor( yR[1] ), ceiling( yR[2] ) )
  }
  
  # plots
  ipar = names(dsim)[str_detect( names(dsim), 'true' )]
  
  if( ipar=='D_true' ){
    
    with(dsim, 
         {
           plot(A, D_true, col=col.alpha(colors[1], a), pch=19,
                xlim=xR, ylim=yR )
           points(A, D_obs, col=col.alpha(colors[2], a), pch=19)
           for(i in 1:nrow(dsim)){
             lines( rep(A[i], 2), c(D_obs[i]-1.96*sX, D_obs[i]+1.96*sX), 
                    col=col.alpha('red', a))
           }
           b = coef( lm(D_obs ~ -1 + A) )
           abline( c(0, b) )
           mtext( paste0('sD: ', sX, ',  bAD: ', round(b, 3)), 
                  3, adj=0, cex=1.5, at=xR[1])
         }
    )
    
  } else if( ipar=='M_true' ){
    
    with(dsim, 
         {
           plot(M_true, D, col=col.alpha(colors[1], a), pch=19,
                xlim=xR, ylim=yR )
           points(M_obs, D, col=col.alpha(colors[2], a), pch=19)
           for(i in 1:nrow(dsim)){
             lines( c(M_obs[i]-1.96*sX, M_obs[i]+1.96*sX), rep(D[i], 2), 
                    col=col.alpha('red', a))
           }
           b = coef( lm(D ~ -1 + A + M_obs) )['M_obs']
           abline( c(0, b) )
           mtext( paste0('sM: ', sX, ',  bMD: ', round(b, 3)), 
                  3, adj=0, cex=1.5, at=xR[1])
         }
    )
    
  }
  
}




# general function plot
# dsim = object generated replicating the appropriate 
#         sim() function
# var = selected variable
# var_range = range of censoring/truncation for var
#
f_plot4 = function(dsim, var='I', var_range){
  
  # # test
  # dsim=d
  # var='all'
  # var_range=out_range[[3]]
  
  with( dsim, 
        {
          plot(dsim$E_full, dsim$I_full, pch=19, col=col.alpha('black',0.3), 
               xlab='E', ylab='I')
          points(dsim$E_cens, dsim$I_cens, pch=1, col='red', cex=1.5)
          b = coef( lm(I_cens ~ E_cens + SES, data=dsim) )
          abline( c(b[1], b[2]) )
          
          if(var=='I'){
            for(j in 1:nrow(d)){
              lines( rep(E_full[j], 2), c(I_full[j], I_cens[j]), 
                     col=col.alpha('red', 0.3), lty=2)
            }
            abline( h=var_range[2], lty=2 )
            
          } else if(var=='E'){
            for(j in 1:nrow(d)){
              lines( c(d$E_full[j], d$E_cens[j]), rep(d$I_full[j], 2), 
                     col=col.alpha('red', 0.3), lty=2)
            }
            abline( v=var_range[1], lty=2 ) 
            
          } else if(var=='all'){
            
            abline( h=var_range[1], lty=2 )
            abline( h=var_range[2], lty=2 )
            abline( v=var_range[1], lty=2 )
            abline( v=var_range[2], lty=2 )
            
            for(j in 1:nrow(dsim)){
              lines( c(E_full[j], E_cens[j]), c(I_full[j], I_cens[j]), 
                     col=col.alpha('red', 0.3), lty=2)
            }
            
          }
          
          mtext( paste0(var, ' observed in (', var_range[1], ', ', var_range[2], ')',
          ',   bEI: ', round(b[2], 3)), 3, adj=0, cex=1.5, at=min(E_full))
        }
        
  )
  
}





# general function plot
# dsim = object generated replicating the appropriate 
#         sim() function
# var = selected variable
# b = breaks for distribution plot
# var_range = range of censoring/truncation for var
#
f_plot5 = function(dsim, var='I', b=30, var_range){
  
  with( dsim, 
        {
          if(var=='I'){
            hist(I_full, col=col.alpha('black', 0.2), breaks=30, 
                 border='white', main='', xlab=var)
            hist(I_cens, col=col.alpha('red', 0.2), breaks=b, 
                 border='white', add=T)
            abline( v=max(I_cens), lty=2 )
            mtext( paste0(var, ' observed in (', var_range[1], ', ', var_range[2], ')'), 
                   3, adj=0, cex=1.5, at=min(I_full)) 
            
          } else if(var=='E'){
            hist(dsim$E_full, col=col.alpha('black', 0.2), breaks=30, 
                 border='white', main='', xlab=var)
            hist(dsim$E_cens, col=col.alpha('red', 0.2), breaks=b, 
                 border='white', add=T)
            abline( v=min(E_cens), lty=2 )
            mtext( paste0(var, ' observed in (', var_range[1], ', ', var_range[2], ')'), 
                   3, adj=0, cex=1.5, at=min(E_full))
          }
          
        }
        
  )
  
}