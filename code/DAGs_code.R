# TO DO ####
# - measurement error (residual counfounding), based on 
#   Causal Diagram course section 4
# - truncation and censoring, with solution (see python examples)
# - add plots for missingsness and it's impact depending on assumptions
# - add example for MAR when you do not get the function right
# - example of post-stratification, sensitivity analysis





# preliminar ####
rm(list=ls())

librerias = c('dagitty','mvtnorm','rethinking','multcomp','tidyverse',
              'lme4','rstan','coda','runjags','rjags','cmdstanr',
              'posterior','bayesplot')
sapply(librerias, require, character.only=T)



wd = 'C:/Users/JRiveraEspejo/Desktop/1. Work/#Classes/PhD Antwerp/presentations/DAGs_PP'
setwd(wd)
source( file.path(wd, 'code', 'figures_code.R') )






# (DAG elementals) ####
# It helps to state you assumptions
# IT is not an statistical model
# It does not tell about the specific functional form: Y = f(X)

# There is one mode type of DAG elemental, but it is related to 
# measurement error (also known as residual confounding). 
# Keep in mind that missing values fall in that category: 
# missing values is just a extreme version of measurement error


# The fork
gen_dag = "dag {
  Z -> {X Y};
  X [exposure];
  Y [outcome]
}"
dag_fork = dagitty( gen_dag )
coordinates( dag_fork ) = list( x=c(Z=0, X=0.2, Y=0.2) ,
                                y=c(Z=0, X=-0.2, Y=0.2) )

# The pipe
gen_dag = "dag {
  X -> Z -> Y;
  X [exposure];
  Y [outcome]
}"
dag_pipe = dagitty( gen_dag )
coordinates( dag_pipe ) = list( x=c(Z=0.2, X=0, Y=0.4) ,
                                y=c(Z=0, X=0, Y=0) )

# The collider
gen_dag = "dag {
  {X Y} -> Z;
  X [exposure];
  Y [outcome]
}"
dag_collider = dagitty( gen_dag )
coordinates( dag_collider ) = list( x=c(Z=0.2, X=0, Y=0) ,
                                    y=c(Z=0, X=-0.2, Y=0.2) )

# The descendant
# descendant on fork
gen_dag = "dag {
  D <- Z -> {X Y};
  X [exposure];
  Y [outcome]
}"
dag_desc_fork = dagitty( gen_dag )
coordinates( dag_desc_fork ) = list( x=c(Z=0.2, X=0.4, Y=0.4, D=0) ,
                                     y=c(Z=0, X=-0.2, Y=0.2, D=0) )

# descendant on pipe
gen_dag = "dag {
  X -> Z -> Y;
  Z -> D;
  X [exposure];
  Y [outcome]
}"
dag_desc_pipe = dagitty( gen_dag )
coordinates( dag_desc_pipe ) = list( x=c(Z=0.2, X=0, Y=0.4, D=0.2) ,
                                     y=c(Z=-0.2, X=0, Y=0, D=0.2) )

# descendant on collider
gen_dag = "dag {
  {X Y} -> Z -> D;
  X [exposure];
  Y [outcome]
}"
dag_desc_coll = dagitty( gen_dag )
coordinates( dag_desc_coll ) = list( x=c(Z=0.2, X=0, Y=0, D=0.4) ,
                                     y=c(Z=0, X=-0.2, Y=0.2, D=0) )



# all DAGS you need
par(mfrow=c(2,3))
drawdag( dag_fork)
drawdag( dag_pipe )
drawdag( dag_collider )

drawdag( dag_desc_fork )
drawdag( dag_desc_pipe )
drawdag( dag_desc_coll )
par(mfrow=c(1,1))













# (fork: spurious) ####
# 
# Location: chapter 05 (p. 125)
# 
# also know as:
#   - spurious association
#   - counfounder
#
# Simulation details 1:
# A = median age at marriage
#   A -> M: negative (more A, less M)
#   A -> D: negative (more A, less D)
# M = marriage rate
#   M -> D: null (to emphasize spurious)
# D = divorce rate
# 
# Hypothesis:
# Does M really cause D?
# in a trivial sense it does, but is it a "true" causal relationship?
#
# DAGs
gen_dag = "dag{ 
  A -> {D M};
  M -> D;
}"
dag_plot1 = dagitty( gen_dag )
coordinates(dag_plot1) = list( x=c(A=0,D=1,M=2) , 
                               y=c(A=0,D=1,M=0) )


gen_dag = 'dag{ 
  A -> {D M};
}'
dag_plot2 = dagitty( gen_dag )
coordinates(dag_plot2) = list( x=c(A=0,D=1,M=2) , 
                               y=c(A=0,D=1,M=0) )

par(mfrow=c(1,2))
drawdag( dag_plot1 )
drawdag( dag_plot2 )
par(mfrow=c(1,1))
# which (causal) process model is the correct?



# EXTRAS
# implied conditional Ind.
impliedConditionalIndependencies( dag_plot1 )
impliedConditionalIndependencies( dag_plot2 )


# Markov Equivalent models
MElist = equivalentDAGs(dag_plot1)
drawdag(MElist)


# adjustments sets
adjustmentSets( dag_plot1 , exposure="M" , outcome="D", 
                type='minimal', effect='direct')
adjustmentSets( dag_plot1 , exposure="A" , outcome="D", 
                type='minimal', effect='direct')


adjustmentSets( dag_plot2 , exposure="M" , outcome="D", 
                type='minimal', effect='direct')
adjustmentSets( dag_plot2 , exposure="A" , outcome="D", 
                type='minimal', effect='direct')







# simulation
# n = simulation sample size
# bAM, bAD, bMD = simulated parameters
# rep = to use in replicatation
#
f_sim = function(n=100, bAM=-1, bAD=-1, bMD=0, rep=F){
  
  # # test
  # n=100; bAM=-1; bAD=-1; bMD=0; rep=F
  
  # sim
  A = rnorm( n ) # sim A
  M = rnorm( n , mean=bAM*A ) # sim A -> M
  D = rnorm( n , mean=bAD*A + bMD*M ) # sim {A M} -> D
  d = data.frame(A=A,M=M,D=D)
  
  # return object
  if(!rep){
    # full data
    return(d)
    
  } else{
    # parameters
    b1 = coef( lm(D ~ A + M, data=d) )['M'] # correct relation
    b2 = coef( lm(D ~ M, data=d) )['M'] # spurious relation
    b = c(b1, b2)
    names(b) = c('M','Ms')
    return( b )
    
  }
  
}



# relationships
d = f_sim(n=100, bAM=-1, bAD=-1, bMD=0, rep=F)

# pdf('fork1_panel.pdf')
psych::pairs.panels(d)
# dev.off()
# notice cor(M,D)>0, when it should cor(M,D)=0


# models
summary(lm(D ~ M, data=d)) # spurious relation
summary(lm(D ~ A + M, data=d)) # controlled relation
summary(lm(D ~ A, data=d)) # true relation
summary(lm(M ~ A, data=d))




# sampling variation
# pdf('fork1_samplesize.pdf')
par(mfrow=c(2,1))
dsim = replicate( 1e4, f_sim(n=20, bAM=-1, bAD=-1, bMD=0, rep=T) )
f_plot1(dsim=dsim, ipar='M', n=20, xR=c(-1.5,1.5), by=0.5, 
        leg=T, legend=c('true','biased'))

dsim = replicate( 1e4, f_sim(n=100, bAM=-1, bAD=-1, bMD=0, rep=T) )
f_plot1(dsim=dsim, ipar='M', n=100, xR=c(-1.5,1.5), by=0.5, 
        leg=F)
par(mfrow=c(1,1))
# dev.off()
# M -> D (spuriously), if not controlled by A 
# equally biased with n=100, but more "confident" of M -> D





# what is going on?
set.seed(12345)
d = f_sim(n=1000, bAM=-1, bAD=-1, bMD=0, rep=F)

# pdf('fork1_triptych.pdf', width=14, height=7)
par(mfrow=c(1,3))

Vlim= round( c( min(d$A), max(d$A) ), 2)
coef_mod = coef(lm(D ~ M, data=d))

plot(d[,c('M','D')], pch=19, col=col.alpha('black',0.2), 
     xlim=c( min(d$M), max(d$M) ), ylim=c( min(d$D), max(d$D) ) )
abline(a=coef_mod[1],b=coef_mod[2], col=col.alpha('black',0.3), lwd=2)
abline(h=0, col=col.alpha('red',0.5), lty=2)
mtext( paste0('bM=', round(coef_mod[2],2), ',  A=[', Vlim[1], ',', Vlim[2], ']'), 
       3, adj=0, cex=1.5)

Vlim=c(-1,1)
idx = d$A>Vlim[1] & d$A<Vlim[2] # stratification
coef_mod = coef(lm(D ~ M, data=d[idx,]))

plot(d[idx,c('M','D')], pch=19, col=col.alpha('black',0.2), 
     xlim=c( min(d$M), max(d$M) ), ylim=c( min(d$D), max(d$D) ) )
abline(a=coef_mod[1],b=coef_mod[2], col=col.alpha('black',0.3), lwd=2)
abline(h=0, col=col.alpha('red',0.5), lty=2)
mtext( paste0('bM=', round(coef_mod[2],2), ',  A=[', Vlim[1], ',', Vlim[2], ']'),
       3, adj=0, cex=1.5)

Vlim=c(-0.1,0.1)
idx = d$A>Vlim[1] & d$A<Vlim[2] # stratification
coef_mod = coef(lm(D ~ M, data=d[idx,]))

plot(d[idx,c('M','D')], pch=19, col=col.alpha('black',0.2), 
     xlim=c( min(d$M), max(d$M) ), ylim=c( min(d$D), max(d$D) ) )
abline(a=coef_mod[1],b=coef_mod[2], col=col.alpha('black',0.3), lwd=2)
abline(h=0, col=col.alpha('red',0.5), lty=2)
mtext( paste0('bM=', round(coef_mod[2],2), ',  A=[', Vlim[1], ',', Vlim[2], ']'),
       3, adj=0, cex=1.5)

par(mfrow=c(1,1))
# dev.off()


















# (fork: masked) ####
# 
# Location: chapter 05 (p. 144)
# 
# also known as:
#   - masked relationships
#   - also related to the omitted variable issue
#
# Simulation details 1: 
# M = mass in kg
#   M -> N: positive (more M, more N, trade off lifespan and learning)
#   M -> K: negative (to emphasize masked, but think about humans)
# N = ratio on neocortex over total brain mass
#   N -> K: positive (more N, more K)
# K = Kilo calories per gram of milk
# 
# Hypothesis:
# Does N cause K?, and how about M?
# Larger brains in mammal, need more energetic milk,
# more mass in mammal, need more energetic milk 
# (but humans are in the middle)
#
# DAGs
gen_dag = "dag{ 
  {M N} -> K;
  M -> N 
}"
dag_plot1 = dagitty( gen_dag )
coordinates(dag_plot1) = list( x=c(M=0,K=1,N=2) , 
                               y=c(M=0,K=1,N=0) )
drawdag( dag_plot1 )



# simulation
# n = simulation sample size
# bMN, bNK, bMK = simulated parameters
# rep = to use in replication
#
f_sim = function(n=100, bMN=1, bNK=1, bMK=-1, rep=F){
  
  # # test
  # n=100; bMN=1; bNK=1; bMK=-1; rep=F
  
  # sim
  M = rnorm( n )
  N = rnorm( n , bMN*M )
  K = rnorm( n , bNK*N + bMK*M ) # negative sign in M is to emphasize the masking
  d = data.frame(N=N,M=M,K=K)
  
  # return object
  if(!rep){
    # full data
    return(d)
    
  } else{
    # parameters
    b1 = coef( lm(K ~ N + M, data=d) )[ c('N','M') ] # correct relation
    b2 = c( coef( lm(K ~ N, data=d) )['N'], # biased relation
            coef( lm(K ~ M, data=d) )['M'] ) # biased relation
    b = c(b1, b2)
    names(b) = c('N','M','Nm','Mm')
    return( b )
    
  }
  
}

# relationships
d = f_sim(n=100, bMN=1, bNK=1, bMK=-1, rep=F)

# pdf('fork2_panel.pdf')
psych::pairs.panels(d)
# dev.off()
# notice cor(M,K)~0, when it should cor(M,K)<0 


# models
summary(lm(K ~ N, data=d)) # biased estimate
summary(lm(K ~ M, data=d)) # biased estimate
summary(lm(K ~ N + M, data=d)) # less biased estimates
summary(lm(N ~ M, data=d))



# sampling variation
# pdf('fork2_samplesize.pdf')
par(mfrow=c(2,1))
dsim = replicate( 1e4, f_sim(n=20, bMN=1, bNK=1, bMK=-1, rep=T) )
f_plot1(dsim=dsim, ipar='M', n=20, xR=c(-2.5,1.5), by=0.5, 
        leg=T, legend=c('true','biased'))

dsim = replicate( 1e4, f_sim(n=100, bMN=1, bNK=1, bMK=-1, rep=T) )
f_plot1(dsim=dsim, ipar='M', n=100, xR=c(-2.5,1.5), by=0.5, leg=F)
par(mfrow=c(1,1))
# dev.off()
# M -/> K (masked), if not controlled by N 
# equally biased with n=100, but more "confident" of M -/> K




# what is going on?
set.seed(12345)
d = f_sim(n=1000, bMN=1, bNK=1, bMK=-1, rep=F)

# pdf('fork2_triptych.pdf', width=14, height=7)
par(mfrow=c(1,3))

Vlim= round( c( min(d$M), max(d$M) ), 2)
coef_mod = coef(lm(K ~ N, data=d))

plot(d[,c('N','K')], pch=19, col=col.alpha('black',0.2), 
     xlim=c( min(d$N), max(d$N) ), ylim=c( min(d$K), max(d$K) ) )
abline(a=coef_mod[1],b=coef_mod[2], col=col.alpha('black',0.3), lwd=2)
abline(h=0, col=col.alpha('red',0.5), lty=2)
mtext( paste0('bN=', round(coef_mod[2],2), ',  M=[', Vlim[1], ',', Vlim[2], ']'), 
       3, adj=0, cex=1.5)

Vlim=c(-1,1)
idx = d$M>Vlim[1] & d$M<Vlim[2] # stratification
coef_mod = coef(lm(K ~ N, data=d[idx,]))

plot(d[idx,c('N','K')], pch=19, col=col.alpha('black',0.2), 
     xlim=c( min(d$N), max(d$N) ), ylim=c( min(d$K), max(d$K) ) )
abline(a=coef_mod[1],b=coef_mod[2], col=col.alpha('black',0.3), lwd=2)
abline(h=0, col=col.alpha('red',0.5), lty=2)
mtext( paste0('bN=', round(coef_mod[2],2), ',  M=[', Vlim[1], ',', Vlim[2], ']'),
       3, adj=0, cex=1.5)

Vlim=c(-0.1,0.1)
idx = d$M>Vlim[1] & d$M<Vlim[2] # stratification
coef_mod = coef(lm(K ~ N, data=d[idx,]))

plot(d[idx,c('N','K')], pch=19, col=col.alpha('black',0.2), 
     xlim=c( min(d$N), max(d$N) ), ylim=c( min(d$K), max(d$K) ) )
abline(a=coef_mod[1],b=coef_mod[2], col=col.alpha('black',0.3), lwd=2)
abline(h=0, col=col.alpha('red',0.5), lty=2)
mtext( paste0('bN=', round(coef_mod[2],2), ',  M=[', Vlim[1], ',', Vlim[2], ']'),
       3, adj=0, cex=1.5)

par(mfrow=c(1,1))
# dev.off()








# Simulation details 3:
# U = unobserved variable 
#   U -> N: positive (more U, more N)
#   U -> M: positive (more U, more M)
# N = ratio on neocortex over total brain mass
#   N -> K: positive (more N, more K)
# M = mass in kg
#   M -> K: negative (to emphasize masked, but think about humans)
# K = Kcalories per gram of milk
# 
# Hypothesis:
# Does N cause K?, and how about M?
# Larger brains in mammal, need more energetic milk,
# more mass in mammal, need more energetic milk 
# (but humans are in the middle)
#
# DAGs
gen_dag = "dag{ 
  {M N} -> K;
  U -> {M N};
  U [unobserved]
}"
dag_plot1 = dagitty( gen_dag )
coordinates(dag_plot1) = list( x=c(M=0,K=1,N=2,U=1) , 
                               y=c(M=0,K=1,N=0,U=-1) )
drawdag( dag_plot1 )


# adjustments sets
adjustmentSets( dag_plot1 , exposure="N" , outcome="K", 
                type='minimal', effect='direct')



# simulation
# n = simulation sample size
# bU, bNK, bMK = simulated parameters
# rep = to use in replication
#
f_sim = function(n=100, bU=1, bNK=1, bMK=-1, rep=F){
  
  # # test
  # n=100; bU=1; bNK=1; bMK=-1; rep=F
  
  # sim
  U = rnorm( n )
  N = rnorm( n , bU*U ) # U affects N and M equally
  M = rnorm( n , bU*U )
  K = rnorm( n , bNK*N + bMK*M ) # negative sign in M is to emphasize the masking
  d = data.frame(U=U,N=N,M=M,K=K)
  
  # return object
  if(!rep){
    # full data
    return(d)
    
  } else{
    # parameters
    b1 = coef( lm(K ~ N + M, data=d) )[ c('N','M') ] # correct relation
    b2 = c( coef( lm(K ~ N, data=d) )['N'], # biased relation
            coef( lm(K ~ M, data=d) )['M'] ) # biased relation
    b = c(b1, b2)
    names(b) = c('N','M','Nm','Mm')
    return( b )
    
  }
  
}

# relationships
d = f_sim(n=100, bU=1, bNK=1, bMK=-1, rep=F)

# pdf('fork3_panel.pdf')
psych::pairs.panels(d[,-1])
# dev.off()
# notice cor(N,K)~0.5 and cor(M,K)~-0.5, 
# when it should cor(N,K)>0.5 and cor(M,K)<-0.5


# models
summary(lm(K ~ N, data=d)) # unobserved path still open
summary(lm(K ~ M, data=d)) # unobserved path still open
summary(lm(K ~ N + M, data=d)) # unobserved path close


# sampling variation
# pdf('fork3_samplesize.pdf')
par(mfrow=c(2,2))
dsim = replicate( 1e4, f_sim(n=20, bU=1, bNK=1, bMK=-1, rep=T) )
f_plot1(dsim=dsim, ipar='N', n=20, xR=c(-1,2), by=0.5, 
        leg=T, legend=c('true','biased'))
f_plot1(dsim=dsim, ipar='M', n=20, xR=c(-2,1), by=0.5, leg=F)

dsim = replicate( 1e4, f_sim(n=100, bU=1, bNK=1, bMK=-1, rep=T) )
f_plot1(dsim=dsim, ipar='N', n=100, xR=c(-1,2), by=0.5, leg=F)
f_plot1(dsim=dsim, ipar='M', n=100, xR=c(-2,1), by=0.5, leg=F)
par(mfrow=c(1,1))
# dev.off()
# N -> K (masked), if not controlled by M 
# M -> K (masked), if not controlled by N 
# equally biased with n=100, but more "confident" of {N M} -> D (masked)





# what is going on?
set.seed(12345)
d = f_sim(n=1000, bU=1, bNK=1, bMK=-1, rep=F)

# pdf('fork3_triptych.pdf', width=14, height=7)
par(mfrow=c(1,3))

Vlim= round( c( min(d$M), max(d$M) ), 2)
coef_mod = coef(lm(K ~ N, data=d))

plot(d[,c('N','K')], pch=19, col=col.alpha('black',0.2), 
     xlim=c( min(d$N), max(d$N) ), ylim=c( min(d$K), max(d$K) ) )
abline(a=coef_mod[1],b=coef_mod[2], col=col.alpha('black',0.3), lwd=2)
abline(h=0, col=col.alpha('red',0.5), lty=2)
mtext( paste0('bN=', round(coef_mod[2],2), ',  M=[', Vlim[1], ',', Vlim[2], ']'), 
       3, adj=0, cex=1.5)

Vlim=c(-1,1)
idx = d$M>Vlim[1] & d$M<Vlim[2] # stratification
coef_mod = coef(lm(K ~ N, data=d[idx,]))

plot(d[idx,c('N','K')], pch=19, col=col.alpha('black',0.2), 
     xlim=c( min(d$N), max(d$N) ), ylim=c( min(d$K), max(d$K) ) )
abline(a=coef_mod[1],b=coef_mod[2], col=col.alpha('black',0.3), lwd=2)
abline(h=0, col=col.alpha('red',0.5), lty=2)
mtext( paste0('bN=', round(coef_mod[2],2), ',  M=[', Vlim[1], ',', Vlim[2], ']'),
       3, adj=0, cex=1.5)

Vlim=c(-0.1,0.1)
idx = d$M>Vlim[1] & d$M<Vlim[2] # stratification
coef_mod = coef(lm(K ~ N, data=d[idx,]))

plot(d[idx,c('N','K')], pch=19, col=col.alpha('black',0.2), 
     xlim=c( min(d$N), max(d$N) ), ylim=c( min(d$K), max(d$K) ) )
abline(a=coef_mod[1],b=coef_mod[2], col=col.alpha('black',0.3), lwd=2)
abline(h=0, col=col.alpha('red',0.5), lty=2)
mtext( paste0('bN=', round(coef_mod[2],2), ',  M=[', Vlim[1], ',', Vlim[2], ']'),
       3, adj=0, cex=1.5)

par(mfrow=c(1,1))
# dev.off()








# (fork: multicollinearity) ####
#
# Location: chapter 6 (p. 163)
#
# also know as:
#   - special case of masked relationship
#
# Simulation details: 
# U = unobserved variable (e.g. genetics)
#   U -> LL: positive (more U, more LL)
#   U -> RL: positive (more U, more RL)
# LL = left leg's longitude
#   LL -> H: positive (more LL, more H, trivial sense)
# RL = right leg's longitude
#   RL -> H: positive (more RL, more H, trivial sense)
# X = observed covariate
#   X -> M: positive (more X, more M)
# M = outcome
# 
# Hypothesis:
# does X impact on M?
#
# DAG
gen_dag = "dag {
  {LL RL} -> H;
  U -> {LL RL};
  U [unobserved]
}"
dag_plot1 = dagitty( gen_dag )
coordinates(dag_plot1) = list( x=c(LL=0,H=1,U=1,RL=2) , 
                               y=c(LL=0,H=1,U=-1,RL=0) )
drawdag( dag_plot1 )



# simulation
# n = simulation sample size
# bEA, bEX, bEM, bXM = simulated parameters
# rep = to use in replication
#
f_sim = function(n=100, m_H=170, pL=0.5, re=1, rep=F){
  
  # # test
  # n=100; m_H=170; pL=0.5; re=1; rep=F
  
  # backward simulation
  H = round( rnorm( n , m_H, 2), 1) # sim total height of each
  Lp = runif( n , pL-0.05, pL+0.05) # leg as proportion of height
  LL = round( Lp*H + rnorm( n , 0, re ), 1) # sim left leg as proportion + error
  RL = round( Lp*H + rnorm( n , 0, re ), 1) # sim right leg as proportion + error
  d = data.frame(LL,RL,H) # combine into data frame
  
  # return object
  if(!rep){
    # full data
    return(d)
    
  } else{
    # parameters
    b1 = coef( lm(H ~ -1 + LL, data=d) )['LL'] # more efficient
    b1 = c(b1, coef( lm(H ~ -1 + RL, data=d) )['RL'] ) # more efficient
    b2 = coef( lm(H ~ -1 + LL + RL, data=d) )[c('LL','RL')] # inefficient
    b = c(b1, b2)
    names(b) = c('LL','RL','LLi','RLi')
    return( b )
    
  }
  
}


# relationships
set.seed(12345)
d = f_sim(n=100, m_H=170, pL=0.5, re=1, rep=F)

# pdf('fork4_panel.pdf')
psych::pairs.panels(d)
# dev.off()
# notice cor(LL,H)~0.2, cor(RL,H)~0.2


# models
summary(lm(H ~ -1 + LL + RL, data=d)) # inefficient (SE large) 
summary(lm(H ~ -1 + LL, data=d)) # unbiased, and efficient (SE lower)
summary(lm(H ~ -1 + RL, data=d)) # unbiased, and efficient (SE lower)



# sampling variation
# pdf('fork4_samplesize.pdf')
par(mfrow=c(2,2))
dsim = replicate( 1e4, f_sim(n=20, m_H=170, pL=0.5, re=1, rep=T) )
f_plot1(dsim=dsim, ipar='LL', n=20, xR=c(-0.5,2.2), by=0.1, 
        leg=T, legend=c('true','biased'))
f_plot1(dsim=dsim, ipar='RL', n=20, xR=c(-0.5,2.2), by=0.1, leg=F)

dsim = replicate( 1e4, f_sim(n=100, m_H=170, pL=0.5, re=1, rep=T) )
f_plot1(dsim=dsim, ipar='LL', n=100, xR=c(-0.5,2.2), by=0.1, leg=F)
f_plot1(dsim=dsim, ipar='RL', n=100, xR=c(-0.5,2.2), by=0.1, leg=F)
par(mfrow=c(1,1))
# dev.off()
# {LL RL} -> H, inefficient if both in model 
# {LL RL} -> H, better efficiency with one in model 
# equally biased with n=100, but less "confident" of {LL RL} -> H




# what is going on?
set.seed(12345)
d = f_sim(n=100, m_H=170, pL=0.5, re=1, rep=F)

par_sample = extract.samples(lm(H ~ -1 + LL + RL, data=d))
par_mean = apply(par_sample,2,mean)

# pdf('fork4_triptych.pdf', width=7, height=7)
plot(par_sample, pch=19, col=col.alpha('black',0.1))
points(par_mean[1],par_mean[2], pch=19, col=2, cex=2)
lines( x=rep(par_mean[1], 2), y=c(-3,par_mean[2]), col=col.alpha('red',0.5), lty=2)
lines( x=c(-3, par_mean[1]), y=rep(par_mean[2], 2), col=col.alpha('red',0.5), lty=2)
mtext( paste0('bLL=', round(par_mean[1],2), ',  bRL=', round(par_mean[2],2) ), 
       3, adj=0, cex=1.5)
# dev.off()










# (fork: neutral control) ####
# 
# Location: Cinelli et al, 2021 (p. 4)
# 
# also known as:
#   - precision booster
#
# Simulation details: 
# G = gender
#   G -> SI: positive (assume it captures other unobservables )
# A = hearing age
#   A -> SI: positive (more A, more SI)
# H = invserve logit of entropy
# 
# Hypothesis:
# Should we include G in our model?
#
# DAGs
gen_dag = "dag{ 
  {G A} -> SI;
}"
dag_plot1 = dagitty( gen_dag )
coordinates(dag_plot1) = list( x=c(G=0,SI=1,A=2) , 
                               y=c(G=0,SI=1,A=0) )
drawdag( dag_plot1 )



# simulation
# n = simulation sample size
# bMN, bNK, bMK = simulated parameters
# rep = to use in replication
#
f_sim = function(n=100, bAS=-1, bGS=-1, rep=F){
  
  # # test
  # n=100; bAS=-1; bGS=-1; rep=F
  
  # sim
  G = sample( 0:1, n, replace=T )
  A = rnorm( n )
  H = rnorm( n , bAS*A + bGS*G )
  d = data.frame(G=G,A=A,H=H)
  
  # return object
  if(!rep){
    # full data
    return(d)
    
  } else{
    # parameters
    b1 = coef( lm(H ~ A, data=d) )['A'] # correct relation
    b2 = coef( lm(H ~ A + G, data=d) )['A'] # more precise relation
    s1 = summary( lm(H ~ A, data=d) )$coefficients[,2][2] # se
    s2 = summary( lm(H ~ A + G, data=d) )$coefficients[,2][2]
    b = c(b1, b2, s1, s2)
    names(b) = c('A','AG','sA','sAG')
    return( b )
    
  }
  
}

# relationships
d = f_sim(n=100, bAS=-1, bGS=-1, rep=F)

# pdf('fork5_panel.pdf')
psych::pairs.panels(d)
# dev.off()
# notice cor(M,K)~0, when it should cor(M,K)<0 


# models
summary(lm(H ~ A, data=d)) # correct estimate
summary(lm(H ~ A + G, data=d)) # correct estimate, more precise



# sampling variation
# pdf('fork5_samplesize.pdf')
par(mfrow=c(2,2))
dsim = replicate( 1e4, f_sim(n=20, bAS=-1, bGS=-1, rep=T) )
f_plot1(dsim=dsim, ipar='A', n=20, xR=c(-2.5,0), by=0.5, 
        leg=T, legend=c('A only','A and G'))
f_plot1(dsim=dsim, ipar='sA', n=20, xR=c(0,0.5), by=0.1, leg=F)

dsim = replicate( 1e4, f_sim(n=100, bAS=-1, bGS=-1, rep=T) )
f_plot1(dsim=dsim, ipar='A', n=100, xR=c(-2.5,0), by=0.5, leg=F)
f_plot1(dsim=dsim, ipar='sA', n=100, xR=c(0,0.5), by=0.1, leg=F)
par(mfrow=c(1,1))
# dev.off()
# M -/> K (masked), if not controlled by N 
# equally biased with n=100, but more "confident" of M -/> K













# (pipe: precision parasite) ####
#
# Location: lecture 06, slides, 2022 course, Cinelli et al, 2021 (p.5)
#
# also an example of:
#   - 
#
# Data details: 
# Z = (e.g.)
#   Z -> X: positive 
# X = (e.g.)
#   X -> Y: positive 
# Y = (e.g.)
# 
# Hypothesis:
# does X affect Y?
#
# DAGs
gen_dag = "dag{ 
  Z -> X;
  X -> Y;
}"
dag_plot1 = dagitty( gen_dag )
coordinates(dag_plot1) = list( x=c(Z=0,X=0,Y=2) , 
                               y=c(Z=-1,X=0,Y=0) )
drawdag( dag_plot1 )


# simulation
# n = simulation sample size
# bZX, bXY = parameter of simulation
# rep = to use in replication
#
f_sim = function(n=100, bZX=1, bXY=1, rep=F){
  
  # # test
  # n=100; bZX=1; bXY=1; rep=F
  
  # sim
  Z = rnorm( n ) 
  X = rnorm( n , bZX*Z ) 
  Y = rnorm( n , bXY*X )  
  d = data.frame(Z,X,Y)
  
  # return object
  if(!rep){
    # full data
    return(d)
    
  } else{
    # parameters
    b1 = coef( lm(Y ~ X, data=d) )['X'] # unbiased effects
    b2 = coef( lm(Y ~ X + Z, data=d) )['X'] # biased effect
    b = c(b1, b2)
    names(b) = c('X','Xp')
    return( b )
    
  }
  
}

# relationships
d = f_sim(n=100, bZX=1, bXY=1, rep=F)

# pdf('pipe1_samplesize.pdf')
psych::pairs.panels(d)
# dev.off()
# no problems with the relationship


# models
summary(lm(Y ~ X, data=d)) # unbiased effect, more precision
summary(lm(Y ~ X + Z, data=d)) # unbiased effects, less precision


# sampling variation
# pdf('pipe1_samplesize.pdf')
par(mfrow=c(2,1))
dsim = replicate( 1e4, f_sim(n=20, bZX=1, bXY=1, rep=T) )
f_plot1(dsim=dsim, ipar='X', n=20, xR=c(0,2), by=0.5, 
        leg=T, legend=c('with X','with X and Z'))

dsim = replicate( 1e4, f_sim(n=100, bZX=1, bXY=1, rep=T) )
f_plot1(dsim=dsim, ipar='X', n=100, xR=c(0,2), by=0.5, leg=F)
par(mfrow=c(1,1))
# dev.off()
# X -> Y (correct), but controlling for Z makes you loose efficiency








# Location: Cinelli et al, 2021 (p. 7)
#
# Data details: 
# Z = (e.g.)
#   Z -> X: positive 
# X = (e.g.)
#   X -> Y: positive 
# Y = (e.g.)
# 
# Hypothesis:
# does X affect Y?
#
# DAGs
gen_dag = "dag{ 
  X -> Z;
  X -> Y;
}"
dag_plot1 = dagitty( gen_dag )
coordinates(dag_plot1) = list( x=c(Z=0,X=0,Y=2) , 
                               y=c(Z=-1,X=0,Y=0) )
drawdag( dag_plot1 )


# simulation 2
# n = simulation sample size
# bXZ, bXY = parameter of simulation
# rep = to use in replication
#
f_sim = function(n=100, bXZ=1, bXY=1, rep=F){
  
  # # test
  # n=100; bZX=1; bXY=1; rep=F
  
  # sim
  X = rnorm( n ) 
  Z = rnorm( n , bXZ*X ) 
  Y = rnorm( n , bXY*X )  
  d = data.frame(Z,X,Y)
  
  # return object
  if(!rep){
    # full data
    return(d)
    
  } else{
    # parameters
    b1 = coef( lm(Y ~ X, data=d) )['X'] # unbiased effects
    b2 = coef( lm(Y ~ X + Z, data=d) )['X'] # biased effect
    b = c(b1, b2)
    names(b) = c('X','Xp')
    return( b )
    
  }
  
}

# relationships
d = f_sim(n=100, bXZ=1, bXY=1, rep=F)

# pdf('pipe1_samplesize.pdf')
psych::pairs.panels(d)
# dev.off()
# no problems with the relationship


# models
summary(lm(Y ~ X, data=d)) # unbiased effect, more precision
summary(lm(Y ~ X + Z, data=d)) # unbiased effects, less precision


# sampling variation
par(mfrow=c(2,1))
dsim = replicate( 1e4, f_sim(n=20, bXZ=1, bXY=1, rep=T) )
f_plot1(dsim=dsim, ipar='X', n=20, xR=c(0,2), by=0.5, 
        leg=T, legend=c('true','less precise'))

dsim = replicate( 1e4, f_sim(n=100, bXZ=1, bXY=1, rep=T) )
f_plot1(dsim=dsim, ipar='X', n=100, xR=c(0,2), by=0.5, leg=F)
par(mfrow=c(1,1))
# X -> Y (correct), but controlling for Z makes you loose efficiency











# (pipe: post-treatment) ####
#
# Location: Chapter 06 (p. 170)
#
# also known as:
#   - included variable bias (bigger grouping)
#
# Simulation details: 
# H_0 = initial plant height
#   H_0 -> H_1: positive (in a trivial sense)
# F = presence of fungus
#   F -> H_1: negative (more F, less H_1)
# T = treatment for fungus
#   T -> F: negative (T=1, less F)
# H_1 = final plant height
# 
# Hypothesis:
# Does the treatment work? (higher H_1 for T=1 than T=0)
# Should we control for F to understand the mechanisms? 
#
# DAGs
gen_dag = "dag{
  {H_0 F} -> H_1;
  T -> F
}"
dag_plot1 = dagitty( gen_dag )
coordinates(dag_plot1) = list( x=c(H_0=0,H_1=1,F=2,T=3) , 
                               y=c(H_0=0,H_1=1,F=0.5,T=0) )
drawdag( dag_plot1 )



# simulation
# n = simulation sample size
# bTF, bFH = simulated parameters
# rep = to use in replication
#
f_sim = function(n=100, bTF=-0.4, bFH=-3, rep=F){
  
  # # test
  # n=20; bTF=-0.4; bFH=-3; rep=F
  
  # sim
  h0 = rnorm( n , 10, 2) # simulate initial heights
  Tr = rep( 0:1 , each=n/2 ) # assign treatments
  Fu = rbinom( n , size=1 , prob=0.5 + bTF*Tr ) # simulate fungus 
  # fungus prob=0.5 (if treatment=0), and prob=0.1 (if treatment=1)
  h1 = h0 + rnorm( n , 5 + bFH*Fu) # simulate growth
  # you loose 3cm in the mean height growth if you have fungus
  
  d = data.frame( h0=h0, h1=h1, Tr=factor(Tr), Fu=factor(Fu) )
  
  
  # return object
  if(!rep){
    # full data
    return(d)
    
  } else{
    # parameters
    b1 = coef( lm(h1-h0 ~ Tr, data=d) )['Tr1'] # only treatment
    #
    # m_res = lm(h1-h0 ~ -1 + T, data=d) # contrast, only treatment
    # K = matrix(c(1, -1), 1) # contrast 
    # t = glht(m_res, linfct = K) # unbiased test, library(multcomp)
    # b1 = coef(t)
    
    
    b2 = coef( lm(h1-h0 ~ Fu, data=d) )['Fu1'] # only fungus
    #
    # m_res = lm(h1-h0 ~ -1 + F, data=d) # contrast, only fungus
    # K = matrix(c(1, -1), 1) # contrast 
    # t = glht(m_res, linfct = K) # unbiased test, library(multcomp)
    # b2 = coef(t)
    
    b3 = coef( lm(h1-h0 ~ Tr + Fu, data=d) )[c('Tr1','Fu1')]
    #
    # m_res = lm(h1-h0 ~ -1 + T + F, data=d) # contrast, both
    # K = matrix(c(-1, 1, 0), 1)  # contrast treatment
    # t = glht(m_res, linfct = K) # unbiased test, library(multcomp)
    # b3 = coef( t )
    
    b = c(b1, b2, b3)
    names(b) = c('Tr','Fu','Trp','Fup')
    return( b )
    
  }
  
}


# relationships
d = f_sim(n=100, bTF=-0.4, bFH=-3, rep=F)
d %>%
  mutate(diff=h1-h0) %>%
  group_by(Tr) %>%
  summarise(mean=mean(diff), sd=sd(diff), n=n() ) %>%
  mutate(se=sd/sqrt(n))
d %>%
  mutate(diff=h1-h0) %>%
  group_by(Fu) %>%
  summarise(mean=mean(diff), sd=sd(diff), n=n() ) %>%
  mutate(se=sd/sqrt(n))
d %>%
  mutate(diff=h1-h0) %>%
  group_by(Tr, Fu) %>%
  summarise(mean=mean(diff), sd=sd(diff), n=n() ) %>%
  mutate(se=sd/sqrt(n))
# notice mean(diff | T) shows an effect
# same mean(diff | F) shows an effect
# same mean(diff | F, T) is not so clear



# models
summary(lm(h1-h0 ~ Tr, data=d)) # only treatment
summary(lm(h1-h0 ~ Fu, data=d)) # only fungus
summary(lm(h1-h0 ~ Tr + Fu, data=d)) # only fungus
#
# m_res = lm(h1-h0 ~ -1 + T + F, data=d) # both treatment and fungus
# K = matrix(c(-1, 1, 0), 1)  # contrast
# t = glht(m_res, linfct = K) # biased test
# summary(t)

m_res = glm(Fu ~ -1 + Tr, data=d, family='binomial')
summary(m_res)
exp(coef(m_res)) # odds
inv_logit(coef(m_res)) # probability



# sampling variation
# pdf('pipe2_samplesize.pdf')
par(mfrow=c(2,2))
dsim = replicate( 1e4, f_sim(n=30, bTF=-0.4, bFH=-3, rep=T) )
f_plot1(dsim=dsim, ipar='T', n=20, xR=c(-1.5,3), by=0.5, 
        leg=T, legend=c('Only T or F','both'))
f_plot1(dsim=dsim, ipar='F', n=20, xR=c(-5,-1), by=0.5, leg=F)

dsim = replicate( 1e4, f_sim(n=100, bTF=-0.4, bFH=-3, rep=T) )
f_plot1(dsim=dsim, ipar='T', n=100, xR=c(-1.5,3), by=0.5, leg=F)
f_plot1(dsim=dsim, ipar='F', n=100, xR=c(-5,-1), by=0.5, leg=F)
par(mfrow=c(1,1))
# dev.off()
# N -> K (masked), if not controlled by M 
# M -> K (masked), if not controlled by N 
# equally biased with n=100, but more "confident" of {N M} -> D (masked)











# (pipe:  simpson's paradox) ####
#
# Location: Chapter 05 (p. 144)
#
# also known as:
#   - masked relationships
#   - mediation
#
# N = ratio on neocortex over total brain mass
#   N -> M: positive (more N, more M, trade off lifespan and learning)
#   N -> K: positive (more N, more K)
# M = mass in kg
#   M -> K: negative (to extreme masked, but think about humans)
# K = Kilo calories per gram of milk
# 
# Hypothesis:
# Does N cause K?, and how about M?
# Larger brains in mammal, need more energetic milk,
# more mass in mammal, need more energetic milk 
# (but humans are in the middle)
#
# DAGs
gen_dag = "dag{ 
  {M N} -> K;
  N -> M 
}"
dag_plot1 = dagitty( gen_dag )
coordinates(dag_plot1) = list( x=c(M=0,K=1,N=2) , 
                               y=c(M=0,K=1,N=0) )
drawdag( dag_plot1 )



# simulation
# n = simulation sample size
# bMN, bNK, bMK = simulated parameters
# rep = to use in replication
#
f_sim = function(n=100, bNM=1, bNK=1, bMK=-1, rep=F){
  
  # # test
  # n=100; bNM=1; bNK=1; bMK=-1; rep=F
  
  # sim
  N = rnorm( n )
  M = rnorm( n , bNM*N )
  K = rnorm( n , bNK*N + bMK*M ) # negative sign in M is to emphasize the masking
  d = data.frame(N=N,M=M,K=K)
  
  # return object
  if(!rep){
    # full data
    return(d)
    
  } else{
    # parameters
    b1 = coef( lm(K ~ N + M, data=d) )[ c('N','M') ] # correct relation
    b2 = c( coef( lm(K ~ N, data=d) )['N'], # biased relation
            coef( lm(K ~ M, data=d) )['M'] ) # biased relation
    b = c(b1, b2)
    names(b) = c('N','M','Nm','Mm')
    return( b )
    
  }
  
}

# relationships
d = f_sim(n=100, bNM=1, bNK=1, bMK=-1, rep=F)

# pdf('pipe3_panel.pdf')
psych::pairs.panels(d)
# dev.off()
# notice cor(N,K)~0, when it should cor(N,K)>0 


# models
summary(lm(K ~ N, data=d)) # biased estimate
summary(lm(K ~ M, data=d)) # biased estimate
summary(lm(K ~ N + M, data=d)) # less biased estimate
summary(lm(M ~ N, data=d))



# sampling variation
# pdf('pipe3_samplesize.pdf')
par(mfrow=c(2,1))
dsim = replicate( 1e4, f_sim(n=20, bNM=1, bNK=1, bMK=-1, rep=T) )
f_plot1(dsim=dsim, ipar='N', n=20, xR=c(-1.5,2.5), by=0.5, 
        leg=T, legend=c('true','biased'))

dsim = replicate( 1e4, f_sim(n=100, bNM=1, bNK=1, bMK=-1, rep=T) )
f_plot1(dsim=dsim, ipar='N', n=100, xR=c(-1.5,2.5), by=0.5, 
        leg=F)
par(mfrow=c(1,1))
# dev.off()
# N -/> K (masked), if not controlled by M 
# equally biased with n=100, but more "confident" of N -/> K




# what is going on?
set.seed(12345)
d = f_sim(n=1000, bNM=1, bNK=1, bMK=-1, rep=F)


# pdf('pipe3_triptych.pdf', width=14, height=7)
par(mfrow=c(1,3))

Vlim= round( c( min(d$M), max(d$M) ), 2)
coef_mod = coef(lm(K ~ N, data=d))

plot(d[,c('N','K')], pch=19, col=col.alpha('black',0.2), 
     xlim=c( min(d$N), max(d$N) ), ylim=c( min(d$K), max(d$K) ) )
abline(a=coef_mod[1],b=coef_mod[2], col=col.alpha('black',0.3), lwd=2)
abline(h=0, col=col.alpha('red',0.5), lty=2)
mtext( paste0('bN=', round(coef_mod[2],2), ',  M=[', Vlim[1], ',', Vlim[2], ']'), 
       3, adj=0, cex=1.5)

Vlim=c(-1,1)
idx = d$M>Vlim[1] & d$M<Vlim[2] # stratification
coef_mod = coef(lm(K ~ N, data=d[idx,]))

plot(d[idx,c('N','K')], pch=19, col=col.alpha('black',0.2), 
     xlim=c( min(d$N), max(d$N) ), ylim=c( min(d$K), max(d$K) ) )
abline(a=coef_mod[1],b=coef_mod[2], col=col.alpha('black',0.3), lwd=2)
abline(h=0, col=col.alpha('red',0.5), lty=2)
mtext( paste0('bN=', round(coef_mod[2],2), ',  M=[', Vlim[1], ',', Vlim[2], ']'),
       3, adj=0, cex=1.5)

Vlim=c(-0.1,0.1)
idx = d$M>Vlim[1] & d$M<Vlim[2] # stratification
coef_mod = coef(lm(K ~ N, data=d[idx,]))

plot(d[idx,c('N','K')], pch=19, col=col.alpha('black',0.2), 
     xlim=c( min(d$N), max(d$N) ), ylim=c( min(d$K), max(d$K) ) )
abline(a=coef_mod[1],b=coef_mod[2], col=col.alpha('black',0.3), lwd=2)
abline(h=0, col=col.alpha('red',0.5), lty=2)
mtext( paste0('bN=', round(coef_mod[2],2), ',  M=[', Vlim[1], ',', Vlim[2], ']'),
       3, adj=0, cex=1.5)

par(mfrow=c(1,1))
# dev.off()






# Location: Chapter 11 (p. 340)
#
# also an example of:
#   - masked relationships
#   - mediation
#
# Data details: 
# G = gender
#   G -> A: null (assumed no effect)
#   G -> D: negative (G=female, specific D's)
# D = department to be admitted
#   D -> A: negative (specific D's, less A)
# A = number of admissions
# 
# Hypothesis:
# G has no relationship with A? 
# D can be a confounding variable
#
# DAGs
gen_dag = "dag{
  {G D} -> A;
  G -> D
}"
dag_plot1 = dagitty( gen_dag )
coordinates(dag_plot1) = list( x=c(G=0,D=0,A=1) , 
                               y=c(G=0,D=-1,A=0) )
drawdag( dag_plot1 )



# simulation
# n = simulation sample size
# pGD, ar = parameter of simulation
# rep = to use in replication
#
f_sim = function(n=100, pGD=c(0.3,0.8), ar=c(0.1,0.3,0.1,0.3), rep=F){
  
  # # test
  # n = 100; pGD=c(0.3,0.8); ar=c(0.1,0.3,0.1,0.3)
  
  G = sample( 1:2, size=n, replace=T ) # even gender distribution 
  D = rbinom( n , size=1, ifelse(G==1,pGD[1],pGD[2]) ) + 1 # gender 1 tends to apply to department 1, 2 to 2 
  ar = matrix(ar, nrow=2)
  
  # simulate acceptance (interaction of G and D)
  A = rbinom( n, size=1, ar[D,G] ) 
  d = data.frame( G=factor(G), D=factor(D), A)
  
  # return object
  if(!rep){
    # full data
    return(d)
    
  } else{
    # parameters
    b1 = coef( glm(A ~ -1 + G + D, data=d, family='binomial') )[c('G1','G2')] # biased effect
    b1 = c( b1[2]-b1[1], inv_logit(b1[2])-inv_logit(b1[1]) )
    b2 = coef( glm(A ~ -1 + G, data=d, family='binomial') )[c('G1','G2')] # unbiased effects
    b2 = c( b2[2]-b2[1], inv_logit(b2[2])-inv_logit(b2[1]) )
    b = c(b1, b2)
    names(b) = c('GC','GP','GCb','GPb')
    return( b )
    
  }
  
}

# no discrimination
ar = c(0.1,0.3,0.1,0.3) # NO discrimination
matrix(ar, nrow=2) # [D,G]

d = f_sim(n=1000, pGD=c(0.3,0.8), ar=ar, rep=F)
with(d, table(G, A)/nrow(d) ) # it seems G=1 gets accepted less
with(d, table(D, A)/nrow(d) ) # notice it is because D=1 accepts less
with(d, table(G, D)/nrow(d) ) # and G=1 applies more to D=1


# with discrimination
ar = c(0.05,0.2,0.1,0.3) # discrimination on G=1
matrix(ar, nrow=2) # [D,G]

d = f_sim(n=1000, pGD=c(0.3,0.8), ar=ar, rep=F)
with(d, table(G, A)/nrow(d) ) # it seems G=1 gets accepted less
with(d, table(D, A)/nrow(d) ) # because D=1 accepts less
with(d, table(G, D)/nrow(d) ) # and G=1 applies more to D=1
# overall we can say in both cases we see the same pattern if 
# we only use G -> A (i.e. not controlling for D)




# models
d = f_sim(n=1000, pGD=c(0.3,0.8), ar=c(0.1,0.3,0.1,0.3), rep=F)

m_res = glm(A ~ -1 + G, data=d, family='binomial')
# summary(m_res) 
K = matrix(c(1, -1), 1) # contrast of interest
t = glht(m_res, linfct = K)
summary(t) # only gender
inv_logit(coef(m_res)) # probability
# notice we have calculated the total effects (not wrong)
# just answers a different research question


m_res = glm(A ~ -1 + G + D, data=d, family='binomial')
# summary(m_res) 
K = matrix(c(1, -1, 0), 1) # contrast of interest
t = glht(m_res, linfct = K)
summary(t) # G and D
inv_logit(coef(m_res)) # probability
# now we observe there is no effect of G -> A



# sampling variation
par(mfrow=c(2,2))
dsim = replicate( 1e4, f_sim(n=100, pGD=c(0.3,0.8), ar=c(0.1,0.3,0.1,0.3), rep=T) )
f_plot1(dsim=dsim, ipar='GC', xR=c(-1,1.5), by=0.5)
f_plot1(dsim=dsim, ipar='GP', xR=c(-0.1,0.3), by=0.1)

dsim = replicate( 1e4, f_sim(n=500, pGD=c(0.3,0.8), ar=c(0.1,0.3,0.1,0.3), rep=T) )
f_plot1(dsim=dsim, ipar='GC', xR=c(-1,1.5), by=0.5)
f_plot1(dsim=dsim, ipar='GP', xR=c(-0.1,0.3), by=0.1)
par(mfrow=c(1,1))
# E -> H (negative), but underestimated
# equally biased with n=100, but more "confident" of E -> H (underestimated)








# real world data
data(UCBadmit)
d = UCBadmit
names(d) = c('D','G','A','R','N')
d$D = factor( with(d, ifelse(D=='A',1, 
                     ifelse(D=='B',2,
                            ifelse(D=='C',3,
                                   ifelse(D=='D',4,
                                          ifelse(D=='E',5,6)))) ) ) )
d$G = factor( with(d, ifelse(G=='female',1,2)) )
# str(d)


# models
m_res = glm(cbind(A, N) ~ -1 + G, data=d, family='binomial')
# summary(m_res) 
K = matrix(c(1, -1), 1) # contrast of interest
t = glht(m_res, linfct = K)
summary(t) # only gender, biased effect
inv_logit(coef(m_res)[1]) - inv_logit(coef(m_res)[2]) # probability
# total effect of G -> A


m_res = glm(cbind(A, N) ~ -1 + G + D, data=d, family='binomial')
# summary(m_res)
K = matrix(c(1, -1, 0, 0, 0, 0, 0), 1) # contrast of interest
t = glht(m_res, linfct = K)
summary(t) # gender controlled by department, unbiased effect
inv_logit(coef(m_res)[1]) - inv_logit(coef(m_res)[2]) # probability
# notice we have calculated the direct effect of G -> A | D
# more specifically the shown difference is G -> A | D=='A' (dept. A)
# in order to know the direct effect across all departments
# we need to marginalize over D






# Simulation details 2: 
#
# Location: chapter 6 (p. 180)
#
# G = grandparent's educational level
#   G -> P: positive (more G, more P)
#   G -> C: null (to emphasize the problem)
# P = parent's educational level
#   P -> C: positive (more P, more C)
# C = child's educational achievement
# 
# Hypothesis:
# G and P impact positively on C?
#
# DAG
gen_dag = "dag {
  G -> {P C};
  P -> C
}"
dag_plot1 = dagitty( gen_dag )
coordinates(dag_plot1) = list( x=c(G=0,P=1,C=1) , 
                               y=c(G=0,P=0,C=1) )
drawdag( dag_plot1 )



# simulation
# n = simulation sample size
# bGP, bPC, bGC = simulated parameters
# rep = to use in replication
#
f_sim = function(n=100, bGP=1, bPC=1, bGC=0, rep=F){
  
  # # test
  # n=100; bGP=1; bPC=1; bGC=0; rep=F
  
  # sim
  G = rnorm( n )
  P = rnorm( n , bGP*G )
  C = rnorm( n , bPC*P + bGC*G )
  d = data.frame(P=P,G=G,C=C)
  
  # return object
  if(!rep){
    # full data
    return(d)
    
  } else{
    # parameters
    b1 = coef( lm(C ~ G + P, data=d) )['G'] # unbiased effect
    b2 = coef( lm(C ~ G, data=d) )['G'] # biased effects
    b = c(b1, b2)
    names(b) = c('G','Gs')
    return( b )
    
  }
  
}

# relationships
d = f_sim(n=100, bGP=1, bPC=1, bGC=0, rep=F)

psych::pairs.panels(d)
# notice cor(G,C)>0, when it should be cor(G,C)=0


# models
summary(lm(C ~ G, data=d)) # biased estimate
summary(lm(C ~ G + P, data=d)) # unbiased estimate
summary(lm(C ~ P, data=d)) # unbiased estimate



# sampling variation
par(mfrow=c(2,1))
dsim = replicate( 1e4, f_sim(n=20, bGP=1, bPC=1, bGC=0, rep=T) )
f_plot1(dsim=dsim, ipar='G', xR=c(-1,2), by=0.5)

dsim = replicate( 1e4, f_sim(n=100, bGP=1, bPC=1, bGC=0, rep=T) )
f_plot1(dsim=dsim, ipar='G', xR=c(-1,2), by=0.5)
par(mfrow=c(1,1))
# G -> C, if we do not control for P 
# equally biased with n=100, but more "confident" of G -> C












# (pipe/fork: good controls) ####
# 
# Location: Cinelli et al, 2021 (p. 3)
# 
# also known as:
#   - 
#
# Z = confounder
# X = exposure
# M = mediator
# Y = outcome
# 
# Hypothesis:
# What is the total effect of X on Y?,
#
# DAGs
gen_dag = "dag{ 
  X -> M;
  Z -> {X M};
  M -> Y;

}"
dag_plot1 = dagitty( gen_dag )
coordinates(dag_plot1) = list( x=c(X=0,Z=0.5,M=1,Y=1.5) , 
                               y=c(X=0,Z=-1,M=0,Y=0) )
drawdag( dag_plot1 )

# adjustments sets
adjustmentSets( dag_plot1 , exposure="X" , outcome="Y", 
                type='minimal', effect='total')



# similar case 1
gen_dag = "dag{ 
  X -> M;
  U -> {Z M}
  Z -> X;
  M -> Y;
  U [unobserved]
}"
dag_plot1 = dagitty( gen_dag )
coordinates(dag_plot1) = list( x=c(X=0,Z=0.3,U=0.6,M=1,Y=1.5) , 
                               y=c(X=0,Z=-0.5,U=-1,M=0,Y=0) )
drawdag( dag_plot1 )

# adjustments sets
adjustmentSets( dag_plot1 , exposure="X" , outcome="Y", 
                type='minimal', effect='total')


# similar case 2
gen_dag = "dag{ 
  X -> M;
  U -> {Z X}
  Z -> M;
  M -> Y;
  U [unobserved]
}"
dag_plot1 = dagitty( gen_dag )
coordinates(dag_plot1) = list( x=c(X=0,Z=0.65,U=0.3,M=1,Y=1.5) , 
                               y=c(X=0,Z=-0.5,U=-1,M=0,Y=0) )
drawdag( dag_plot1 )

# adjustments sets
adjustmentSets( dag_plot1 , exposure="X" , outcome="Y", 
                type='minimal', effect='total')





# simulation
# n = simulation sample size
# bZX, bXM, bZM, bMY = simulated parameters
# rep = to use in replication
#
f_sim = function(n=100, bZX=1, bXM=0, bZM=1, bMY=0.5, rep=F){
  
  # # test
  # n=100; bZX=1; bXM=0; bZM=1; bMY=0; rep=F
  # # special case: bXM=0 or bMY=0
  
  # sim
  Z = rnorm( 100 )
  X = rnorm( 100 , 0*Z )
  M = rnorm( 100 , 0.5*X + 1*Z ) 
  Y = rnorm( 100 , 0.5*M )
  d = data.frame(Z=Z,X=X,M=M,Y=Y)
  
  # return object
  if(!rep){
    # full data
    return(d)
    
  } else{
    # parameters
    b1 = coef( lm(Y ~ X + Z, data=d) )['X'] # correct relation
    b2 = coef( lm(Y ~ X , data=d) )['X'] # biased relation
    s1 = summary( lm(Y ~ X + Z, data=d) )$coefficients[,2][2] # se
    s2 = summary( lm(Y ~ X, data=d) )$coefficients[,2][2]
    b = c(b1, b2, s1, s2)
    names(b) = c('X','Xb','sX','sXb')
    return( b )
    
  }
  
}

# relationships
d = f_sim(n=100, bZX=1, bXM=0, bZM=1, bMY=0.5, rep=F)

# pdf('pipefork1_panel.pdf')
psych::pairs.panels(d)
# dev.off()


# models
summary(lm(Y ~ X, data=d)) # biased
summary(lm(Y ~ X + Z, data=d)) # unbiased



# sampling variation
# pdf('pipefork1_samplesize.pdf')
par(mfrow=c(2,1))
dsim = replicate( 1e4, f_sim(n=20, bZX=1, bXM=0, bZM=1, bMY=0.5, rep=T) )
f_plot1(dsim=dsim, ipar='X', n=20, xR=c(-1,1), by=0.1, 
        leg=T, legend=c('Z stratification','no Z stratification'))

dsim = replicate( 1e4, f_sim(n=100, bZX=1, bXM=0, bZM=1, bMY=0.5, rep=T) )
f_plot1(dsim=dsim, ipar='X', n=100, xR=c(-1,1), by=0.1, leg=F)
par(mfrow=c(1,1))
# dev.off()




# special cases

# precision parasite
# pdf('pipefork1_samplesize2.pdf')
par(mfrow=c(2,1))
dsim = replicate( 1e4, f_sim(n=20, bZX=1, bXM=0.5, bZM=1, bMY=0, rep=T) )
f_plot1(dsim=dsim, ipar='X', n=20, xR=c(-1,1), by=0.1, 
        leg=T, legend=c('Z stratification','no Z stratification'))

dsim = replicate( 1e4, f_sim(n=100, bZX=1, bXM=0.5, bZM=1, bMY=0, rep=T) )
f_plot1(dsim=dsim, ipar='X', n=100, xR=c(-1,1), by=0.1, leg=F)
par(mfrow=c(1,1))
# dev.off()



# precision booster
# pdf('pipefork1_samplesize3.pdf')
par(mfrow=c(2,2))
dsim = replicate( 1e4, f_sim(n=20, bZX=0, bXM=0.5, bZM=1, bMY=0.5, rep=T) )
f_plot1(dsim=dsim, ipar='X', n=20, xR=c(-1,1), by=0.1, 
        leg=T, legend=c('Z stratification','no Z stratification'))
f_plot1(dsim=dsim, ipar='sX', n=20, xR=c(0.1,0.5), by=0.1, leg=F)

dsim = replicate( 1e4, f_sim(n=100, bZX=0, bXM=0.5, bZM=1, bMY=0.5, rep=T) )
f_plot1(dsim=dsim, ipar='X', n=100, xR=c(-1,1), by=0.1, leg=F)
f_plot1(dsim=dsim, ipar='sX', n=100, xR=c(0.1,0.5), by=0.1, leg=F)
par(mfrow=c(1,1))
# dev.off()







# (pipe/fork: bias amplification) ####
# 
# Location: chapter 14 (p. 455), Cinelli et al, 2021 (p. 5)
# 
# also known as:
#   - related to instrumental variables
#   - also related to the omitted variable issue
#
# Simulation details 1:
# U = unobserved variable (e.g. ability)
#   U -> E: positive (more U, more E)
#   U -> W: positive (more U, more W)
# E = education level
#   E -> W: positive (more E, more W)
# Q = instrumental variable (e.g.)
#   M -> N: positive 
# W = (future) wages
# 
# Hypothesis:
# Does E cause W?,
# Education affects future wages, but by how much?
#
# DAGs
gen_dag = "dag{ 
  E -> W;
  Q -> E;
  U -> {E W};
  U [unobserved]
}"
dag_plot1 = dagitty( gen_dag )
coordinates(dag_plot1) = list( x=c(Q=0,E=0,U=0.5,W=1) , 
                               y=c(Q=-1,E=0,U=-1,W=0) )
drawdag( dag_plot1 )


# adjustments sets
adjustmentSets( dag_plot1 , exposure="E" , outcome="W", 
                type='minimal', effect='direct')




# simulation
# n = simulation sample size
# bU, bQE, bEW = simulated parameters
# rep = to use in replication
#
f_sim = function(n=100, bU=1, bQE=1, bEW=0, rep=F){
  
  # # test
  # n=100; bU=1; bQE=1; bEW=0; rep=F
  
  # sim
  U = rnorm( n )
  Q = sample( 1:4, n, replace=T )
  E = rnorm( n , bQE*Q + bU*U ) # U affects E and W equally
  W = rnorm( n , bEW*E + bU*U )
  d = data.frame(U=U,Q=Q,E=E,W=W)
  
  # return object
  if(!rep){
    # full data
    return(d)
    
  } else{
    # parameters
    b1 = coef( lm(W ~ E + U, data=d) )['E'] # correct relation
    b2 = c( coef( lm(W ~ E + Q, data=d) )['E'], # more biased relation
            coef( lm(W ~ E, data=d) )['E'] ) # biased relation
    b = c(b1, b2)
    names(b) = c('E','Emb','Eb')
    return( b )
    
  }
  
}

# relationships
d = f_sim(n=100, bU=1, bQE=1, bEW=0, rep=F)

# pdf('pipefork2_panel.pdf')
psych::pairs.panels(d[,-1])
# dev.off()
# notice cor(M,K)>0, when it should be cor(M,K)~0 (need control)


# models
summary(lm(W ~ E + U, data=d)) # unbiased (not possible)
summary(lm(W ~ E, data=d)) # biased
summary(lm(W ~ E + Q, data=d)) # more biased


# sampling variation
# pdf('pipefork2a_samplesize.pdf')
par(mfrow=c(2,1))
dsim = replicate( 1e4, f_sim(n=20, bU=1, bQE=1, bEW=0, rep=T) )
f_plot1(dsim=dsim, ipar='E', n=20, xR=c(-1,1.5), by=0.1, 
        leg=T, legend=c('true','more biased','biased'))

dsim = replicate( 1e4, f_sim(n=100, bU=1, bQE=1, bEW=0, rep=T) )
f_plot1(dsim=dsim, ipar='E', n=100, xR=c(-1,1.5), by=0.1, 
        leg=F)
par(mfrow=c(1,1))
# dev.off()
# N -> K (biased), if not controlled by M 
# N -> K (more biased), if controlled by M 
# equally biased with n=100, but more "confident" of N -> K (overestimated)




# how to solve it??

# bayesian way
m = ulam(
  alist(
    c(W,E) ~ multi_normal( c(muW,muE), R , S ),
    muW <- aW + bEW*E,
    muE <- aE + bQE*Q,
    c(aW,aE) ~ normal( 0 , 0.2 ),
    c(bEW,bQE) ~ normal( 0 , 0.5 ),
    R ~ lkj_corr( 2 ),
    S ~ exponential( 1 ) ), 
  data=d , chains=4 , cores=4 )

precis( m , depth=3 )


# frequentist way
s1 = lm( E ~ Q, data=d)
Ehat = s1$fitted.values
s2 = lm( W ~ Ehat, data=d)
# se not corrected

require(AER)
tsls = ivreg( W ~ E | Q, data=d)
# se corrected


summary(s2)
summary(tsls)












# Location: lecture 06, slides, 2022 course
#
# also an example of:
#   - contextual confounds
#
# Simulation details 2: 
# W = win the lottery
#   W -> H: positive (W=1, more H)
# H = happiness
#   H -> L: null (to emphasize problem)
# U = contextual confound
#   U -> H: positive
#   U -> L: positive
# L = lifespan
# 
# Hypothesis:
# does wining the lottery (W) affects you lifespan (L)?
#
# DAGs
gen_dag = "dag{ 
  W -> H;
  H -> L;
  U -> {H L}
}"
dag_plot1 = dagitty( gen_dag )
coordinates(dag_plot1) = list( x=c(W=0,H=1,U=1.5,L=2) , 
                               y=c(W=0,H=0,U=-1,L=0) )
drawdag( dag_plot1 )


# simulation
# n = simulation sample size
# bWH, bWH, bHL, bU = parameter of simulation
# rep = to use in replication
#
f_sim = function(n=100, bWH=1, bHL=0, bU=1, rep=F){
  
  # # test
  # n=100; bWH=1; bHL=0; bU=1; rep=F
  
  # sim
  W = rbinom(n, size=1, p=0.3) 
  U = rnorm(n) 
  H = rnorm(n, bWH*W + bU*U) 
  L = rpois(n, lambda=exp(bHL*H + bU*U) ) 
  d = data.frame(U,W,H,L)
  
  # return object
  if(!rep){
    # full data
    return(d)
    
  } else{
    # parameters
    b1 = coef( glm(L ~ W, data=d, family='poisson') )['W'] # total (unbiased) effects
    b2 = coef( glm(L ~ W + H, data=d, family='poisson') )['W'] # direct (biased) effect
    b = c(b1, b2)
    names(b) = c('W','Wb')
    return( b )
    
  }
  
}

# relationships
d = f_sim(n=100, bU=1, bWH=1, bHL=0, rep=F)
d %>%
  group_by(W) %>%
  summarise(mean=mean(L), sd=sd(L), n=n() ) %>%
  mutate(se=sd/sqrt(n))
# notice equal mean(L) for either group in W
# backdoor already closed


# models
summary(glm(L ~ W, data=d, family='poisson')) # total (unbiased) effects 
summary(glm(L ~ W + H, data=d, family='poisson')) # direct (biased) effects


# sampling variation
# pdf('pipefork2b_samplesize.pdf')
par(mfrow=c(2,1))
dsim = replicate( 1e4, f_sim(n=20, bWH=1, bHL=0, bU=1, rep=T) )
f_plot1(dsim=dsim, ipar='W', n=20, xR=c(-2.5,2), by=0.5,
        leg=T, legend=c('total','direct'))

dsim = replicate( 1e4, f_sim(n=100, bWH=1, bHL=0, bU=1, rep=T) )
f_plot1(dsim=dsim, ipar='W', n=100, xR=c(-2.5,2), by=0.5, leg=F)
par(mfrow=c(1,1))
# dev.off()
# W -\> L, in "total" wining the lottery does not reduce your lifespan
# but the direct effect is biased









# (collider: Berkson's paradox) ####
#
# Location: chapter 6 (p. 161)
#
# also known as:
#   - selection-distortion effect
#   - selection bias
#   - "convenience sample" bias
#
# Simulation details: 
# NW = research's news worthiness
#   NW -> S: positive (more NW, more prob. S)
# TW = research's trust worthiness
#   TW -> S: positive (more TW, more prob. S)
# S = binary research selection
# 
# Hypothesis:
# TW and NW has no relationship after we select by S?
#
# DAGs
gen_dag = "dag{ 
  {NW TW} -> S
}"
dag_plot1 = dagitty( gen_dag )
coordinates(dag_plot1) = list( x=c(NW=0,S=1,TW=2) , 
                               y=c(NW=0,S=1,TW=0) )
drawdag( dag_plot1 )



# simulation
# n = simulation sample size
# p = proportion of selection
# rep = to use in replication
#
f_sim = function(n=100, p=0.1, rep=F){
  
  # # test
  # n=100; p=0.1; rep=F
  
  # sim
  NW = rnorm( n ) # uncorrelated
  TW = rnorm( n )
  Sc = NW + TW  # total score
  q = quantile( Sc , 1-p ) # top 10% threshold
  S = ifelse( Sc >= q , 1 , 0 ) # select top 10%
  d = data.frame( NW=NW, TW=TW, S=S )
  
  # return object
  if(!rep){
    # full data
    return(d)
    
  } else{
    # parameters
    b1 = coef( lm(NW ~ TW, data=d) )['TW'] # unbiased effects
    b2 = coef( lm(NW ~ TW, data=d[d$S==1,]) )['TW'] # biased effect
    b = c(b1, b2)
    names(b) = c('TW','TWs')
    return( b )
    
  }
  
}

# relationships
d = f_sim(n=200, p=0.1, rep=F)

# stratified relationship
# jpeg('collider1_panel2.jpg')
psych::pairs.panels(d[d$S==1,-3])
# dev.off()


# (ideally) un-stratified relationship
# jpeg('collider1_panel1.jpg')
psych::pairs.panels(d[,-3])
# dev.off()


# models
summary(lm(NW ~ TW, data=d[d$S==1,])) # biased effects
summary(lm(NW ~ TW, data=d)) # unbiased effects (ideal)


# sampling variation
# pdf('collider1_samplesize.pdf')
par(mfrow=c(2,1))
dsim = replicate( 1e4, f_sim(n=50, p=0.1, rep=T) )
f_plot1(dsim=dsim, ipar='TW', n=50, xR=c(-2,1), by=0.5, 
        leg=T, legend=c('true','biased'))

dsim = replicate( 1e4, f_sim(n=200, p=0.1, rep=T) )
f_plot1(dsim=dsim, ipar='TW', n=200, xR=c(-2,1), by=0.5, leg=F)
par(mfrow=c(1,1))
# dev.off()
# NW -> TW, if we regress on group S==1 
# equally biased with n=100, but more "confident" of NW -> TW



# what is going on?
# pdf('collider1_triptych.pdf')
with(d, 
     { plot(NW, TW, col=col.alpha('black',0.2), pch=19,
            xlab='news worthiness', ylab='trust worthiness')
       points(NW[S==1], TW[S==1] , 
              col=col.alpha('blue',0.3), pch=19,)
       abline( lm(TW[S==1] ~ NW[S==1]) , 
               col=col.alpha('blue',0.3))
     } )
# dev.off()








# what can I do?
# DAGs
gen_dag = "dag{ 
  {NW Z} -> S
  TW -> Z
}"
dag_plot1 = dagitty( gen_dag )
coordinates(dag_plot1) = list( x=c(NW=0,S=1,TW=2,Z=1.5) , 
                               y=c(NW=0,S=1,TW=0,Z=0.5) )
drawdag( dag_plot1 )



# simulation
# n = simulation sample size
# p = proportion of selection
# rep = to use in replication
#
f_sim = function(n=100, bTZ=1, p=0.1, rep=F){
  
  # # test
  # n=100; p=0.1; rep=F
  
  # sim
  NW = rnorm( n ) # uncorrelated
  TW = rnorm( n )
  Z = rnorm( n, bTZ*TW )
  Sc = NW + Z  # total score
  q = quantile( Sc , 1-p ) # top 10% threshold
  S = ifelse( Sc >= q , 1 , 0 ) # select top 10%
  d = data.frame( NW=NW, TW=TW, Z=Z, S=S )
  
  # return object
  if(!rep){
    # full data
    return(d)
    
  } else{
    # parameters
    b1 = coef( lm(NW ~ TW, data=d) )['TW'] # unbiased effects
    b2 = coef( lm(NW ~ TW, data=d[d$S==1,]) )['TW'] # biased effect
    b3 = coef( lm(NW ~ TW + Z, data=d[d$S==1,]) )['TW'] # corrected effect
    b = c(b1, b2, b3)
    names(b) = c('TW','TWs','TWc')
    return( b )
    
  }
  
}

# relationships
d = f_sim(n=200, bTZ=1, p=0.1, rep=F)


# models
summary(lm(NW ~ TW, data=d[d$S==1,])) # biased effects
summary(lm(NW ~ TW + Z, data=d[d$S==1,])) # corrected effects
summary(lm(NW ~ TW, data=d)) # unbiased effects (ideal)


# sampling variation
# pdf('collider1_samplesize2.pdf')
par(mfrow=c(2,1))
dsim = replicate( 1e4, f_sim(n=50, bTZ=1, p=0.1, rep=T) )
f_plot1(dsim=dsim, ipar='TW', n=50, xR=c(-2,1), by=0.5, 
        leg=T, legend=c('true','biased','corrected'))

dsim = replicate( 1e4, f_sim(n=200, bTZ=1, p=0.1, rep=T) )
f_plot1(dsim=dsim, ipar='TW', n=200, xR=c(-2,1), by=0.5, leg=F)
par(mfrow=c(1,1))
# dev.off()
# NW -> TW, if we regress on group S==1 
# equally biased with n=100, but more "confident" of NW -> TW


















# (collider: M-bias) ####
#
# Location: lecture 06, slides, 2022 course, Cinelli et al, 2021 (p. 4)
#
# also known as:
#   - Collider bias: pre-treatment
#
# Ux = hobbies person 1
#   Ux -> X: positive (better Ux, better X)
#   Ux -> Z: positive (better Ux, better Z)
# Uy = hobbies person 2
#   Uy -> Y: positive (better Uy, better X)
#   Uy -> Z: positive (better Uy, better Z)
# Z = type of friends (defined as continuum)
#     more Z, better type of friends
# X = health of person 1
#   X -> Y: null (to emphasize error)
# Y = health of person 2
# 
# Hypothesis:
# does X impact Y of another?
#
# DAGs
gen_dag = "dag{ 
  Ux -> {Z X};
  Uy -> {Z Y};
  X -> Y
}"
dag_plot1 = dagitty( gen_dag )
coordinates(dag_plot1) = list( x=c(Ux=0,X=0,Z=1,Y=2,Uy=2) , 
                               y=c(Ux=-1,X=1,Z=0,Y=1,Uy=-1) )
drawdag( dag_plot1 )



# simulation
# n = simulation sample size
# bU, bUX, bUY, bXY = parameter of simulation
# rep = to use in replication
#
f_sim = function(n=100, bUZ=0.5, bU=1, bXY=0, bZ=0, rep=F){
  
  # # test
  # n=100; bUZ=1; bU=1; bXY=0; bZ=0; rep=F
  
  # sim
  U1 = sample(1:5, n , replace=T)
  U2 = sample(1:5, n , replace=T)
  Z = rnorm( n , bUZ*U1 + bUZ*U2)
  X = rnorm( n , bU*U1 + bZ*Z)
  Y = rnorm( n , bU*U2 + bXY*X + bZ*Z)
  d = data.frame(U1,U2,Z,X,Y)
  
  # return object
  if(!rep){
    # full data
    return(d)
    
  } else{
    # parameters
    b1 = coef( lm(Y ~ X, data=d) )['X'] # unbiased effects
    b2 = coef( lm(Y ~ X + Z, data=d) )['X'] # biased effect
    b = c(b1, b2)
    names(b) = c('X','Xp')
    return( b )
    
  }
  
}

# relationships
d = f_sim(n=100, bUZ=0.5, bU=1, bXY=0, bZ=0, rep=F)

# pdf('collider4_panel.pdf')
psych::pairs.panels(d[,3:5])
# dev.off()
# cor(X,Y)~0, when not controlled by Z (already good)


# models
summary(lm(Y ~ X, data=d)) # unbiased effects (efficient enough to reject)
summary(lm(Y ~ X + Z, data=d)) # biased effects (efficient enough to reject)


# sampling variation
# pdf('collider4_samplesize.pdf')
par(mfrow=c(2,1))
dsim = replicate( 1e4, f_sim(n=20, bUZ=0.5, bU=1, bXY=0, bZ=0, rep=T) )
f_plot1(dsim=dsim, ipar='X', n=20, xR=c(-0.7,0.5), by=0.2, 
        leg=T, legend=c('true','biased'))

dsim = replicate( 1e4, f_sim(n=100, bUZ=0.5, bU=1, bXY=0, bZ=0, rep=T) )
f_plot1(dsim=dsim, ipar='X', n=100, xR=c(-0.7,0.5), by=0.2, leg=F)
par(mfrow=c(1,1))
# dev.off()
# X -/> Y, when n=20 (still some bias)
# X -> Y, when n=100 (same bias), 
# but now is more "confident" of X -> Y











# (descendant: proxies) ####
#
# also know as:
#   - a solution for unobserved confounding
#
# Simulation details 1: 
#
# Location: 
#
# A = proxy variable (e.g. age)
# E = unobserved variable (e.g. instruction type)
#   E -> A: negative (specific E, less A)
#   E -> X: negative (specific E, less prob. X)
#   E -> M: positive (specific E, more M)
# X = observed covariate (e.g. teaching experience)
#   X -> M: null (to emphasize problem)
# M = outcome (e.g. attainment in mathematics)
# 
# Hypothesis:
# does X impact on M?
#
# DAG
gen_dag = "dag {
  E -> {X M A};
  X -> M;
  E [unobserved]
}"
dag_plot1 = dagitty( gen_dag )
coordinates(dag_plot1) = list( x=c(X=0,E=1,A=1,M=2) , 
                               y=c(X=0,E=-1,A=-2,M=0) )
drawdag( dag_plot1 )



# simulation
# n = simulation sample size
# bEA, bEX, bEM, bXM = simulated parameters
# rep = to use in replication
#
f_sim = function(n=100, bEA=-8, bEX=-0.2, bEM=2, bXM=0, rep=F){
  
  # # test
  # n=100; bEA=-8; bEX=-0.2; bEM=2; bXM=0; rep=F
  
  # sim
  E = sample( 1:2, n, replace=T) # two types of education:
  A = round( rnorm( n , 35 + bEA*E)) # proxy
  X = rbinom( n , size=A , prob=0.5 + bEX*E ) # count variable
  M = rnorm( n , bXM*X + bEM*E )
  d = data.frame(X=X,A=A,E=E,M=M)
  
  # return object
  if(!rep){
    # full data
    return(d)
    
  } else{
    # parameters
    b1 = coef( lm(M ~ X + E, data=d) )['X'] # unbiased effect (not possible)
    b2 = coef( lm(M ~ X, data=d) )['X'] # biased effects
    b3 = coef( lm(M ~ X + A, data=d) )['X'] # unbiased effect by proxy
    b = c(b1, b2, b3)
    names(b) = c('X','Xb','Xc')
    return( b )
    
  }
  
}

# relationships
d = f_sim(n=100, bEA=-8, bEX=-0.2, bEM=2, bXM=0, rep=F)

# pdf('descendant1_panel.pdf')
psych::pairs.panels(d[,-3])
# dev.off()
# notice cor(X,M)<0, when it should be cor(X,M)=0


# models
summary(lm(M ~ X, data=d)) # biased effect (X)
summary(lm(M ~ X + A, data=d))  # unbiased effect by proxy (A)
summary(lm(M ~ X + E, data=d)) # unbiased effect (not possible)



# sampling variation
# pdf('descendant1a_samplesize.pdf')
par(mfrow=c(2,1))
dsim = replicate( 1e4, f_sim(n=20, bEA=-8, bEX=-0.2, bEM=2, bXM=0, rep=T) )
f_plot1(dsim=dsim, ipar='X', n=20, xR=c(-0.5,0.5), by=0.1, 
        leg=T, legend=c('true','without A','with A'))

dsim = replicate( 1e4, f_sim(n=100, bEA=-8, bEX=-0.2, bEM=2, bXM=0, rep=T) )
f_plot1(dsim=dsim, ipar='X', n=100, xR=c(-0.5,0.5), by=0.1, leg=F)
par(mfrow=c(1,1))
# dev.off()
# X -> M, if we do not control for E (not possible) 
# X -> M, fixed if controlled by A (proxy) 
# equally biased with n=100, but more "confident" of X -> M









# Simulation details 2:
# A = cause variable (e.g. age)
# E = unobserved variable (e.g. instruction type)
#   E -> X: negative (specific E, less prob. X)
#   E -> M: positive (specific E, more M)
# X = observed covariate (e.g. teaching experience)
#   X -> M: null (to emphasize problem)
# M = outcome (e.g. attainment in mathematics)
# 
# Hypothesis:
# does X impact on M?
#
# DAG
gen_dag = "dag {
  E -> {X M};
  A -> E;
  X -> M;
  E [unobserved]
}"
dag_plot1 = dagitty( gen_dag )
coordinates(dag_plot1) = list( x=c(X=0,E=1,A=1,M=2) , 
                               y=c(X=0,E=-1,A=-2,M=0) )
drawdag( dag_plot1 )


# simulation
# n = simulation sample size
# bEA, bEX, bEM, bXM = simulated parameters
# rep = to use in replication
#
f_sim = function(n=100, bAE=1, bEX=-0.2, bEM=2, bXM=0, rep=F){
  
  # # test
  # n=100; bAE=1; bEX=-0.2; bEM=2; bXM=0; rep=F
  
  # sim
  A = round( c( rnorm( n/2 , 35 ), rnorm( n/2 , 45 ) ) ) # cause of a cause
  As = c(standardize(A))
  E = rbinom( n, size=1, prob=inv_logit( bAE*As ) ) # two types of education:
  X = rbinom( n , size=A , prob=0.5 + bEX*E ) # count variable
  M = rnorm( n , bXM*X + bEM*E )
  d = data.frame(X=X,A=A,E=E,M=M)
  
  # return object
  if(!rep){
    # full data
    return(d)
    
  } else{
    # parameters
    b1 = coef( lm(M ~ X + E, data=d) )['X'] # unbiased effect (not possible)
    b2 = coef( lm(M ~ X, data=d) )['X'] # biased effects
    b3 = coef( lm(M ~ X + A, data=d) )['X'] # unbiased effect by proxy
    b = c(b1, b2, b3)
    names(b) = c('X','Xb','Xp')
    return( b )
    
  }
  
}


# sampling variation
# pdf('descendant1b_samplesize.pdf')
par(mfrow=c(2,1))
dsim = replicate( 1e4, f_sim(n=20, bAE=1, bEX=-0.2, bEM=2, bXM=0, rep=T) )
f_plot1(dsim=dsim, ipar='X', n=20, xR=c(-0.5,0.5), by=0.1, 
        leg=T, legend=c('true','without A','with A'))

dsim = replicate( 1e4, f_sim(n=100, bAE=1, bEX=-0.2, bEM=2, bXM=0, rep=T) )
f_plot1(dsim=dsim, ipar='X', n=100, xR=c(-0.5,0.5), by=0.1, leg=F)
par(mfrow=c(1,1))
# dev.off()
# X -> M, if we do not control for E (not possible) 
# X -> M, fixed if controlled by A (proxy) 
# equally biased with n=100, but more "confident" of X -> M









# TO DO ####
# (descendant: proxies complex) ####

# Simulation details 3:
# Location: slides lecture 10, 2022 course
#
# Data details: 
# G = gender
#   G -> A: null (assumed no effect)
#   G -> D: negative (G=female, specific D's)
# D = department to be admitted
#   D -> A: negative (specific D's, less A)
# U = unobserved confound (e.g. ability, qualification)
#   U -> D: positive (more U, specific D's, the hardest to get in)
#   U -> A: positive (more U, more A)
#   U -> {T1, T2, T3}: positive (more U, higher Ti)
# T1, T2, T3: observed exam scores
# A = admission {0,1}
# 
# Hypothesis:
# G has no relationship with A? 
# D can be a confounding variable
#
# DAGs
gen_dag = "dag{
  {G D} -> A;
  G -> D
  U -> {D A};
  U [unobserved]
}"
dag_plot1 = dagitty( gen_dag )
coordinates(dag_plot1) = list( x=c(G=0,D=0,U=1,A=1) , 
                               y=c(G=0,D=-1,U=-1,A=0) )
drawdag( dag_plot1 )
# clearly you cannot control for U (unobserved)
# opens an unwanted backdoor path


# DAGs
gen_dag = "dag{
  {G D} -> A;
  G -> D
  U -> {D A T1 T2 T3};
  U [unobserved]
}"
dag_plot1 = dagitty( gen_dag )
coordinates(dag_plot1) = list( x=c(G=0,D=0,U=0.5,A=1,T1=1,T2=1,T3=1) , 
                               y=c(G=0,D=-1,U=-1,A=0,T1=-1,T2=-0.8,T3=-0.6) )
drawdag( dag_plot1 )
# but now you have observed values for scores (proxies)




# simulation
# n = simulation sample size
# pGD, pUD, ar = parameter of simulation
# rep = to use in replication
#
f_sim = function(n=1000, pGD=c(0.3,0.8), pUD=0.5, bU=1,
                 ar=list( c(0.1,0.1,0.1,0.3),
                          c(0.2,0.3,0.2,0.5)), rep=F){
  
  # test
  n = 100; pGD=c(0.3,0.8); pUD=0.2; bU=1
  ar=list( c(0.1,0.3,0.1,0.3),
           c(0.2,0.5,0.2,0.5))
  
  # sim
  G = sample( 1:2, size=n, replace=T ) # even gender distribution 
  U = rbern(n, 0.1) # ability, high (1) to average (0)
  # gender 1 tends to apply to department 1, 2 to 2
  # and G=1 with greater ability tend to apply to 2 as well
  D = rbern( n, ifelse( G==1, U*pUD, pGD[2] ) ) + 1 # D=2 discriminatory
  # matrix of acceptance rates [dept,gender]
  ar[[1]] = matrix( ar[[1]], nrow=2 ) # different acceptance per U
  ar[[2]] = matrix( ar[[2]], nrow=2 )
  
  # now we have proxies (different levels of reliability)
  T1 = rnorm(n, bU*U, 0.1)
  T2 = rnorm(n, bU*U, 0.5)
  T3 = rnorm(n, bU*U, 0.25)
  
  # simulate acceptance (interaction of G and D)
  p = sapply( 1:n , function(i){
    ifelse( U[i]==0, ar[[1]][D[i],G[i]], ar[[2]][D[i],G[i]] ) } )
  A = rbern( n , p ) 
  
  d = data.frame(U,T1,T2,T3,G=factor(G),D=factor(D),A)
  
  
  # return object
  if(!rep){
    # full data
    return(d)
    
  } else{
    # parameters
    b1 = coef( glm(A ~ -1 + G + D + U, data=d, family='binomial') )[c('G1','G2')] # unbiased effects
    b1 = inv_logit(b1[2]-b1[1])
    b2 = coef( glm(A ~ -1 + G + D, data=d, family='binomial') )[c('G1','G2')] # biased effects
    b2 = inv_logit(b2[2]-b2[1])
    b3 = coef( glm(A ~ -1 + G, data=d, family='binomial') )[c('G1','G2')] # biased effect
    b3 = inv_logit(b3[2]-b3[1])
    b = c(b1, b2, b3)
    names(b) = c('GC','GCb1','GCb2')
    return( b )
    
  }
  
}


# models
d = f_sim(n=1000, pGD=c(0.1,0.8), pUD=0.5, bU=1,
          ar=list( c(0.1,0.1,0.1,0.3), # discrimination
                   c(0.2,0.3,0.2,0.5)), rep=F) 

# pdf('descendant3_panel.pdf')
psych::pairs.panels(d[-1])
# dev.off()


m_res = glm(A ~ -1 + G, data=d, family='binomial') # regression without intercept
# summary(m_res) 
K = matrix(c(1, -1), 1) # contrast of interest
t = glht(m_res, linfct = K)
summary(t) # only gender, biased effect
inv_logit(coef(m_res)[1]) - inv_logit(coef(m_res)[2]) # probability
# we can infer discrimination, but is it the whole story?
# besides this is the total effect


m_res = glm(A ~ -1 + G + D, data=d, family='binomial')
# summary(m_res) 
K = matrix(c(1, -1, 0), 1) # contrast of interest
t = glht(m_res, linfct = K)
summary(t) # G and D, less biased effect
inv_logit(coef(m_res)[1]) - inv_logit(coef(m_res)[2]) # probability
# controlling by D tells you there is no discrimination
# it would be wrong (we know there is)


m_res = glm(A ~ -1 + G + D + T2, data=d, family='binomial')
# summary(m_res) 
K = matrix(c(1, -1, 0, 0), 1) # contrast of interest
t = glht(m_res, linfct = K)
summary(t) # G and D, less biased effect
inv_logit(coef(m_res)[1]) - inv_logit(coef(m_res)[2]) # probability
# we can control by unobserved confounder (Ti)
# but not always fix the problem (sometimes worsen it)
# why?, measurement error


sumT = with(d, T1+T2+T3) # assuming we can sum the scores
m_res = glm(A ~ -1 + G + D + sumT, data=d, family='binomial')
# summary(m_res) 
K = matrix(c(1, -1, 0, 0), 1) # contrast of interest
t = glht(m_res, linfct = K)
summary(t) # G and D, less biased effect
inv_logit(coef(m_res)[1]) - inv_logit(coef(m_res)[2]) # probability
# we can control by unobserved confounder (sumT)
# but is does not fix the problem (sometimes worsen it)
# why?, measurement error





# sampling variation
# pdf('descendant2_samplesize.pdf')
par(mfrow=c(2,2))
dsim = replicate( 1e4, f_sim(n=300, pGD=c(0.1,0.8), pUD=0.5, bU=1,
                             ar=list( c(0.1,0.1,0.1,0.3), # discrimination
                                      c(0.2,0.3,0.2,0.5)), rep=T) )
f_plot1(dsim=dsim, ipar='GC', n=300, xR=c(-0.5,2.2), by=0.1, 
        leg=T, legend=c('true','with G and D','with G only'))

dsim = replicate( 1e4, f_sim(n=500, pGD=c(0.1,0.8), pUD=0.5, bU=1,
                             ar=list( c(0.1,0.1,0.1,0.3), # discrimination
                                      c(0.2,0.3,0.2,0.5)), rep=T) )
f_plot1(dsim=dsim, ipar='GC', n=500, xR=c(-0.5,2.2), by=0.1, leg=F)
par(mfrow=c(1,1))
# dev.off()





# how should we fix it?
# putting data into list
dlist = with(d, list(N=nrow(d),G=G,D=D,T1=T1,T2=T2,T3=T3,A=A))

# total effects
m = ulam(
  alist(
    # A model
    A ~ bernoulli(p),
    logit(p) <- a[G,D],
    matrix[G,D]:a ~ normal(0,1) ), 
  data=dlist, chains=4, cores=4 ) 


# proxy effects
mp = ulam(
  alist(
    # A model
    A ~ bernoulli(p),
    logit(p) <- a[G,D] + b*u[i], # notice is generated from Ti
    matrix[G,D]:a ~ normal(0,1), # priors
    b ~ normal(0,1),
    
    # u and T model
    vector[N]:u ~ normal(0,1), # we estimate U from Ti
    T1 ~ normal(u,tau[1]), # three different reliabilities
    T2 ~ normal(u,tau[2]),
    T3 ~ normal(u,tau[3]),
    vector[3]:tau ~ exponential(1) ), 
  data=dlist, chains=4, cores=4, 
  constraints=list(b="lower=0") ) # we know direction of U->Ti


# # results
# precis(m, 3, pars=c("a"))
# precis(mp, 3, pars=c("a","b","tau"))
# # relative effects (hard to eyeball)


# contrasts
post_m = extract.samples(m)
post_m$cont_D1 = post_m$a[,1,1] - post_m$a[,2,1]
post_m$cont_D2 = post_m$a[,1,2] - post_m$a[,2,2]

post_mp = extract.samples(mp)
post_mp$cont_D1 = post_mp$a[,1,1] - post_mp$a[,2,1]
post_mp$cont_D2 = post_mp$a[,1,2] - post_mp$a[,2,2]


# plot 
dens( post_m$cont_D1, lwd=1, col=3, xlim=c(-2,1.5), 
      xlab="F-M contrast in each department" )
abline(v=0,lwd=1, col=1, lty=2)
abline(v=mean(post_m$cont_D1),lwd=1, col=3, lty=2)
dens( post_m$cont_D2, lwd=1, col=4, add=T )
abline(v=mean(post_m$cont_D2),lwd=1, col=4, lty=2)
dens( post_mp$cont_D1, lwd=3, col=3, add=T )
abline(v=mean(post_mp$cont_D1),lwd=3, col=3, lty=2)
dens( post_mp$cont_D2, lwd=3, col=4, add=T )
abline(v=mean(post_mp$cont_D2),lwd=3, col=4, lty=2)
# notice for D=2 the effect is even less favorable for women
# because U mask the discrimination, but we know is shouldn't


with(d, plot( jitter(U), apply(post_mp$u,2,mean), 
              xlab="u (true)", ylab="posterior mean u", 
              xaxt="n", lwd=2, 
              col=ifelse(G==1,col.alpha(2,0.6), col.alpha(4,0.6)) ) )
axis(1,at=0:1,labels=0:1)
# we managed to identify the different U by G













# (descendant: case control) ####
#
# Location: lecture 06, slides, 2022 course, Cinelli et al, 2021 (p. 8, 19)
#
# also an example of:
#   - Virtual collider
#
# Data details: 
# E = education
#   E -> H: negative (more E, less H)
# H = hours in occupation (continuum)
#   H -> I: negative (more H, less I)
# I = Income
# 
# Hypothesis:
# does E affect O?
#
# DAGs
gen_dag = "dag{ 
  E -> H;
  H -> I;
}"
dag_plot1 = dagitty( gen_dag )
coordinates(dag_plot1) = list( x=c(E=0,H=1,I=1) , 
                               y=c(E=0,H=0,I=-1) )
drawdag( dag_plot1 )



# simulation
# n = simulation sample size
# bEO, bHI = parameter of simulation
# rep = to use in replication
#
f_sim = function(n=100, bEO=-1, bHI=-1, rep=F){
  
  # # test
  # n=100; bEO=-1; bHI=1; rep=F
  
  # sim
  E = rnorm( n ) 
  H = rnorm( n , bEO*E ) 
  I = rnorm( n , bHI*H ) 
  d = data.frame(E,H,I)
  
  # return object
  if(!rep){
    # full data
    return(d)
    
  } else{
    # parameters
    b1 = coef( lm(H ~ E, data=d) )['E'] # unbiased effects
    b2 = coef( lm(H ~ E + I, data=d) )['E'] # biased effect
    b = c(b1, b2)
    names(b) = c('E','Eb')
    return( b )
    
  }
  
}

# relationships
d = f_sim(n=100, bEO=-1, bHI=1, rep=F)

# pdf('descendant3_panel.pdf')
psych::pairs.panels(d)
# dev.off()
# notice equal E -> O


# models
summary(lm(H ~ E, data=d)) # unbiased effects
summary(lm(H ~ E + I, data=d)) # biased effects


# sampling variation
# pdf('descendant3_samplesize.pdf')
par(mfrow=c(2,1))
dsim = replicate( 1e4, f_sim(n=20, bEO=-1, bHI=-1, rep=T) )
f_plot1(dsim=dsim, ipar='E', n=20, xR=c(-2,0.5), by=0.5, 
        leg=T, legend=c('true','biased'))

dsim = replicate( 1e4, f_sim(n=100, bEO=-1, bHI=-1, rep=T) )
f_plot1(dsim=dsim, ipar='E', n=100, xR=c(-2,0.5), by=0.5, leg=F)
par(mfrow=c(1,1))
# dev.off()
# E -> H (negative), but underestimated
# equally biased with n=100, but more "confident" of E -> H (underestimated)







# simulation 2
#
# Location: lecture 06, slides, 2022 course, Cinelli et al, 2021 (p. 8, 19)
#
# DAGs
gen_dag = "dag{ 
  X -> {Y Z};
}"
dag_plot1 = dagitty( gen_dag )
coordinates(dag_plot1) = list( x=c(Z=0,X=0,Y=1) , 
                               y=c(Z=-1,X=0,Y=0) )
drawdag( dag_plot1 )



# n = simulation sample size
# bEO, bHI = parameter of simulation
# rep = to use in replication
#
f_sim = function(n=100, bXY=1, bXZ=1, rep=F){
  
  # # test
  # n=100; bXY=1; bXZ=1; rep=F
  
  # sim
  X = rnorm( n ) 
  Z = rnorm( n , bXZ*X ) 
  Y = rnorm( n , bXY*X ) 
  d = data.frame(Z,X,Y)
  
  # return object
  if(!rep){
    # full data
    return(d)
    
  } else{
    # parameters
    b1 = coef( lm(Y ~ X, data=d) )['X'] # unbiased effects
    b2 = coef( lm(Y ~ X + Z, data=d) )['X'] # less precise effect
    s1 = summary( lm(Y ~ X, data=d) )$coefficients[,2][2] # se
    s2 = summary( lm(Y ~ X + Z, data=d) )$coefficients[,2][2]
    b = c(b1, b2, s1, s2)
    names(b) = c('X','Xb','sX','sXb')
    return( b )
    
  }
  
}

# relationships
d = f_sim(n=100, bXY=1, bXZ=1, rep=F)



# sampling variation
# pdf('descendant3b_samplesize.pdf')
par(mfrow=c(2,2))
dsim = replicate( 1e4, f_sim(n=20, bXY=1, bXZ=1, rep=T) )
f_plot1(dsim=dsim, ipar='X', n=20, xR=c(-0.5,2), by=0.5, 
        leg=T, legend=c('only X','X and Z'))
f_plot1(dsim=dsim, ipar='sX', n=20, xR=c(0,0.5), by=0.1, leg=F)

dsim = replicate( 1e4, f_sim(n=100, bXY=1, bXZ=1, rep=T) )
f_plot1(dsim=dsim, ipar='X', n=100, xR=c(-0.5,2), by=0.5, leg=F)
f_plot1(dsim=dsim, ipar='sX', n=100, xR=c(0,0.5), by=0.1, leg=F)
par(mfrow=c(1,1))
# dev.off()









# (descendant: rc, outcome only) ####
#
# Location: chapter 15 (p. 491)
# see also: https://sites.google.com/view/robertostling/home/teaching
#
# also known as:
#   - measurement error (ME)
#   - residual confounding (RC)
#
# data details: 
# A = median age at marriage
#   A -> M: negative (more A, less M)
#   A -> D: negative (more A, less D)
# M = marriage rate
#   M -> D: null (expected)
# D_obs = observed divorce rate
# e_D = (unobserved) measurement error in divorce rate
# D = (unobserved) 'true' divorce rate
# 
# Hypothesis:
# does M or A impact on D (true)?
#
# DAG
gen_dag = "dag{ 
  A -> {M D};
  M -> D;
  {D e_D} -> D_obs;
  
  D[latent];
  e_D[latent];
}"

dagME = dagitty( gen_dag )

coordinates( dagME ) = list(
  x=c(M=-0.2, A=0, D=0, D_obs=0.2, e_D=0.4),
  y=c(M=0, A=-0.2, D=0, D_obs=0, e_D=0) )

drawdag(dagME)




# n = simulation sample size
# bAM, bAD, bMD, sD = parameter of simulation
# rep = to use in replication
#
f_sim = function(n=100, bAM=-1, bAD=-1, bMD=0, sD=1, var='A', rep=F){
  
  # # test
  # n=100; bAM=-1; bAD=-1; bMD=0; sD=0.5; var='A'; rep=F

  # sim
  A = rnorm( n )
  M = rnorm( n , bAM*A )
  D_true = rnorm( n, bAD*A + bMD*M )
  D_obs = rnorm( n, mean=D_true, sd=sD )
  d = data.frame( A, M, D_true, D_obs )
  
  # plot
  if(!rep){
    return(d)
  } else{
    res = summary(lm(D_true ~ A + M, data=d))
    idx = str_detect( rownames(res$coefficients), var)
    b1 = res$coefficients[idx, c('Estimate','Std. Error')]
    
    res = summary(lm(D_obs ~ A + M, data=d))
    idx = str_detect( rownames(res$coefficients), var)
    b2 = res$coefficients[idx, c('Estimate','Std. Error')]
    b = c(b1, b2)
    names(b) = c('bAt','sAt','bAo','sAo')
    return(b)
  }
  
}


# what's going on
s = c(0.2,1,2)


# pdf('descendant4_me1.pdf', width=14, height=7)
par(mfrow=c(1,3))
for(i in 1:length(s)){
  set.seed(123456)
  d = f_sim(n=100, bAM=-1, bAD=-1, bMD=0, sD=s[i], rep=F)
  f_plot3(dsim=d, sX=s[i], a=0.2, xR=c(-4,4), yR=c(-7,7))
  
  if(i==1){
    legend('bottomleft', legend=c('latent','observed'), pch=19, bty='n',
           col=c(col.alpha('black', 0.3), col.alpha('red', 0.3)))
  }
}
par(mfrow=c(1,1))
# dev.off()
# not so pervasive


# pdf('descendant4_me2.pdf', width=14, height=7)
par(mfrow=c(1,3))
for(i in 1:length(s)){
  set.seed(123456)
  d = f_sim(n=100, bAM=-1, bAD=-1, bMD=0, sD=s[i], rep=F)
  f_plot2(dsim=d, sX=s[i], a=0.2, xR=c(-4,4), yR=c(-7,7))
  
  if(i==1){
    legend('bottomleft', legend=c('latent','observed'), pch=19, bty='n',
           col=c(col.alpha('black', 0.3), col.alpha('red', 0.3)))
  }
}
par(mfrow=c(1,1))
# dev.off()
# not so pervasive




d = f_sim(n=100, bAM=-1, bAD=-1, bMD=0, sD=1, rep=F)


# pdf('descendant4_panel.pdf')
psych::pairs.panels(d[,c('A','M','D_obs')])
# dev.off()


# models
summary( lm(D_true ~ A + M, data=d) )
summary( lm(D_obs ~ A + M, data=d) )



# sampling variation
# pdf('descendant4_samplesize.pdf')
par(mfrow=c(3,2))
# i=0.2
for(i in c(0.2, 1, 2)){
  # dsim = replicate( 1e4, f_sim(n=20, bAD=-1, sD=i, rep=T) )
  # f_plot1(dsim=dsim, ipar='bA', n=20, xR=c(-2,2), by=0.2, 
  #         leg=T, legend=c('true','observed'))
  # f_plot1(dsim=dsim, ipar='sA', n=20, xR=c(0,0.5), by=0.1, leg=F)
  
  dsim = replicate( 1e4, f_sim(n=100, bAD=-1, sD=i, rep=T) )
  f_plot1(dsim=dsim, ipar='bA', n=100, xR=c(-2,2), by=0.2, 
          leg=T, legend=c('true', paste0('sD = ', i)), loc='topright')
  f_plot1(dsim=dsim, ipar='sA', n=100, xR=c(0,0.5), by=0.1, leg=F)
}
par(mfrow=c(1,1))
# dev.off()



# frequentist fix
require(metafor)
me_model1 = rma(yi = D_obs,
                sei = rep(1, 100),
                data = d,
                method = "REML",
                mods = ~ A + M,
                test = "t")
summary(me_model1)
# different estimates and SE, but the correction is done
# notice we do not estimate a latent variable



# Bayesian fix
dlist = list(
  N = nrow(d),
  D_obs = d$D_obs,
  D_sd = rep(1, nrow(d)),
  A = d$A,
  M = d$M
)

me_model2 = ulam(
  alist(
    D_obs ~ dnorm( D_true , D_sd ),
    vector[N]:D_true ~ dnorm( mu , sigma ),
    mu <- a + bA*A + bM*M,
    a ~ dnorm(0,0.2),
    bA ~ dnorm(0,0.5),
    bM ~ dnorm(0,0.5),
    sigma ~ dexp(1)) , 
  data=dlist , chains=4 , cores=4 )
precis( me_model , depth=1 )









# data (true example)
# library(rethinking)
data(WaffleDivorce)
d = WaffleDivorce

d = list(
  D_obs = standardize( d$Divorce ),
  D_sd = d$Divorce.SE / sd( d$Divorce ),
  M = standardize( d$Marriage ),
  A = standardize( d$MedianAgeMarriage ),
  N = nrow(d) 
)
# str(d)


# chapter 5 model # unbiased effects, less efficient
m5.3 = quap(
  alist(
    D_obs ~ dnorm( mu , sigma ) ,
    mu <- a + bA*A + bM*M ,
    a ~ dnorm( 0 , 0.2 ) ,
    bM ~ dnorm( 0 , 0.5 ) ,
    bA ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )) , 
  data = d)
precis( m5.3 )


# ulam # unbiased effects, a bit more efficient
m15.1b = ulam(
  alist(
    D_obs ~ dnorm( D_true , D_sd ),
    vector[N]:D_true ~ dnorm( mu , sigma ),
    mu <- a + bA*A + bM*M,
    a ~ dnorm(0,0.2),
    bA ~ dnorm(0,0.5),
    bM ~ dnorm(0,0.5),
    sigma ~ dexp(1)) , 
  data=d , chains=4 , cores=4 )
precis( m15.1b , depth=1 )












# (descendant: rc, predictor only) ####
#
# Location: chapter 15 (p. 495)
#
# also known as:
#   - measurement error (ME)
#   - residual confounding (RC)
#
# data details: 
# A_obs = observed median age at marriage
# e_A = (unobserved) measurement error in median age at marriage
# A = (unobserved) 'true' median age at marriage
#   A -> M: negative (more A, less M)
#   A -> D: negative (more A, less D)
# M_obs = observed marriage rate
# e_M = (unobserved) measurement error in marriage rate
# M = (unobserved) 'true' marriage rate
#   M -> D: null (expected)
# D_obs = observed divorce rate
# e_D = (unobserved) measurement error in divorce rate
# D = (unobserved) 'true' divorce rate
# 
# Hypothesis:
# does A (true) or M impact on D_obs?
#
# DAG
gen_dag = "dag{ 
  A -> {M D};
  M -> {D M_obs};
  e_M -> M_obs;

  M[latent];
  e_M[latent];
}"
dagME = dagitty( gen_dag )

coordinates( dagME ) = list(
  x=c(e_M=-0.6, M_obs=-0.4, M=-0.2, A=0, D=0),
  y=c(e_M=0, M_obs=0, M=0, A=-0.2, D=0) )

drawdag(dagME)



# n = simulation sample size
# bAM, bAD, bMD, sA = parameter of simulation
# rep = to use in replication
#
f_sim = function(n=100, bAM=-1, bAD=-1, bMD=1, sM=0.5, var='M', rep=F){
  
  # # test
  # n=100; bAM=-1; bAD=-1; bMD=0; sM=0.5; var='A'; rep=F
  
  # sim
  A = rnorm( n )
  M_true = rnorm( n , bAM*A )
  M_obs = rnorm( n , mean=M_true, sd=sM )
  D = rnorm( n , bAD*A + bMD*M_true)
  d = data.frame(A, M_true, M_obs, D)
  
  # plot
  if(!rep){
    return(d)
  } else{
    res = summary(lm(D ~ A + M_true, data=d))
    idx = str_detect( rownames(res$coefficients), var)
    b1 = res$coefficients[idx, c('Estimate','Std. Error')]
    
    res = summary(lm(D ~ A + M_obs, data=d))
    idx = str_detect( rownames(res$coefficients), var)
    b2 = res$coefficients[idx, c('Estimate','Std. Error')]
    b = c(b1, b2)
    names(b) = c('bMt','sMt','bMo','sMo')
    return(b)
  }
  
}


# what's going on
s = c(0.2,1,2)


# pdf('descendant5_me1.pdf', width=14, height=7)
par(mfrow=c(1,3))
for(i in 1:length(s)){
  set.seed(123456)
  d = f_sim(n=100, bAM=-1, bAD=-1, bMD=1, sM=s[i], rep=F)
  f_plot3(dsim=d, sX=s[i], a=0.2, xR=c(-7,7), yR=c(-4,4))
  
  if(i==1){
    legend('bottomleft', legend=c('latent','observed'), pch=19, bty='n',
           col=c(col.alpha('black', 0.3), col.alpha('red', 0.3)))
  }
}
par(mfrow=c(1,1))
# dev.off()
# not so pervasive


# pdf('descendant5_me2.pdf', width=14, height=7)
par(mfrow=c(1,3))
for(i in 1:length(s)){
  set.seed(123456)
  d = f_sim(n=100, bAM=-1, bAD=-1, bMD=1, sM=s[i], rep=F)
  f_plot2(dsim=d, sX=s[i], a=0.2, xR=c(-7,7), yR=c(-4,4))
  
  if(i==1){
    legend('bottomleft', legend=c('latent','observed'), pch=19, bty='n',
           col=c(col.alpha('black', 0.3), col.alpha('red', 0.3)))
  }
}
par(mfrow=c(1,1))
# dev.off()
# not so pervasive



d = f_sim(n=100, bAM=-1, bAD=-1, bMD=1, sM=1, rep=F)

# pdf('descendant5_panel.pdf')
psych::pairs.panels(d[,c('A','M_obs','D')])
# dev.off()


# models
summary( lm(D ~ A + M_true, data=d) )
summary( lm(D ~ A + M_obs, data=d) )



# sampling variation
# pdf('descendant5_samplesize.pdf')
par(mfrow=c(3,2))
# i=0.2
for(i in c(0.2, 1, 2)){
  # dsim = replicate( 1e4, f_sim(n=20, bAM=-1, bAD=-1, bMD=1, sM=i, rep=T) )
  # f_plot1(dsim=dsim, ipar='bM', n=20, xR=c(-2,2), by=0.2,
  #         leg=T, legend=c('true','observed'))
  # f_plot1(dsim=dsim, ipar='sM', n=20, xR=c(0,0.5), by=0.1, leg=F)
  
  dsim = replicate( 1e4, f_sim(n=100, bAM=-1, bAD=-1, bMD=1, sM=i, rep=T) )
  f_plot1(dsim=dsim, ipar='bM', n=100, xR=c(-2,2), by=0.2, 
          leg=T, legend=c('true', paste0('sM = ', i)), loc='topleft')
  f_plot1(dsim=dsim, ipar='sM', n=100, xR=c(0,0.5), by=0.1, leg=F)
}
par(mfrow=c(1,1))
# dev.off()




# frequentist fix? (also uses bayesian model)
require(eivtools)
serror = diag(c(0, 1))
dimnames(serror) = list(c('A','M_obs'), c('A','M_obs'))
me_model = eivreg(D ~ A + M_obs, data=d, Sigma_error=serror)
summary(me_model)



# bayesian fix
dlist = list(
  N = nrow(d),
  D = d$D,
  A = d$A,
  M_obs = d$M_obs,
  M_sd = rep(1, nrow(d))
)

me_model = ulam(
  alist(
    D ~ dnorm( mu , sigma ),
    mu <- a + bA*A + bM*M_true[i],
    M_obs ~ dnorm( M_true , M_sd ),
    vector[N]:M_true ~ dnorm( 0 , sigma ),
    a ~ dnorm(0,0.2),
    bA ~ dnorm(0,0.5),
    bM ~ dnorm(0,0.5),
    sigma ~ dexp(1)) , 
  data=dlist , chains=4 , cores=4 )
precis( me_model , depth=1 )









# (descendant: rc, out and pred) ####
#
# Location: chapter 15 (p. 495)
#
# also known as:
#   - measurement error (ME)
#   - residual confounding (RC)
#
# data details: 
# A = median age at marriage
#   A -> M: negative (more A, less M)
#   A -> D: negative (more A, less D)
# M_obs = observed marriage rate
# e_M = (unobserved) measurement error in marriage rate
# M = (unobserved) 'true' marriage rate
#   M -> D: null (expected)
# D_obs = observed divorce rate
# e_D = (unobserved) measurement error in divorce rate
# D = (unobserved) 'true' divorce rate
# 
# Hypothesis:
# does M (true) or A impact on D (true)?
#
# DAG
gen_dag = "dag{ 
  A -> {M D};
  M -> D;
  {M e_M} -> M_obs;
  {D e_D} -> D_obs;
  
  D[latent];
  e_D[latent];
  M[latent];
  e_M[latent];
}"
dagME = dagitty( gen_dag )
coordinates( dagME ) = list(
  x=c(A=-0.2, M=0, D=0, D_obs=0.2, M_obs=0.2, e_D=0.4, e_M=0.4),
  y=c(A=0, M=-0.2, D=0.2, D_obs=0.2, M_obs=-0.2, e_D=0.2, e_M=-0.2) )
drawdag(dagME)





# n = simulation sample size
# bAM, bAD, bMD, sA = parameter of simulation
# rep = to use in replication
#
f_sim = function(n=100, bAM=-1, bAD=-1, bMD=1, s=0.5, rep=F){
  
  # # test
  # n=100; bAM=-1; bAD=-1; bMD=0; sM=0.5; var='A'; rep=F
  
  # sim
  A = rnorm( n )
  M_true = rnorm( n , bAM*A )
  M_obs = rnorm( n , mean=M_true, sd=s )
  D_true = rnorm( n, bAD*A + bMD*M_true )
  D_obs = rnorm( n, mean=D_true, sd=s )
  d = data.frame(A, M_true, M_obs, D_true, D_obs)
  
  # plot
  if(!rep){
    return(d)
  } else{
    res = summary(lm(D_true ~ A + M_true, data=d))
    idx = str_detect( rownames(res$coefficients), 'M_true')
    b1 = res$coefficients[idx, c('Estimate','Std. Error')]
    
    res = summary(lm(D_obs ~ A + M_obs, data=d))
    idx = str_detect( rownames(res$coefficients), 'M_obs')
    b2 = res$coefficients[idx, c('Estimate','Std. Error')]
    b = c(b1, b2)
    names(b) = c('bMt','sMt','bMo','sMo')
    return(b)
  }
  
}


d = f_sim(n=100, bAM=-1, bAD=-1, bMD=1, s=1, rep=F)


# models
summary( lm(D_true ~ A + M_true, data=d) )
summary( lm(D_obs ~ A + M_obs, data=d) )


# sampling variation
# pdf('descendant6_samplesize.pdf')
par(mfrow=c(3,2))
# i=0.2
for(i in c(0.2, 1, 2)){
  # dsim = replicate( 1e4, f_sim(n=100, bAM=-1, bAD=-1, bMD=0, sA=i, var='A', rep=T) )
  # f_plot1(dsim=dsim, ipar='bA', n=20, xR=c(-2,2), by=0.2,
  #         leg=T, legend=c('true','observed'))
  # f_plot1(dsim=dsim, ipar='sA', n=20, xR=c(0,0.5), by=0.1, leg=F)
  
  dsim = replicate( 1e4, f_sim(n=100, bAM=-1, bAD=-1, bMD=1, s=i, rep=T) )
  f_plot1(dsim=dsim, ipar='bM', n=100, xR=c(-2,2), by=0.2, 
          leg=T, legend=c('true', paste0('sD = sM = ', i)), loc='topleft')
  f_plot1(dsim=dsim, ipar='sM', n=100, xR=c(0,0.5), by=0.1, leg=F)
}
par(mfrow=c(1,1))
# dev.off()





# only bayesian solution
d = list(
  D_obs = d$D_obs,
  D_sd = rep(1, 100),
  M_obs = d$M_obs,
  M_sd = rep(1, 100),
  A = d$A,
  N = nrow(d)
)

me_model = ulam(
  alist(
    D_obs ~ dnorm( D_true , D_sd ),
    vector[N]:D_true ~ dnorm( mu , sigma ),
    mu <- a + bA*A + bM*M_true[i],
    M_obs ~ dnorm( M_true , M_sd ),
    vector[N]:M_true ~ dnorm( 0 , 1 ),
    a ~ dnorm(0,0.2),
    bA ~ dnorm(0,0.5),
    bM ~ dnorm(0,0.5),
    sigma ~ dexp( 1 )) , 
  data=d , chains=4 , cores=4 )
precis(me_model)









# true data
data(WaffleDivorce)
d = WaffleDivorce

d = list(
  D_obs = standardize( d$Divorce ),
  D_sd = d$Divorce.SE / sd( d$Divorce ),
  M_obs = standardize( d$Marriage ),
  M_sd = d$Marriage.SE / sd( d$Marriage ),
  A = standardize( d$MedianAgeMarriage ),
  N = nrow(d)
)
# str(d)


# chapter 5 model # unbiased effects, less efficient
simple_model = quap(
  alist(
    D_obs ~ dnorm( mu , sigma ) ,
    mu <- a + bA*A + bM*M_obs ,
    a ~ dnorm( 0 , 0.2 ) ,
    bM ~ dnorm( 0 , 0.5 ) ,
    bA ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )) , 
  data = d)
precis( simple_model )


# ulam # unbiased effects, a bit more efficient
me_model = ulam(
  alist(
    D_obs ~ dnorm( D_true , D_sd ),
    vector[N]:D_true ~ dnorm( mu , sigma ),
    mu <- a + bA*A + bM*M_true[i],
    M_obs ~ dnorm( M_true , M_sd ),
    vector[N]:M_true ~ dnorm( 0 , 1 ),
    a ~ dnorm(0,0.2),
    bA ~ dnorm(0,0.5),
    bM ~ dnorm(0,0.5),
    sigma ~ dexp( 1 )) , 
  data=d , chains=4 , cores=4 )
precis(me_model)














# (descendant: censoring, outcome only) ####
#
# Location: https://sites.google.com/view/robertostling/home/teaching
#
# also known as:
#   - special case of measurement error
#
# data details (to make a point): 
# SES = socio-economical status
#   SES -> If: positive (more SES, more I)
#   SES -> E: positive (more SES, more E)
# E = year of education (standardized)
#   E -> If: positive (more E, more If)
# If = income (fully observed)
#   If -> It: positive (one is an observation of the other)
# It = truncated observation of If
# 
# Hypothesis:
# How E explains I?
#
# DAG
gen_dag = "dag{ 
  SES -> {E I_full};
  E -> I_full;
  {I_full e_C} -> I_cens;
  
  I_full[latent];
  e_C[latent];
}"

dagME = dagitty( gen_dag )

coordinates( dagME ) = list(
  x=c(E=-0.2, SES=0, I_full=0, I_cens=0.2, e_C=0.4),
  y=c(E=0, SES=-0.2, I_full=0, I_cens=0, e_C=0) )

drawdag(dagME)



# simulation
# n = simulation sample size
# bSE, bEI, bSI, xT, yT = simulation parameters (last two are thresholds)
# rep = to use in replication
#
f_sim = function(n=100, bSE=1, bEI=1, bSI=0.3, xT=c(-3,3), yT=c(-3,3), rep=F){
  
  # # test
  # n = 100; bSE=1; bEI=1; bSI=0.3; xT=c(-2,2); yT=c(-2,2); rep=F
  
  # sim
  SES = rnorm( n )
  E = rnorm( n , bSE*SES)
  I = rnorm( n , bEI*E + bSI*SES)
  d = data.frame( SES=SES, E_full=E, E_cens=E, I_full=I, I_cens=I)
  
  if( !is.null(xT) ){
    d$E_cens[E <= xT[1]] = xT[1]
    d$E_cens[E >= xT[2]] = xT[2]
  }
  
  if( !is.null(yT) ){
    d$I_cens[I <= yT[1]] = yT[1]
    d$I_cens[I >= yT[2]] = yT[2]  
  }
  
  
  # plot
  if(!rep){
    return(d)
  } else{
    res = summary(lm(I_full ~ E_full + SES, data=d))
    idx = str_detect( rownames(res$coefficients), 'E_full')
    b1 = res$coefficients[idx, c('Estimate','Std. Error')]
    
    res = summary(lm(I_cens ~ E_cens + SES, data=d))
    idx = str_detect( rownames(res$coefficients), 'E_cens')
    b2 = res$coefficients[idx, c('Estimate','Std. Error')]
    b = c(b1, b2)
    names(b) = c('bEf','sEf','bEt','sEt')
    return(b)
  }
  
}



# what's going on
out_range = list(a=c(-5,5), b=c(-5,1), c=c(-5,-1))


# pdf('descendant7_cens1.pdf', width=14, height=7)
par(mfrow=c(1,3))
# i=1
for(i in 1:length(out_range)){
  set.seed(123456)
  d = f_sim(100, bSE=1, bEI=1, bSI=0.3, xT=c(-10,10), yT=out_range[[i]], rep=F)
  
  f_plot4(dsim=d, var='I', var_range =out_range[[i]] )
  if(i==1){
    legend('topleft', legend=c('true data','censored data'), pch=c(19, 1), 
           bty='n', col=c(col.alpha('black', 0.3), 'red'))
  }
}
par(mfrow=c(1,1))
# dev.off()




# pdf('descendant7_cens2.pdf', width=14, height=7)
par(mfrow=c(1,3))
# i=1
for(i in 1:length(out_range)){
  set.seed(123456)
  d = f_sim(100, bSE=1, bEI=1, bSI=0.3, xT=c(-10,10), yT=out_range[[i]], rep=F)
  
  f_plot5(dsim=d, var='I', b=ifelse(i==1, 30, 10), var_range=out_range[[i]])
  if(i==1){
    legend('topleft', legend=c('true data','censored data'), pch=c(19, 1), 
           bty='n', col=c(col.alpha('black', 0.3), 'red'))
  }
}
par(mfrow=c(1,1))
# dev.off()





d = f_sim(100, bSE=1, bEI=1, bSI=0.3, xT=c(-10,10), yT=c(-5,0), rep=F)


# pdf('descendant7_panel.pdf')
psych::pairs.panels(d[,c('SES','E_cens','I_cens')])
# dev.off()


# models
summary( lm(I_full ~ E_full + SES, data=d) )
summary( lm(I_cens ~ E_cens + SES, data=d) )



# sampling variation
out_range = list(a=c(-5,5), b=c(-5,1), c=c(-5,-1))

# pdf('descendant7_samplesize.pdf')
# i=1
par(mfrow=c(3,2))
for(i in 1:length(out_range)){
  # dsim = replicate( 1e4, f_sim(20, bSE=1, bEI=1, bSI=0.3, xT=c(-10,10), yT=c(-5,0), rep=T) )
  # f_plot1(dsim=dsim, ipar='bM', n=20, xR=c(-2,2), by=0.2,
  #         leg=T, legend=c('true','observed'))
  # f_plot1(dsim=dsim, ipar='sM', n=20, xR=c(0,0.5), by=0.1, leg=F)
  
  dsim = replicate( 1e4, f_sim(100, bSE=1, bEI=1, bSI=0.3, xT=c(-10,10), yT=out_range[[i]], rep=T) )
  f_plot1(dsim=dsim, ipar='bE', n=100, xR=c(-2,2), by=0.2, 
          leg=T, legend=c('true', paste0('I in (', out_range[[i]][1], ', ', out_range[[i]][2], ')')), loc='topleft')
  f_plot1(dsim=dsim, ipar='sE', n=100, xR=c(0,0.5), by=0.1, leg=F)
}
par(mfrow=c(1,1))
# dev.off()





# bayesian fix
idx_rc = d$I_cens >= max(d$I_cens)
idx_nc = !idx_rc

dlist = list(
  N = nrow(d),
  E_rc = d$E_cens[idx_rc],
  SES_rc = d$SES[idx_rc],
  I_rc = d$I_cens[idx_rc],
  E_nc = d$E_cens[idx_nc],
  SES_nc = d$SES[idx_nc],
  I_nc = d$I_cens[idx_nc]
)


cens_model = ulam(
  alist(
    I_nc ~ dnorm( mu_nc , sigma ),
    I_rc ~ custom( normal_lccdf( I_rc | mu_rc , sigma) ),
    mu_nc <- a + bE*E_nc + bS*SES_nc,
    mu_rc <- a + bE*E_rc + bS*SES_rc,
    a ~ dnorm(0,0.2),
    bE ~ dnorm(0,0.5),
    bS ~ dnorm(0,0.5),
    sigma ~ dexp(1)) , 
  data=dlist , chains=4 , cores=4 )
precis( cens_model , depth=1 )

# stancode(cens_model)
# # what are you doing











# (descendant: censoring, predictor only) ####
#
# Location: https://sites.google.com/view/robertostling/home/teaching
#
# also known as:
#   - special case of measuring error
#
# data details (to make a point): 
# SES = socio-economical status
#   SES -> If: positive (more SES, more I)
#   SES -> E: positive (more SES, more E)
# E = year of education (standardized)
#   E -> If: positive (more E, more If)
# If = income (fully observed)
#   If -> It: positive (one is an observation of the other)
# It = truncated observation of If
# 
# Hypothesis:
# How E explains I?
#
# DAG
gen_dag = "dag{ 
  SES -> {E_full I};
  E_full -> I;
  {E_full e_C} -> E_cens;
  
  E_full[latent];
  e_C[latent];
}"

dagME = dagitty( gen_dag )

coordinates( dagME ) = list(
  x=c(e_C=-0.6, E_cens=-0.4, E_full=-0.2, SES=0, I=0),
  y=c(e_C=0, E_cens=0, E_full=0, SES=-0.2, I=0) )

drawdag(dagME)



# simulation
# n = simulation sample size
# bSE, bEI, bSI, xT, yT = simulation parameters (last two are thresholds)
# rep = to use in replication
#
f_sim = function(n=100, bSE=1, bEI=1, bSI=0.3, xT=c(-3,3), yT=c(-3,3), rep=F){
  
  # # test
  # n = 100; bSE=1; bEI=1; bSI=0.3; xT=c(-2,2); yT=c(-2,2); rep=F
  
  # sim
  SES = rnorm( n )
  E = rnorm( n , bSE*SES)
  I = rnorm( n , bEI*E + bSI*SES)
  d = data.frame( SES=SES, E_full=E, E_cens=E, I_full=I, I_cens=I)
  
  if( !is.null(xT) ){
    d$E_cens[E <= xT[1]] = xT[1]
    d$E_cens[E >= xT[2]] = xT[2]
  }
  
  if( !is.null(yT) ){
    d$I_cens[I <= yT[1]] = yT[1]
    d$I_cens[I >= yT[2]] = yT[2]  
  }
  
  
  # plot
  if(!rep){
    return(d)
  } else{
    res = summary(lm(I_full ~ E_full + SES, data=d))
    idx = str_detect( rownames(res$coefficients), 'E_full')
    b1 = res$coefficients[idx, c('Estimate','Std. Error')]
    
    res = summary(lm(I_cens ~ E_cens + SES, data=d))
    idx = str_detect( rownames(res$coefficients), 'E_cens')
    b2 = res$coefficients[idx, c('Estimate','Std. Error')]
    b = c(b1, b2)
    names(b) = c('bEf','sEf','bEt','sEt')
    return(b)
  }
  
}



# what's going on
out_range = list(a=c(-5,5), b=c(-1,5), c=c(1,5))


# pdf('descendant8_cens1.pdf', width=14, height=7)
par(mfrow=c(1,3))
# i=1
for(i in 1:length(out_range)){
  set.seed(123456)
  d = f_sim(100, bSE=1, bEI=1, bSI=0.3, xT=out_range[[i]], yT=c(-10,10), rep=F)
  
  f_plot4(d=d, var='E', var_range =out_range[[i]] )
  if(i==1){
    legend('topleft', legend=c('true data','censored data'), pch=c(19, 1), 
           bty='n', col=c(col.alpha('black', 0.3), 'red'))
  }
}
par(mfrow=c(1,1))
# dev.off()




# pdf('descendant8_cens2.pdf', width=14, height=7)
par(mfrow=c(1,3))
# i=3
for(i in 1:length(out_range)){
  set.seed(123456)
  d = f_sim(100, bSE=1, bEI=1, bSI=0.3, xT=out_range[[i]], yT=c(-10,10), rep=F)
  
  f_plot5(dsim=d, var='E', b=c(30,20,10)[i], var_range=out_range[[i]])
  if(i==1){
    legend('topleft', legend=c('true data','censored data'), pch=c(19, 1), 
           bty='n', col=c(col.alpha('black', 0.3), 'red'))
  }
}
par(mfrow=c(1,1))
# dev.off()






d = f_sim(100, bSE=1, bEI=1, bSI=0.3, xT=c(0,6), yT=c(-10,10), rep=F)


# pdf('descendant8_panel.pdf')
psych::pairs.panels(d[,c('SES','E_cens','I_cens')])
# dev.off()


# models
summary( lm(I_full ~ E_full + SES, data=d) )
summary( lm(I_cens ~ E_cens + SES, data=d) )



# sampling variation
out_range = list(a=c(-5,5), b=c(-1,5), c=c(1,5))

# pdf('descendant8_samplesize.pdf')
# i=1
par(mfrow=c(3,2))
for(i in 1:length(out_range)){
  # dsim = replicate( 1e4, f_sim(20, bSE=1, bEI=1, bSI=0.3, xT=c(-10,10), yT=c(-5,0), rep=T) )
  # f_plot1(dsim=dsim, ipar='bM', n=20, xR=c(-2,2), by=0.2,
  #         leg=T, legend=c('true','observed'))
  # f_plot1(dsim=dsim, ipar='sM', n=20, xR=c(0,0.5), by=0.1, leg=F)
  
  dsim = replicate( 1e4, f_sim(100, bSE=1, bEI=1, bSI=0.3, xT=out_range[[i]], yT=c(-10,10), rep=T) )
  f_plot1(dsim=dsim, ipar='bE', n=100, xR=c(-2,2), by=0.2, 
          leg=T, legend=c('true', paste0('E in (', out_range[[i]][1], ', ', out_range[[i]][2], ')')), loc='topleft')
  f_plot1(dsim=dsim, ipar='sE', n=100, xR=c(0,0.5), by=0.1, leg=F)
}
par(mfrow=c(1,1))
# dev.off()






# TO DO ####
# bayesian fix
idx_lc = d$E_cens <= min(d$E_cens)
idx_nc = !idx_lc


dlist = list(
  N = nrow(d),
  #E_lc = d$E_cens[idx_lc],
  SES_lc = d$SES[idx_lc],
  I_lc = d$I_cens[idx_lc],
  E_nc = d$E_cens[idx_nc],
  SES_nc = d$SES[idx_nc],
  I_nc = d$I_cens[idx_nc]
)

cens_model = ulam(
  alist(
    I_nc ~ dnorm( mu_nc , sigma ),
    I_lc ~ dnorm( mu_lc , sigma ),
    mu_nc <- a + bE*E_nc[i] + bS*SES_nc,
    mu_lc <- a + bE*E_lc[i] + bS*SES_lc,
    E_lc ~ custom( normal_lcdf( E_lc | mu_rc , sigma) ),
    a ~ dnorm(0,0.2),
    bE ~ dnorm(0,0.5),
    bS ~ dnorm(0,0.5),
    sigma ~ dexp(1)) , 
  data=dlist , chains=4 , cores=4 )
precis( cens_model , depth=1 )

# stancode(cens_model)
# # what are you doing









# (descendant: censoring, both) ####
#
# Location: https://sites.google.com/view/robertostling/home/teaching
#
# also known as:
#   - special case of measuring error
#
# data details (to make a point): 
# SES = socio-economical status
#   SES -> If: positive (more SES, more I)
#   SES -> E: positive (more SES, more E)
# E = year of education (standardized)
#   E -> If: positive (more E, more If)
# If = income (fully observed)
#   If -> It: positive (one is an observation of the other)
# It = truncated observation of If
# 
# Hypothesis:
# How E explains I?
#
# DAG
gen_dag = "dag{ 
  SES -> {E_full I_full};
  E_full -> I_full;
  {E_full e_C1} -> E_cens;
  {I_full e_C2} -> I_cens;
  
  I_full[latent];
  e_C1[latent];
  E_full[latent];
  e_C2[latent];
}"

dagME = dagitty( gen_dag )

coordinates( dagME ) = list(
  x=c(e_C1=-0.6, E_cens=-0.4, E_full=-0.2, SES=0, I_full=0, I_cens=0.2, e_C2=0.4),
  y=c(e_C1=0, E_cens=0, E_full=0, SES=-0.2, I_full=0, I_cens=0, e_C2=0) )

drawdag(dagME)



# simulation
# n = simulation sample size
# bSE, bEI, bSI, xT, yT = simulation parameters (last two are thresholds)
# rep = to use in replication
#
f_sim = function(n=100, bSE=1, bEI=1, bSI=0.3, xT=c(-3,3), yT=c(-3,3), rep=F){
  
  # # test
  # n = 100; bSE=1; bEI=1; bSI=0.3; xT=c(-2,2); yT=c(-2,2); rep=F
  
  # sim
  SES = rnorm( n )
  E = rnorm( n , bSE*SES)
  I = rnorm( n , bEI*E + bSI*SES)
  d = data.frame( SES=SES, E_full=E, E_cens=E, I_full=I, I_cens=I)
  
  if( !is.null(xT) ){
    d$E_cens[E <= xT[1]] = xT[1]
    d$E_cens[E >= xT[2]] = xT[2]
  }
  
  if( !is.null(yT) ){
    d$I_cens[I <= yT[1]] = yT[1]
    d$I_cens[I >= yT[2]] = yT[2]  
  }
  
  
  # plot
  if(!rep){
    return(d)
  } else{
    res = summary(lm(I_full ~ E_full + SES, data=d))
    idx = str_detect( rownames(res$coefficients), 'E_full')
    b1 = res$coefficients[idx, c('Estimate','Std. Error')]
    
    res = summary(lm(I_cens ~ E_cens + SES, data=d))
    idx = str_detect( rownames(res$coefficients), 'E_cens')
    b2 = res$coefficients[idx, c('Estimate','Std. Error')]
    b = c(b1, b2)
    names(b) = c('bEf','sEf','bEt','sEt')
    return(b)
  }
  
}



# what's going on
out_range = list(a=c(-3,3), b=c(-2,2), c=c(-1,1))


# pdf('descendant9_cens1.pdf', width=14, height=7)
par(mfrow=c(1,3))
# i=1
for(i in 1:length(out_range)){
  set.seed(123456)
  d = f_sim(100, bSE=1, bEI=1, bSI=0.3, xT=out_range[[i]], yT=out_range[[i]], rep=F)
  
  f_plot4(d=d, var='all', var_range=out_range[[i]] )
  if(i==1){
    legend('topleft', legend=c('true data','censored data'), pch=c(19, 1), 
           bty='n', col=c(col.alpha('black', 0.3), 'red'))
  }
}
par(mfrow=c(1,1))
# dev.off()




# pdf('descendant9_cens2.pdf', width=14, height=7)
par(mfrow=c(2,3))
# i=3
for(i in 1:length(out_range)){
  set.seed(123456)
  d = f_sim(100, bSE=1, bEI=1, bSI=0.3, xT=out_range[[i]], yT=out_range[[i]], rep=F)
  
  f_plot5(dsim=d, var='I', b=c(20,10,5)[i], var_range=out_range[[i]])
  if(i==1){
    legend('topleft', legend=c('true data','censored data'), pch=c(19, 1), 
           bty='n', col=c(col.alpha('black', 0.3), 'red'))
  }
  
}
for(i in 1:length(out_range)){
  set.seed(123456)
  d = f_sim(100, bSE=1, bEI=1, bSI=0.3, xT=out_range[[i]], yT=out_range[[i]], rep=F)
  
  f_plot5(dsim=d, var='E', b=c(30,20,10)[i], var_range=out_range[[i]])
  if(i==1){
    legend('topleft', legend=c('true data','censored data'), pch=c(19, 1), 
           bty='n', col=c(col.alpha('black', 0.3), 'red'))
  }
}
par(mfrow=c(1,1))
# dev.off()






d = f_sim(100, bSE=1, bEI=1, bSI=0.3, xT=c(0,6), yT=c(-5,0), rep=F)


# pdf('descendant9_panel.pdf')
psych::pairs.panels(d[,c('SES','E_cens','I_cens')])
# dev.off()


# models
summary( lm(I_full ~ E_full + SES, data=d) )
summary( lm(I_cens ~ E_cens + SES, data=d) )



# sampling variation
out_range = list(a=c(-3,3), b=c(-2,2), c=c(-1,1))

# pdf('descendant9_samplesize.pdf')
# i=1
par(mfrow=c(3,2))
for(i in 1:length(out_range)){
  # dsim = replicate( 1e4, f_sim(20, bSE=1, bEI=1, bSI=0.3, xT=c(-10,10), yT=c(-5,0), rep=T) )
  # f_plot1(dsim=dsim, ipar='bM', n=20, xR=c(-2,2), by=0.2,
  #         leg=T, legend=c('true','observed'))
  # f_plot1(dsim=dsim, ipar='sM', n=20, xR=c(0,0.5), by=0.1, leg=F)
  
  dsim = replicate( 1e4, f_sim(100, bSE=1, bEI=1, bSI=0.3, xT=out_range[[i]], yT=out_range[[4-i]], rep=T) )
  f_plot1(dsim=dsim, ipar='bE', n=100, xR=c(-2,2), by=0.2, 
          leg=T, loc='topleft', 
          legend=c('true', 
                   paste0('E in (', out_range[[i]][1], ', ', out_range[[i]][2], ')',
                   ' and I in (', out_range[[4-i]][1], ', ', out_range[[4-i]][2], ')')))
  f_plot1(dsim=dsim, ipar='sE', n=100, xR=c(0,0.5), by=0.1, leg=F)
}
par(mfrow=c(1,1))
# dev.off()

















# (descendant: miss, MCAR) ####
#
# Location: chapter 15 (p. 500)
#
# also known as:
#   - extreme case of measurement error
#
# data details: 
# S = level of studying
#   S -> H: positive (more S, more value in H)
# D = dog either eats homework (D=1) or not (D=0)
#   D -> H_m: positive (D=1, H_m missing)
# H = complete homework standardized score (NO missing)
#   H -> H_m: positive (one is an observation of the other)
# H_m = observed homework standardized score (SOME missing)
# 
# Hypothesis:
# How the missingness affect the estimates?
#
# DAG
gen_dag = "dag{ 
  S -> H;
  {H D} -> H_m;
  H[latent];
}"
dagMCAR = dagitty( gen_dag )
coordinates( dagMCAR ) = list(
  x=c( S=-0.2, D=-0.2, H=0, H_m=0),
  y=c( S=-0.2, D=0, H=-0.2, H_m=0) )
drawdag(dagMCAR)




# simulation
# n = simulation sample size
# bS = parameter of simulation
# rep = to use in replication
#
f_sim = function(n=100, bS=1, pMiss=0.1, rep=F){
  
  # # test
  # n = 100; bS=1; pMiss=0.1; rep=F
  
  # sim
  S = rnorm( n )
  H = rnorm( n , bS*S )
  D = rbern( n , pMiss)
  Hm = H
  Hm[D==1] = NA
  d = data.frame(S=S, H=H, Hm=Hm)
  
  if(!rep){
    # full data
    return(d)
    
  } else {
    # parameters
    b1 = coef( lm(H ~ S, data=d) )['S']
    b2 = coef( lm(H ~ S, data=d[!is.na(d$Hm),]) )['S']
    b = c(b1, b2)
    names(b) = c('St','So')
    return( b )
  }
  
}




# TO DO ####
# what's going on
p = c(0.1, 0.3, 0.6)


# pdf('descendant10_cens1.pdf', width=14, height=7)
par(mfrow=c(1,3))
# i=1
for(i in 1:length(out_range)){
  set.seed(123456)
  d = f_sim(100, bSE=1, bEI=1, bSI=0.3, xT=c(-10,10), yT=out_range[[i]], rep=F)
  
  f_plot4(dsim=d, var='I', var_range =out_range[[i]] )
  if(i==1){
    legend('topleft', legend=c('true data','censored data'), pch=c(19, 1), 
           bty='n', col=c(col.alpha('black', 0.3), 'red'))
  }
}
par(mfrow=c(1,1))
# dev.off()




# pdf('descendant10_cens2.pdf', width=14, height=7)
par(mfrow=c(1,3))
# i=1
for(i in 1:length(out_range)){
  set.seed(123456)
  d = f_sim(100, bSE=1, bEI=1, bSI=0.3, xT=c(-10,10), yT=out_range[[i]], rep=F)
  
  f_plot5(dsim=d, var='I', b=ifelse(i==1, 30, 10), var_range=out_range[[i]])
  if(i==1){
    legend('topleft', legend=c('true data','censored data'), pch=c(19, 1), 
           bty='n', col=c(col.alpha('black', 0.3), 'red'))
  }
}
par(mfrow=c(1,1))
# dev.off()





d = f_sim(100, bSE=1, bEI=1, bSI=0.3, xT=c(-10,10), yT=c(-5,0), rep=F)


# pdf('descendant10_panel.pdf')
psych::pairs.panels(d[,c('SES','E_cens','I_cens')])
# dev.off()


# models
summary( lm(I_full ~ E_full + SES, data=d) )
summary( lm(I_cens ~ E_cens + SES, data=d) )



# sampling variation
out_range = list(a=c(-5,5), b=c(-5,1), c=c(-5,-1))

# pdf('descendant10_samplesize.pdf')
# i=1
par(mfrow=c(3,2))
for(i in 1:length(out_range)){
  # dsim = replicate( 1e4, f_sim(20, bSE=1, bEI=1, bSI=0.3, xT=c(-10,10), yT=c(-5,0), rep=T) )
  # f_plot1(dsim=dsim, ipar='bM', n=20, xR=c(-2,2), by=0.2,
  #         leg=T, legend=c('true','observed'))
  # f_plot1(dsim=dsim, ipar='sM', n=20, xR=c(0,0.5), by=0.1, leg=F)
  
  dsim = replicate( 1e4, f_sim(100, bSE=1, bEI=1, bSI=0.3, xT=c(-10,10), yT=out_range[[i]], rep=T) )
  f_plot1(dsim=dsim, ipar='bE', n=100, xR=c(-2,2), by=0.2, 
          leg=T, legend=c('true', paste0('I in (', out_range[[i]][1], ', ', out_range[[i]][2], ')')), loc='topleft')
  f_plot1(dsim=dsim, ipar='sE', n=100, xR=c(0,0.5), by=0.1, leg=F)
}
par(mfrow=c(1,1))
# dev.off()








par(mfrow=c(3,2))
for(i in c(0.1,0.2,0.4,0.8,1.6,2)){
  d = f_sim(n=100, bS=1, seed=NULL, rep=F)
  f_plot3(dsim=d, sX)
}
par(mfrow=c(1,1))
# not so pervasive



# sampling variability
par(mfrow=c(3,2))
for(i in c(0.1,0.2,0.4,0.8,1.6,2)){
  
  bAD = -1
  dsim = replicate( 1e4, f_sim(n=100, bAD=bAD, sD=i, plot=F) )
  
  dens( dsim, lwd=3, xaxt='n', col=2,
        xlab="posterior mean", xlim=c(-2,0) )
  axis(side=1, at=seq(-2, 0, by=0.2))
  mtext( paste0('sD: ', i), 3, adj=0, cex=2, at=-2)
  mtext( paste0('mean slope: ', round(mean(dsim), 3)), 3, adj=0, cex=2, at=-0.6)
  abline( v=mean(dsim), lty=2, lwd=3, col=2)
  abline( v=bAD, lty=2, lwd=3)
  
}
par(mfrow=c(1,1))



# relationships
d = f_sim(n=100, bS=1, seed=NULL, rep=F)
psych::pairs.panels(d)
# cor(X,Y)~0, when not controlled by Z (already good)


# model
model_res = glm(cbind(H, 10) ~ -1 + S, data=d, family='binomial')
coef(model_res) # 'true' log odd effect (not known)
inv_logit(coef(model_res)) # 'true' probability

model_res = glm(cbind(Hm, 10) ~ -1 + S, data=d, family='binomial')
coef(model_res) # unbiased log odd effect, less efficient
inv_logit(coef(model_res)) # unbiased probability


# sampling variation
par(mfrow=c(2,2))
dsim = replicate( 1e4, f_sim(n=20, bS=1, seed=NULL, rep=T) )
f_plot1(dsim=dsim, ipar='SR', xR=c(-0.3,1), by=0.2)
f_plot1(dsim=dsim, ipar='SP', xR=c(0.4,0.8), by=0.1)

dsim = replicate( 1e4, f_sim(n=100, bS=1, seed=NULL, rep=T) )
f_plot1(dsim=dsim, ipar='SR', xR=c(-0.3,1), by=0.2)
f_plot1(dsim=dsim, ipar='SP', xR=c(0.4,0.8), by=0.1)
par(mfrow=c(1,1))
# S -> H unbiased when n=20 
# we only loose precision 










# (descendant: miss, MAR easy) ####
#
# Location: chapter 15 (p. 501)
#
# also known as:
#   - extreme case of measurement error
#
# data details: 
# S = level of studying
#   S -> D: negative (more S, less prob. D)
#   S -> H: positive (more S, more H)
# D = dog either eats homework (D=1) or not (D=0)
#   D -> H_m: positive (D=1, H_m missing)
# H = complete homework score 0-10 (NO missing)
#   H -> H_m: positive (one is an observation of the other)
# H_m = observed homework score 0-10 (SOME missing)
# 
# Hypothesis:
# How the missingness affect the estimates?
#
# DAG
gen_dag = "dag{ 
  S -> {H D};
  {H D} -> H_m;
  H[latent];
}"
dagMCAR = dagitty( gen_dag )
coordinates( dagMCAR ) = list(
  x=c( S=-0.2, D=-0.2, H=0, H_m=0),
  y=c( S=-0.2, D=0, H=-0.2, H_m=0) )
drawdag(dagMCAR)


# simulation
# n = simulation sample size
# pS = parameter of simulation
# rep = to use in replication
#
f_sim = function(n=100, bS=1, seed=NULL, rep=F){
  
  # # test
  # n = 100; pS=1, seed=NULL; rep=F
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  # sim
  S = rnorm( n )
  H = rbinom( n , size=10 , inv_logit(bS*S) )
  D = ifelse( S < 0 , 1 , 0 ) # 1 not observed
  Hm = H
  Hm[D==1] = NA
  d = data.frame(S=S, H=H, Hm=Hm)
  
  if(!rep){
    # full data
    return(d)
    
  } else {
    # parameters
    b1 = coef( glm(cbind(H, 10) ~ -1 + S, data=d, family='binomial') )['S'] # unbiased effect
    b1 = c( b1, inv_logit(b1) )
    b2 = coef( glm(cbind(Hm, 10) ~ -1 + S, data=d[!is.na(d$Hm),], family='binomial') )['S'] # biased effects
    b2 = c( b2, inv_logit(b2) )
    b = c(b1, b2)
    names(b) = c('SR','SP','SRb','SPb')
    return( b )
  }

}

# relationships
d = f_sim(n=100, bS=1, seed=NULL, rep=F)
psych::pairs.panels(d)
# cor(X,Y)~0, when not controlled by Z (already good)


# model
d = f_sim(n=100, bS=1, seed=NULL, rep=F)

m_res = glm(cbind(H, 10) ~ -1 + S, data=d, family='binomial')
# summary(m_res) # 'true' log odd effect (not known)
coef(m_res) # 'true' odds
inv_logit(coef(m_res)) # 'true' probability

m_res = glm(cbind(Hm, 10) ~ -1 + S, data=d[!is.na(d$Hm),], family='binomial')
# summary(m_res) # biased log odd effect, less efficient
coef(m_res) # biased odds
inv_logit(coef(m_res)) # biased probability



# sampling variation
par(mfrow=c(2,2))
dsim = replicate( 1e4, f_sim(n=100, pS=1, seed=NULL, rep=T) )
f_plot1(dsim=dsim, ipar='SR', xR=c(-0.4,0.6), by=0.1)
f_plot1(dsim=dsim, ipar='SP', xR=c(0.40,0.65), by=0.02)

dsim = replicate( 1e4, f_sim(n=200, pS=1, seed=NULL, rep=T) )
f_plot1(dsim=dsim, ipar='SR', xR=c(-0.4,0.6), by=0.1)
f_plot1(dsim=dsim, ipar='SP', xR=c(0.40,0.65), by=0.02)
par(mfrow=c(1,1))
# S -> H is smaller in observed data 
# equally biased with n=100, but more "confident" of S -> H







# (descendant: miss, imputation) ####
#
# Location: slides lecture 10, 2022 course
#
# also known as:
#   - 
#
# data details: 
#   same as previous
# 
# Hypothesis:
# How can we minimize the effect of missingness in the estimates?
#
f_imp = function(d, K=10, to_imp='Hm', ipar='S'){
  
  # # test
  # d = f_sim(n=100, pS=1, seed=NULL, rep=F)
  # K=10; to_imp='Hm'; ipar='S'
  
  # imputation model
  m_res = glm(cbind(Hm, 10) ~ 1, data=d[!is.na(d$Hm),], family='binomial')
  pH = inv_logit(coef(m_res)) # probability of observing H
  # pH
  
  # imputation
  coef_imp = data.frame(matrix(NA, nrow=K, ncol=2)) # saving parameters
  names(coef_imp) = c(paste0('b',ipar),'se')
  
  # k=1
  for(k in 1:K){
    nas = is.na(d$Hm)
    
    # imputation process
    d2 = d
    d2[nas, to_imp] = rbinom( sum(nas), size=10 , p=pH )
    
    # fitting imputed model
    m_res = glm(cbind(Hm, 10) ~ -1 + S, data=d2, family='binomial')
    m_res = summary(m_res)
    
    # saving parameters
    coef_imp[k, ] = m_res$coefficients[ipar, 1:2]
    
  }
  
  # combining parameters
  b = mean(coef_imp[,1]) # imputed parameter
  
  # combining SEs (to get CIs)
  v_W = mean(coef_imp[,2]^2)
  v_B = (K-1)^(-1) * sum( (coef_imp[,1] - b)^2 )
  v_T = v_W + (1 + K^(-1)) * v_B
  se_T = sqrt(v_T)
  
  # return object
  return(c(b, se_T))
  
}


# generate data (with missing values)
d = f_sim(n=100, pS=1, seed=NULL, rep=F)

m_res = glm(cbind(H, 10) ~ -1 + S, data=d, family='binomial')
# summary(m_res) # 'true' log odd effect (not known)
coef(m_res) # 'true' odds
inv_logit(coef(m_res)) # 'true' probability


m_res = glm(cbind(Hm, 10) ~ -1 + S, data=d[!is.na(d$Hm),], family='binomial')
# summary(m_res) # biased log odd effect, less efficient
coef(m_res) # biased odds
inv_logit(coef(m_res)) # biased probability


m_res = f_imp(d=d, K=10, to_imp='Hm', ipar='S')
m_res[1]
inv_logit(m_res[1])
# way less biased effects




# sampling variation
f_SimImp = function(n=100, pS=1, ipar='S'){
  
  # # test
  # n=100; pS=1; ipar='S'
  
  # generate data
  d = f_sim(n=n, pS=pS, seed=NULL, rep=F)

  # values
  b1 = coef( glm(cbind(H, 10) ~ -1 + S, data=d, family='binomial') )[ipar]
  b1 = c(b1, inv_logit(b1) ) # true
  b2 = coef( glm(cbind(Hm, 10) ~ -1 + S, data=d[!is.na(d$Hm),], family='binomial') )[ipar]
  b2 = c(b2, inv_logit(b2) ) # biased 
  b3 = f_imp(d=d, K=10, to_imp='Hm', ipar='S')[1]
  b3 = c(b3, inv_logit(b3) ) # imputed
  
  b = c(b1, b2, b3)
  names(b) = c('SR','SP','SRb','SPb','SRi','SPi')
  
  return(b)
}


par(mfrow=c(2,2))
dsim = replicate( 1000, f_SimImp(n=100, pS=1, ipar='S') )
f_plot1(dsim=dsim, ipar='SR', xR=c(-0.4,0.6), by=0.1)
f_plot1(dsim=dsim, ipar='SP', xR=c(0.4,0.65), by=0.02)

dsim = replicate( 1000, f_SimImp(n=200, pS=1, ipar='S') )
f_plot1(dsim=dsim, ipar='SR', xR=c(-0.4,0.6), by=0.1)
f_plot1(dsim=dsim, ipar='SP', xR=c(0.4,0.65), by=0.02)
par(mfrow=c(1,1))
# S -> H is less biased with imputation 
# equally biased with n=100, but more closer to true S -> H














# (descendant: miss, MAR hard) ####
#
# Location: chapter 15 (p. 501)
#
# also known as:
#   - extreme case of measurement error
#
# data details: 
# S = level of studying
#   S -> H: positive (more S, more H)
# X = noise level in child's house
#   X -> D: negative (more S, less prob. D)
#   X -> H: negative (more X, less S)
# D = dog either eats homework (D=1) or not (D=0)
#   D -> H_m: positive (D=1, H_m missing)
# H = complete homework score 0-10 (NO missing)
#   H -> H_m: positive (one is an observation of the other)
# H_m = observed homework score 0-10 (SOME missing)
# 
# Hypothesis:
# How the missingness affect the estimates?
#
gen_dag = "dag{ 
  S -> {H};
  X -> {D H};
  {H D} -> H_m;
  X [latent]
  H[latent];
}"
dagMCAR = dagitty( gen_dag )
coordinates( dagMCAR ) = list(
  x=c( S=-0.2, D=-0.2, X=-0.1, H=0, H_m=0),
  y=c( S=-0.2, D=0, X=-0.1, H=-0.2, H_m=0) )
drawdag(dagMCAR)


# data
N <- 1000
set.seed(501)
X <- rnorm(N)
S <- rnorm(N)
H <- rbinom( N , size=10 , inv_logit( 2 + S - 2*X ) )
D <- ifelse( X > 1 , 1 , 0 )
Hm <- H
Hm[D==1] <- NA

d = data.frame(S=S[D==0], Hm=Hm[D==0]) # only observed data
# nrow(d)


# model
model_res = glm(cbind(H, 10) ~ -1 + S, family='binomial')
coef(model_res) # 'true' log odd effect (not known)
exp(coef(model_res)) # 'true' odds
inv_logit(coef(model_res)) # 'true' probability

model_res = glm(cbind(Hm, 10) ~ -1 + S, data=d, family='binomial')
coef(model_res) # biased log odd effect, less efficient
exp(coef(model_res)) # biased odds
inv_logit(coef(model_res)) # biased probability



# solution 1 (if you have X)
d = data.frame(S=S[D==0], Hm=Hm[D==0], X=X[D==0]) # only observed data
model_res = glm(cbind(Hm, 10) ~ -1 + S + X, data=d, family='binomial')
coef(model_res)[1] # less biased log odd effect, less efficient
exp(coef(model_res))[1] # less biased odds
inv_logit(coef(model_res))[1] # less biased probability



# solution 2 (if you do not have X, not quite good)
# imputation model
model_res = glm(cbind(Hm, 10) ~ 1, family='binomial')
pH = inv_logit(coef(model_res)) # probability of observing H
pH

d = data.frame(S=S, Hm=Hm) # full data
d2 = d # copy for modification

K = 10 # number of imputations
coef_imp = data.frame(matrix(NA, nrow=K, ncol=2)) # saving parameters
names(coef_imp) = c('bS','se')

for(k in 1:10){
  nas = is.na(d$Hm)
  
  # imputation process
  d2$Hm[nas] = rbinom( sum(nas) , size=10 , p=pH )
  
  # fitting imputed model
  model_res = glm(cbind(Hm, 10) ~ -1 + S, data=d2, family='binomial')
  model_res = summary(model_res)
  
  # saving parameters
  coef_imp[k, ] = model_res$coefficients[1:2]
  
}

# combining parameters
bS = mean(coef_imp$bS) # imputed parameter
bS # less biased (than just use observed data)
exp(bS) # less biased odds
inv_logit(bS) # less biased probability






# (descendant: miss, MNAR) ####
#
# Location: chapter 15 (p. 503)
#
# also known as:
#   - extreme case of measurement error
#
# data details: 
# S = level of studying
#   S -> H: positive (more S, more H)
# D = dog either eats homework (D=1) or not (D=0)
#   D -> H_m: positive (D=1, H_m missing)
# H = complete homework score 0-10 (NO missing)
#   H -> D: negative (more H, less prob. D)
#   H -> H_m: positive (one is an observation of the other)
# H_m = observed homework score 0-10 (SOME missing)
# 
# Hypothesis:
# How the missingness affect the estimates?
#
# DAG
gen_dag = "dag{ 
  S -> H;
  H -> {H_m D};
  D -> H_m;
  H[latent];
}"
dagMCAR = dagitty( gen_dag )
coordinates( dagMCAR ) = list(
  x=c( S=-0.2, D=-0.2, H=0.2, H_m=0.2),
  y=c( S=-0.2, D=0.2, H=-0.2, H_m=0.2) )
drawdag(dagMCAR)


# simulations
N = 1000
set.seed(501)
S = rnorm(N)
H = rbinom( N , size=10 , inv_logit(S) )
D = ifelse( H < 5 , 1 , 0 )
Hm <- H
Hm[D==1] <- NA

d = data.frame(S=S[D==0], Hm=Hm[D==0]) # incomplete data
# nrow(d)


# model
model_res = glm(cbind(H, 10) ~ -1 + S, family='binomial')
coef(model_res) # 'true' log odd effect (not known)
exp(coef(model_res)) # 'true' odds
inv_logit(coef(model_res)) # 'true' probability

model_res = glm(cbind(Hm, 10) ~ -1 + S, data=d, family='binomial')
coef(model_res) # completely biased log odd effect, less efficient
exp(coef(model_res)) # completely biased odds
inv_logit(coef(model_res)) # completely biased probability


# solution
# None, you are f...
# but you can do sensitivity analysis with imputation









# (descendant: truncation, outcome only) ####
#
# Location: https://sites.google.com/view/robertostling/home/teaching
#
# also known as:
#   - special case of missing data
#
# data details (to make a point): 
# SES = socio-economical status
#   SES -> If: positive (more SES, more I)
#   SES -> E: positive (more SES, more E)
# E = year of education (standardized)
#   E -> If: positive (more E, more If)
# If = income (fully observed)
#   If -> It: positive (one is an observation of the other)
# It = truncated observation of If
# 
# Hypothesis:
# How E explains I?
#
# DAG
gen_dag = "dag{ 
  SES -> {E I_full};
  E -> I_full;
  {I_full e_T} -> I_trunc;
  
  I_full[latent];
  e_T[latent];
}"

dagME = dagitty( gen_dag )

coordinates( dagME ) = list(
  x=c(E=-0.2, SES=0, I_full=0, I_trunc=0.2, e_T=0.4),
  y=c(E=0, SES=-0.2, I_full=0, I_trunc=0, e_T=0) )

drawdag(dagME)



# simulation
# n = simulation sample size
# bXY, xT, yT = simulation parameters
# rep = to use in replication
#
f_sim = function(n=100, bSE=1, bEI=1, bSI=0.3, xT=c(-3,3), yT=c(-3,3), rep=F){
  
  # # test
  # n = 100; bSE=1; bEI=1; bSI=0.3; xT=c(-2,2); yT=c(-2,2); rep=F
  
  # sim
  SES = rnorm( n )
  E = rnorm( n , bSE*SES)
  I = rnorm( n , bEI*E + bSI*SES)
  d = data.frame( SES=SES, E_full=E, E_trunc=E, I_full=I, I_trunc=I)
  
  if( !is.null(xT) ){
    keep = (E >= xT[1]) & (E <= xT[2])
    d$E_trunc[!keep] = NA
  }
  
  if( !is.null(yT) ){
    keep = (I >= yT[1]) & (I <= yT[2])
    d$I_trunc[!keep] = NA  
  }
  
  
  # plot
  if(!rep){
    return(d)
  } else{
    res = summary(lm(I_full ~ E_full + SES, data=d))
    idx = str_detect( rownames(res$coefficients), 'E_full')
    b1 = res$coefficients[idx, c('Estimate','Std. Error')]
    
    res = summary(lm(I_trunc ~ E_trunc + SES, data=d))
    idx = str_detect( rownames(res$coefficients), 'E_trunc')
    b2 = res$coefficients[idx, c('Estimate','Std. Error')]
    b = c(b1, b2)
    names(b) = c('bEf','sEf','bEt','sEt')
    return(b)
  }
  
}



# what's going on
out_range = list(a=c(-5,5), b=c(-5,1), c=c(-5,-1))


# pdf('descendant7_trunc.pdf', width=14, height=7)
par(mfrow=c(1,3))
# i=1
for(i in 1:length(out_range)){
  set.seed(123456)
  d = f_sim(100, bSE=1, bEI=1, bSI=0.3, xT=c(-10,10), yT=out_range[[i]], rep=F)
  
  plot(d$E_full, d$I_full, pch=19, col=col.alpha('black',0.3))
  points(d$E_trunc, d$I_trunc, pch=1, col='red', cex=1.5)
  b = coef( lm(I_trunc ~ E_trunc + SES, data=d) )
  abline( c(b[1], b[2]) )
  abline( h=out_range[[i]][2], lty=2 )
  mtext( paste0('I observed in (', out_range[[i]][1], ', ', out_range[[i]][2], '),  bEI: ', round(b[2], 3)), 
         3, adj=0, cex=1.5, at=min(d$E_full))
  
  if(i==1){
    legend('topleft', legend=c('true data','truncated data'), pch=c(19, 1), 
           bty='n', col=c(col.alpha('black', 0.3), 'red'))
  }
}
par(mfrow=c(1,1))
# dev.off()




d = f_sim(100, bSE=1, bEI=1, bSI=0.3, xT=c(-10,10), yT=c(-5,0), rep=F)


# pdf('descendant7_panel.pdf')
psych::pairs.panels(d[,c('SES','E_trunc','I_trunc')])
# dev.off()


# models
summary( lm(I_full ~ E_full + SES, data=d) )
summary( lm(I_trunc ~ E_trunc + SES, data=d) )



# sampling variation
# pdf('descendant7_samplesize.pdf')
par(mfrow=c(3,2))
# i=0.2
for(i in c(0.2, 1, 2)){
  # dsim = replicate( 1e4, f_sim(20, bSE=1, bEI=1, bSI=0.3, xT=c(-10,10), yT=c(-5,0), rep=T) )
  # f_plot1(dsim=dsim, ipar='bM', n=20, xR=c(-2,2), by=0.2,
  #         leg=T, legend=c('true','observed'))
  # f_plot1(dsim=dsim, ipar='sM', n=20, xR=c(0,0.5), by=0.1, leg=F)
  
  dsim = replicate( 1e4, f_sim(100, bSE=1, bEI=1, bSI=0.3, xT=c(-10,10), yT=c(-5,0), rep=F) )
  f_plot1(dsim=dsim, ipar='bM', n=100, xR=c(-2,2), by=0.2, 
          leg=T, legend=c('true', paste0('sM = ', i)), loc='topleft')
  f_plot1(dsim=dsim, ipar='sM', n=100, xR=c(0,0.5), by=0.1, leg=F)
}
par(mfrow=c(1,1))
# dev.off()



# TO DO ####
# bayesian fix
dlist = list(
  N = nrow(d),
  D = d$D,
  A = d$A,
  M_obs = d$M_obs,
  M_sd = rep(1, nrow(d))
)

me_model = ulam(
  alist(
    D ~ dnorm( mu , sigma ),
    mu <- a + bA*A + bM*M_true[i],
    M_obs ~ dnorm( M_true , M_sd ),
    vector[N]:M_true ~ dnorm( 0 , sigma ),
    a ~ dnorm(0,0.2),
    bA ~ dnorm(0,0.5),
    bM ~ dnorm(0,0.5),
    sigma ~ dexp(1)) , 
  data=dlist , chains=4 , cores=4 )
precis( me_model , depth=1 )











# truncated predictor
gen_dag = "dag {
  X -> Y;
  T -> X;
  T [unobserved]
}"
dag_plot1 = dagitty( gen_dag )
coordinates(dag_plot1) = list( x=c(X=0,T=0,Y=1) , 
                               y=c(X=0,T=-1,Y=0) )
drawdag( dag_plot1 )



# truncating predictor
par(mfrow=c(3,2))
for(i in c(-2,-1,-0.5,0,0.5,1)){
  f_sim(n=1000, bXY=1, xT=i, yT=NULL, rep=F)
}
par(mfrow=c(1,1))
# kind of pervasive












































# TO DO ####
# # (pipe: post-stratification) ####
# #
# # Location: lecture 09, 2022 course
# #
# # also an example of:
# #   - marginalization
# #
# # Data details: 
# # G = gender
# #   G -> A: null (assumed no effect)
# #   G -> D: negative (G=female, specific D's)
# # D = department to be admitted
# #   D -> A: negative (specific D's, less A)
# # A = number of admissions
# # 
# # Hypothesis:
# # what is the marginal effect of G -> A? 
# # 
# # parameter posterior
# b = extract.samples(m_res) 
# b = b[sample(1:nrow(b), 500),] # smaller posterior sample
# bC = data.frame(b[1] - b[2],
#                 b[1] - b[2] + b[3],
#                 b[1] - b[2] + b[4],
#                 b[1] - b[2] + b[5],
#                 b[1] - b[2] + b[6],
#                 b[1] - b[2] + b[7])
# names(bC) = c('A','B','C','D','E','F')
# 
# for(i in 1:6){
#   if(i==1){
#     dens( bC[,i], lwd=2, col=i+1, xlim=c(-3,0.5),
#           xlab="effect of gender perception" )
#   } else{
#     dens( bC[,i], lwd=2, col=i+1, add=T )
#   }
# }
# abline( v=0, lty=2)
# 
# 
# 
# # prediction function
# f_pred = function(beta, d){ # log-odds function 
#   
#   # # test
#   # beta=b; d=dsim; prob=F
#   
#   # storage
#   res = matrix(NA, ncol=nrow(beta), nrow=nrow(d))
#   # dim(res)
#   
#   # i=1
#   for(i in 1:nrow(beta)){
#     D = as.integer(d$D)
#     bG = unlist( beta[i, d$G] )
#     bD = unlist( ifelse( D==1, 0, beta[i,D+1] ) )
#     res[,i] = bG + bD # P(A|G,D), log-odds
#   }
#   
#   res = colMeans(res) # P(A|G) = sum[ P(A|G,D).P(D) ]
#   # notice marginalization is converted into a simple difference
#   # of conditional parameters, based on the law of total probability
#   # P(A) = sum[ P(A|B).P(B) ], where P(B)=1
#   
#   # return object
#   return( data.frame(RE=res, P=inv_logit(res)) )
#   
# }
# 
# # simulation
# total_apps = sum(d$N) # number of applications per department 
# apps_per_dept = sapply( 1:6 , function(i) sum(d$N[d$D==i]) ) 
# 
# 
# # simulate as if all apps from women 
# dsim = data.frame( 
#   D=rep(1:6, times=apps_per_dept), 
#   G=rep(1,total_apps))
# p_G1 = f_pred(beta=b, d=dsim) 
# 
# # simulate as if all apps from men
# dsim = data.frame( 
#   D=rep(1:6, times=apps_per_dept), 
#   G=rep(2,total_apps))
# p_G2 = f_pred(beta=b, d=dsim)
# # notice the post-stratification is done when we create 
# # the data. We can create a data that reflects the 
# # appropriate population sizes or weights.
# 
# 
# 
# par(mfrow=c(1,2))
# # in relative terms
# dens( p_G1[,'RE'] - p_G2[,'RE'] , lwd=4 , col=2 , 
#       xlab="effect of gender perception" )
# abline( v=mean(p_G1[,'RE'] - p_G2[,'RE']), lty=2)
# abline( v=coef(t), lty=2, col=2)
# 
# # in probability terms
# dens( p_G1[,'P'] - p_G2[,'P'] , lwd=4 , col=2 , 
#       xlab="effect of gender perception" )
# abline( v=mean(p_G1[,'P'] - p_G2[,'P']), lty=2)
# abline( v=inv_logit(coef(m_res)[1]) - inv_logit(coef(m_res)[2]), lty=2, col=2)
# par(mfrow=c(1,1))
# # notice both vertical lines are equal because after controlling
# # by D the G effect is no longer biased. However, the correct 
# # form of calculating the marginal effects is the long one
# # NOTICE: by modal distribution







# # (pipe/collider: sensitivity analysis) ####
# #
# # Location: slides lecture 10, 2022 course
# #
# # also know as:
# #   - 
# #
# # Simulation details 1: 
# # G = gender
# #   G -> A: null (assumed no effect)
# #   G -> D: negative (G=female, specific D's)
# # D = department to be admitted
# #   D -> A: negative (specific D's, less A)
# # U = unobserved confound (e.g. ability, qualification)
# #   U -> D: positive (more U, specific D's, the hardest to get in)
# #   U -> A: positive (more U, more A)
# #   U -> {T1, T2, T3}: positive (more U, higher Ti)
# # T1, T2, T3: observed exam scores
# # A = admission {0,1}
# # 
# # Hypothesis:
# # G has no relationship with A? 
# # D can be a confounding variable
# #
# # DAGs
# gen_dag = "dag{
#   {G D} -> A;
#   G -> D
#   U -> {D A};
#   U [unobserved]
# }"
# dag_plot1 = dagitty( gen_dag )
# coordinates(dag_plot1) = list( x=c(G=0,D=0,U=1,A=1) , 
#                                y=c(G=0,D=-1,U=-1,A=0) )
# drawdag( dag_plot1 )
# # clearly you cannot control for U (unobserved)
# # opens an unwanted backdoor path
# 
# 
# 
# # simulation
# # n = simulation sample size
# # pGD, pUD, ar = parameter of simulation
# # rep = to use in replication
# #
# f_sim = function(n=1000, pGD=c(0.3,0.8), pUD=0.5,
#                  ar=list( c(0.1,0.1,0.1,0.3),
#                           c(0.2,0.3,0.2,0.5)), rep=F){
#   
#   # # test
#   # n = 100; pGD=c(0.3,0.8); pUD=0.2; ar=list( c(0.1,0.3,0.1,0.3),
#   #                                            c(0.2,0.5,0.2,0.5))
#   
#   # set.seed(17)
#   G = sample( 1:2, size=n, replace=T ) # even gender distribution 
#   U = rbern(n, 0.1) # ability, high (1) to average (0)
#   # gender 1 tends to apply to department 1, 2 to 2
#   # and G=1 with greater ability tend to apply to 2 as well
#   D = rbern( n, ifelse( G==1, U*pUD, pGD[2] ) ) + 1 # D=2 discriminatory
#   # matrix of acceptance rates [dept,gender]
#   ar[[1]] = matrix( ar[[1]], nrow=2 ) # different acceptance per U
#   ar[[2]] = matrix( ar[[2]], nrow=2 )
#   
#   # simulate acceptance
#   p = sapply( 1:n , function(i){
#     ifelse( U[i]==0, ar[[1]][D[i],G[i]], ar[[2]][D[i],G[i]] ) } )
#   A = rbern( n , p )
#   
#   d = data.frame(U,G=factor(G),D=factor(D),A)
#   
#   
#   # return object
#   if(!rep){
#     # full data
#     return(d)
#     
#   } else{
#     # parameters
#     b1 = coef( glm(A ~ -1 + G + D + U, data=d, family='binomial') )[c('G1','G2')] # unbiased effects
#     b1 = c( b1[2]-b1[1], inv_logit(b1[2])-inv_logit(b1[1]) )
#     b2 = coef( glm(A ~ -1 + G + D, data=d, family='binomial') )[c('G1','G2')] # biased effects
#     b2 = c( b1[2]-b1[1], inv_logit(b1[2])-inv_logit(b1[1]) )
#     b3 = coef( glm(A ~ -1 + G, data=d, family='binomial') )[c('G1','G2')] # biased effect
#     b3 = c( b2[2]-b2[1], inv_logit(b2[2])-inv_logit(b2[1]) )
#     b = c(b1, b2, b3)
#     names(b) = c('GC','GP','GCb1','GPb1','GCb2','GPb2')
#     return( b )
#     
#   }
#   
# }
# 
# 
# # models
# d = f_sim(n=1000, pGD=c(0.3,0.8), pUD=0.5,
#           ar=list( c(0.1,0.1,0.1,0.3), # discrimination
#                    c(0.2,0.3,0.2,0.5)), rep=F) 
# 
# m_res = glm(A ~ -1 + G, data=d, family='binomial')
# # summary(m_res) 
# K = matrix(c(1, -1), 1) # contrast of interest
# t = glht(m_res, linfct = K)
# summary(t) # only gender, biased effect
# inv_logit(coef(m_res)[1]) - inv_logit(coef(m_res)[2]) # probability
# # we can infer discrimination total effect terms (no problem)
# 
# 
# m_res = glm(A ~ -1 + G + D, data=d, family='binomial')
# # summary(m_res) 
# K = matrix(c(1, -1, 0), 1) # contrast of interest
# t = glht(m_res, linfct = K)
# summary(t) # G and D, less biased effect
# inv_logit(coef(m_res)[1]) - inv_logit(coef(m_res)[2]) # probability
# # but controlling by D tells you there is no discrimination
# # but we know there is, what to do?
# # U masks the discrimination
# 
# 
# 
# # sampling variation
# par(mfrow=c(2,2))
# dsim = replicate( 1e4, f_sim(n=200, pGD=c(0.3,0.8), pUD=0.5,
#                              ar=list( c(0.1,0.1,0.1,0.3), # discrimination
#                                       c(0.2,0.3,0.2,0.5)), rep=T) )
# f_plot1(dsim=dsim, ipar='GC', xR=c(-2,2.5), by=0.5)
# f_plot1(dsim=dsim, ipar='GP', xR=c(-0.3,0.4), by=0.05)
# 
# dsim = replicate( 1e4, f_sim(n=500, pGD=c(0.3,0.8), pUD=0.5,
#                              ar=list( c(0.1,0.1,0.1,0.3), # discrimination
#                                       c(0.2,0.3,0.2,0.5)), rep=T) )
# f_plot1(dsim=dsim, ipar='GC', xR=c(-2,2.5), by=0.5)
# f_plot1(dsim=dsim, ipar='GP', xR=c(-0.3,0.4), by=0.05)
# par(mfrow=c(1,1))
# # E -> H (negative), but underestimated
# # equally biased with n=100, but more "confident" of E -> H (underestimated)
# 
# 
# 
# 
# 
# # putting data into list
# dlist = with(d, list(N=nrow(d),G=G,D=D,A=A))
# 
# # model application to D2
# dlist$D2 = ifelse( dlist$D==2 , 1 , 0 ) 
# 
# # similar effect U[G] -> A
# dlist$b = c(1,1) 
# # logic: high ability (U), more prob. of being admitted (A)
# #         equal increment for both G
# 
# # different effect U -> D | G
# dlist$g = c(1,0) 
# # logic: high ability (U=1), more prob. to apply to D=2
# # How large g[1]=1 has to be to change our conclusions?
# 
# # normally in sensitivity analysis you try a range of values
# # for b and g.
# 
# 
# # controlled effect
# m = ulam(
#   alist(
#     # A model
#     A ~ bernoulli(p),
#     logit(p) <- a[G,D],
#     matrix[G,D]:a ~ normal(0,1) ), 
#   data=dlist, chains=4, cores=4 ) 
# 
# 
# # sensitivity effect
# mp = ulam(
#   alist( 
#     # A model
#     A ~ bernoulli(p),
#     logit(p) <- a[G,D] + b[G]*u[i], # assume u (sample)
#     matrix[G,D]:a ~ normal(0,1),
#     
#     # D model
#     D2 ~ bernoulli(q),
#     logit(q) <- delta[G] + g[G]*u[i], # assume u (sample)
#     delta[G] ~ normal(0,1),
#     
#     # declare unobserved u
#     vector[N]:u ~ normal(0,1) ), # distribution of u
#   data=dlist , chains=4 , cores=4 )
# # notice the change on paradigm: 
# # rather than estimating the b and g, we provide them, 
# # and we find a U that fits with the assumptions of b and g.
# 
# 
# # # results
# # precis(m, 3, pars=c("a"))
# # precis(mp, 3, pars=c("a"))
# # # relative effects (hard to eyeball)
# 
# 
# # contrasts
# post_m = extract.samples(m)
# post_m$cont_D1 = post_m$a[,1,1] - post_m$a[,2,1]
# post_m$cont_D2 = post_m$a[,1,2] - post_m$a[,2,2]
# 
# post_mp = extract.samples(mp)
# post_mp$cont_D1 = post_mp$a[,1,1] - post_mp$a[,2,1]
# post_mp$cont_D2 = post_mp$a[,1,2] - post_mp$a[,2,2]
# 
# 
# # plot 
# dens( post_m$cont_D1, lwd=1, col=3, xlim=c(-1.5,1.5), 
#       xlab="F-M contrast in each department" )
# abline(v=mean(post_m$cont_D1),lwd=1, col=3, lty=2)
# dens( post_m$cont_D2, lwd=1, col=4, add=T )
# abline(v=mean(post_m$cont_D2),lwd=1, col=4, lty=2)
# dens( post_mp$cont_D1, lwd=3, col=3, add=T )
# abline(v=mean(post_mp$cont_D1),lwd=3, col=3, lty=2)
# dens( post_mp$cont_D2, lwd=3, col=4, add=T )
# abline(v=mean(post_mp$cont_D2),lwd=3, col=4, lty=2)
# # notice for D=2 the effect is even favorable for women
# # because U mask the discrimination, but we know is shouldn't
# # NOW: if b[U->D2] = 1, the b[G->A] ~ -0.27
# #     i.e. if the effect of U in application and admission
# #     is large, we will be masking a discrimination as 
# #     large as -0.27 (against women)
# 
# 
# with(d, plot( jitter(U), apply(post_mp$u,2,mean), 
#               xlab="u (true)", ylab="posterior mean u", 
#               xaxt="n", lwd=2, 
#               col=ifelse(G==1,col.alpha(2,0.6), col.alpha(4,0.6)) ) )
# axis(1,at=0:1,labels=0:1)
# # we managed to identify the different U by G
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # real world data
# data(UCBadmit)
# d = UCBadmit
# 
# # to long format
# dat_long1 <- uncount( d , admit )
# dat_long1$Y <- 1
# dat_long1$reject <- NULL
# dat_long0 <- uncount( d , reject )
# dat_long0$Y <- 0
# dat_long0$admit <- NULL
# dat_long01 <- rbind( dat_long1 , dat_long0 )
# dat_long01$applications <- NULL
# dat_long01$admit <- dat_long01$Y
# dat_long01$Y <- NULL
# 
# # to list
# dlist <- list(
#   A = dat_long01$admit,
#   G = ifelse(dat_long01$applicant.gender=="female",1,2),
#   D = as.integer(dat_long01$dept) )
# 
# # application to D=1 (stands out from the data/modelling)
# dlist$D1 = ifelse(dlist$D==1,1,0) 
# dlist$N = length(dlist$D)
# 
# # similar effect U[G] -> A
# dlist$b = c(1,1)
# # logic: high ability (U), more prob. of being admitted (A)
# #         equal increment for both G
# 
# # different effect U -> D | G
# dlist$g = c(1,0)
# # logic: high ability (U=1), more prob. to apply to D=2
# # How large g[1]=1 has to be to change our conclusions?
# 
# # controlled effect
# m = ulam(
#   alist(
#     # A model
#     A ~ bernoulli(p),
#     logit(p) <- a[G,D],
#     matrix[G,D]:a ~ normal(0,1) ), 
#   data=dlist, chains=4, cores=4 ) 
# 
# 
# # sensitivity effects
# mp = ulam(
#   alist( 
#     # A model
#     A ~ bernoulli(p),
#     logit(p) <- a[G,D] + b[G]*u[i],
#     matrix[G,D]:a ~ normal(0,1),
#     
#     # D model
#     D1 ~ bernoulli(q),
#     logit(q) <- delta[G] + g[G]*u[i],
#     delta[G] ~ normal(0,1),
#     
#     # declare unobserved u
#     vector[N]:u ~ normal(0,1) ), 
#   data=dlist, chains=4, cores=4 )
# 
# 
# # # results
# # precis(m, 3, pars=c("a"))
# # precis(mp, 3, pars=c("a"))
# # # relative effects (hard to eyeball)
# 
# 
# # contrasts
# post_m <- extract.samples(m)
# post_m$cont_D1 <- post_m$a[,1,1] - post_m$a[,2,1]
# 
# post_mp <- extract.samples(mp)
# post_mp$cont_D1 <- post_mp$a[,1,1] - post_mp$a[,2,1]
# 
# 
# # plot 
# dens( post_m$cont_D1, lwd=1, col=2, xlim=c(-1,2), 
#       xlab="F-M contrast in each department" )
# abline(v=0,lwd=1, col=1, lty=2)
# abline(v=mean(post_m$cont_D1),lwd=1, col=2, lty=2)
# dens( post_mp$cont_D1, lwd=3, col=2, add=T )
# abline(v=mean(post_mp$cont_D1),lwd=3, col=2, lty=2)
# # notice for D=1 the effect is favorable for women
# # because U might be masking the effects
# # NOW: if b[U->D1] = 1, the b[G->A] ~ 0.28
# #     i.e. if the effect of U in application and admission
# #     is large, we will be biasing an effect (0.98) 
# #     as small as 0.27 (in favor of women)
# 
# 
# 
# 
# 
# 
# 
# 
# # Simulation details 2: 
# #
# # Location: chapter 6 (p. 180)
# #
# # U = unobserved variable (e.g. neighborhood)
# #   U -> P: positive (different U's, more P)
# #   U -> C: positive (different U's, more C)
# # G = grandparent's educational level
# #   G -> P: positive (more G, more P)
# #   G -> C: null (to emphasize the problem)
# # P = parent's educational level
# #   P -> C: positive (more P, more C)
# # C = child's educational achievement
# # 
# # Hypothesis:
# # G and P impact positively on C?
# #
# # DAG
# gen_dag = "dag {
#   G -> {P C};
#   P -> C;
#   U -> {P C};
#   U [unobserved]
# }"
# dag_plot1 = dagitty( gen_dag )
# coordinates(dag_plot1) = list( x=c(G=0,P=1,C=1,U=2) , 
#                                y=c(G=0,P=0,C=1,U=0.5) )
# drawdag( dag_plot1 )
# 
# 
# 
# 
# # simulation
# # n = simulation sample size
# # bU, bGP, bPC, bGC = simulated parameters
# # rep = to use in replication
# #
# f_sim = function(n=100, bU=2, bGP=1, bPC=1, bGC=0, rep=F){
#   
#   # # test
#   # n=100; bU=2; bGP=1; bPC=1; bGC=0; rep=F
#   
#   # sim
#   U = 2*rbern( n , 0.5 ) - 1
#   G = rnorm( n )
#   P = rnorm( n , bGP*G + bU*U )
#   C = rnorm( n , bPC*P + bGC*G + bU*U )
#   d = data.frame(U=U,P=P,G=G,C=C)
#   
#   # return object
#   if(!rep){
#     # full data
#     return(d)
#     
#   } else{
#     # parameters
#     b1 = coef( lm(C ~ G + P + U, data=d) )['G'] # unbiased effect
#     b2 = coef( lm(C ~ G, data=d) )['G'] # more biased effects
#     b3 = coef( lm(C ~ G + P, data=d) )['G'] # biased effects
#     b = c(b1, b2, b3)
#     names(b) = c('G','Gb','Gs')
#     return( b )
#     
#   }
#   
# }
# 
# # relationships
# d = f_sim(n=100, bU=2, bGP=1, bPC=1, bGC=0, rep=F)
# psych::pairs.panels(d)
# # notice cor(G,C)>0, when it should be cor(G,C)=0
# 
# 
# 
# # models
# summary(lm(C ~ G, data=d)) # more biased effect
# summary(lm(C ~ G + P, data=d))  # less biased effect (change sign)
# summary(lm(C ~ G + P + U, data=d))  # unbiased effect (not possible)
# 
# # summary(lm(P ~ G, data=d)) # biased effect
# # summary(lm(P ~ G + U, data=d)) # unbiased effect (not possible)
# # summary(lm(C ~ P, data=d))  # biased effect
# 
# 
# 
# 
# # sampling variation
# par(mfrow=c(2,1))
# dsim = replicate( 1e4, f_sim(n=20, bU=2, bGP=1, bPC=1, bGC=0, rep=T) )
# f_plot1(dsim=dsim, ipar='G', xR=c(-2,3), by=0.5)
# 
# dsim = replicate( 1e4, f_sim(n=100, bU=2, bGP=1, bPC=1, bGC=0, rep=T) )
# f_plot1(dsim=dsim, ipar='G', xR=c(-2,3), by=0.5)
# par(mfrow=c(1,1))
# # G -> C, if we do not control for P and U (but it is not possible) 
# # equally biased with n=100, but more "confident" of G -> C
# # notice how relationship changes between models
# 
# 
# 
# 
# # What is going on?
# P_lim = quantile(d$P, c(0.45, 0.60))
# P_index = d$P>P_lim[1] & d$P<P_lim[2] 
# # parents at specific levels of education
# 
# U_index = d$U==-1
# 
# plot(d$G[U_index], d$C[U_index], col='black', 
#      xlim=range(d$G), ylim=range(d$C),
#      xlab='grandparent education (G)', ylab='granchild education (C)')
# points(d$G[!U_index], d$C[!U_index], col='blue')
# 
# # plotting parents
# points(d$G[P_index & U_index], d$C[P_index & U_index], 
#        col=col.alpha('black', 0.8), pch=19)
# points(d$G[P_index & !U_index], d$C[P_index & !U_index], 
#        col=col.alpha('blue', 0.8), pch=19)
# abline( lm(d$C[P_index] ~ d$G[P_index]) )
# # here we can show the negative association that we observe in m6.11





