library(numDeriv)
library(goftest)
#type-I basis function####
varphi_D=function(tau){
  re=tau^2*(alpha*tau*trigamma(alpha*tau)-1)
  return(re)
}
l_varphi_D=function(tau){
  return(log(varphi_D(tau)))
}
d_l_varphi_D=function(tau){
  return(1/2*grad(l_varphi_D,tau))
}

#grid search####
grid_search_D=function(){
  re=list()
  objold=0
  index=1
  for(n in 1:floor((1-co*dt)/(cm+ci))){
    for(m in 1:floor((1-n*ci)/(n*cm+co*dt))){
      obj=n^2*m^2*varphi_D((1-n*ci-m*n*cm)/(co*m))
      if(obj>objold){
        re=list(n=n,m=m,tau=(1-n*ci-m*n*cm)/(co*m),Te=(1-n*ci-m*n*cm)/(co*m)*m,value=1/obj)
        objold=obj
      }
      index=index+1
    }
  }
  return(re)
}

#exact D-optimal design####
exact_design_D=function(){
  #case 8
  if(ci+cm+co*dt==1){
    obj=varphi_D(dt)
    re=list(n=1,m=1,tau=dt,Te=dt,value=1/obj,case="case8")
    return(re)
  }
  #case 4
  tau_upper=(1-ci-cm)/co
  ccost=max(co/(ci+cm),co/(1-ci))
  if(d_l_varphi_D(tau_upper)>ccost){
    obj=varphi_D(tau_upper)
    re=list(n=1,m=1,tau=tau_upper,Te=tau_upper,value=1/obj,case="case4")
  }
  #case 1
  case2u1=ci*cm/(co*(1-2*ci))
  case2u2=(1-ci-cm)/co
  case2root=function(tau){
    d_l_varphi_D(tau)-co/(cm+co*tau)
  }
  if(1-2*ci-cm>=0&case2u1>dt){
    tau_star=tryCatch(uniroot(case2root,c(dt,case2u1))$root,error=function(e)e)
    if(is.list(tau_star)){
      
    }else{
      m=(1-ci)/(cm+co*tau_star)
      obj=m^2*varphi_D(tau_star)
      return(list(n=1,m=m,tau=tau_star,Te=tau_star*m,value=1/obj,case="case1"))
    }
  }else if(1-2*ci-cm<0&case2u2>dt){
    tau_star=tryCatch(uniroot(case2root,c(dt,case2u2))$root,error=function(e)e)
    if(is.list(tau_star)){
      
    }else{
      m=(1-ci)/(cm+co*tau_star)
      obj=m^2*varphi_D(tau_star)
      return(list(n=1,m=m,tau=tau_star,Te=tau_star*m,value=1/obj,case="case1"))
    }
  }
  #case 2
  case3u=case2u2
  case3l=ci/(co*(cm+2*ci))
  case3root=function(tau){
    d_l_varphi_D(tau)-co/(1-co*tau)
  }
  if(1-2*ci-cm>0&case3u>max(dt,case3l)){
    tau_star=tryCatch(uniroot(case3root,c(max(dt,case3l),case3u))$root,error=function(e)e)
    if(is.list(tau_star)){
      
    }else{
      n=(1-co*tau_star)/(ci+cm)
      obj=n^2*varphi_D(tau_star)
      return(list(n=n,m=1,tau=tau_star,Te=tau_star,value=1/obj,case="case2"))
    }
  }
  #case 3
  case4l=case2u1
  case4u=case3l
  case4root=function(tau){
    n=(-ci+sqrt(ci^2+cm*ci/(co*tau)))/(cm*ci/(co*tau))
    d_l_varphi_D(tau)-co/(cm*n+co*tau)
  }
  if(1-2*ci-cm>0&case4u>max(dt,case4l)){
    tau_star=tryCatch(uniroot(case4root,c(max(dt,case4l),case4u))$root,error=function(e)e)
    if(is.list(tau_star)){
      
    }else{
      n=(-ci+sqrt(ci^2+cm*ci/(co*tau_star)))/(cm*ci/(co*tau_star))
      m=(-ci+sqrt(ci^2+cm*ci/(co*tau_star)))/(cm)
      obj=n^2*m^2*varphi_D(tau_star)
      return(list(n=n,m=m,tau=tau_star,Te=tau_star*m,value=1/obj,case="case3"))
    }
  }
  #case 5
  case2u1=ci*cm/(co*(1-2*ci))
  if(1-2*ci-cm>=0&case2u1>dt&d_l_varphi_D(dt)<co/(cm+co*dt)){
    tau_star=tryCatch(uniroot(case2root,c(dt,case2u1))$root,error=function(e)e)
    m=(1-ci)/(cm+co*dt)
    obj=m^2*varphi_D(dt)
    return(list(n=1,m=m,tau=dt,value=1/obj,case="case5"))
  }else if(1-2*ci-cm<0&&d_l_varphi_D(dt)<co/(cm+co*dt)){
    tau_star=tryCatch(uniroot(case2root,c(dt,case2u1))$root,error=function(e)e)
    m=(1-ci)/(cm+co*dt)
    obj=m^2*varphi_D(dt)
    return(list(n=1,m=m,tau=dt,Te=dt*m,value=1/obj,case="case5"))
  }
  #case 6
  case3l=ci/(co*(cm+2*ci))
  if(1-2*ci-cm>0&dt>case3l&d_l_varphi_D(dt)<co/(1-co*dt)){
    n=(1-co*dt)/(ci+cm)
    obj=n^2*varphi_D(dt)
    return(list(n=n,m=1,tau=dt,Te=dt,value=1/obj,case="case6"))
  }
  #case 7
  case4l=case2u1
  case4u=case3l
  n=(-ci+sqrt(ci^2+cm*ci/(co*dt)))/(cm*ci/(co*dt))
  m=(-ci+sqrt(ci^2+cm*ci/(co*dt)))/(cm)
  if(1-2*ci-cm>0&case4u>dt&dt>case4l&d_l_varphi_D(dt)<co/(cm*n+co*dt)){
    tau_star=tryCatch(uniroot(case4root,c(max(dt,case4l),case4u))$root,error=function(e)e)
    obj=n^2*m^2*varphi_D(dt)
    return(list(n=n,m=m,tau=dt,Te=dt*m,value=1/obj,case="case7"))
  }
}

#basis A function####
varphi_A=function(tau){
  re=(1/(tau^2*trigamma(alpha*tau)-tau/alpha)+1/(alpha*tau))^{-1}
  return(re)
}

l_varphi_A=function(tau){
  return(log(varphi_A(tau)))
}
d_l_varphi_A=function(tau){
  return(grad(l_varphi_A,tau))
}

#exact A-optimal design####
grid_search_A=function(){
  re=list()
  objold=0
  index=1
  for(n in 1:floor((1-co*dt)/(cm+ci))){
    for(m in 1:floor((1-n*ci)/(n*cm+co*dt))){
      obj=n*m*varphi_A((1-n*ci-m*n*cm)/(co*m))
      # print(obj)
      # if(obj<0){
      #   browser()
      # }
      if(obj>objold){
        re=list(n=n,m=m,tau=(1-n*ci-m*n*cm)/(co*m),Te=(1-n*ci-m*n*cm)/(co*m)*m,value=1/obj)
        objold=obj
      }
      index=index+1
    }
  }
  return(re)
}
exact_design_A=function(){
  #case 8
  if(ci+cm+co*dt==1){
    obj=varphi_A(dt)
    re=list(n=1,m=1,tau=dt,Te=dt,value=1/obj,case="case8")
    return(re)
  }
  #case 4
  tau_upper=(1-ci-cm)/co
  ccost=max(co/(ci+cm),co/(1-ci))
  if(d_l_varphi_A(tau_upper)>ccost){
    obj=varphi_A(tau_upper)
    re=list(n=1,m=1,tau=tau_upper,Te=tau_upper,value=1/obj,case="case4")
  }
  #case 1
  case2u1=ci*cm/(co*(1-2*ci))
  case2u2=(1-ci-cm)/co
  case2root=function(tau){
    d_l_varphi_A(tau)-co/(cm+co*tau)
  }
  if(1-2*ci-cm>=0&case2u1>dt){
    tau_star=tryCatch(uniroot(case2root,c(dt,case2u1))$root,error=function(e)e)
    if(is.list(tau_star)){
      
    }else{
      m=(1-ci)/(cm+co*tau_star)
      obj=m*varphi_A(tau_star)
      return(list(n=1,m=m,tau=tau_star,Te=tau_star*m,value=1/obj,case="case1"))
    }
  }else if(1-2*ci-cm<0&case2u2>dt){
    tau_star=tryCatch(uniroot(case2root,c(dt,case2u2))$root,error=function(e)e)
    if(is.list(tau_star)){
      
    }else{
      m=(1-ci)/(cm+co*tau_star)
      obj=m*varphi_A(tau_star)
      return(list(n=1,m=m,tau=tau_star,Te=tau_star*m,value=1/obj,case="case1"))
    }
  }
  #case 2
  case3u=case2u2
  case3l=ci/(co*(cm+2*ci))
  case3root=function(tau){
    d_l_varphi_A(tau)-co/(1-co*tau)
  }
  if(1-2*ci-cm>0&case3u>max(dt,case3l)){
    tau_star=tryCatch(uniroot(case3root,c(max(dt,case3l),case3u))$root,error=function(e)e)
    if(is.list(tau_star)){
      
    }else{
      n=(1-co*tau_star)/(ci+cm)
      obj=n*varphi_A(tau_star)
      return(list(n=n,m=1,tau=tau_star,Te=tau_star,value=1/obj,case="case2"))
    }
  }
  #case 3
  case4l=case2u1
  case4u=case3l
  case4root=function(tau){
    n=(-ci+sqrt(ci^2+cm*ci/(co*tau)))/(cm*ci/(co*tau))
    d_l_varphi_A(tau)-co/(cm*n+co*tau)
  }
  if(1-2*ci-cm>0&case4u>max(dt,case4l)){
    tau_star=tryCatch(uniroot(case4root,c(max(dt,case4l),case4u))$root,error=function(e)e)
    if(is.list(tau_star)){
      
    }else{
      n=(-ci+sqrt(ci^2+cm*ci/(co*tau_star)))/(cm*ci/(co*tau_star))
      m=(-ci+sqrt(ci^2+cm*ci/(co*tau_star)))/(cm)
      obj=n*m*varphi_A(tau_star)
      return(list(n=n,m=m,tau=tau_star,Te=tau_star*m,value=1/obj,case="case3"))
    }
  }
  #case 5
  case2u1=ci*cm/(co*(1-2*ci))
  if(1-2*ci-cm>=0&case2u1>dt&d_l_varphi_A(dt)<co/(cm+co*dt)){
    tau_star=tryCatch(uniroot(case2root,c(dt,case2u1))$root,error=function(e)e)
    m=(1-ci)/(cm+co*dt)
    obj=m*varphi_A(dt)
    return(list(n=1,m=m,tau=dt,Te=dt*m,value=1/obj,case="case5"))
  }else if(1-2*ci-cm<0&&d_l_varphi_A(dt)<co/(cm+co*dt)){
    tau_star=tryCatch(uniroot(case2root,c(dt,case2u1))$root,error=function(e)e)
    m=(1-ci)/(cm+co*dt)
    obj=m*varphi_A(dt)
    return(list(n=1,m=m,tau=dt,Te=dt*m,value=1/obj,case="case5"))
  }
  #case 6
  case3l=ci/(co*(cm+2*ci))
  if(1-2*ci-cm>0&dt>case3l&d_l_varphi_A(dt)<co/(1-co*dt)){
    n=(1-co*dt)/(ci+cm)
    obj=n*varphi_A(dt)
    return(list(n=n,m=1,tau=dt,Te=dt,value=1/obj,case="case6"))
  }
  #case 7
  case4l=case2u1
  case4u=case3l
  n=(-ci+sqrt(ci^2+cm*ci/(co*dt)))/(cm*ci/(co*dt))
  m=(-ci+sqrt(ci^2+cm*ci/(co*dt)))/(cm)
  if(1-2*ci-cm>0&case4u>dt&dt>case4l&d_l_varphi_A(dt)<co/(cm*n+co*dt)){
    tau_star=tryCatch(uniroot(case4root,c(max(dt,case4l),case4u))$root,error=function(e)e)
    obj=n*m*varphi_A(dt)
    return(list(n=n,m=m,tau=dt,Te=dt*m,value=1/obj,case="case7"))
  }
}


#basis V function####
Qcdf=function(par){
  t=par[1]
  alpha=par[2]
  gamma=par[3]
  re=1-pgamma(eta,shape=alpha*t,rate=alpha*exp(-gamma))
  return(re)
}
xip=function(par){
  alpha=par[1]
  gamma=par[2]
  Qcdfroot=function(par){
    xi=par
    if(xi==0){
      return(-p)
    }else{
      re=1-pgamma(eta,shape=alpha*xi,rate=alpha*exp(-gamma))-p
    }
    return(re)
  }
  re=uniroot(Qcdfroot,c(0,exp(-gamma)*eta))
  return(re$root)
}
DQcdf=function(par){
  re=grad(Qcdf,c(xip(par),par))
  return(re)
}
varphi_V=function(tau){
  gradQcdf=DQcdf(c(alpha,gamma))
  h1=gradQcdf[2]/gradQcdf[1]
  h2=gradQcdf[3]/gradQcdf[1]
  re=((h1^2)/(tau^2*trigamma(alpha*tau)-tau/alpha)+(h2^2)/(alpha*tau))^{-1}
  return(re)
}

l_varphi_V=function(tau){
  return(log(varphi_V(tau)))
}
d_l_varphi_V=function(tau){
  return(grad(l_varphi_V,tau))
}
#exact V-optimal design####
exact_design_V=function(){
  #case 8
  if(ci+cm+co*dt==1){
    obj=varphi_V(dt)
    re=list(n=1,m=1,tau=dt,Te=dt,avar=1/obj,case="case8")
    return(re)
  }
  #case 4
  tau_upper=(1-ci-cm)/co
  ccost=max(co/(ci+cm),co/(1-ci))
  if(d_l_varphi_V(tau_upper)>ccost){
    obj=varphi_V(tau_upper)
    re=list(n=1,m=1,tau=tau_upper,Te=tau_upper,avar=1/obj,case="case4")
  }
  #case 1
  case2u1=ci*cm/(co*(1-2*ci))
  case2u2=(1-ci-cm)/co
  case2root=function(tau){
    d_l_varphi_V(tau)-co/(cm+co*tau)
  }
  if(1-2*ci-cm>=0&case2u1>dt){
    tau_star=tryCatch(uniroot(case2root,c(dt,case2u1))$root,error=function(e)e)
    if(is.list(tau_star)){
      
    }else{
      m=(1-ci)/(cm+co*tau_star)
      obj=m*varphi_V(tau_star)
      return(list(n=1,m=m,tau=tau_star,Te=tau_star*m,avar=1/obj,case="case1"))
    }
  }else if(1-2*ci-cm<0&case2u2>dt){
    tau_star=tryCatch(uniroot(case2root,c(dt,case2u2))$root,error=function(e)e)
    if(is.list(tau_star)){
      
    }else{
      m=(1-ci)/(cm+co*tau_star)
      obj=m*varphi_V(tau_star)
      return(list(n=1,m=m,tau=tau_star,Te=tau_star*m,avar=1/obj,case="case1"))
    }
  }
  #case 2
  case3u=case2u2
  case3l=ci/(co*(cm+2*ci))
  case3root=function(tau){
    d_l_varphi_V(tau)-co/(1-co*tau)
  }
  if(1-2*ci-cm>0&case3u>max(dt,case3l)){
    tau_star=tryCatch(uniroot(case3root,c(max(dt,case3l),case3u))$root,error=function(e)e)
    if(is.list(tau_star)){
      
    }else{
      n=(1-co*tau_star)/(ci+cm)
      obj=n*varphi_V(tau_star)
      return(list(n=n,m=1,tau=tau_star,Te=tau_star,avar=1/obj,case="case2"))
    }
  }
  #case 3
  case4l=case2u1
  case4u=case3l
  case4root=function(tau){
    n=(-ci+sqrt(ci^2+cm*ci/(co*tau)))/(cm*ci/(co*tau))
    d_l_varphi_V(tau)-co/(cm*n+co*tau)
  }
  if(1-2*ci-cm>0&case4u>max(dt,case4l)){
    tau_star=tryCatch(uniroot(case4root,c(max(dt,case4l),case4u))$root,error=function(e)e)
    if(is.list(tau_star)){
      
    }else{
      n=(-ci+sqrt(ci^2+cm*ci/(co*tau_star)))/(cm*ci/(co*tau_star))
      m=(-ci+sqrt(ci^2+cm*ci/(co*tau_star)))/(cm)
      obj=n*m*varphi_V(tau_star)
      return(list(n=n,m=m,tau=tau_star,Te=tau_star*m,avar=1/obj,case="case3"))
    }
  }
  #case 6
  case2u1=ci*cm/(co*(1-2*ci))
  if(1-2*ci-cm>=0&case2u1>dt&d_l_varphi_V(dt)<co/(cm+co*dt)){
    tau_star=tryCatch(uniroot(case2root,c(dt,case2u1))$root,error=function(e)e)
    m=(1-ci)/(cm+co*dt)
    obj=m*varphi_V(dt)
    return(list(n=1,m=m,tau=dt,avar=1/obj,case="case5"))
  }else if(1-2*ci-cm<0&&d_l_varphi_V(dt)<co/(cm+co*dt)){
    tau_star=tryCatch(uniroot(case2root,c(dt,case2u1))$root,error=function(e)e)
    m=(1-ci)/(cm+co*dt)
    obj=m*varphi_V(dt)
    return(list(n=1,m=m,tau=dt,Te=dt*m,avar=1/obj,case="case5"))
  }
  #case 7
  case3l=ci/(co*(cm+2*ci))
  if(1-2*ci-cm>0&dt>case3l&d_l_varphi_V(dt)<co/(1-co*dt)){
    n=(1-co*dt)/(ci+cm)
    obj=n*varphi_V(dt)
    return(list(n=n,m=1,tau=dt,Te=dt*1,avar=1/obj,case="case6"))
  }
  #case 8
  case4l=case2u1
  case4u=case3l
  n=(-ci+sqrt(ci^2+cm*ci/(co*dt)))/(cm*ci/(co*dt))
  m=(-ci+sqrt(ci^2+cm*ci/(co*dt)))/(cm)
  if(1-2*ci-cm>0&case4u>dt&dt>case4l&d_l_varphi_V(dt)<co/(cm*n+co*dt)){
    tau_star=tryCatch(uniroot(case4root,c(max(dt,case4l),case4u))$root,error=function(e)e)
    obj=n*m*varphi_V(dt)
    return(list(n=n,m=m,tau=dt,Te=dt*m,avar=1/obj,case="case7"))
  }
}
grid_search_V=function(){
  re=list()
  objold=0
  index=1
  for(n in 1:floor((1-co*dt)/(cm+ci))){
    for(m in 1:floor((1-n*ci)/(n*cm+co*dt))){
      obj=n*m*varphi_V((1-n*ci-m*n*cm)/(co*m))
      if(obj>objold){
        re=list(n=n,m=m,tau=(1-n*ci-m*n*cm)/(co*m),Te=(1-n*ci-m*n*cm)/(co*m)*m,value=1/obj)
        objold=obj
      }
      index=index+1
    }
  }
  return(re)
}
#type-II-basis function####
varrho_D=function(gm,te,alpha,gamma=NULL,dt,eta=NULL,p=NULL,log=T){
  m=gm
  omega=te-(m-1)*dt
  re=alpha*te*((m-1)*dt^2*psigamma(alpha*dt,1)+(omega)^2*psigamma(alpha*omega,1))-te^2
  if(log==T){
    return(log(re))
  }else{
    return(re)
  }
}
phi_D=function(n,m,te,alpha,gamma=NULL,dt,eta=NULL,p=NULL){
  re=varrho_D(m,te,alpha,gamma=NULL,dt,eta=NULL,p=NULL,F)
  return(1/(n^2*re))
}
varphi_Dm=function(m,te,alpha,gamma=NULL,dt,eta=NULL,p=NULL){
  re=grad(varrho_D,x=m,te=te,alpha=alpha,dt=dt,log=T)/2
  return(re)
}
varphi_DT=function(m,te,alpha,gamma=NULL,dt,eta=NULL,p=NULL){
  re=grad(varrho_D,x=te,gm=m,alpha=alpha,dt=dt,log=T)/2
  return(re)
}
varrho_A=function(gm,te,alpha,gamma=NULL,dt,eta=NULL,p=NULL,log=T){
  m=gm
  omega=te-(m-1)*dt
  re1=(m-1)*dt^2*psigamma(alpha*dt,1)+(omega)^2*psigamma(alpha*omega,1)-te/alpha
  re2=alpha*te
  re=(1/re1+1/re2)^(-1)
  if(log==T){
    return(log(re))
  }else{
    return(re)
  }
}
phi_A=function(n,m,te,alpha,gamma=NULL,dt,eta=NULL,p=NULL){
  re=varrho_A(m,te,alpha,gamma=NULL,dt,eta=NULL,p=NULL,F)
  return(1/(n*re))
}
varphi_Am=function(m,te,alpha,gamma=NULL,dt,eta=NULL,p=NULL){
  re=grad(varrho_A,x=m,te=te,alpha=alpha,dt=dt,log=T)
  return(re)
}
varphi_AT=function(m,te,alpha,gamma=NULL,dt,eta=NULL,p=NULL){
  re=grad(varrho_A,x=te,gm=m,alpha=alpha,dt=dt,log=T)
  return(re)
}
varrho_V=function(gm,te,alpha,gamma,dt,eta,p,log=T){
  m=gm
  Qcdf=function(par){
    t=par[1]
    alpha=par[2]
    gamma=par[3]
    re=1-pgamma(eta,shape=alpha*t,rate=alpha*exp(-gamma))
    return(re)
  }
  xip=function(par){
    alpha=par[1]
    gamma=par[2]
    Qcdfroot=function(par){
      xi=par
      if(xi==0){
        return(-p)
      }else{
        re=1-pgamma(eta,shape=alpha*xi,rate=alpha*exp(-gamma))-p
      }
      return(re)
    }
    re=uniroot(Qcdfroot,c(0,2*exp(-gamma)*eta))
    return(re$root)
  }
  DQcdf=function(par){
    re=grad(Qcdf,c(xip(par),par))
    return(re)
  }
  gradQcdf=DQcdf(c(alpha,gamma))
  h1=gradQcdf[2]/gradQcdf[1]
  h2=gradQcdf[3]/gradQcdf[1]
  omega=te-(m-1)*dt
  re1=(m-1)*dt^2*psigamma(alpha*dt,1)+(omega)^2*psigamma(alpha*omega,1)-te/alpha
  re2=alpha*te
  re=(h1^2/re1+h2^2/re2)^(-1)
  if(log==T){
    return(log(re))
  }else{
    return(re)
  }
}
phi_V=function(n,m,te,alpha,gamma,dt,eta,p){
  re=varrho_V(m,te,alpha,gamma,dt,eta,p,F)
  return(1/(n*re))
}
varphi_Vm=function(m,te,alpha,gamma,dt,eta,p){
  re=grad(varrho_V,x=m,te=te,alpha=alpha,gamma=gamma,dt=dt,eta=eta,p=p,log=T)
  return(re)
}
varphi_VT=function(m,te,alpha,gamma,dt,eta,p){
  re=grad(varrho_V,x=te,gm=m,alpha=alpha,gamma=gamma,dt=dt,eta=eta,p=p,log=T)
  return(re)
}

Qcdf=function(par){
  t=par[1]
  alpha=par[2]
  gamma=par[3]
  re=1-pgamma(eta,shape=alpha*t,rate=alpha*exp(-gamma))
  return(re)
}
xip=function(par){
  alpha=par[1]
  gamma=par[2]
  Qcdfroot=function(par){
    xi=par
    if(xi==0){
      return(-p)
    }else{
      re=1-pgamma(eta,shape=alpha*xi,rate=alpha*exp(-gamma))-p
    }
    return(re)
  }
  re=uniroot(Qcdfroot,c(0,exp(-gamma)*eta))
  return(re$root)
}
DQcdf=function(par){
  re=grad(Qcdf,c(xip(par),par))
  return(re)
}
varphi_V=function(tau){
  gradQcdf=DQcdf(c(alpha,gamma))
  h1=gradQcdf[2]/gradQcdf[1]
  h2=gradQcdf[3]/gradQcdf[1]
  re=((h1^2)/(tau^2*trigamma(alpha*tau)-tau/alpha)+(h2^2)/(alpha*tau))^{-1}
  return(re)
}

l_varphi_V=function(tau){
  return(log(varphi_V(tau)))
}
d_l_varphi_V=function(tau){
  return(grad(l_varphi_V,tau))
}
#optimal type-II design####
exact_typII=function(ci,cm,co,alpha,gamma=NULL,dt,eta=NULL,p=NULL,varphim,varphiT,phi){
  # browser()
  TC=function(n,m,te){
    ci*n+cm*n*m+co*te
  }
  #m,te,alpha,gamma=NULL,dt,eta=NULL,p=NULL
  #case1
  Tm=function(m){
    (1-cm*m-ci)/co
  }
  f=function(m){
    varphim(m,Tm(m),alpha,gamma,dt,eta,p)/varphiT(m,Tm(m),alpha,gamma,dt,eta,p)-cm/co
  }
  m1=try(uniroot(f,c(1,(1-ci-co*dt)/cm)))
  if(is.list(m1)){
    m1=m1$root
    if(co/(ci+cm*m1)<varphiT(m1,Tm(m1),alpha,gamma,dt,eta,p)){
      ns=1
      ms=m1
      Ts=Tm(m1)
      return(list(n=ns,m=ms,Te=Ts,obj=phi(ns,ms,Ts,alpha,gamma,dt,eta,p),cost=TC(ns,ms,Ts),case="case1"))
    }
  }
  
  #case2
  Tn=function(n){
    (1-n*(cm+ci))/co
  }
  f=function(n){
    n*varphiT(1,Tn(n),alpha,gamma,dt,eta,p)-co/(ci+cm)
  }
  n2=try(uniroot(f,c(1,(1-co*dt)/(ci+cm))))
  if(is.list(n2)){
    n2=n2$root
    if(varphim(1,Tn(n2),alpha,gamma,dt,eta,p)<cm/(ci+cm)){
      ns=n2
      ms=1
      Ts=Tn(n2)
      return(list(n=ns,m=ms,Te=Ts,obj=phi(ns,ms,Ts,alpha,gamma,dt,eta,p),cost=TC(ns,ms,Ts),case="case2"))
    }
  }
  #case3
  # browser()
  Tnm=function(n,m){
    (1-ci*n-cm*n*m)/co
  }
  f=function(par){
    # browser()
    n=par[1]
    m=par[2]
    if(n<1|m<1|m>(1-ci*n)/(cm*n+co*dt)|n>(1-co*dt)/(ci+cm)){
      return(Inf)
    }
    re1=abs(varphim(m,Tnm(n,m),alpha,gamma,dt,eta,p)-cm/(ci+cm*m))
    re2=abs(n*varphiT(m,Tnm(n,m),alpha,gamma,dt,eta,p)-co/(ci+cm*m))
    return(re1+re2)
  }
  # f1=function(n,m){
  #   # browser()
  #   re1=varphim(m,Tnm(n,m),alpha,gamma,dt,eta,p)-cm/(ci+cm*m)
  #   return(re1)
  # }
  # f2=function(n,m){
  #   # browser()
  #   re2=n*varphiT(m,Tnm(n,m),alpha,gamma,dt,eta,p)-co/(ci+cm*m)
  #   return(re2)
  # }
  n3=m3=2
  # browser()
  # for(iter in 1:1){
  #   n3=uniroot(f1,c(1,floor((1-ci)/(cm+co*dt))/2),n=m3)$root
  #   m3=uniroot(f2,c(1,floor((1-co*dt)/(ci+cm))/2),m=n3)$root
  # }
  # nini=mean(c(1,(1-co*dt)/(ci+cm)))
  # mini=mean(c(1,(1-ci*nini)/(cm*nini+co*dt)))
  
  n3m3=list(value=Inf)
  for(nini in 1:floor((1-co*dt)/(ci+cm))){
    mini=floor((1-ci*nini)/(cm*nini+co*dt))
    n3m3new=optim(c(nini,mini),f,list(maxit=10000))
    if(n3m3new$value<n3m3$value){
      n3m3=n3m3new
      # print(n3m3$value)
      if(n3m3$value<10^-10){
        break
      }
    }
  }
  
  n3=n3m3$par[1]
  m3=n3m3$par[2]
  if(n3>1&m3>1&Tnm(n3,m3)>m3*dt+10^-4){
    ns=n3
    ms=m3
    Ts=Tnm(n3,m3)
    return(list(n=ns,m=ms,Te=Ts,obj=phi(ns,ms,Ts,alpha,gamma,dt,eta,p),cost=TC(ns,ms,Ts),case="case3"))
  }
  #case4
  T4=(1-cm-ci)/co
  if(co/(ci+cm)<varphiT(1,T4,alpha,gamma,dt,eta,p)&varphim(1,T4,alpha,gamma,dt,eta,p)/varphiT(1,T4,alpha,gamma,dt,eta,p)<cm/co){
    ns=1
    ms=1
    Ts=T4
    return(list(n=ns,m=ms,Te=Ts,obj=phi(ns,ms,Ts,alpha,gamma,dt,eta,p),cost=TC(ns,ms,Ts),case="case4"))
  }
  #case5
  m5=(1-ci)/(cm+co*dt)
  if((cm+dt*co)/(ci+cm*m5)<varphim(m5,m5*dt,alpha,gamma,dt,eta,p)+dt*varphiT(m5,m5*dt,alpha,gamma,dt,eta,p)){
    ns=1
    ms=m5
    Ts=m5*dt
    return(list(n=ns,m=ms,Te=Ts,obj=phi(ns,ms,Ts,alpha,gamma,dt,eta,p),cost=TC(ns,ms,Ts),case="case5"))
  }
  #case6
  n6=(1-co*dt)/(ci+cm)
  if(n6*(varphim(1,dt,alpha,gamma,dt,eta,p)+dt*varphiT(1,dt,alpha,gamma,dt,eta,p))<(co*dt+cm*n6)/(ci+cm)&n6*varphiT(1,dt,alpha,gamma,dt,eta,p)<co/(ci+cm)){
    ns=n6
    ms=1
    Ts=dt
    return(list(n=ns,m=ms,Te=Ts,obj=phi(ns,ms,Ts,alpha,gamma,dt,eta,p),cost=TC(ns,ms,Ts),case="case6"))
  }
  #case7
  nm=function(m){
    (1-co*m*dt)/(ci+cm*m)
  }
  f=function(m){
    nm(m)*(varphim(m,m*dt,alpha,gamma,dt,eta,p)+dt*varphiT(m,m*dt,alpha,gamma,dt,eta,p))-(co*dt+cm*nm(m))/(ci+cm*m)
  }
  m7=try(uniroot(f,c(1,(1-ci)/(co*dt+cm))))
  if(is.list(m7)){
    m7=m7$root
    if(nm(m7)*varphiT(m7,m7*dt,alpha,gamma,dt,eta,p)<co/(ci+cm*m7)){
      ns=nm(m7)
      ms=m7
      Ts=m7*dt
      return(list(n=ns,m=ms,Te=Ts,obj=phi(ns,ms,Ts,alpha,gamma,dt,eta,p),cost=TC(ns,ms,Ts),case="case7"))
    }
  }
  #case8
  if(ci+cm+co*dt==1){
    ns=1
    ms=1
    Ts=dt
    return(list(n=ns,m=ms,Te=Ts,obj=phi(ns,ms,Ts,alpha,gamma,dt,eta,p),cost=TC(ns,ms,Ts),case="case8"))
  }
}
grid_search_II=function(phi){
  re=list()
  objold=Inf
  index=1
  for(n in 1:floor((1-co*dt)/(cm+ci))){
    for(m in 1:floor((1-n*ci-co*dt)/(n*cm))){
      if((1-n*ci-m*n*cm)/(co)<m*dt){
        
      }else{
        obj=phi(n,m,(1-n*ci-m*n*cm)/(co),alpha,gamma,dt,eta,p)
        if(obj<objold){
          re=list(n=n,m=m,te=(1-n*ci-m*n*cm)/(co),value=obj,cost=ci*n+cm*n*m+co*(1-n*ci-m*n*cm)/(co))
          objold=obj
        }
        index=index+1
      }
    }
  }
  return(re)
}
#figures####
#Lim (2015)
pdf("Lim_2015.pdf",width = 7,height=7,pointsize = 18)
alpha=0.065;gamma=-0.77;eta=0.5;p=0.1
curve(log(1/varphi_V(x)),20,2000,xlab=expression(tau),ylab=expression(phi[V](tau)),main="Lim (2015)",n=10000)
dev.off()
#Tseng(2009)
pdf("Tseng_2009.pdf",width = 7,height=7,pointsize = 18)
alpha=exp(4.17-4058.79/323);beta=1/0.0625;gamma=log(alpha/beta);eta=5;p=0.05
curve(log(1/varphi_V(x)),0.01,2000,xlab=expression(tau),ylab=expression(phi[V](tau)),main="Tseng (2009)")
dev.off()


#plot of K(tau) lim2015####
library(numDeriv)
alpha=0.065;gamma=-0.77
cb=1000;ci=30/cb;cm=1.9/cb;co=2.7/cb;eta=0.5;p=0.1;dt=5
K1=function(tau){
  (cm/co+tau)^(-1)
}
K2=function(tau){
  (1/co-tau)^(-1)
}
K3=function(tau){
  (tau*sqrt(1+cm/(co*ci*tau)))^(-1)
}
pdf("k_tau_ex1.pdf",width=7,height=7,pointsize=18)
par(mar=c(4,1,2.5,1))
plot(10,10,xlim=c(-9,log((1-ci-cm)/co)+1),ylim=log(c(co*(cm+2*ci)/(ci+cm),co/cm)),xlab=expression(tau),ylab="",xaxt="n",yaxt="n")
axis(1,at=log(c(ci*cm/(co*(1-2*ci)),ci/(co*(cm+2*ci)))),labels = expression(frac(C[it]*C[mea],C[op]*(1-2*C[it])),frac(C[it],C[op]*(C[mea]+2*C[it]))),padj=0.5,cex.axis=1.24)
axis(3,at=log(c((1-ci-cm)/co)),labels = expression(frac(1-C[it]-C[mea],C[op])),padj=0.5,cex.axis=1.24)
abline(v=log(c(ci*cm/(co*(1-2*ci)),ci/(co*(cm+2*ci)),(1-ci-cm)/co)),col="gray")
#abline(h=c(co/cm,K1(ci*cm/(co*(1-2*ci))),K3(ci/(co*(cm+2*ci))),K2((1-ci-cm)/co)),col="gray")
K1x=seq(-9,log(ci*cm/(co*(1-2*ci))),length.out=100)
lines((K1x),log(K1(exp(K1x))))
K3x=seq(log(ci*cm/(co*(1-2*ci))),log(ci/(co*(cm+2*ci))),length.out=100)
lines((K3x),log(K3(exp(K3x))))
K2x=seq(log(ci/(co*(cm+2*ci))),log((1-ci-cm)/co),length.out=100)
lines((K2x),log(K2(exp(K2x))))
abline(v=log(dt),lty=2)
axis(3,at=log(dt),labels = expression(Delta~t),padj=0.5,cex.axis=1.24)
curve(log(Vectorize(d_l_varphi_D)(exp(x))),log(ci*cm/(co*(1-2*ci))*0.1),log((1-ci-cm)/co),add = T,col="red",n=100)
case4root=function(tau){
  n=(-ci+sqrt(ci^2+cm*ci/(co*tau)))/(cm*ci/(co*tau))
  d_l_varphi_D(tau)-2*co/(cm*n+co*tau)
}
tau_star=tryCatch(uniroot(case4root,c(ci*cm/(co*(1-2*ci)),ci/(co*(cm+2*ci))))$root,error=function(e)e)
curve(log(Vectorize(d_l_varphi_A)(exp(x))),log(ci*cm/(co*(1-2*ci))*0.1),log((1-ci-cm)/co),add = T,col="green",n=100)
curve(log(Vectorize(d_l_varphi_V)(exp(x))),log(ci*cm/(co*(1-2*ci))*0.1),log((1-ci-cm)/co),add = T,col="blue",n=100)
dev.off()


#example 2
test_time=c(0,50,100,150,200,250)
LED_deg=matrix(c(13.4,17.9,17.3,20.2,24.9,16.3,27,13.8,18.8,33.2,33.9,23.5,
                 21.3,28.6,29.7,31.7,33.3,26.0,35.0,32.4,35.0,36.7,35.8,38.3,
                 24,34.6,36,37.7,37.2,32.6,39.3,37.3,39.4,40.7,40.6,38.7,
                 28.4,38.3,38.7,40,41,37,41.7,40,40.7,42.7,42.2,40.3,
                 32.0,42.0,40.7,41,46,38.7,42,40.3,42.7,43.5,44.7,44.0),ncol=5)
LED_deg=cbind(rep(10,12),LED_deg)
LED_deg=100-LED_deg

pdf("ex2data.pdf",width=7,height=7,pointsize = 18)
par(mar=c(4,4,1,1))
plot(-100,xlim=c(0,250),ylim=c(50,90),xlab="Time (hr)",ylab="Light intensity")
for(i in 1:nrow(LED_deg)){
  points(test_time,LED_deg[i,],pch=16)
  lines(test_time,LED_deg[i,])
}
dev.off()

likelihood=function(par){
  # browser()
  alpha=par[1]
  gam=par[2]
  time_tran=test_time
  inc=sapply(1:nrow(LED_deg),function(x){-diff(LED_deg[x,])})
  dtime_tran=diff(time_tran)
  re=sapply(1:ncol(inc),function(x){
    dt=dtime_tran
    dz=inc[,x]
    re2=alpha*dt*log(alpha)-lgamma(alpha*dt)+(alpha*dt-1)*log(dz)-alpha*(dz*exp(-gam)+gam*dt)
    return(re2)
  })
  res=sum(re)
  return(res)
}
gammaopt=optim(c(0.02,2.11),likelihood,control = list(fnscale=-1))
alpha_hat=gammaopt$par[1]
gamma_hat=gammaopt$par[2]

pdf("ex2dataprobplot.pdf",width=7,height=7,pointsize = 18)
par(mar=c(4,4,1,1))
inc=sapply(1:nrow(LED_deg),function(x){-diff(LED_deg[x,])})
py=((1:length(inc))-0.5)/length(inc)
qy=qgamma(py,shape=alpha_hat*50,rate=alpha_hat*exp(-gamma_hat))
qqplot(inc,qy,xlab="Sample quantiles",ylab="Probability",pch=16,yaxt="n")
# atp=seq(qgamma(c(0.001),shape=alpha_hat*50,rate=alpha_hat*exp(-gamma_hat)),qgamma(c(0.99),shape=alpha_hat*50,rate=alpha_hat*exp(-gamma_hat)),length.out=10)
atp=c(0.001,0.3,0.6,0.8,0.9,0.95,0.98,0.99)
# axis(2,at=atp,label=round(pgamma(atp,shape=alpha_hat*50,rate=alpha_hat*exp(-gamma_hat)),3))
axis(2,at=qgamma(atp,shape=alpha_hat*50,rate=alpha_hat*exp(-gamma_hat)),label=atp)
axissmallt=function(atp){
  n=length(atp)
  satp=c()
  for(i in 2:n){
    satp=c(satp,seq(atp[i-1],atp[i],length.out=11))
  }
  axis(2,at=qgamma(satp,shape=alpha_hat*50,rate=alpha_hat*exp(-gamma_hat)),tck=-0.015,labels=F)
}
axissmallt(atp)
abline(a=0,b=1)
dev.off()
ad.test(sort(inc),"pgamma",shape=alpha_hat*50,rate=alpha_hat*exp(-gamma_hat))


#tables####
#example 1
alpha=0.065;gamma=-0.77
cb=1000;ci=30/cb;cm=1.9/cb;co=2.7/cb;eta=0.5;p=0.1;dt=5
summarizeoptIandII=function(){
  opiiD=as.numeric(unlist(exact_typII(ci=ci,cm=cm,co=co,alpha=alpha,gamma=gamma,dt=dt,eta=eta,p=p,varphi_Dm,varphi_DT,phi_D)))
  # browser()
  opiiA=as.numeric(unlist(exact_typII(ci=ci,cm=cm,co=co,alpha=alpha,gamma=gamma,dt=dt,eta=eta,p=p,varphi_Am,varphi_AT,phi_A)))
  opiiV=as.numeric(unlist(exact_typII(ci=ci,cm=cm,co=co,alpha=alpha,gamma=gamma,dt=dt,eta=eta,p=p,varphi_Vm,varphi_VT,phi_V)))
  opiD=as.numeric(unlist(exact_design_D()))
  opiA=as.numeric(unlist(exact_design_A()))
  opiV=as.numeric(unlist(exact_design_V()))
  re=matrix(c(unlist(opiD)[c(5,5,1,2,4)]/c(1,unlist(opiiD)[4],1,1,1),
              unlist(opiiD)[c(4,4,1,2,3)]/c(1,unlist(opiiD)[4],1,1,1),
              unlist(opiA)[c(5,5,1,2,4)]/c(1,unlist(opiiA)[4],1,1,1),
              unlist(opiiA)[c(4,4,1,2,3)]/c(1,unlist(opiiA)[4],1,1,1),
              unlist(opiV)[c(5,5,1,2,4)]/c(1,unlist(opiiV)[4],1,1,1),
              unlist(opiiV)[c(4,4,1,2,3)]/c(1,unlist(opiiV)[4],1,1,1)),ncol=5,byrow=T)
  re[,2]=1/re[,2]
  re=data.frame(re)
  colnames(re)=c("obj","RE","n","m","Te")
  return(re)
}
options(digits = 3)
table1=summarizeoptIandII()
#example 2
alpha=0.02825;gamma=-2.0725
cb=1000;ci=7.56*10^(-2);cm=1.06*10^(-3);co=1.17*10^(-4);eta=50;p=0.05;dt=5
summarizeoptIandII=function(){
  opiiD=as.numeric(unlist(exact_typII(ci=ci,cm=cm,co=co,alpha=alpha,gamma=gamma,dt=dt,eta=eta,p=p,varphi_Dm,varphi_DT,phi_D)))
  # browser()
  opiiA=as.numeric(unlist(exact_typII(ci=ci,cm=cm,co=co,alpha=alpha,gamma=gamma,dt=dt,eta=eta,p=p,varphi_Am,varphi_AT,phi_A)))
  opiiV=as.numeric(unlist(exact_typII(ci=ci,cm=cm,co=co,alpha=alpha,gamma=gamma,dt=dt,eta=eta,p=p,varphi_Vm,varphi_VT,phi_V)))
  opiiDgrid=as.numeric(unlist(grid_search_II(phi_D)))
  opiiAgrid=as.numeric(unlist(grid_search_II(phi_A)))
  opiiVgrid=as.numeric(unlist(grid_search_II(phi_V)))
  opiD=as.numeric(unlist(exact_design_D()))
  opiDgrid=as.numeric(unlist(grid_search_D()))
  opiA=as.numeric(unlist(exact_design_A()))
  opiAgrid=as.numeric(unlist(grid_search_A()))
  opiV=as.numeric(unlist(exact_design_V()))
  opiVgrid=as.numeric(unlist(grid_search_V()))
  re=matrix(c(unlist(opiD)[c(5,5,1,2,4)]/c(1,unlist(opiiD)[4],1,1,1),
              unlist(opiDgrid)[c(5,5,1,2,4)]/c(1,unlist(opiiD)[4],1,1,1),
              unlist(opiiD)[c(4,4,1,2,3)]/c(1,unlist(opiiD)[4],1,1,1),
              unlist(opiiDgrid)[c(4,4,1,2,3)]/c(1,unlist(opiiD)[4],1,1,1),
              c(1/3600/varphi_D(250),1/3600/varphi_D(250)/unlist(opiiD)[4],12,5,250),
              unlist(opiA)[c(5,5,1,2,4)]/c(1,unlist(opiiA)[4],1,1,1),
              unlist(opiAgrid)[c(5,5,1,2,4)]/c(1,unlist(opiiA)[4],1,1,1),
              unlist(opiiA)[c(4,4,1,2,3)]/c(1,unlist(opiiA)[4],1,1,1),
              unlist(opiiAgrid)[c(4,4,1,2,3)]/c(1,unlist(opiiA)[4],1,1,1),
              c(1/60/varphi_A(250),1/60/varphi_A(250)/unlist(opiiA)[4],12,5,250),
              unlist(opiV)[c(5,5,1,2,4)]/c(1,unlist(opiiV)[4],1,1,1),
              unlist(opiVgrid)[c(5,5,1,2,4)]/c(1,unlist(opiiV)[4],1,1,1),
              unlist(opiiV)[c(4,4,1,2,3)]/c(1,unlist(opiiV)[4],1,1,1),
              unlist(opiiVgrid)[c(4,4,1,2,3)]/c(1,unlist(opiiV)[4],1,1,1),
              c(1/60/varphi_V(250),1/60/varphi_V(250)/unlist(opiiV)[4],12,5,250)),ncol=5,byrow=T)
  re[,2]=1/re[,2]
  re=data.frame(re)
  colnames(re)=c("obj","RE","n","m","Te")
  return(re)
}
options(digits = 3)
table2=summarizeoptIandII()