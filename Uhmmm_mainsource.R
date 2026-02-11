

Uhmmm.MLfit<-function(table,modelobs=NULL,filemodelobs=NULL,modellat=NULL,filemodellat=NULL,modelcomp2=NULL,start=NULL,startvec=NULL,startmth="SANN",method="NM",method.inv="Broyden",print.level=0,
                      finalHessian=TRUE,tol=1e-09,iterlim=1000,NMiter=0,numeric=FALSE,
                      fixed=NULL,compare=FALSE,recalc=FALSE,
                      control=list(allowSingular=FALSE)){
                  if(!is.null(filemodelobs)){load(filemodelobs)}
                  if(!is.null(filemodellat)){load(filemodellat)}
  contr<-list(tol=tol,gradtol=tol,reltol=tol)
 
  
                        lev<-prod(modelobs$modello$livelli)
                       # global_p<-NULL
                       
                        strata<-modelobs$modello$strata
                        nlat<-length(modelobs$modello$livelli)
                        stratalat<-2^nlat
                        stratax<-strata/stratalat
        
                        lev2<-lev*stratalat



if(is.null(modelcomp2)){
marglist<-NULL
for(i in 1:nlat){
  a<-combinations(nlat,i)
  m<-matrix("m",nlat,dim(a)[1])
  for(j in 1: dim(a)[1]){
    m[a[j,],j]<-"b"
  }
  marglist<- c(marglist,apply(m,MARGIN=2,FUN=paste,collapse="-")) # tolto c() ad apply
}

modelcomp2<-marg.list(paste(paste(rep("m-",nlat),collapse=
""),marglist,sep=""),mflag="m")
}

                        modelcomp2<-hmmm.model(marg=modelcomp2,lev=c(modelobs$modello$livelli,modellat$model$livelli),strata=stratax,
                                               names=c(modelobs$names,modellat$names))
                        
                        
                   

sparse.dir.sum <-  function(...){
  nmat <- nargs()
  allmat<- list(...)
  C<-allmat[[1]]
  
  for(i in 2:nmat) {
    B<-allmat[[i]]
    #A<-C
    C <- rbind(cbind(C, Matrix(0, nrow = nrow(C), ncol = ncol(B))), 
               cbind(Matrix(0, nrow = nrow(B), ncol = ncol(C)), B))}
  return(C)
}
                    
C<-sparse.dir.sum(Matrix(modelcomp2$matrici$CMAT),Matrix(modelobs$matrici$CMAT))
                      
                     M<- rbind2(Matrix(modelcomp2$matrici$MMAT),Matrix(modelobs$matrici$MMAT))
remove(modelcomp2)

if(nlat==1){
modellat1<-hmmm.model(marg=modellat[[1]][1],lev=modellat$modello$livelli,strata=1,names=modellat$names) }
else{modellat1<-hmmm.model(lev=modellat$modello$livelli,strata=1,names=modellat$names) }

                       modelcomp1<-hmmm.model(marg=modelobs$marglist,lev=modelobs$modello$livelli,strata=1,names=modelobs$names)
                    
                       # print("computing arrays for model given lat and obs")
                     Xl<-array(c( modellat$matrici$X),dim=c(stratalat-1,stratax,dim(modellat$matrici$X)[2]))
                        Xo<-array(c(modelobs$matrici$X),dim=c(lev-1,stratax,stratalat,dim(modelobs$matrici$X)[2]))
                        Xo<-aperm(Xo,c(1,3,2,4))
                       
                      nbeta<-dim(modellat$matrici$X)[2]+dim(modelobs$matrici$X)[2]
                       ZF<-Matrix(kronecker(diag(strata/stratalat),matrix(1,lev2,1)))
                      
                       H<-diag(1,lev)
                       h<-matrix(1,lev)
                       H<-cbind(-h, H[,-1])
                       H<-H[-1,]
                       H<-Matrix(diag(1,strata/stratalat)%x%H)
                       L<-Matrix(diag(1,stratax)%x%matrix(1,1,stratalat)%x%diag(1,lev))
                    
                    
                       if(is.null(start)) start<-rep(0,nbeta)
                      # print("starting values")
                      # print(start)
                        


                      
                       inv_GMI<-function(etpar,mod
                                         ,start,interf=NULL,interfder=NULL,Zlist=NULL){
                      
                        if(is.null(interf)){
                         interf<-mod$functions$L.fct 
                         interfder<-mod$functions$derLt.fct}
                        if(is.null(Zlist)){
                         Zlist<-cocadise(matrix(1,length(start)+1))} 
                         
                         
                         myfun<-function(x){
                        
                           x<-as.matrix(exp(Zlist$DMAT%*%x))#
                          
                           x<-x/sum(x)          
                           interf(x)-etpar}
                         myder<-function(x){
                           x<-as.matrix(exp(Zlist$DMAT%*%x))#
                           x<-pmax(0.001*sum(x),x)
                           x<-x/sum(x)
                           D<-t( t(Zlist$DMAT)%*%
                           Diagonal(x=c(x)) %*%interfder(x))
                        
                         as.matrix(D)}
                         r<-try(nleqslv(x=start,myfun,method=method.inv,jac=myder,
                                    control=control,global="qline",jacobian=FALSE)$x,silent=TRUE)#,silent=TRUE
                         if(class(r)=="try-error"){
                         myfun2<-function(x){ sum(myfun(x)*myfun(x))}
                         r<-optim(par=start,fn=myfun2,method="Nelder-Mead",
                                  control=list(abstol=10^-8,maxit=1000))
                         
                        r<-r$par}
                        
                         r<-as.matrix(exp(Zlist$DMAT%*%r))
                         r<-r/sum(r) 
                       
                         r
                       }
                       
                       
                        
                        f_loglik<-function(beta){
                          beta_l<-beta[1:dim( modellat$matrici$X)[2]]
                          beta_o<-beta[(dim( modellat$matrici$X)[2]+1):(dim( modellat$matrici$X)[2]+dim(modelobs$matrici$X)[2])]
                          
                          myf<-function(Xlocal,beta,mod){
                          eta<-Xlocal%*%beta
                          
                          startheta<-rep(0,length(eta))
                        p<-inv_GMI(eta,mod,start=startheta)
                        
                        
                        

                       }
                       
                       p<-apply(X=Xo,MARGIN=c(2,3),FUN=myf,beta=beta_o,mod=modelcomp1)#MARGIN=c(2,3)
                      
                       plat<-apply(X=Xl,MARGIN=2,FUN=myf,beta=beta_l,mod=modellat1)
                       p<-matrix(c(p),lev,stratalat*stratax)
                       p<-p%*%diag(c(plat))
                       p<-matrix(c(p),lev*stratalat,stratax)
                       
                       P<-matrix(1,1,stratalat)%x%diag(1,lev)%*%p
                       plik<-pmax(0.0001,vec(P))
                       fl<-t(vec(table))%*%log(plik)
                       if(!recalc){
                         if(!is.na(fl)&!is.infinite(fl) ){
                       if(fl>=global_loglik){
                         global_loglik<<-fl
                       global_p<<-p}}
                       }
                       fl
                        }
                       f_grad<-function(beta){
                         p<-global_p
                         if(recalc||compare){
                           
                           beta_l<-beta[1:dim( modellat$matrici$X)[2]]
                           beta_o<-beta[(dim( modellat$matrici$X)[2]+1):(dim( modellat$matrici$X)[2]+dim(modelobs$matrici$X)[2])]
                           
                           myf<-function(Xlocal,beta,mod){
                             eta<-Xlocal%*%beta
                             #
                             startheta<-rep(0,length(eta))
                             p<-inv_GMI(eta,mod,start=startheta)
                             
                             
                             
                           }
                           
                           p<-apply(X=Xo,MARGIN=c(2,3),FUN=myf,beta=beta_o,mod=modelcomp1)
                           
                           plat<-apply(X=Xl,MARGIN=2,FUN=myf,beta=beta_l,mod=modellat1)
                           p<-matrix(c(p),lev,stratalat*stratax)
                           p<-p%*%diag(c(plat))
                           p<-matrix(c(p),lev*stratalat,stratax)
                          
                          
                         }
                          
                        
                         m<-p%*%diag(c(colSums((table))),dim(table)[2])
                        
                         x<-vec(p)
                         
                         q<-matrix(L%*%x)
                         m<-vec(m)
                         
              
                         omega  <- Diagonal(x=c(x))-((ZF*c(x))%*%t(ZF*c(x)))
                        X<-modellat$matrici$X%s%matrix(c(Xo),(lev-1)*stratalat*stratax,dim(modelobs$matrici$X)[2])
                        
                        xx<-matrix(M%*%x)
                        
                        R<-try(solve(C%*%Diagonal(x=c(1/(xx)))%*%M%*%omega%*%(diag(1,stratax)%x%diag(1,lev2)[,-1]))%*%X,silent=TRUE)###
                        if(class(R)=="try-error"){
                         # print("try-error I'll use ginv")
                        R<-ginv(as.matrix(C%*%Diagonal(x=c(1/(xx)))%*%M%*%omega%*%(diag(1,stratax)%x%diag(1,lev2)[,-1])))%*%X###
                        } ###
                         Q<-H%*%Diagonal(x=1/c(q))%*%L%*%omega%*%(diag(1,strata/stratalat)%x%diag(1,lev2)[,-1])
                         
                         D<-Q%*%R
                         
                        
                         grad<-colSums(t(vec(table)-L%*%vec(m))%*%(diag(1,strata/stratalat)%x%diag(1,lev)[,-1])%*%D)
                              

                        
                       }#grad

                       f_hess<-function(beta){
                         beta_l<-beta[1:dim( modellat$matrici$X)[2]]
                         beta_o<-beta[(dim( modellat$matrici$X)[2]+1):(dim( modellat$matrici$X)[2]+dim(modelobs$matrici$X)[2])]
                         
                         myf<-function(Xlocal,beta,mod){
                           eta<-Xlocal%*%beta
                           
                           startheta<-rep(0,length(eta))
                           p<-inv_GMI(eta,mod,start=startheta)
                           
                           
                           
                           
                         }
                         
                         p<-apply(X=Xo,MARGIN=c(2,3),FUN=myf,beta=beta_o,mod=modelcomp1)#MARGIN=c(2,3)
                         
                         plat<-apply(X=Xl,MARGIN=2,FUN=myf,beta=beta_l,mod=modellat1)
                         p<-matrix(c(p),lev,stratalat*stratax)
                         p_array<-array(c(p),dim=c(lev,stratalat,stratax))
                         p<-p%*%diag(c(plat))
                         pp<-matrix(c(p),lev*stratalat,stratax)
                         
                         PQ<-matrix(1,1,stratalat)%x%diag(1,lev)%*%pp
                         
                         q<-vec(PQ)
                         X<-modellat$matrici$X%s%matrix(c(Xo),(lev-1)*stratalat*stratax,dim(modelobs$matrici$X)[2])
                         #removed staff
                         
                         
                         ZFM<-kronecker(diag(strata/stratalat),matrix(1,lev,1))
                         
                         m<-pp%*%diag(c(colSums(matrix(table))),dim(table)[2])
                         x<-p<-vec(pp)
                         q<-vec(PQ)
                         m<-vec(m)
                         taum<-vec(table)
                         tau<-vec(prop.table(table,2))
                         taumu<-rep(colSums(table),each=lev)
                         omega  <- Diagonal(x=c(x))-((ZF*c(x))%*%t(ZF*c(x)))
                         xx<-matrix(M%*%x)
                         R<-try(solve(C%*%Diagonal(x=c(1/(xx)))%*%M%*%omega%*%(diag(1,stratax)%x%diag(1,lev2)[,-1])),silent=TRUE)
                         if(class(R)=="try-error"){
                           print("try-error I'll use ginv")
                           R<-ginv(C%*%Diagonal(x=c(1/(xx)))%*%M%*%omega%*%(diag(1,stratax)%x%diag(1,lev2)[,-1]))
                         }
                         R<-R%*%X
                      
                         
                         Q<-H%*%Diagonal(x=1/c(q))%*%L%*%omega%*%(diag(1,strata/stratalat)%x%diag(1,lev2)[,-1])
                         D<-Q%*%R
                         q<-vec(PQ)
                         qm<-matrix(L%*%m)
                        
                         omegabarp<-Diagonal(x=c(q))-((ZFM*c(q))%*%t(ZFM*c(q)))
                         omegabarm<-diag(c(taum))-(ZFM*c(taum))%*%t(ZFM*c(tau)) 
                         
                         
                         B<-Diagonal(x=1/c(q))%*%omegabarp%*%(diag(1,strata/stratalat)%x%diag(1,lev)[,-1])%*%D
                         
                         
                         -t(B)%*%(omegabarm+(ZFM*c(taumu*(tau-q)))%*%(t(ZFM*c(tau-q))))%*%B
                         
                       }#hess
                       
                       
                       

                      
beta<-start
beta_l<-beta[1:dim( modellat$matrici$X)[2]]
beta_o<-beta[(dim( modellat$matrici$X)[2]+1):(dim( modellat$matrici$X)[2]+dim(modelobs$matrici$X)[2])]

myf<-function(Xlocal,beta,mod){
  eta<-Xlocal%*%beta
  
  startheta<-rep(0,length(eta))
  p<-inv_GMI(eta,mod,start=startheta)
  
 
  
}

p<-apply(X=Xo,MARGIN=c(2,3),FUN=myf,beta=beta_o,mod=modelcomp1)

plat<-apply(X=Xl,MARGIN=2,FUN=myf,beta=beta_l,mod=modellat1)
p<-matrix(c(p),lev,stratalat*stratax)
p<-p%*%diag(c(plat))
p<-matrix(c(p),lev*stratalat,stratax)


  global_p<-p
global_loglik<--9999999

                       if(compare){
                         
                         return(comp<-compareDerivatives(f=f_loglik,grad=f_grad,t0=start,print=FALSE))}
                       
                   if(numeric){
                     
                     
                     fit<-maxLik(
                       logLik=f_loglik,grad=NULL,hess=NULL,start=start,method=method,print.level=print.level,
                       iterlim=iterlim,tol=tol,finalHessian=FALSE,fixed=fixed)
                   }
                   else{
                       
                      
                       if(NMiter>0){
                         #cat ("computing", startmth," -starting values","\n")
                         fit<-maxLik(logLik=f_loglik,
                                     ,grad=f_grad,
                                     start=start,method=startmth,print.level=print.level,
                                     iterlim=NMiter,finalHessian=FALSE,fixed=fixed,tol=tol)
                         
                         start<-fit$estimate
                         #cat(startmth,"- starting values","\n")
                         #print(start)
                       }
                       save(file="MyDayStart",start)
                       #print("starting values saved in file  MyDayStart")
                       contr<-list(tol=tol,gradtol=tol,reltol=tol)
                       
                       fit<-maxLik(logLik=f_loglik,grad=f_grad,start=start,method=method,print.level=print.level,
                                   hess= f_hess,
                                   iterlim=iterlim,finalHessian=FALSE,fixed=fixed,
                                   control=contr
                                   )
                  
                      
                      
                      ###             #####
                      
                      beta<-fit$estimate
                      beta_l<-beta[1:dim( modellat$matrici$X)[2]]
                      beta_o<-beta[(dim( modellat$matrici$X)[2]+1):(dim( modellat$matrici$X)[2]+dim(modelobs$matrici$X)[2])]
                      
                      myf<-function(Xlocal,beta,mod){
                        eta<-Xlocal%*%beta
                        
                        startheta<-rep(0,length(eta))
                        p<-inv_GMI(eta,mod,start=startheta)
                        
                        
                        
                        
                      }
                      
                      p<-apply(X=Xo,MARGIN=c(2,3),FUN=myf,beta=beta_o,mod=modelcomp1)#MARGIN=c(2,3)
                      
                      plat<-apply(X=Xl,MARGIN=2,FUN=myf,beta=beta_l,mod=modellat1)
                      p<-matrix(c(p),lev,stratalat*stratax)
                      p_array<-array(c(p),dim=c(lev,stratalat,stratax))
                      p<-p%*%diag(c(plat))
                      pp<-matrix(c(p),lev*stratalat,stratax)
                      
                      PQ<-matrix(1,1,stratalat)%x%diag(1,lev)%*%pp
                      
                      q<-vec(PQ)
                   }
                   
                       if(!finalHessian){
                         X<-modellat$matrici$X%s%matrix(c(Xo),(lev-1)*stratalat*stratax,dim(modelobs$matrici$X)[2])
                         vecpar<-as.matrix(fit$estimate)
                         rownames(vecpar)<-c(colnames( modellat$matrici$X),colnames(modelobs$matrici$X))
                         veceta=X%*%fit$estimate
                         etami<-(C%*%log(M%*%vec(p)))
                        
                         controllo<-cbind(etami,veceta)
                         controllo<-controllo[-(1:((stratalat-1)*stratax)),]
                         #lev<-c(3,3,3)
                         #levlat<-c(2,2,2)
                         controllo<-cbind(controllo, controllo[,1]-controllo[,2])
                         controllo<-array(t(controllo),dim=c(3,lev-1,stratalat,stratax))
                         controllo<-aperm(controllo,c(2,1,3,4))
                         
                         
                         vecetalat<-matrix(veceta[(1:((stratalat-1)*stratax)),1])
                         vecetalatmat<-matrix(veceta[(1:((stratalat-1)*stratax))],stratalat-1,stratax)
                         aa<-array(veceta[-(1:((stratalat-1)*stratax))],dim=c(lev-1,stratalat,stratax))
                         
                         a<-aa[,stratalat,]
                         veceta<-matrix(c(a),((lev-1)*stratax),1)
                         
                         veceta<-rbind(vecetalat,veceta)
                         rownames(veceta)<-c( paste("Linklat",(1:((stratalat-1)*stratax)),sep=""),paste("Linkobs",1:((lev-1)*stratax),sep=""))
   #controllo
   p_array<-array(c(p),dim=c(lev,stratalat,stratax))
   myf.eta<-function(plocal,mod){
     eta<-mod$matrici$C%*%log(mod$matrici$M%*%plocal)
     }
   
   etaobs.contr<-apply(X=p_array,MARGIN=c(2,3),FUN=myf.eta,mod=modelcomp1)#MARGIN=c(2,3)
   
   
   
   etalat.contr<-apply(X=plat,MARGIN=2,FUN=myf.eta,mod=modellat1)
                      
                         
                         
                         
                         list(fit=fit,vecpar=vecpar,veceta=veceta,vecPtr=vec(plat),vecPobs=q,
                              etarray=aa,etalatarray=vecetalatmat,contretaobs=etaobs.contr,
                              contretalat=etalat.contr,
                              vecjoint=p,controllo=controllo)}
                       
                      






else{ 





                       
                        X<-modellat$matrici$X%s%matrix(c(Xo),(lev-1)*stratalat*stratax,dim(modelobs$matrici$X)[2])
                   
L_pi<-diag(1,stratax)%x%(diag(1,stratalat)%x%matrix(1,1,lev))

ZFM<-kronecker(diag(strata/stratalat),matrix(1,lev,1))

            f_Hess<-function(beta){

m<-pp%*%diag(c(colSums(matrix(table))),dim(table)[2])
x<-p<-vec(pp)
q<-vec(PQ)
m<-vec(m)

omega  <- Diagonal(x=c(x))-((ZF*c(x))%*%t(ZF*c(x)))
xx<-matrix(M%*%x)
R<-try(solve(C%*%Diagonal(x=c(1/(xx)))%*%M%*%omega%*%(diag(1,stratax)%x%diag(1,lev2)[,-1])),silent=TRUE)
if(class(R)=="try-error"){
  print("try-error I'll use ginv")
  R<-ginv(C%*%Diagonal(x=c(1/(xx)))%*%M%*%omega%*%(diag(1,stratax)%x%diag(1,lev2)[,-1]))
}
R<-R%*%X


Q<-H%*%Diagonal(x=1/c(q))%*%L%*%omega%*%(diag(1,strata/stratalat)%x%diag(1,lev2)[,-1])
D<-Q%*%R
q<-vec(PQ)
qm<-matrix(L%*%m)
omegabarm<-(Diagonal(x=c(qm))-((ZFM*c(qm))%*%t(ZFM*c(q))))
omegabarp<-Diagonal(x=c(q))-((ZFM*c(q))%*%t(ZFM*c(q)))





B<-Diagonal(x=1/c(q))%*%omegabarp%*%(diag(1,strata/stratalat)%x%diag(1,lev)[,-1])%*%D

-t(B)%*%omegabarm%*%B


}#hess
 
                  
m<-pp%*%diag(c(colSums(matrix(table))),dim(table)[2])
x<-p<-vec(pp)
q<-vec(PQ)
m<-vec(m)
                   
                          
                           
                     
                       
                       Fisher<--f_Hess(beta)
                       

omega  <- Diagonal(x=c(x))-((ZF*c(x))%*%t(ZF*c(x)))
xx<-matrix(M%*%x)

R<-try(solve(C%*%Diagonal(x=c(1/(xx)))%*%M%*%omega%*%(diag(1,stratax)%x%diag(1,lev2)[,-1])),silent=TRUE)
if(class(R)=="try-error"){
  print("try-error I'll use ginv")
  R<-svd.inverse(as.matrix(C%*%Diagonal(x=c(1/(xx)))%*%M%*%omega%*%(diag(1,stratax)%x%diag(1,lev2)[,-1])))
}
R<-R%*%X


qm<-matrix(L%*%m)
Q<-H%*%Diagonal(x=1/c(q))%*%L%*%omega%*%(diag(1,strata/stratalat)%x%diag(1,lev2)[,-1])
D<-Q%*%R


omegabarm<-(Diagonal(x=c(qm))-((ZFM*c(qm))%*%t(ZFM*c(q))))/sum(table)
omegabarp<-Diagonal(x=c(q))-((ZFM*c(q))%*%t(ZFM*c(q)))
B<-Diagonal(x=1/c(q),length(c(q)))%*%omegabarp%*%(diag(1,strata/stratalat)%x%diag(1,lev)[,-1])%*%D

Fisher<-t(B)%*%omegabarm%*%B
Z<-Matrix(diag(1,strata/stratalat)%x%diag(1,lev2)[,-1])
###################################################################
V_beta=try(solve(Fisher)/sum(table),silent=TRUE)
if(class(V_beta)=="try-error"){
  print("try-error singular Fisher") 
  print("Fisher rank")
  print(matrix.rank(as.matrix(Fisher)))
  V_beta<-svd.inverse(as.matrix(Fisher)/sum(table))}
V_eta<-X%*%V_beta%*%t(X)

V_bj<-L%*%omega%*%Z%*%R%*%V_beta%*%t(L%*%omega%*%Z%*%R)
V_pi<-L_pi%*%omega%*%Z%*%R%*%V_beta%*%t(L_pi%*%omega%*%Z%*%R)



vecpar<-cbind(beta,sqrt(diag(V_beta)))
zeta<-vecpar[,1]/vecpar[,2]
vecpar<-cbind(vecpar,zeta,2-2*pnorm(abs(zeta)))
vecpar<-cbind(1:dim(vecpar)[1],vecpar)

rownames(vecpar)<-c(colnames( modellat$matrici$X),colnames(modelobs$matrici$X))
colnames(vecpar)<-c("#","MLE","S.E.","Z","P")
vecPtr<-cbind(vec(plat),sqrt(diag(V_pi)))
zeta<-(vecPtr[,1]-1/stratalat)/vecPtr[,2]
vecPtr<-cbind(vecPtr,zeta, 2-2*pnorm(abs(zeta)))
colnames(vecPtr)<-c("MLE","S.E.","Z","P")
veceta<-cbind(X%*%beta,sqrt(diag(V_eta)))
zeta<-veceta[,1]/veceta[,2]
veceta<-cbind(veceta,zeta, 2-2*pnorm(abs(zeta)))
veceta[veceta[,3]==Inf|veceta[,3]==-Inf,3:4]<-NaN
vecetalat<-veceta[(1:((stratalat-1)*stratax)),]
a<-array(veceta[-(1:((stratalat-1)*stratax)),],dim=c(lev-1,stratalat,stratax,4))
etarray<-a[,,,1]
a<-a[,stratalat,,]
veceta<-matrix(c(a),((lev-1)*stratax),4)
veceta<-rbind(vecetalat,veceta)
rownames(veceta)<-c( paste("Linklat",(1:((stratalat-1)*stratax)),sep=""),paste("Linkobs",1:((lev-1)*stratax),sep=""))


colnames(veceta)<-c("MLE","S.E.","Z","P")
vecbj<-cbind(q,sqrt(diag(V_bj)))
#print(vecbj)

colnames(vecbj)<-c("MLE","S.E.")
vecres<-vec(prop.table(table,2))-q
V_res<-diag(1/rowSums(ZFM%*%t(ZFM)%*%(ZFM*c(qm))))%*%omegabarp-V_bj
vecres<-cbind(vecres,sqrt(diag(V_res)))
zeta<-vecres[,1]/vecres[,2]
vecres<-cbind(vecres, zeta,2-2*pnorm(abs(zeta)))
colnames(vecres)<-c("MLE","S.E.","Z","P")
lista<-list(fit=fit,
            VarMatrices=list(Fisher=Fisher,V_beta=V_beta,V_eta=V_eta,V_bj=V_bj,V_pi=V_pi,V_res=V_res),
            vecpar=vecpar,veceta=veceta,etarray=etarray, vecjoint=p,vecPtr=vecPtr,vecPobs=vecbj,vecres=vecres)

}
                       
                       
                       
                      }




hmmm.model.U<-function(marg,lev,names,Formula=NULL,strata=1,fnames=NULL,
                       cocacontr=NULL,ncocacontr=NULL,replace=TRUE,INDIP="UNR",UNC="uniform"){
  #alternative a "uniform" per UNC sono "g.uniform","b.parabolic""g.parabolic","l.parabolic" 
  #"sg.uniform","sb.parabolic""sg.parabolic","l.parabolic"
  #da usare solo
  #con i corrispondenti logits g b l per le marginali delle osservabili
  #con logit global usare  g.triangular e g.uniform e NON USARE  uniform
 if(length(UNC)==1){UNC<-rep(UNC,length(lev)) }
 
   labelfac<-fnames
  labellat<-paste("L",LETTERS[1:length(names)],sep="")
  model<-hmmm.model(marg=marg,lev=lev,names=names)
  
  
  descr<-summary(model,printflag=FALSE)
  label<-names
  intcode<-descr[,1]
  intnames<-descr[,2]
  f<-Formula
  if(is.null(Formula)){
    
    
    
    
    f<-list()
    for (i in 1:length(label)){
      f[[i]]<-as.formula((paste("~",paste(paste(labellat[i],labelfac,sep="*"),collapse="+"),"+",labellat[i],"*",label[i],sep="")))
      #paste(intnames[i],"=",
    }
  }
  num.lat<-length(labellat)
  ll<-paste(rep(labellat,each=num.lat),rep(labellat,num.lat),sep="*")
  LL<-matrix(ll,num.lat,num.lat)
  names.lint<-vech(as.matrix(LL[-1,-num.lat]))
  
  length.lint<-length(label)+choose(length(label),2)
  if(length(Formula)==length(label)){
  for (i in (1+length(label)):length.lint){
    ii<-i-length(label)
    if(INDIP[[ii]]=="IND"){f[[i]]="zero"}
    
    else{f[[i]]<-as.formula((paste("~",intnames[i],":","(",names.lint[ii],")",sep="")))}
    
  }
  }
  if(length.lint<length(intnames)){
    for (i in (1+length.lint):length(intnames)){f[[i]]="zero"}}
  names(f)<-intnames
  
  
  fnames<-c(fnames,labellat)
  os<-sum(strata-1)
  strataold<-strata
  strata<-c(strata,rep(2,length(names)))
  str<-prod(strata)
  if(strata[1]==1){strata=strata[-1]}
  
  X<-create.XMAT(model,Formula=f,strata=strata,fnames=fnames,
                 cocacontr=cocacontr,ncocacontr=ncocacontr,replace=FALSE)
  XFULL<-X
  
  nvar<-length(lev)
  l<-length(X)
  c<-dim(X[[l]])[2]
  ########
  px<-0
  npar<-0
  Xnames="zero"
  
  npar<-as.numeric(descr[,"npar"])
  ######
  for(i in 1:nvar){
 
    ###################################################################################################
    if(UNC[[i]]=="POLY"){
      k<-lev[i]
    llg<--(contr.poly(1:(k-1))[,1])
    idn<-grep(labellat[i],colnames(X[[i]]))
    ncol<- dim(X[[i]])[2]-length(idn)+1
    X[[i]]<-X[[i]][,-((ncol-lev[i]+3):ncol)]
    X[[i]][,1]<-(1-X[[i]][,2])*llg
    X[[i]][,3:(ncol-lev[i]+2)]<-X[[i]][,3:(ncol-lev[i]+2)]*X[[i]][,1]
    colnames(X[[i]])<-c("phi1",colnames(X[[i]])[-1])
    }
    if(UNC[[i]]=="sPOLY"){
      k<-lev[i]
      llg<--(contr.poly(1:(k-1))[,1])
      idn<-grep(labellat[i],colnames(X[[i]]))
      ncol<- dim(X[[i]])[2]-length(idn)+1
      
      xxx<-(1-X[[i]][,2]) 
      
      X[[i]]<-X[[i]][,-((ncol-lev[i]+3):ncol)]
      X[[i]][,1]<-(1-X[[i]][,2])*llg
      X[[i]][,3:(ncol-lev[i]+2)]<-X[[i]][,3:(ncol-lev[i]+2)]*X[[i]][,1]
      X[[i]]<-cbind(xxx,X[[i]])
      colnames(X[[i]])<-c("phi0","phi1",colnames(X[[i]])[-(1:2)])
    }
    
    
     if(UNC[[i]]=="TUTZ0"){
      k<-lev[i]
      llg<-rep(-1,k-1)
      if(floor(k/2)==k/2){llg[1:(k/2-1)]<-1
      llg[k/2]=0}
      if(floor(k/2)!=k/2){
        M<-(k+1)/2
        llg[1:(M-1)]<-1
      }
     
      idn<-grep(labellat[i],colnames(X[[i]]))
      ncol<- dim(X[[i]])[2]-length(idn)+1
      X[[i]]<-X[[i]][,-((ncol-lev[i]+3):ncol)]
      X[[i]][,1]<-(1-X[[i]][,2])*llg
      X[[i]][,3:(ncol-lev[i]+2)]<-X[[i]][,3:(ncol-lev[i]+2)]*X[[i]][,1]
      colnames(X[[i]])<-c("phi1",colnames(X[[i]])[-1])
    }
    
    if(UNC[[i]]=="sTUTZ0"){
      k<-lev[i]
      llg<-rep(-1,k-1)
      if(floor(k/2)==k/2){llg[1:(k/2-1)]<-1
      llg[k/2]=0}
      if(floor(k/2)!=k/2){
        M<-(k+1)/2
        llg[1:(M-1)]<-1
      }
      
      idn<-grep(labellat[i],colnames(X[[i]]))
      ncol<- dim(X[[i]])[2]-length(idn)+1
      
      xxx<-(1-X[[i]][,2]) 
      
      X[[i]]<-X[[i]][,-((ncol-lev[i]+3):ncol)]
      X[[i]][,1]<-(1-X[[i]][,2])*llg
      X[[i]][,3:(ncol-lev[i]+2)]<-X[[i]][,3:(ncol-lev[i]+2)]*X[[i]][,1]
      X[[i]]<-cbind(xxx,X[[i]])
      colnames(X[[i]])<-c("phi0","phi1",colnames(X[[i]])[-(1:2)])
    }
    
    
    
    
    
    if(UNC[[i]]=="TUTZ1"){
      k<-lev[i]
      if(floor(k/2)==k){
        m<-k/2
        llg<--((1:(k-1))-m )}
      if(floor(k/2)!=k){ m<-floor(k/2)+1
      llg<--((1:(k-1))-m+0.5) 
      }
      idn<-grep(labellat[i],colnames(X[[i]]))
     ncol<- dim(X[[i]])[2]-length(idn)+1
     X[[i]]<-X[[i]][,-((ncol-lev[i]+3):ncol)]
     X[[i]][,1]<-(1-X[[i]][,2])*llg
     X[[i]][,3:(ncol-lev[i]+2)]<-X[[i]][,3:(ncol-lev[i]+2)]*X[[i]][,1]
     colnames(X[[i]])<-c("phi1",colnames(X[[i]])[-1])
    }
    
    if(UNC[[i]]=="sTUTZ1"){
      k<-lev[i]
      if(floor(k/2)==k){
        m<-k/2
        llg<--((1:(k-1))-m )}
      if(floor(k/2)!=k){ m<-floor(k/2)+1
      llg<--((1:(k-1))-m+0.5) 
      }
      
      idn<-grep(labellat[i],colnames(X[[i]]))
      ncol<- dim(X[[i]])[2]-length(idn)+1
      
      xxx<-(1-X[[i]][,2]) 
      
  
      
      X[[i]]<-X[[i]][,-((ncol-lev[i]+3):ncol)]
      X[[i]][,1]<-(1-X[[i]][,2])*llg
      X[[i]][,3:(ncol-lev[i]+2)]<-X[[i]][,3:(ncol-lev[i]+2)]*X[[i]][,1]
      X[[i]]<-cbind(xxx,X[[i]])
      colnames(X[[i]])<-c("phi0","phi1",colnames(X[[i]])[-(1:2)])
    }
    
    ###################################################################################################
    
    
    
    
    
    
    
    
    
    
    
    
    if(UNC[[i]]=="PARAB"|UNC[[i]]=="TUTZ2"){
      m<-lev[i]
      x<-1:m
      
      p<-6*(x)*(m+1-x)/(m+2)/(m+1)/(m)
      
      llg<-log(p[-1]/p[-length(p)])
      
      idn<-grep(labellat[i],colnames(X[[i]]))
      ncol<- dim(X[[i]])[2]-length(idn)+1
      X[[i]]<-X[[i]][,-((ncol-lev[i]+3):ncol)]
      X[[i]][,1]<-(1-X[[i]][,2])*llg
      X[[i]][,3:(ncol-lev[i]+2)]<-X[[i]][,3:(ncol-lev[i]+2)]*X[[i]][,1]
      colnames(X[[i]])<-c("phi1",colnames(X[[i]])[-1])
    }
    if(UNC[[i]]=="sPARAB"|UNC[[i]]=="sTUTZ2"){
      m<-lev[i]
      x<-1:m
      
      p<-6*(x)*(m+1-x)/(m+2)/(m+1)/(m)
      
      llg<-log(p[-1]/p[-length(p)])
      
      idn<-grep(labellat[i],colnames(X[[i]]))
      ncol<- dim(X[[i]])[2]-length(idn)+1
      
      xxx<-(1-X[[i]][,2]) 
       
  
      
      X[[i]]<-X[[i]][,-((ncol-lev[i]+3):ncol)]
      X[[i]][,1]<-(1-X[[i]][,2])*llg
      X[[i]][,3:(ncol-lev[i]+2)]<-X[[i]][,3:(ncol-lev[i]+2)]*X[[i]][,1]
      X[[i]]<-cbind(xxx,X[[i]])
      colnames(X[[i]])<-c("phi0","phi1",colnames(X[[i]])[-(1:2)])
    }
    
   # ]##################################################################################################
    if(UNC[[i]]=="l.poly"){
      idn<-grep(labellat[i],colnames(X[[i]]))
      X[[i]]<-cbind(X[[i]][,1],X[[i]][,idn])
      k<-lev[i]
    llg<--(contr.poly(1:(k-1))[,1])
    X[[i]][,1]<-(1-X[[i]][,2])*llg
    colnames(X[[i]])<-c("phi1",colnames(X[[i]])[-1])
    }
    
    
    
     if(UNC[[i]]=="l.tutz0"){
      idn<-grep(labellat[i],colnames(X[[i]]))
      X[[i]]<-cbind(X[[i]][,1],X[[i]][,idn])
      k<-lev[i]
      llg<-rep(-1,k-1)
      if(floor(k/2)==k/2){llg[1:(k/2-1)]<-1
      llg[k/2]=0}
      if(floor(k/2)!=k/2){
        M<-(k+1)/2
        llg[1:(M-1)]<-1
      }
      
      X[[i]][,1]<-(1-X[[i]][,2])*llg
      colnames(X[[i]])<-c("phi1",colnames(X[[i]])[-1])
    }
    
    
    
    
    if(UNC[[i]]=="l.tutz1"){
      idn<-grep(labellat[i],colnames(X[[i]]))
      X[[i]]<-cbind(X[[i]][,1],X[[i]][,idn])
      k<-lev[i]
      if(floor(k/2)==k/2){
        m<-k/2
        gloa<--((1:(k-1))-m )}
      if(floor(k/2)!=k/2){ m<-floor(k/2)+1
      gloa<--((1:(k-1))-m+0.5) 
      }
      
      X[[i]][,1]<-(1-X[[i]][,2])*gloa 
      colnames(X[[i]])<-c("phi1",colnames(X[[i]])[-1])
    }
   
    
    
    
     
   
    if(UNC[[i]]=="l.parabolic"){
      idn<-grep(labellat[i],colnames(X[[i]]))
      X[[i]]<-cbind(X[[i]][,1],X[[i]][,idn])
      m<-lev[i]
      x<-1:m
      
      p<-6*(x)*(m+1-x)/(m+2)/(m+1)/(m)
      
      gloa<-log(p[-1]/p[-length(p)])
      
      X[[i]][,1]<-(1-X[[i]][,2])*gloa 
      colnames(X[[i]])<-c("phi1",colnames(X[[i]])[-1])
    }
    if(UNC[[i]]=="b.parabolic"){
      idn<-grep(labellat[i],colnames(X[[i]]))
      X[[i]]<-cbind(X[[i]][,1],X[[i]][,idn])
      m<-lev[i]
      x<-1:m
      
      p<-6*(x)*(m+1-x)/(m+2)/(m+1)/(m)
      
      gloa<-log(p[-1]/p[1])
      
      X[[i]][,1]<-(1-X[[i]][,2])*gloa 
      colnames(X[[i]])<-c("phi1",colnames(X[[i]])[-1])
    } 
    if(UNC[[i]]=="g.parabolic"){
      idn<-grep(labellat[i],colnames(X[[i]]))
      X[[i]]<-cbind(X[[i]][,1],X[[i]][,idn])
      m<-lev[i]
      x<-1:m
      
      P<-x*(x+1)*((m+1)*3-2*x-1)/(m+2)/(m+1)/(m)
      gloa<--log(P[-length(P)]/(1-P[-length(P)]))
      X[[i]][,1]<-(1-X[[i]][,2])*gloa
      colnames(X[[i]])<-c("phi1",colnames(X[[i]])[-1])
      }
    
    if(UNC[[i]]=="g.uniform"){
      idn<-grep(labellat[i],colnames(X[[i]]))
      X[[i]]<-cbind(X[[i]][,1],X[[i]][,idn])
      gloa<-rep(log((1:oa)/rev(1:oa)),str)  
      X[[i]][,1]<-(1-X[[i]][,2])*gloa
      colnames(X[[i]])<-c("phi1",colnames(X[[i]])[-1])
    }
    if(UNC[[i]]=="sl.poly"){
      idn<-grep(labellat[i],colnames(X[[i]]))
      X[[i]]<-cbind(X[[i]][,1],X[[i]][,idn])
      k<-lev[i] 
    llg<--(contr.poly(1:(k-1))[,1])
    xxx<-(1-X[[i]][,2]) 
    X[[i]][,1]<-(1-X[[i]][,2])*llg
    X[[i]]<-cbind(xxx,X[[i]])
    colnames(X[[i]])<-c("phi0","phi1",colnames(X[[i]])[-(1:2)])
    }
    
    
    
    
    if(UNC[[i]]=="sl.tutz0"){
      idn<-grep(labellat[i],colnames(X[[i]]))
      X[[i]]<-cbind(X[[i]][,1],X[[i]][,idn])
      k<-lev[i]
      llg<-rep(-1,k-1)
      if(floor(k/2)==k/2){llg[1:(k/2-1)]<-1
      llg[k/2]=0}
      if(floor(k/2)!=k/2){
        M<-(k+1)/2
        llg[1:(M-1)]<-1
      }
      
      xxx<-(1-X[[i]][,2]) 
      X[[i]][,1]<-(1-X[[i]][,2])*llg
      X[[i]]<-cbind(xxx,X[[i]])
      colnames(X[[i]])<-c("phi0","phi1",colnames(X[[i]])[-(1:2)])
    }
    
    
    
    
    if(UNC[[i]]=="sl.tutz1"){
      idn<-grep(labellat[i],colnames(X[[i]]))
      X[[i]]<-cbind(X[[i]][,1],X[[i]][,idn])
      k<-lev[i]
      if(floor(k/2)==k/2){
        m<-k/2
        gloa<--((1:(k-1))-m )}
      if(floor(k/2)!=k/2){ m<-floor(k/2)+1
      gloa<--((1:(k-1))-m+0.5) 
      }
      
      xxx<-(1-X[[i]][,2]) 
      X[[i]][,1]<-(1-X[[i]][,2])*gloa 
      X[[i]]<-cbind(xxx,X[[i]])
      colnames(X[[i]])<-c("phi0","phi1",colnames(X[[i]])[-(1:2)])
    }
    
    
    
    if(UNC[[i]]=="sl.parabolic"){
      idn<-grep(labellat[i],colnames(X[[i]]))
      X[[i]]<-cbind(X[[i]][,1],X[[i]][,idn])
      m<-lev[i]
      x<-1:m
      
      p<-6*(x)*(m+1-x)/(m+2)/(m+1)/(m)
     
      gloa<-log(p[-1]/p[-length(p)])
      xxx<-(1-X[[i]][,2]) 
      X[[i]][,1]<-(1-X[[i]][,2])*gloa 
     X[[i]]<-cbind(xxx,X[[i]])
     colnames(X[[i]])<-c("phi0","phi1",colnames(X[[i]])[-(1:2)])
    }
    if(UNC[[i]]=="sb.parabolic"){
      idn<-grep(labellat[i],colnames(X[[i]]))
      X[[i]]<-cbind(X[[i]][,1],X[[i]][,idn])
      m<-lev[i]
      x<-1:m
      
      p<-6*(x)*(m+1-x)/(m+2)/(m+1)/(m)
     
      gloa<-log(p[-1]/p[1])
      
      xxx<-(1-X[[i]][,2]) 
      X[[i]][,1]<-(1-X[[i]][,2])*gloa 
      X[[i]]<-cbind(xxx,X[[i]])
      colnames(X[[i]])<-c("phi0","phi1",colnames(X[[i]])[-(1:2)])
    } 
    if(UNC[[i]]=="sg.parabolic"){
      idn<-grep(labellat[i],colnames(X[[i]]))
      X[[i]]<-cbind(X[[i]][,1],X[[i]][,idn])
      m<-lev[i]
      x<-1:m
     
      P<-x*(x+1)*((m+1)*3-2*x-1)/(m+2)/(m+1)/(m)
      gloa<--log(P[-length(P)]/(1-P[-length(P)]))
      xxx<-(1-X[[i]][,2]) 
      X[[i]][,1]<-(1-X[[i]][,2])*gloa 
      X[[i]]<-cbind(xxx,X[[i]])
      colnames(X[[i]])<-c("phi0","phi1",colnames(X[[i]])[-(1:2)])
    }
    if(UNC[[i]]=="sg.uniform"){
      idn<-grep(labellat[i],colnames(X[[i]]))
      X[[i]]<-cbind(X[[i]][,1],X[[i]][,idn])
      gloa<-rep(log((1:oa)/rev(1:oa)),str)  
      xxx<-(1-X[[i]][,2]) 
      X[[i]][,1]<-(1-X[[i]][,2])*gloa 
      X[[i]]<-cbind(xxx,X[[i]])
      colnames(X[[i]])<-c("phi0","phi1",colnames(X[[i]])[-(1:2)])
    }
    
    if(UNC[[i]]=="sl.uniform"){
      idn<-grep(labellat[i],colnames(X[[i]]))
      X[[i]]<-cbind(X[[i]][,1],X[[i]][,idn])
      X[[i]][,1]<-(1-X[[i]][,2])
      colnames(X[[i]])<-c("phi0",colnames(X[[i]])[-1])
      
    }
      if(UNC[[i]]=="uniform") {
        idn<-grep(labellat[i],colnames(X[[i]]))
        X[[i]]<-cbind(X[[i]][,1],X[[i]][,idn])
        X[[i]]<-X[[i]][,-1]
      
      }
    #
   
      px<-c(px,dim(X[[i]])[2])
    Xnames<-c(Xnames,colnames(X[[i]]))
    #
   
  }
  if(length(Formula)==length(label)){  
  for(i in (nvar+1):length.lint){
    ii<-i-length(label)
    if (f[[i]]!="zero"){
      c<-dim(X[[i]])[1]
      c1<-c/str
      
      X[[i]]<-X[[i]][,(dim(X[[i]])[2]-c1+1):dim(X[[i]])[2]]
      if(INDIP[[ii]]=="CST"){
        na<-colnames(X[[i]])[1]
        X[[i]]<-matrix(rowSums(X[[i]]))
        colnames(X[[i]])<-na}
      #
      
      
      if(INDIP[[ii]]=="SYM"){
       
        na<-colnames(X[[i]])
        sdim<-dim(X[[i]])[2]
        X[[i]]<-X[[i]]%*%D.matrix(sqrt(sdim))
        colnames(X[[i]])<-vech(matrix(na,sqrt(sdim),sqrt(sdim)))
       
      }
      
      if(class(INDIP[[ii]])=="formula"){
        print("associazion formula used")
        va<-intcode[i]
        va<-as.numeric(unlist(strsplit(va,split=character(0))))
        lev1<-lev[va[1]]-1
        lev2<-lev[va[2]]-1
        x<-factor(rep(1:lev1,lev2))
        xx<-factor(rep(1:lev2,each=lev1))
        da<-data.frame(x,xx)
        colnames(da)<-c("R","C")
       disegno<- model.matrix(INDIP[[ii]],da)
        X[[i]]<-
          matrix(rep(c(0,1),c(str-prod(strataold),prod(strataold))))%x%
          disegno
        colnames(X[[i]])<-colnames(disegno)
      }
      px<-c(px,dim(X[[i]])[2])
      
      
      
      
      Xnames<-c(Xnames,colnames(X[[i]]))
      #
    }
  }
  }
  else{
    lbl<-paste(labellat,"1",sep="")
    ll<-paste(rep(lbl,each=num.lat),rep(lbl,num.lat),sep=":")
    LL<-matrix(ll,num.lat,num.lat)
    names.lint2<-vech(as.matrix(LL[-1,-num.lat]))
    for(i in (nvar+1):length.lint){
      ii<-i-length(label)
     if(X[[i]]!="zero"){
      idn<-grep(names.lint2[ii],colnames(X[[i]]))
      npro<-colnames(X[[i]])[idn]
      X[[i]]<-X[[i]][,idn]
      if(!is.matrix(X[[i]])) { X[[i]]<-matrix(X[[i]])
       colnames(X[[i]])<-npro}
      px<-c(px,dim(X[[i]])[2])
      Xnames<-c(Xnames,colnames(X[[i]]))}
    }
  }
  
  px<-px[-1]
  Xnames<-Xnames[-1]
  XX<-matrix(0,sum(npar)*prod(strata),sum(px))
  or<-rep(rep((1:dim(descr)[1]),times=npar),prod(strata))
  pstart<-0
  for(i in 1:dim(descr)[1]){
    if (f[[i]]!="zero"){
      XX[or==i,(pstart+1):(pstart+px[i])]<-X[[i]]
      pstart<-pstart+px[i]}
  }
  
  colnames(XX)<-Xnames
  
  if(replace){
    model<-hmmm.model(marg=marg,lev=lev,names=names,strata=str,X=XX)
    model$levcov<-strataold
    model$fnames<-fnames
    model$marglist<-marg
    model
  }
  
  else{list(X=X,XFULL=XFULL,XX=XX,Formula=Formula)}
}



######################################

#########
cocadise<-function(Z=NULL,names=NULL,formula=NULL,lev)
{
  if (is.null(formula))
  {
    ncz<-dim(Z)[2]
    zi<-Z[,1]
    zi<-zi[zi>0]
    DD<-diag(c(zi))[,2:length(zi)] #######
    DD<-as.matrix(DD)
    zi<-zi[-1]
    DI<-cbind(-zi,diag(c(zi)))
    if( ncz>1){
      for(i in 2:ncz){
        zi<-Z[,i]
        zi<-zi[zi>0]
        DMATI<-diag(c(zi))[,2:length(zi)]
        DMATI<-as.matrix(DMATI)########
        zi<-zi[-1]
        CMATI<-cbind(-zi,diag(c(zi)))
        DD<-direct.sum(DD,DMATI)
        DI<-direct.sum(DI,CMATI)
      }
    }
    matrici<-list(IMAT=DI,DMAT=DD)
  }
  else
    LDMatrix(lev,formula,names)
  
}

  hmmm_equiv<-function(f,marg,names,modelobs){
                        lev<-prod(modelobs$modello$livelli)
                        strata<-modelobs$modello$strata
                        nlat<-length(modelobs$modello$livelli)
                        stratalat<-2^nlat
                        stratax<-strata/stratalat
        
 XX<-array(c(modelobs$matrici$X),dim=c(lev-1,stratax,stratalat,dim(modelobs$matrici$X)[2]))
XX<-XX[,,stratalat,]
XX<-matrix(c(XX),(lev-1)*stratax,dim(modelobs$matrici$X)[2])
model<-hmmm.model(marg=marg,lev=modelobs$modello$livelli,names=names,strata=stratax,X=XX)
hmmmfit<-hmmm.mlfit(f,
                       model=model, 
                       maxit =60,   
                       chscore.criterion=1,  
                       y.eps=0.1,
                       norm.diff.conv = 1e-06, 
                       norm.score.conv = 1e-06        ) 
  }  
  

  #################################################
  ###############################################################
  ###############################################################
  ###############################################################
  ##################################################################à
  ############################################  
  hmmm.model.T<-function(marg,lev,names,Formula=NULL,strata=1,fnames=NULL,
                         cocacontr=NULL,ncocacontr=NULL,replace=TRUE,INDIP="CST",UNC="TUTZ2",MOD="A"){
    #alternative a "uniform" per UNC sono "g.uniform","b.parabolic""g.parabolic","l.parabolic" 
    #"sg.uniform","sb.parabolic""sg.parabolic","l.parabolic"
    #da usare solo
    #con i corrispondenti logits g b l per le marginali delle osservabili
    #con logit global usare  g.triangular e g.uniform e NON USARE  uniform
    if(length(UNC)==1){UNC<-rep(UNC,length(lev)) }
    
    labelfac<-fnames
    labellat<-paste("L",LETTERS[1:length(names)],sep="")
    #model<-hmmm.model(marg=marg,lev=lev,names=names)
    model<-hmmm.model(marg=marg,lev=lev,names=names)
    
    descr<-summary(model,printflag=FALSE)
    label<-names
    intcode<-descr[,1]
    intnames<-descr[,2]
    f<-Formula
    
    
    
    num.lat<-length(labellat)
    label1<-paste(labellat,"1",sep="")
    llcst<-paste(rep(label1,each=num.lat),rep(label1,num.lat),sep=":")
    llcst<-matrix(llcst,num.lat,num.lat)
    llcst<-vech(as.matrix(llcst[-1,-num.lat]))
    
    
    
    
    
    if(is.null(Formula)){
      
      
      
      
      f<-list()
      for (i in 1:length(label)){
        f[[i]]<-as.formula((paste("~",paste(paste(labellat[i],labelfac,sep="*"),collapse="+"),"+",labellat[i],"*",label[i],sep="")))
        #paste(intnames[i],"=",
      }
    }
    num.lat<-length(labellat)
    ll<-paste(rep(labellat,each=num.lat),rep(labellat,num.lat),sep="*")
    LL<-matrix(ll,num.lat,num.lat)
    names.lint<-vech(as.matrix(LL[-1,-num.lat]))
    
    length.lint<-length(label)+choose(length(label),2)
    if(  (length(Formula)==length(label))&(length(label)<length.lint) ){#new
      for (i in (1+length(label)):length.lint){
        ii<-i-length(label)
        if(INDIP[[ii]]=="IND"){f[[i]]="zero"}
        
        else{f[[i]]<-as.formula((paste("~",intnames[i],":","(",names.lint[ii],")",sep="")))}
        
      }
    }
    if(length.lint<length(intnames)){
      for (i in (1+length.lint):length(intnames)){f[[i]]="zero"}}
    names(f)<-intnames
    #Formula<-f
    
    fnames<-c(fnames,labellat)
    os<-sum(strata-1)
    strataold<-strata
    strata<-c(strata,rep(2,length(names)))
    str<-prod(strata)
    if(strata[1]==1){strata=strata[-1]}
    
    X<-create.XMAT(model,Formula=f,strata=strata,fnames=fnames,
                   cocacontr=cocacontr,ncocacontr=ncocacontr,replace=FALSE)
    XFULL<-X
    
    nvar<-length(lev)
    l<-length(X)
    c<-dim(X[[l]])[2]
    ########
    px<-0
    npar<-0
    Xnames="zero"
    
    npar<-as.numeric(descr[,"npar"])
    ######
    for(i in 1:nvar){
    
      ###################################################################################################
      if(UNC[[i]]=="TUTZ0"){
        k<-lev[i]
        llg<-rep(-1,k-1)
        if(floor(k/2)==k/2){llg[1:(k/2-1)]<-1
        llg[k/2]=0}
        if(floor(k/2)!=k/2){
          M<-(k+1)/2
          llg[1:(M-1)]<-1
        }
        
        idn<-grep(labellat[i],colnames(X[[i]]))
        ncol<- dim(X[[i]])[2]-length(idn)+1
        X[[i]]<-X[[i]][,-((ncol-lev[i]+3):ncol)]
        X[[i]][,1]<-(1-X[[i]][,2])*llg
        X[[i]][,3:(ncol-lev[i]+2)]<-X[[i]][,3:(ncol-lev[i]+2)]*X[[i]][,1]
        LA<-X[[i]][,2]
        X[[i]][LA==0,]<-X[[i]][LA==0,]+X[[i]][LA==1,]
        #LA->X[[i]][,2]
        if(MOD=="B"){
          X[[i]][,(dim(X[[i]])[2]-k+2):dim(X[[i]])[2]]<-
            X[[i]][,(dim(X[[i]])[2]-k+2):dim(X[[i]])[2]]*LA}
        colnames(X[[i]])<-c("phi1",colnames(X[[i]])[-1])
      }
      
      
      
      
      
      
      
      if(UNC[[i]]=="TUTZ1"){
        k<-lev[i]
        if(floor(k/2)==k/2){
          m<-k/2
          llg<--((1:(k-1))-m )}
        if(floor(k/2)!=k/2){ m<-floor(k/2)+1
        llg<--((1:(k-1))-m+0.5) 
        }
        idn<-grep(labellat[i],colnames(X[[i]]))
        ncol<- dim(X[[i]])[2]-length(idn)+1
        X[[i]]<-X[[i]][,-((ncol-lev[i]+3):ncol)]
        X[[i]][,1]<-(1-X[[i]][,2])*llg
        X[[i]][,3:(ncol-lev[i]+2)]<-X[[i]][,3:(ncol-lev[i]+2)]*X[[i]][,1]
        LA<-X[[i]][,2]
        X[[i]][LA==0,]<-X[[i]][LA==0,]+X[[i]][LA==1,]
        
        #####
        if(MOD=="B"){
        X[[i]][,(dim(X[[i]])[2]-k+2):dim(X[[i]])[2]]<-
          X[[i]][,(dim(X[[i]])[2]-k+2):dim(X[[i]])[2]]*LA}
        #LA->X[[i]][,2]
        colnames(X[[i]])<-c("phi1",colnames(X[[i]])[-1])
      }
      
      if(UNC[[i]]=="TUTZ2"){
        m<-lev[i]
        x<-1:m
        
        p<-6*(x)*(m+1-x)/(m+2)/(m+1)/(m)
        
        llg<-log(p[-1]/p[-length(p)])
        
        idn<-grep(labellat[i],colnames(X[[i]]))
        ncol<- dim(X[[i]])[2]-length(idn)+1
        X[[i]]<-X[[i]][,-((ncol-lev[i]+3):ncol)]
        X[[i]][,1]<-(1-X[[i]][,2])*llg
        X[[i]][,3:(ncol-lev[i]+2)]<-X[[i]][,3:(ncol-lev[i]+2)]*X[[i]][,1]
        LA<-X[[i]][,2]
        X[[i]][LA==0,]<-X[[i]][LA==0,]+X[[i]][LA==1,]
        #LA->X[[i]][,2]
        if(MOD=="B"){
          X[[i]][,(dim(X[[i]])[2]-k+2):dim(X[[i]])[2]]<-
            X[[i]][,(dim(X[[i]])[2]-k+2):dim(X[[i]])[2]]*LA}
        colnames(X[[i]])<-c("phi1",colnames(X[[i]])[-1])
      }
     
      # ]##################################################################################################
      
     
      
      px<-c(px,dim(X[[i]])[2])
      Xnames<-c(Xnames,colnames(X[[i]]))
      #
      
    }
    if(length(names)>1){
    if(length(Formula)==length(label)){  
      for(i in (nvar+1):length.lint){
        ii<-i-length(label)
        if (f[[i]]!="zero"){
          c<-dim(X[[i]])[1]
          c1<-c/str
          
   ######################################################################################################
          #####################################################################
         
          if(INDIP[[ii]]=="CST4"){
           # print("association given uncertainty")
           # na<-colnames(X[[i]])[1]
            idnLA1LB1<-grep(llcst[ii],colnames(X[[i]]))
            splitlab<-strsplit(llcst[ii],":")
            idnLA1<-grep(splitlab[[1]][1],colnames(X[[i]]))
            idnLB1<-grep(splitlab[[1]][2],colnames(X[[i]]))
            #idnLA1LB1<-grep(llcst[ii],colnames(X[[i]]))
            idnLA2<-grep(gsub("1","2",splitlab[[1]][1]),colnames(X[[i]]))
           
            X[[i]]<-cbind(
              rowSums(X[[i]][,idnLA1])-rowSums(X[[i]][,idnLB1])+rowSums(X[[i]][,idnLA1LB1]),
              rowSums(X[[i]][,idnLA2])-rowSums(X[[i]][,idnLA1LB1]),
              rowSums(X[[i]][,idnLB1])-2*rowSums(X[[i]][,idnLA1LB1]),
              rowSums(X[[i]][,idnLA1LB1]))
              
            colnames(X[[i]])<-c("ass00","ass10","ass01","ass11")}
          #
          
          if(INDIP[[ii]]=="CST3"){
            
            idnLA1LB1<-grep(llcst[ii],colnames(X[[i]]))
            splitlab<-strsplit(llcst[ii],":")
            idnLA1<-grep(splitlab[[1]][1],colnames(X[[i]]))
            idnLB1<-grep(splitlab[[1]][2],colnames(X[[i]]))
            #idnLA1LB1<-grep(llcst[ii],colnames(X[[i]]))
            idnLA2<-grep(gsub("1","2",splitlab[[1]][1]),colnames(X[[i]]))
            
            
            X[[i]]<-cbind(
              rowSums(X[[i]][,idnLA1])-rowSums(X[[i]][,idnLB1])+rowSums(X[[i]][,idnLA1LB1]),
              (rowSums(X[[i]][,idnLA2])-rowSums(X[[i]][,idnLA1LB1])+
              rowSums(X[[i]][,idnLB1])-2*rowSums(X[[i]][,idnLA1LB1])),
              rowSums(X[[i]][,idnLA1LB1]))
            #,
            #rowSums(X[[i]][,2:(c1+1)]),
            #rowSums(X[[i]][,(c1+2):(2*c1+1)]),
            #rowSums(X[[i]][,(2*c1+2):(3*c1+1)]),
            #rowSums(X[[i]][,(dim(X[[i]])[2]-c1+1):dim(X[[i]])[2]])
            
            #X[[i]]<-cbind(matrix(rowSums(X[[i]])),1-matrix(rowSums(X[[i]])))
            colnames(X[[i]])<-c("ass00","ass10-01","ass11")}
          #
          
          
          
          
          if(INDIP[[ii]]=="CST2"){
            #print("association given uncertainty")
           
            idnLA1LB1<-grep(llcst[ii],colnames(X[[i]]))
           
            
            
              X[[i]]<-cbind(
                X[[i]][,1]-rowSums(X[[i]][,idnLA1LB1]),
                rowSums(X[[i]][,idnLA1LB1]))
            #,
            #rowSums(X[[i]][,2:(c1+1)]),
            #rowSums(X[[i]][,(c1+2):(2*c1+1)]),
            #rowSums(X[[i]][,(2*c1+2):(3*c1+1)]),
            #rowSums(X[[i]][,(dim(X[[i]])[2]-c1+1):dim(X[[i]])[2]])
            
            #X[[i]]<-cbind(matrix(rowSums(X[[i]])),1-matrix(rowSums(X[[i]])))
            colnames(X[[i]])<-c("assU","assA")}
          #
          
          
          if(INDIP[[ii]]=="CST1"){
            
            idnLA1LB1<-grep("LA1:LB1",colnames(X[[i]]))
            
            
            X[[i]]<-
              matrix(rowSums(X[[i]][,idnLA1LB1]))
            
            colnames(X[[i]])<-"assA"}
          #
          
          
         
          
          if(class(INDIP[[ii]])=="formula"){
            print("association formula used")
            va<-intcode[i]
            va<-as.numeric(unlist(strsplit(va,split=character(0))))
            lev1<-lev[va[1]]-1
            lev2<-lev[va[2]]-1
            x<-factor(rep(1:lev1,lev2))
            xx<-factor(rep(1:lev2,each=lev1))
            da<-data.frame(x,xx)
            colnames(da)<-c("R","C")
            disegno<- model.matrix(INDIP[[ii]],da)
            X[[i]]<-
              matrix(rep(c(0,1),c(str-prod(strataold),prod(strataold))))%x%
              disegno
            colnames(X[[i]])<-colnames(disegno)
          }
          px<-c(px,dim(X[[i]])[2])
          
          
          
          
          Xnames<-c(Xnames,colnames(X[[i]]))
          #
        }
      }
    }
    else{
      lbl<-paste(labellat,"1",sep="")
      ll<-paste(rep(lbl,each=num.lat),rep(lbl,num.lat),sep=":")
      LL<-matrix(ll,num.lat,num.lat)
      names.lint2<-vech(as.matrix(LL[-1,-num.lat]))
      for(i in (nvar+1):length.lint){
        ii<-i-length(label)
        if(any(X[[i]]!="zero")){  
          idn<-grep(names.lint2[ii],colnames(X[[i]]))
          npro<-colnames(X[[i]])[idn]
          X[[i]]<-X[[i]][,idn]
          if(!is.matrix(X[[i]])) { X[[i]]<-matrix(X[[i]])
          colnames(X[[i]])<-npro}
          px<-c(px,dim(X[[i]])[2])
          Xnames<-c(Xnames,colnames(X[[i]]))}
      }
    }#fineasselse
    }#new fine parte ass
    px<-px[-1]
    Xnames<-Xnames[-1]
    XX<-matrix(0,sum(npar)*prod(strata),sum(px))
    or<-rep(rep((1:dim(descr)[1]),times=npar),prod(strata))
    pstart<-0
    for(i in 1:dim(descr)[1]){
      if (f[[i]]!="zero"){
        XX[or==i,(pstart+1):(pstart+px[i])]<-X[[i]]
        pstart<-pstart+px[i]}
    }
    
    colnames(XX)<-Xnames
    
    if(replace){
      model<-hmmm.model(marg=marg,lev=lev,names=names,strata=str,X=XX)
      model$levcov<-strataold
      model$fnames<-fnames
      model$marglist<-marg
      model
    }
    
    else{list(X=X,XFULL=XFULL,XX=XX,Formula=Formula)}
  }
  
  
  
  ######################################
  
  
