#!/usr/bin/env Rscript
# Shared numeric covariates and participant selection for every external tool.
suppressPackageStartupMessages({library(data.table); library(jsonlite)})
arg <- function(k, default=NULL) {a<-commandArgs(TRUE); i<-match(k,a); if(is.na(i)) default else a[i+1L]}
repo <- '${ADSP_ROOT}'
out <- arg('--out'); arm <- arg('--arm'); yfile <- arg('--ystar')
stopifnot(arm %in% c('retained_additive','retained_joint','full_joint_diagnostic'))
pool <- fread(file.path(repo,'outputs/plasmode/plasmode_pool.csv'), colClasses=list(character='IID'))
y <- fread(yfile, colClasses=list(character='IID'))
stopifnot(identical(pool$IID,y$IID), all(y$Y %in% c(0,1)), !anyNA(pool$Study), !anyNA(pool$center))
pool[, Y:=as.integer(y$Y)]
pool[, unit:=paste(Study,center,sep='|')]
units <- pool[,.(n=.N, cases=sum(Y), controls=sum(Y==0)),by=unit]
units[, informative:=cases>0 & controls>0]
an <- if(arm=='full_joint_diagnostic') rep(TRUE,nrow(pool)) else pool$unit %in% units[informative==TRUE,unit]
meta <- list(arm=arm, n_pool=nrow(pool), n_analyzed=sum(an), n_cases=sum(pool$Y[an]),
             n_units=nrow(units), n_informative_units=sum(units$informative),
             n_single_outcome_units=sum(!units$informative), n_single_outcome_people=sum(units[!informative,n]),
             project_collapse='none', phenotype_source=yfile, accepted_design=FALSE)
fwrite(units,paste0(out,'_units.tsv'),sep='\t')
fwrite(data.table(IID=pool$IID,keep=as.integer(an)),paste0(out,'_selection.tsv'),sep='\t')
if(arm=='full_joint_diagnostic' && any(!units$informative)) {
  meta$status <- 'structurally_separated'
  meta$reason <- 'Full project-by-center fixed effects include units with only cases or only controls. A finite unpenalized logistic maximum-likelihood solution does not exist. No association fits or scans are attempted.'
  writeLines(toJSON(meta,auto_unbox=TRUE,pretty=TRUE),paste0(out,'_meta.json'))
  quit(status=0)
}
pcs <- paste0('PC',1:20)
zs <- function(v) {v<-as.numeric(v); v[is.na(v)]<-mean(v[an],na.rm=TRUE); s<-sd(v[an]); stopifnot(is.finite(s),s>0); (v-mean(v[an]))/s}
pool[,Age_f:=zs(Age)]; pool[,Sex_f:=as.numeric(Sex)]
stopifnot(!anyNA(pool$Sex_f))
for(p in pcs) set(pool,j=p,value=zs(pool[[p]]))
pool[,center_f:=factor(center)]; pool[,Study_f:=factor(Study)]; pool[,unit_f:=factor(unit)]
rhs <- as.formula(paste('~ Age_f + Sex_f +',paste(pcs,collapse=' + '),
                       if(arm=='retained_additive') '+ center_f + Study_f' else '+ unit_f'))
# Drop unused levels first. Keep only QR-independent columns, without changing
# the span of the specified fixed-effects model.
dat <- droplevels(as.data.frame(pool[an]))
Xall <- model.matrix(rhs,dat)
q <- qr(Xall); keep <- sort(q$pivot[seq_len(q$rank)])
stopifnot(1L %in% keep)
X <- Xall[,keep,drop=FALSE]
cn <- make.names(colnames(X)[-1L],unique=TRUE)
colnames(X) <- c('(Intercept)',cn)
fitwarn <- character()
fit <- withCallingHandlers(glm.fit(X,dat$Y,family=binomial(),control=glm.control(maxit=100,epsilon=1e-8)),
                           warning=function(w){fitwarn<<-c(fitwarn,conditionMessage(w));invokeRestart('muffleWarning')})
mu<-fit$fitted.values
valid<-isTRUE(fit$converged) && all(is.finite(fit$coefficients)) && all(is.finite(mu)) && all(mu>0 & mu<1)
meta$n_design_cols<-ncol(Xall);meta$design_rank<-q$rank
meta$dropped_columns<-setdiff(colnames(Xall),colnames(Xall)[keep])
meta$covariate_names<-cn;meta$glm_preflight_converged<-isTRUE(fit$converged)
meta$glm_preflight_warnings<-fitwarn;meta$glm_preflight_max_abs_coef<-max(abs(fit$coefficients))
meta$glm_preflight_mu_range<-range(mu)
meta$glm_preflight_information_rcond<-rcond(crossprod(X,X*as.numeric(mu*(1-mu))))
meta$accepted_design<-valid;meta$status<-if(valid) 'design_ready' else 'invalid_fixed_effects_preflight'
writeLines(toJSON(meta,auto_unbox=TRUE,pretty=TRUE,digits=NA),paste0(out,'_meta.json'))
if(!valid) quit(status=0)
covdf<-as.data.frame(X[,-1L,drop=FALSE],check.names=FALSE)
ph<-cbind(data.frame(IID=dat$IID,Y=dat$Y),covdf)
fullcov<-as.data.frame(matrix(0,nrow(pool),length(cn)));colnames(fullcov)<-cn
fullcov[an,]<-covdf
fullY<-rep(NA_real_,nrow(pool));fullY[an]<-pool$Y[an]
saveRDS(list(ph=ph,covdf=fullcov,cn=cn,an=an,Y=fullY,meta=meta),paste0(out,'.rds'))
fwrite(ph,paste0(out,'_phenocovar.tsv'),sep='\t')
