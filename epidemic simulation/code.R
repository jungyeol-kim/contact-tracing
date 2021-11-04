## Required packages
list.of.packages <- c("data.table","matrixStats", "reshape2","dplyr","igraph","fastmatch")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)

# define fast search function
`%fin%` <- function(x, table) {
  fmatch(x, table, nomatch = 0L) > 0L
}

#### Model Parameters ####
latency.period<-2
f.Ip_pre.to.Ip<-function(Ip_pre.list,t){
  sel<-rbinom(length(Ip_pre.list),1,1./latency.period)
  status<-ifelse(sel==1,'Ip','Ip-L')
  return(status)
}
latency.period<-2
f.Ia_pre.to.Ia<-function(Ia_pre.list,t){
  sel<-rbinom(length(Ia_pre.list),1,1./latency.period)
  status<-ifelse(sel==1,'Ia','Ia-L')
  return(status)
}
incubation.period<-5
post.latency.period<-incubation.period-latency.period
f.Ip.to.Is<-function(Ip.list,t){
  sel<-rbinom(length(Ip.list),1,1./post.latency.period)
  status<-ifelse(sel==1,'Is','Ip')
  return(status)
}
wait.see.period<-4
f.Is.to.RTs<-function(Is.list,t){
  sel<-rbinom(length(Is.list),1,1./wait.see.period)
  status<-ifelse(sel==1,'RT','Is')
  return(status)
}
rec.death.period<-14-wait.see.period
death.rate<-0.0065
f.RTs.to.D<-function(RT.list,t){
  status<-base::sample(c('R', 'D', 'RT'), size = length(RT.list), replace = TRUE, prob = c((1./rec.death.period)*(1.-death.rate), (1./rec.death.period)*death.rate, (1.-1./rec.death.period)))
  return(status)
}
infectious.Ia.period<-7
f.Ia.to.R<-function(Ia.list,t){
  sel<-rbinom(length(Ia.list),1,1./infectious.Ia.period)
  status<-ifelse(sel==1,'R','Ia')
  return(status)
}
#### Contact Tracing Function ####
contact.tracing<-function(hops,TP.list.for.CT,edges.ct.all){
  multihop.temp<-vector(); multihop<-vector();
  for(i in 1:hops){
    temp<-if(i==1) TP.list.for.CT else multihop.temp
    multihop.temp<-unique(edges.ct.all[edges.ct.all$node1 %fin% temp]$node2)
    multihop<-unique(c(multihop,multihop.temp))
  }
  multihop<-multihop[!multihop %fin% TP.list.for.CT]
  return(multihop)
}


#### Simulations ####
intervention.start.inf.rate<-0.1 # percent (%)
turnaround<-1  # Turnaround times for coronavirus test results (days)
FP.rate<-0 # percent (%)
FN.rate<-11 # percent (%)
cooperativity<-100  # cooperativity (%)
#infness.of.asymp<-0.5 # infectiousness of asymptomatic individuals relative to symptomatic

mainDir.0<-"/Users/Shared/Dropbox/1. PENN - Doctorate/Research/COVID-19/Covid19-2 Non Shared/Multihop Contact Tracing/Updates/08-2021/Code & Results/WS"
# mainDir.1<-file.path(mainDir.0, paste0('Inf asym ',infness.of.asymp))
# dir.create(mainDir.1)
mainDir.2<-file.path(mainDir.0, paste0('intervention from ',intervention.start.inf.rate,'% infections'))
dir.create(mainDir.2)
mainDir.3<-file.path(mainDir.2, paste0('Turnaround ',turnaround,' day'))
dir.create(mainDir.3)
mainDir.4<-file.path(mainDir.3, paste0('FP ',FP.rate,'%'))
dir.create(mainDir.4)
mainDir.5<-file.path(mainDir.4, paste0('FN ',FN.rate,'%'))
dir.create(mainDir.5)
mainDir.6<-file.path(mainDir.5, paste0('Cooperativity (Scenario 1) ',cooperativity,'%'))
dir.create(mainDir.6)





print(system.time({
  for(aaaa in c(4)){
    for(bbbb in c(0.01,0.1,1)){
      for(cccc in c(0.1,0.15,0.2,0.25,0.3)){
        topology<-paste0(mainDir.0,"/ws_network_p_",bbbb,"_k",aaaa,".RData")
        edges.all.undir.final<-get(load(topology))
        
        subDir<-paste0(aaaa," neighbors : p=",bbbb," : PoI ",cccc)
        dir.create(file.path(mainDir.6, subDir))
        setwd(file.path(mainDir.6, subDir))
        
        N<-edges.all.undir.final[[3]] # Total number of nodes in erdos-renyi and scale-free networks
        T.max = edges.all.undir.final[[2]] # total number of days
        test.x.days.later<-3
        
        prob.inf.from.symp<-cccc # Probability with which a symptomatic individual infects a susceptible in an interaction
        prob.inf.from.presymp<-prob.inf.from.symp # Probability with which a presymptomatic individual infects a susceptible in an interaction
        infness.of.asymp<-1 # infectiousness of asymptomatic individuals relative to symptomatic
        prob.inf.from.asymp<-prob.inf.from.symp*infness.of.asymp # Probability with which a asymptomatic individual infects a susceptible in an interaction
        prop.asymp<-0.4 # proportion of asymptomatic individuals among infected people.

        I0<-3 # the number of initially infected individuals
        symptom<-c('Is','RT') # compartments with symptoms
        nosymptom<-c('S', 'Ia-L', 'Ip-L', 'Ip','Ia','R') # compartments without symptoms
        
        edge.list.undir<-vector();
        edge.list.undir<-copy(edges.all.undir.final[[1]][,c(1,2)])
        edge.list.undir$edge_num<-1:nrow(edge.list.undir)
        
        epidemic<-function(iteration,hhh){
          hops<-hhh
          
          wdata.new <- rep(NA,N)
          wdata.ct.new <- data.table('1'=rep(NA,N))
          wdata.tp.new <- data.table('1'=rep(NA,N))
          
          wdata.qt.temp <- rep(Inf,N)
          wdata.qt.actual.temp <- rep(Inf,N)
          wdata.tp.temp <- rep(0,N)
          wdata.ct.day <- rep(0,T.max) 
          wdata.ct.day.cum <- 0
          
          wdata.inf.count<-rep(0,T.max);
          wdata.test.count<-rep(0,T.max);
          wdata.qt.count<-rep(0,T.max);
          
          ### initial condition
          I0.node<-base::sample.int(N, I0)
          symptomatic<-I0.node[1]; presymptomatic<-I0.node[2]; asymptomatic<-I0.node[3];
          susceptible<-c(1:N)[-I0.node];
          wdata.inf.count[1]<-I0; 
          
          
          if(length(susceptible)>0){wdata.new[susceptible]<-'S'}
          if(length(presymptomatic)>0){wdata.new[presymptomatic]<-'Ip'}
          if(length(symptomatic)>0){wdata.new[symptomatic]<-'Is'}
          if(length(asymptomatic)>0){wdata.new[asymptomatic]<-'Ia'}
          
          temp.num<-round(N*(1.-cooperativity/100.))
          non.cooperative.nodes<-if(temp.num==0){vector()}else{sample(1:N, temp.num, replace = F)}
          cooperative.nodes<-if(length(non.cooperative.nodes)==0){1:N}else{c(1:N)[-non.cooperative.nodes]}
          
          edges.all.final.transmission<-vector(mode='list',length=T.max)
          remove.edges.HA<-vector(mode='list',length=T.max)
          remove.edges.real<-vector(mode='list',length=T.max)
          edges.all.final.transmission[[1]]<-rbindlist(list(edge.list.undir[,c(1,2)],edge.list.undir[,c(2,1)]),use.names=FALSE)
          
          test.rt.all<-vector()
          
          set.seed(iteration)
          for(t in 2:T.max){
            wdata.old<-wdata.new
            #### Epidemic process; status update from day t-1 to day t ####
            S.list<-vector();
            temp1<-which(wdata.old %fin% c('S'))
            if(length(temp1)>0){
              wdata.new[temp1]<-'S'
              S.list<-data.table(node1=as.integer(temp1), status='S')
            }
            
            I.list<-vector();
            temp2<-which(wdata.old %fin% c('Ip','Is','Ia'))
            if(length(temp1)>0 & length(temp2)>0){
              I.list<-wdata.old[temp2]
              I.list<-data.table(node2=as.integer(temp2), status=I.list)
              
              neighbor.of.I<-vector(); inf.list<-vector();
              neighbor.of.I<-edges.all.final.transmission[[t-1]][node1 %fin% S.list$node1 & node2 %fin% I.list$node2]
              if(nrow(neighbor.of.I)>0){
                neighbor.of.I<-I.list[neighbor.of.I, on="node2"]  
                neighbor.of.I<-neighbor.of.I[,.(count = .N), by=.(node1,status)]
                neighbor.of.I$status<-factor(neighbor.of.I$status, levels=c("Ip", "Is", "Ia"))
                neighbor.of.I<-setDT(dcast(neighbor.of.I, node1 ~ status, fun = sum, value.var = "count", drop=F))
                neighbor.of.I[,prob.inf:=1-((1-prob.inf.from.presymp)^Ip)*((1-prob.inf.from.symp)^Is)*((1-prob.inf.from.asymp)^Ia)]
                neighbor.of.I[,new.inf:=rbinom(length(neighbor.of.I$prob.inf), size = 1, prob=neighbor.of.I$prob.inf)]
                
                temp<-rbinom(sum(neighbor.of.I$new.inf==1), size = 1, prob=prop.asymp)
                neighbor.of.I[new.inf==1,new.inf.Ia:=temp]
                neighbor.of.I$new.inf.Ia[is.na(neighbor.of.I$new.inf.Ia)]<-0
                neighbor.of.I[new.inf==1 & new.inf.Ia==1,stt:='Ia-L']
                neighbor.of.I[new.inf==1 & new.inf.Ia!=1,stt:='Ip-L']
                neighbor.of.I<-neighbor.of.I[!is.na(stt)]
                
                new.inf<-nrow(neighbor.of.I)
                if(new.inf>0){
                  wdata.new[neighbor.of.I$node1]<-neighbor.of.I$stt
                  wdata.inf.count[t]<-new.inf
                }
              }
            }
            
            Ip_pre.list<-which(wdata.old=='Ip-L')
            if(length(Ip_pre.list)>=1){
              temp<-f.Ip_pre.to.Ip(Ip_pre.list,t-1)
              wdata.new[Ip_pre.list]<-temp
            }
            
            Ia_pre.list<-which(wdata.old=='Ia-L')
            if(length(Ia_pre.list)>=1){
              temp<-f.Ia_pre.to.Ia(Ia_pre.list,t-1)
              wdata.new[Ia_pre.list]<-temp
            }
            
            Ip.list<-which(wdata.old=='Ip')
            if(length(Ip.list)>=1){
              temp<-f.Ip.to.Is(Ip.list,t-1)
              wdata.new[Ip.list]<-temp
            }
            
            Is.list<-which(wdata.old=='Is')
            RT.new<-vector()
            if(length(Is.list)>=1){
              temp<-f.Is.to.RTs(Is.list,t-1)
              wdata.new[Is.list]<-temp
              RT.new<-Is.list[which(temp=='RT')]
            }
            
            RT.list<-which(wdata.old=='RT')
            if(length(RT.list)>=1){
              temp<-f.RTs.to.D(RT.list,t-1)
              wdata.new[RT.list]<-temp
            }
            
            Ia.list<-which(wdata.old=='Ia')
            if(length(Ia.list)>=1){
              temp<-f.Ia.to.R(Ia.list,t-1)
              wdata.new[Ia.list]<-temp
            }
            
            R.list<-which(wdata.old=='R')
            if(length(R.list)>=1){
              wdata.new[R.list]<-'R'
            }
            
            D.list<-which(wdata.old=='D')
            if(length(D.list)>=1){
              wdata.new[D.list]<-'D'
            }
            
            D.list<-which(wdata.new=='D')
            if(length(D.list)>0){
              wdata.qt.temp[D.list] <- -Inf
              wdata.qt.actual.temp[D.list] <- Inf
            }
            
            
            #### Tests & Quarantine ####
            Test.available<-vector()
            Test.available<-which(wdata.qt.temp!=-Inf & wdata.tp.temp==0) # people 1) who are alive & 2) who have not tested positive
            #Test.available<-Test.available[!Test.available %fin% non.cooperative.nodes]
            
            test.ct<-vector(); test.rt<-vector(); test<-vector();
            if((sum(wdata.inf.count[1:t])/N)>=(intervention.start.inf.rate/100.)){ 
              if(t>test.x.days.later){
                test.ct<-which(wdata.ct.new[[as.character(t-test.x.days.later)]] %fin% c('CT'))
                test.ct<-test.ct[test.ct %fin% Test.available]
              }
            }
            
            test.rt<-RT.new
            test.rt<-test.rt[test.rt %fin% Test.available]
            test.rt<-test.rt[!test.rt %fin% test.ct]
            test.rt.all<-unique(c(test.rt.all,test.rt))
            
            test<-unique(c(test.ct,test.rt))
            wdata.test.count[t]<-length(test)
            
            temp<-vector(); P.list<-vector(); N.list<-vector(); FP.list<-vector(); FN.list<-vector(); TP.list<-vector(); TN.list<-vector(); test.rt.tp<-vector(); test.ct.tp<-vector(); test.rt.tn<-vector(); test.ct.tn<-vector(); test.rt.tn.non.coop<-vector(); test.ct.tn.non.coop<-vector();
            temp<-which(wdata.new %fin% c('Ia-L','Ip-L','Ip','Is','Ia','RT'))
            P.list<-test[test %fin% temp]      
            P.list<-data.table(node=as.integer(P.list), status=wdata.new[P.list], rate=ifelse(wdata.new[P.list] %fin% c('Ia-L','Ip-L'),1,(FN.rate/100.)))
            P.list[,fn:=rbinom(length(P.list$rate), size = 1, prob=P.list$rate)]
            
            N.list<-test[!test %fin% temp] 
            N.list<-data.table(node=as.integer(N.list), fp=rbinom(length(N.list), size=1, prob=(FP.rate/100.)))
            FP.list<-N.list[fp==1][['node']]
            TN.list<-N.list[fp==0][['node']]
            
            FN.list<-P.list[fn==1][['node']]
            TP.list<-P.list[fn==0][['node']]
            
            test.rt.tp<-test.rt[test.rt %fin% TP.list]
            test.rt.tn<-test.rt[!test.rt %fin% TP.list]
            test.rt.tn.non.coop<-test.rt.tn[test.rt.tn %in% non.cooperative.nodes]

            test.ct.tp<-test.ct[test.ct %fin% TP.list]
            test.ct.tn<-test.ct[!test.ct %fin% TP.list]
            test.ct.tn.non.coop<-test.ct.tn[test.ct.tn %in% non.cooperative.nodes]
            
            if(length(test.ct)>0){
              if(length(test.ct.tp)>0){
                if(t<=(T.max-turnaround)){
                  set(wdata.tp.new, i=test.ct.tp, j=as.character(t+turnaround), value='TP')
                }
              }
              if(length(test.ct.tn.non.coop)>0){
                wdata.qt.actual.temp[test.ct.tn.non.coop] <- -turnaround
              }
            }
            
            if(length(test.rt)>0){
              wdata.qt.temp[test.rt] <- -14
              wdata.qt.actual.temp[test.rt] <- -14
              if(length(test.rt.tp)>0){
                if(t<=(T.max-turnaround)){
                  set(wdata.tp.new, i=test.rt.tp, j=as.character(t+turnaround), value='TP')
                }
              }
              if(length(test.rt.tn.non.coop)>0){
                wdata.qt.actual.temp[test.rt.tn.non.coop] <- -turnaround
              }
            }
            
            rm.col<-vector();
            rm.col<-colnames(wdata.tp.new)[colnames(wdata.tp.new)!='1'  & (as.integer(colnames(wdata.tp.new)) < (t))]
            if(length(rm.col)>0){wdata.tp.new[, (rm.col):=NULL]}
            
            qt.nodes<-vector(); qt.nodes.real<-vector();
            #qt.nodes<-which(wdata.qt.temp>=-14 & wdata.qt.temp<0)
            #qt.nodes.real<-which((wdata.qt.temp>=-14 & wdata.qt.temp<0) & (wdata.qt.actual.temp!=-Inf & wdata.qt.actual.temp<0))
            qt.nodes.real<-which(wdata.qt.actual.temp>=-14 & wdata.qt.actual.temp<0)
            if(length(qt.nodes.real)>0){
              wdata.qt.count[t]<-length(qt.nodes.real)
            }
            
            
            #### Contact Tracing ####
            ## Contact network on day t
            temp.rm<-vector(); temp.rm.real<-vector(); temp.edge<-vector(); 
            temp.rm<-which(wdata.qt.temp<0) # people who have died or are in quarantine from HA's perspective
            temp.rm.real<-which(wdata.qt.actual.temp<0) # people who have died or are actually in quarantine

            if(length(temp.rm)>0){
              remove.edges.real[[t]]<-edge.list.undir[(node1 %fin% temp.rm.real | node2 %fin% temp.rm.real)][['edge_num']]
              temp.edge<-edge.list.undir[!(edge_num %fin% remove.edges.real[[t]]),c(1,2)]
              remove.edges.HA[[t]]<-edge.list.undir[(node1 %fin% temp.rm | node2 %fin% temp.rm)][['edge_num']]
            }else{
              temp.edge<-edge.list.undir[,c(1,2)]
            }    
            edges.all.final.transmission[[t]]<-rbindlist(list(temp.edge,temp.edge[,c(2,1)]), use.names=FALSE)  
            
            ## tracable network on day t
            TP.list.for.CT<-vector();
            TP.list.for.CT<-which(wdata.tp.new[[as.character(t)]]=='TP')
            if(length(TP.list.for.CT)>0){
              wdata.tp.temp[TP.list.for.CT]<-1
              if((sum(wdata.inf.count[1:t])/N)>=(intervention.start.inf.rate/100.)){  
                wdata.ct.day[t]<-1
                
                if(hops>0){
                  tt<-vector(); temp.rm.all.HA<-vector(); temp.rm.all.real<-vector(); edges.ct.undir<-vector(); edges.ct.all<-vector();
                  tt<-ifelse(t-13>1,t-13,1):(t)            
                  
                  temp.rm.all.real<-unlist(remove.edges.real[tt])
                  temp.rm.all.real<-which(tabulate(temp.rm.all.real)==length(tt)) # edges which have not appeared for entire last 14 days. 
                  
                  temp.rm.all.HA<-unlist(remove.edges.HA[tt])
                  temp.rm.all.HA<-which(tabulate(temp.rm.all.HA)==length(tt)) # edges which have not been identified by health authority for entire last 14 days.
                  
                  if(length(temp.rm.all.HA)>0){
                    edges.ct.undir<-edge.list.undir[!(edge_num %fin% temp.rm.all.HA),c(1,2)] 
                  }else{
                    edges.ct.undir<-edge.list.undir[,c(1,2)]
                  }
                  edges.ct.all<-rbindlist(list(edges.ct.undir,edges.ct.undir[,c(2,1)]), use.names=FALSE)
                  
                  #edges.ct.all<-edges.ct.all[!(node1 %fin% non.cooperative.nodes[non.cooperative.nodes %fin% test.rt.all])]
                  
                  
                  multihop<-vector(); QT.available<-vector()
                  multihop<-contact.tracing(hops,TP.list.for.CT,edges.ct.all)
                  QT.available<-which(wdata.qt.temp>=0 & wdata.tp.temp==0) # people who can potentially be informed that they are in close contact with an infected person
                                                                            # people 1) who are alive & 2) who are not in quarantine from HA's perspective & 3) who have not tested positive
                  #QT.available<-QT.available[!QT.available %fin% non.cooperative.nodes]
                  multihop<-multihop[multihop %fin% QT.available] # people who will comply with the health authority's request among the traced contacts
                  if(length(multihop)>0){
                    set(wdata.ct.new, i=as.integer(multihop), j=as.character(t), value='CT')
                    
                    rm.col<-vector();
                    rm.col<-colnames(wdata.ct.new)[colnames(wdata.ct.new)!='1'  & (as.integer(colnames(wdata.ct.new)) < (t-test.x.days.later))]
                    if(length(rm.col)>0){wdata.ct.new[, (rm.col):=NULL]}
                    wdata.qt.temp[multihop] <- -(14+1)
                    wdata.qt.actual.temp[multihop] <- -(14+1)
                  }
                }
              }
            }
            wdata.qt.temp<-wdata.qt.temp+1
            wdata.qt.actual.temp<-wdata.qt.actual.temp+1
            
            if(hops==0){
              wdata.ct.day.cum<-cumsum(wdata.ct.day[1:t]);
              if((t==90 & sum(wdata.inf.count[1:t])<40) | (t==365 & wdata.ct.day.cum[t]==0)){
                wdata.inf.count[1:T.max]<-NA
                wdata.test.count[1:T.max]<-NA
                wdata.qt.count[1:T.max]<-NA
                wdata.ct.day[1:T.max]<-NA
                iteration.num<-NA
                break
              }else{
                contact.tracing.day<-vector();
                contact.tracing.day<-min(which(wdata.ct.day.cum==1))
                if((contact.tracing.day+201)==t){
                  wdata.inf.count[t:T.max]<-NA
                  wdata.test.count[t:T.max]<-NA
                  wdata.qt.count[t:T.max]<-NA
                  wdata.ct.day[t:T.max]<-NA
                  iteration.num<-iteration
                  break
                }
              }
            }else{
              wdata.ct.day.cum<-cumsum(wdata.ct.day[1:t]);
              contact.tracing.day<-vector();
              contact.tracing.day<-min(which(wdata.ct.day.cum==1))
              if((contact.tracing.day+201)==t){
                wdata.inf.count[t:T.max]<-NA
                wdata.test.count[t:T.max]<-NA
                wdata.qt.count[t:T.max]<-NA
                wdata.ct.day[t:T.max]<-NA
                iteration.num<-iteration
                break
              }
            }
            
            
          }
          
          list(dailynewinf=wdata.inf.count, dailytests=wdata.test.count, dailyquarantines=wdata.qt.count, contacttracingday=wdata.ct.day, iterationnum=iteration.num)
        }
        
        num_iteration<-150 # the number of runs
        max.hops<-5
        strTmp.1<-as.character(1:num_iteration)
        for(k in 0:max.hops){ 
          daily.new.infections.hop<-setNames(data.table(matrix(,nrow = T.max, ncol = num_iteration)), strTmp.1)
          daily.tests.hop<-setNames(data.table(matrix(,nrow = T.max, ncol = num_iteration)), strTmp.1)
          daily.quarantines.hop<-setNames(data.table(matrix(,nrow = T.max, ncol = num_iteration)), strTmp.1)
          contact.tracing.day.hop<-setNames(data.table(matrix(,nrow = T.max, ncol = num_iteration)), strTmp.1)
          if(k==0){
            iteration.index<-rep(NA,num_iteration)
            pp<-1; qq<-1;
            while(pp<=num_iteration){
              if(qq==(num_iteration+1) & pp==1){break}
              set.seed(qq)
              temp<-epidemic(qq,as.integer(k))
              if(is.na(temp$iterationnum) == TRUE){
                print(paste('<k>=',aaaa,' : p=',bbbb,' : PoI=', cccc, ' hop ', k, ' iteration ', qq,  ' | ', temp$iterationnum, sep = ''))
                qq<-qq+1
              }else{
                daily.new.infections.hop[[pp]] <- (temp$dailynewinf)
                daily.tests.hop[[pp]] <- (temp$dailytests)
                daily.quarantines.hop[[pp]] <- (temp$dailyquarantines)
                contact.tracing.day.hop[[pp]]<-(temp$contacttracingday)
                iteration.index[pp]<-(temp$iterationnum)
                print(paste('<k>=',aaaa,' : p=',bbbb,' : PoI=', cccc, ' hop ', k, ' iteration ', qq,  ' | ', temp$iterationnum, sep = ''))
                qq<-qq+1
                pp<-pp+1
              }
              rm('temp','wdata.old','wdata.new','wdata.ct.new','wdata.tp.new','wdata.qt.new','wdata.qt.temp','wdata.tp.temp','wdata.inf.count','wdata.test.count','wdata.qt.count','wdata.ct.day')
              gc()
            }
          }
          else{
            if(sum(is.na(iteration.index))==num_iteration){break}
            for(pp in 1:num_iteration){
              set.seed(iteration.index[pp])
              temp<-epidemic(iteration.index[pp],as.integer(k))
              
              daily.new.infections.hop[[pp]] <- (temp$dailynewinf)
              daily.tests.hop[[pp]] <- (temp$dailytests)
              daily.quarantines.hop[[pp]] <- (temp$dailyquarantines)
              contact.tracing.day.hop[[pp]]<-(temp$contacttracingday)
              print(paste('<k>=',aaaa,' : p=',bbbb,' : PoI=', cccc, ' hop ', k, ' iteration ', iteration.index[pp],  ' | ', temp$iterationnum, sep = ''))
              
              rm('temp','wdata.old','wdata.new','wdata.ct.new','wdata.tp.new','wdata.qt.new','wdata.qt.temp','wdata.tp.temp','wdata.inf.count','wdata.test.count','wdata.qt.count','wdata.ct.day')
              gc()
            }
          }
          
          file.name.1<-paste("daily_new_infections_",k,"-hop.RData",sep="")
          file.name.2<-paste("daily_tests_",k,"-hop.RData",sep="")
          file.name.3<-paste("daily_quarantines_",k,"-hop.RData",sep="")
          file.name.4<-paste("contact_tracing_day_",k,"-hop.RData",sep="")
          
          save(daily.new.infections.hop, file=file.name.1, compress = T)
          save(daily.tests.hop, file=file.name.2, compress = T)
          save(daily.quarantines.hop, file=file.name.3, compress = T)
          save(contact.tracing.day.hop, file=file.name.4, compress = T)
          
          rm('daily.new.infections.hop','daily.tests.hop','daily.quarantines.hop','contact.tracing.day.hop')
          gc()
        }
        file.name.5<-"iteration_index.RData"
        save(iteration.index, file=file.name.5, compress = T)
        
        print(paste('<k>=',aaaa,' : p=',bbbb,' : PoI=',cccc, sep = ''))
        rm(list=ls()[! ls() %fin% c("aaaa","bbbb","cccc","mainDir.0","mainDir.6","intervention.start.inf.rate","FN.rate","FP.rate","turnaround",'cooperativity','infness.of.asymp',
                                    'latency.period','f.Ip_pre.to.Ip','f.Ia_pre.to.Ia','incubation.period',
                                    'post.latency.period','f.Ip.to.Is','wait.see.period','f.Is.to.RTs',
                                    'rec.death.period','death.rate','f.RTs.to.D','infectious.Ia.period','f.Ia.to.R','contact.tracing')])
        gc() #free up memrory and report the memory usage.
      }
    }
  }
}))





#### Generation of Networks ####
#mainDir.1<-"~/Dropbox/WS"
mainDir.1<-"/Users/Shared/Dropbox/1. PENN - Doctorate/Research/COVID-19/Covid19-2 Non Shared/Impact of contact tracing/Updates/04-2021/WS"
topology<-vector();
N<-100000
T.max<-600
for(k in c(4,8)){
  for(j in c(0.01,0.1,1)){
    setwd(mainDir.1)
    topology<-paste0(mainDir.1,"/ws_network_p_",j,"_k",k,".RData")
    
    set.seed(2020)
    g <- watts.strogatz.game(dim=1,size=N,nei=(k/2), p=j)
    
    V(g)$name<-seq(1:N)
    edges<-get.edgelist(g, names=TRUE)
    edges.all.temp.undir<-data.frame('node1'=edges[,1],'node2'=edges[,2], stringsAsFactors = FALSE)
    edges.all.temp.undir <- data.table(edges.all.temp.undir)
    
    
    edges.all.undir.final<-vector(mode='list',length=3)
    edges.all.undir.final[[1]]<-edges.all.temp.undir
    edges.all.undir.final[[2]]<-T.max
    edges.all.undir.final[[3]]<-N
    
    save(edges.all.undir.final, file=topology)
  }
}







#### ETC ####