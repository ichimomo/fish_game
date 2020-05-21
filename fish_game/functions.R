# define functions ----
get.SPR0 <- function(M,maa,waa,output="simple"){
  nage <- length(M)
  S <- exp(-M)
  N <- numeric()
  N[1] <- 1
  for(i in 2:(nage-1)) N[i] <- N[i-1]*S[i-1]
  N[nage] <- N[nage-1] * S[nage]/(1-S[nage])
  SPR0 <- sum(N * maa * waa)
  if(output=="simple") return(SPR0) else return(listN2(N,SPR0))
}

get.ab.bh <- function(h,R0,biopars){
  SPR0 <- get.SPR0(biopars$M,biopars$maa,biopars$waa)
  S0 <- R0*SPR0    
  beta <- (5*h-1)/(4*h*R0)
  alpha <- SPR0*(1-h)/(4*h)
  a <- 1/alpha
  b <- beta/alpha
  return(tibble::lst(SPR0,R0,h,S0,a,b))
}
SRF_BH <- function (x, a, b)  a * x/(1 + b * x)

caa.est.mat <- function(naa,saa,waa,M,catch.obs,Pope,set_max1=TRUE){
  if(set_max1==TRUE) saa <- saa/max(saa)
  tmpfunc <- function(logx,catch.obs=catch.obs,naa=naa,saa=saa,waa=waa,M=M,out=FALSE,Pope=Pope){
    x <- exp(logx)
    if(isTRUE(Pope)){
      caa <- naa*(1-exp(-saa*x))*exp(-M/2)
    }
    else{
      caa <- naa*(1-exp(-saa*x-M))*saa*x/(saa*x+M)
    }
    wcaa <- caa*waa
    if(out==FALSE){
      return(log((sum(wcaa,na.rm=T)-catch.obs)^2))
    }
    else{
      return(caa)
    }
  }
  tmp <- optimize(tmpfunc,log(c(0.000001,10)),catch.obs=catch.obs,naa=naa,saa=saa,waa=waa,M=M,Pope=Pope,out=FALSE)#,tol=.Machine$double.eps)
  tmp2 <- tmpfunc(logx=tmp$minimum,catch.obs=catch.obs,naa=naa,saa=saa,waa=waa,M=M,Pope=Pope,out=TRUE)
  return(list(x=exp(tmp$minimum),caa=tmp2))
}

# popluation dynamics ----
# determine dimension of population
pop_dyn <- function(nyear=100, nage=5, sd=0, seed=1, M=0.1, Fc=rep(0.3,nyear), h=0.7, R0=100,
                    init_number=10,
                    interactive=FALSE, # when true, Fc is used for population dynamics                    
                    xngrid=12, yngrid=10,
                    debug_mode=FALSE,
                    fish_image="fish.png",
                    catch_number_given=rep(10,nyear),
                    shiny_mode=FALSE, 
                    label_all=FALSE,
                    background_color=rgb(0.1,0.3,0.7,0.4)){
  nfish <- array(0,dim=c(nage,nyear),
                 dimnames=list(age=1:nage,year=1:nyear))
  rec_error <- array(0,dim=c(nyear),
                     dimnames=list(year=1:nyear))
  catch_number <- rep(0,nyear)
  
  # fish position
  fish_position <- expand.grid(x=(1:xngrid)*10, y=(1:yngrid)*10,
                               year=str_c("Y",1:nyear)) %>% 
    as_tibble() %>%
    mutate(value=0) %>% 
    pivot_wider(names_from=year, values_from=value) %>%
    mutate(cell_number=1:(xngrid*yngrid),image=fish_image) 
  
  sd_vector <- sd
  set.seed(sd)
  rec_error[] <- rnorm(length(rec_error)/length(sd_vector), -0.5*sd_vector^2, sd_vector)
  
  # definition of SR parameter 
  SR_pars_true <- list(BH=get.ab.bh(h=h,R0=R0,tibble(M=M,maa=c(0,0,1,1,1),waa=rep(1,5)))[c("a","b")],
                       HS=list(a=2, b=100))
  
  # population dynamics
  nfish[,1] <- init_number
  for(i in 1:nyear){
    if(i>1){
      if(isTRUE(interactive)){
        if(!isTRUE(debug_mode)&&!isTRUE(shiny_mode)){
          tmp <- NA
          while(is.na(tmp)){
            tmp <- readline(str_c(i-1,"Year: how many fish?")) %>%
              as.numeric()
            if(is.na(tmp)) print("Enter number again\n")
          }
          catch_number[i-1] <- tmp
        }
        else{ catch_number[i-1] <- catch_number_given[i-1] }
        if(is.na(catch_number[i-1])) break()
        if(catch_number[i-1]>sum(nfish[,i-1])){
          g1 <- g2 <- ggplot() + 
            theme_solid()+scale_colour_pander() + 
            theme(plot.background = element_rect(fill = "gray")) 
          g1 <- g1 + geom_image(data=tibble(x=1,y=1,img="extinction.png"),aes(x=x, y=y, image=img),size=1,asp=3)
          break         
         }
        Fc[i-1] <- caa.est.mat(nfish[,i-1],rep(1,nage),rep(1,nage),rep(M,nage),
                               catch_number[i-1],Pope=FALSE,set_max1=TRUE)$x
      }
      if(!isTRUE(interactive)){
        catch_number[i-1] <- sum(Fc[i-1]/(Fc[i-1]+M) * nfish[,i-1] * (1-exp(-Fc[i-1]-M)))
      }
      # forward cal
      for(a in 2:nage){
        nfish[a,i] <- nfish[a-1,i-1] * exp(-Fc[i-1]-M) 
      }
      nfish[nage,i] <- nfish[nage,i] + nfish[a,i-1] * exp(-Fc[i-1]-M)
      SRF <- SRF_BH
      nfish[1,i] <- SRF(x=sum(nfish[-1:-2,i]),
                        a=SR_pars_true$BH$a, 
                        b=SR_pars_true$BH$b) * exp(rec_error[i])
    }
    
    # plot
    if(isTRUE(interactive)){
      # fish distribution
      dat <- tibble(age=1:nage,fish=nfish[,i])
      fish_dist <- purrr::map(1:nrow(dat),
                              function(x) rep(dat$age[x], round(dat$fish[x]))) %>% unlist()
      
      fish_position[,i+2] <- 0
      fish_position[sample(fish_position$cell_number, length(fish_dist)),i+2] <- fish_dist
      
      g1 <- ggplot() +
        geom_point(data=fish_position[fish_position[,i+2]==1,],
                   aes(x=x, y=y),color="white", size=8) +
        geom_image(data=fish_position[fish_position[,i+2]==1,],
                   aes(x=x, y=y, image=image),asp=3, size=0.04) +
        geom_image(data=fish_position[fish_position[,i+2]>1 & fish_position[,i+2]<4,],
                   aes(x=x, y=y, image=image),asp=3, size=0.08) +
        geom_image(data=fish_position[fish_position[,i+2]>3,],
                   aes(x=x, y=y, image=image),asp=3, size=0.13) +
        coord_cartesian(xlim=c(0,xngrid+1)*10, ylim=c(0,yngrid+1)*10) +
        theme_solid()+scale_colour_pander() +              
        theme(plot.background = element_rect(fill = background_color)) 
      
      # dotplot of catch
      ymax <- 40
      dat <- tibble(year=1:nyear,catch=catch_number)
      dat2 <- map_dfr(1:nrow(dat), function(x) tibble(year=dat$year[x],catch=0:dat$catch[x])) %>%
        mutate(fyear=as.factor(year), img="fish_isaki2.png") %>%
        dplyr::filter(catch>0) %>%
        mutate(year=floor(catch/ymax)*0.1+year) %>%
        mutate(catch = catch%%ymax)
      
      g2 <- dat2 %>% ggplot() +
        geom_point(aes(x=year,y=catch-0.2), size=3, col="white", shape =17) +
        geom_point(aes(x=year,y=catch), size=6, col="white", shape =19) +
        geom_point(aes(x=year,y=catch+0.1), size=1, col="black", shape =19) +
        coord_cartesian(ylim=c(-0.5,ymax+0.5),xlim=c(1.5,nyear+0.5))+
        theme_solid()+scale_colour_pander() +              
        theme(plot.background = element_rect(fill = background_color),
              legend.position="none")
      
      #          g1 <- g1 + geom_label(data=tibble(x=xngrid*10/2, y=yngrid*10*0.9, label=str_c(i,"年め: 魚の数 ", round(sum(nfish[,i])), "匹")),
      g1 <- g1 + geom_label(data=tibble(x=xngrid*10/2, y=yngrid*10*0.9, label=str_c(i," year: ", round(sum(nfish[,i])), "fish")),          
                            mapping=aes(x=x,y=y,label=label), size=7,
                            fill=rgb(0.1,0.3,0.7,0.4), col="white")
      
      g2 <- g2+ geom_label(data=tibble(x=nyear/2+1, y=ymax, label=str_c("Total catch: ", sum(catch_number))),
                           mapping=aes(x=x,y=y,label=label), size=7,
                           fill=rgb(0.1,0.3,0.7,0.4), col="white")
      
      if(isTRUE(label_all)){
        g2 <- g2 +  geom_label(data=tibble(1:nyear, y=catch_number),
                               mapping=aes(x=x,y=y,label=y), size=7,
                               fill=rgb(0.1,0.3,0.7,0.4), col="white")
      }
      
      if(!isTRUE(shiny_mode)) gridExtra::grid.arrange(g1,g2)
      
    }
  }      
  # if(isTRUE(shiny_mode)) gridExtra::grid.arrange(g1,g2)
  if(!isTRUE(interactive)){
    g1 <- ggplot(data=tibble(x=1:nyear,y=apply(nfish,2,sum)))+
      geom_point(aes(x=x,y=y))
  }
  ssb <- colSums(nfish[-1:-2,])
  biom <- colSums(nfish)
  res <- tibble::lst(nfish,ssb,g1,catch_number,Fc,biom)
  if(isTRUE(shiny_mode)) res$g2 <- g2
  return(res)
  
}

nage <- 5
nyear <- 15
M <- 0.3
K <- 120
pop0 <- pop_dyn(Fc=rep(0,100),R0=10,M=M,nage=nage,nyear=100,init_number=30)

R0 <- 10*120/colSums(pop0$nfish)[100] # determin R0
pop0 <- pop_dyn(Fc=rep(0,100),R0=R0,M=M,nage=nage,nyear=100,init_number=30)

# calc Fmsy
x <- optimize(function(x)  pop_dyn(Fc=rep(x,100),R0=R0,M=M,nage=nage,nyear=100,init_number=30)$catch_number[99],
         c(0,1), maximum=TRUE)
pop_msy <- pop_dyn(Fc=rep(x$maximum,100),R0=R0,M=M,nage=nage,nyear=100,init_number=30)
MSY_15year <- sum(pop_msy$catch_number[1:14])
Fmsy <- x$maximum
Bmsy <- sum(pop_msy$nfish[,100])

init_number <- pop0$nfish[,100]

#res <- pop_dyn(Fc=rep(0,100),R0=R0,M=M,nage=nage,nyear=nyear,init_number=init_number,
               #interactive=TRUE,debug_mode=TRUE,shiny_mode=TRUE,
#               fish_image="fish_isaki.png",xngrid=15,yngrid=8,
               #catch_number_given=rep(10,15))

             