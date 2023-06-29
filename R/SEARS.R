SEARS <- function(p.p, p.d, p.tox, k1, k2, pi_t, pi_e, pT, eff_a = 0.5, eff_b = 0.5, plac_a = 0.5,
                  plac_b = 0.5, tox_a = 1, tox_b = 1, csize, csize2, p.star = 0.2, q.star = 0.6,
                  f.star = 0.06, p.star2 = 0.2, q.star2 = 0.98, d.cs, p.cs, phase1_size,
                  n_earlystop, extrasafe_BOIN = FALSE, offset_BOIN = 0.05, Nsim, n_catchup,
                  control_arm = "fixed", power_c = .5, lower_bound = .05, weight1, weight2, seed = 100){
  set.seed(seed)
  d = length(p.d) # dose levels

  ### BOIN design settings
  p_saf_BOIN <- 0.6 * pT
  p_tox_BOIN <- 1.4 * pT
  if(csize > 1){
    temp <- get.boundary(target = pT, ncohort = 100, cohortsize = csize,
                         n.earlystop = 1000, p.saf = p_saf_BOIN,
                         p.tox = p_tox_BOIN, cutoff.eli = k1,
                         extrasafe = extrasafe_BOIN, offset = offset_BOIN)$full_boundary_tab
  }
  else{
    temp <- get.boundary(target = pT, ncohort = 100, cohortsize = csize,
                         n.earlystop = 1000, p.saf = p_saf_BOIN,
                         p.tox = p_tox_BOIN, cutoff.eli = k1,
                         extrasafe = extrasafe_BOIN, offset = offset_BOIN)$boundary_tab
  }
  b.e = temp[2, ]
  b.d = temp[3, ]
  b.elim = temp[4, ]

  final.tol.sampsize=rep(0,Nsim)
  final.dose.sampsize=matrix(0,Nsim,d)
  final.plac.sampsize=rep(0,Nsim)
  final.count = matrix(0,Nsim,d) ## number of counts of the each doses being selected
  dose.select = matrix(0,Nsim,d)
  final.d.cum.esca = matrix(0,Nsim,d) ## accumulating sample size of dose arms at dose escalation stage
  count.typeI = rep(0,Nsim) ## for computing overall Type I error.
  tox.count = matrix(0,Nsim,d) ## matrix for containing toxic responses every simulation 2-2-2013
  for(j in 1:Nsim) {
    d.a = rep(eff_a,d); d.b = rep(eff_b,d)
    p.a = plac_a; p.b = plac_b
    d.resp = rep(0,d); p.resp = 0 # for storing responses of doses and placebo
    q = c() # storing the posterior unit probabilities
    ## sample size constraints for dose, placebo and total
    #     d.cs = 36; p.cs = 36 ;  t.cs = d*36
    t.cs = d*d.cs

    t.a = rep(tox_a,d); t.b = rep(tox_b,d) # # priors for the toxicity parameters
    tox.resp = rep(0,d) # storing the toxic responses


    select = rep(0,d) ## accumulate sample size for dose.select matrix.
    p.cum.size = 0; d.cum.size = rep(0,d); t.cum.size = 0 ## t.cum.size for accumulating the total sampsize.
    d.cum.esca = rep(0,d) ## accumulating sample size of dose arms at dose escalation stage
    phase1.cum.size = 0 ## for computing total sample size in phase I
    start.dose=1
    D=start.dose

    set1 = c(1:d); set2 = c(); ## set1 is the dose set for toxicity trial; set2 is the dose set for the part2 adaptive randomization efficacy trial.

    count = rep(0,d) ## count the selection counts of dose arms


    stop.now = 0 # stop.now =1 indicates we meet some stopping rules

    #     over_toxic = 0 ## count the over toxic number of patients.date: 12-24-2012
    #     over_toxic_vector = rep(0,Nsim)
    #     over_toxic_rate = rep(0,Nsim)



    while( stop.now==0 ){


      if( length(set1)!=0 ){

        ###### counting toxic numbers ######
        #         if( (p.tox[D]==0.3)||(p.tox[D]==0.5) ) over_toxic = over_toxic + 1  ## date:12-24-2012
        ####################################

        tox.resp[ D ] = rbinom(1,csize,p.tox[ D ])
        ## update the posterior beta parameters for toxicity study
        t.a[ D ] = tox.resp[ D ] + t.a[ D ]; t.b[ D ] = csize - tox.resp[ D ] + t.b[ D ]
        tox.resp <- t.a - tox_a
        # accumulate each dose's sampsize
        d.cum.size[ D ] = d.cum.size[ D ] + csize

        if(1 %in% set1){
          if (extrasafe_BOIN == FALSE){
            if((!is.na(b.elim[d.cum.size[ 1 ]])) && (tox.resp[ 1 ] >= b.elim[d.cum.size[ 1 ]])){
              stop.now = 1
            }
          }
          else{
            if( (d.cum.size[1] >= 3) && (1-pbeta(pT, t.a[ 1 ], t.b[ 1 ])>k1 - offset_BOIN)){
              stop.now=1}
          }
        }

        ##  Meanwhile we have the efficacy responses information and update the posterior parameters of the current dose for the efficacy
        d.resp[ D ] = rbinom(1,csize,p.d[ D ])
        ## update the posterior parameters
        d.a[ D ] = d.resp[ D ] + d.a[ D ]; d.b[ D ] = csize- d.resp[ D ] +d.b[ D ]
        d.resp <- d.a - eff_a


        ### the graduation rule
        n = length(set1)
        DD = D; ii = which(set1==D)
        aa = ((1-pbeta(pi_e,d.a[ D ],d.b[ D ])) >q.star)&&(pbeta(pi_t,t.a[D],t.b[D])>p.star) ## to see if the current dose can be graduated
        bb <- d.cum.size[D] >= n_earlystop
        ##cat("aa = ", aa, "\n")
        ##cat(t.a, "\n", t.b, "\n")
        ##readline()

        ## count the sample accumulation at dose escalation stage
        d.cum.esca[ D ] = d.cum.size[D]
        final.d.cum.esca[j,D]= d.cum.esca[D]
        elim_criterion <- b.elim[d.cum.size[D]]
        escalate_criterion <- b.e[d.cum.size[D]]
        deescalate_criterion <- b.d[d.cum.size[D]]

        ##     if( length(set1)!=0 )
        ##       {

        if( D==min(set1) ) {

          ## if the first dose in set1 is too toxic, the trial will be terminated
          if( (!is.na(elim_criterion)) && (tox.resp[ D ] >= elim_criterion) ){
            if( n>1 ) {
              set1 <- NULL
            } else { set1 = NULL }

            ##stop.now=1
          } else {

            if( (escalate_criterion < tox.resp[ D ]) && (tox.resp[ D ] < deescalate_criterion)){
              ## this is the condition of 'stay'
              if( (aa==TRUE) || (bb == TRUE) ){

                if( n>1 ){
                  D = set1[ii+1]; set1 = set1[-ii]; set2 = unique(c(set2,DD));
                } else { set2 = unique(c(set2,DD)); set1 = NULL;}
              } else { D = D }

            } else {# Ends the if(stay).

            if( tox.resp[ D ] <= escalate_criterion){
              ## this is the condition of 'escalcate'

              if( (aa==TRUE) || (bb == TRUE) ){
                if( n>1 ){
                  D = set1[ii+1]; set1 = set1[-ii]; set2 = unique(c(set2,DD));
                } else { set2 = unique(c(set2,DD)); set1 = NULL;  }
              } else {
                if( n>1 ) { D = set1[ii+1] } else { D = D }
              }

            } # Ends the if(escalate).

            }
          } # Ends the if(tox)... else.

        } else {

          if ( D==max(set1) ){

            if( (!is.na(elim_criterion)) && (tox.resp[ D ] >= elim_criterion) ){

              ## too toxic, we should remove the dose out of the trial now.And if possible, deescalate to the adjacent lower dose
              if( n>1 ) {
                D = set1[ii-1]; set1 = set1[-ii]
              } else { set1 = NULL }

            } else {

              if(tox.resp[ D ] >= deescalate_criterion){

                ## this is the condition of 'descalation'
                if( (aa ==TRUE) || (bb == TRUE) ) {
                  D <- set1[ii-1]
                } else {
                  D = set1[ii-1]
                }


              } else {# Ends the if(descalate).

              if((escalate_criterion < tox.resp[ D ]) && (tox.resp[ D ] < deescalate_criterion)){

                ## this is the condition of 'stay'
                if( (aa==TRUE) || (bb == TRUE) ){
                  if( n>1 ){
                    D = set1[ii-1]; set1 = set1[-ii]; set2 = unique(c(set2,DD))
                  } else { set1 = NULL; set2 = unique(c(set2,DD)) }

                } else {
                  D = D
                }


              } # Ends the if(stay).
              }
            } # Ends the if(tox)...else.

          } else {

            if( (D>min(set1))&&(D<max(set1)) ){

              if( (!is.na(elim_criterion)) && (tox.resp[ D ] >= elim_criterion) ){

                ## here the current dose is too toxic, we should remove this dose and higher doses out of the trial
                set1 = set1[ -c(ii:n) ]; D = set1[ ii-1 ]

              } else {

                if( tox.resp[ D ] >= deescalate_criterion ){

                  ## this is the condition of 'descalate'
                  if( (aa==TRUE) || (bb == TRUE) ){
                    D <- set1[ii-1]
                  } else { D = set1[ii-1] }

                } else {# Ends the if(descalcate).

                if(  (escalate_criterion < tox.resp[ D ]) && (tox.resp[ D ] < deescalate_criterion )){

                  ## this the the condition of 'stay'
                  if( (aa==TRUE) || (bb == TRUE) ){
                    set1 = set1[-ii]; set2 = unique(c(set2,D)); D = set1[ii] ## because we should move to the next higer dose level,but the current dose at the position has been graduated,so the position of next higher dose in the updated set2 is still ii.
                  } else {
                    D = D
                  }

                } else {# Ends the if(stay).

                if( tox.resp[ D ] <= escalate_criterion ){

                  # this is the condition of 'escalate'
                  if( (aa==TRUE) || (bb == TRUE) ){
                    set1 = set1[-ii]; set2 = unique(c(set2,D)); D = set1[ii] ## because we should move to the next higer dose level,but the current dose at the position has been graduated,so the position of next higher dose in the updated set2 is still ii.
                  } else {
                    D = set1[ii+1]
                  }

                } # Ends the if(escalate).
                  }
                }
              } # Ends the if(tox)...else.


            }
          }
        }

        ##if(DD %in% set2) d.cum.size[ DD ] = d.cum.size[ DD ] - 3 ## we should avoid adding sample size twice to this dose arm after they are graduated to phase II at the first time.

        ##   }  # Ends the if( length(set1)!=0 ).

        #   ######  Exclude very unefficious dose arm from phase I trial though it maybe safe!
        ##          if(length(set1)!=0)
        ##            {
        cc = which((1-pbeta(p.p, d.a[set1], d.b[set1])) < f.star)
        if(length(cc) > 0)
        {
          ##                {
          ##                  set1 = set1[-cc]
          ##                }
          ##            }
          ##cat(cc, "\n")


          del.D <- (D %in% set1[cc])

          ##        cc = ( (1-pbeta(p.p,d.a[DD],d.b[DD])) < f.star )
          if(!del.D){
            set1 = set1[-cc]
          }
          else{
            ii = which(set1==D)

            if(length(set1==1)){
              set1=NULL
            } else{
              if(D==min(set1)){
                D = set1[ii+1]
                set1 = set1[-ii]
              }
              if(D==max(set1)){
                D = set1[ii-1]
                set1 = set1[-ii]
              }
              if((D>min(set1))&&(D<max(set1))){
                D= set1[ii-1]
                set1 = set1[-ii]
              }
            }
          }
        }

        ## } ## Ends if( length(set1)!=0 ).

        phase1.cum.size = phase1.cum.size + csize  ## current accumulated total sample size in phase I.
        if( phase1.cum.size>=phase1_size ) { set2=unique(c(set2,set1)); set1=NULL; }  ## stop phase I when the max number of patients in phase I reaches 30, and all doses in the
        ##cat(set1,"\n", set2)
        ##readline()

      } # # Ends if(length(set1!=0)).-- phase I loop.


      #    if( (length(set1)==0)&&(length(set2)==0) ) { stop.now==1; }

      ##################################### Part 2: Adaptive Randomization #############################################################
      d.ind = c()
      p.ind = 0

      ## Exclude toxic arms in set2--over-toxic rule in phase II.
      exclude = NULL
      if( length(set2)!=0 ){

        ########^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        ########^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        max.s2 <- which.max(set2)

        for(kk in 1:length(set2)){
          if( ((1-pbeta(pT,t.a[set2[kk]],t.b[set2[kk]]))> k2) ){
            exclude = c(exclude,kk:max.s2) ## exclude set2[kk] and all the higher doses in set2
            ##cat("kk=", kk, "\n")
            ##cat(set2, "\n",t.a, "\n", t.b,"\n")
            ##readline()
          }

          ########^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          ########^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

        }
        if(length(exclude)!=0){
          set2 = set2[-exclude]
        }
      }
      ## over-toxic rule in phase II loop.


      if( length(set2)!=0 ) {
        nn = length(set2)
        rand=1000
        success_prob <- c(p.p, p.d[set2])
        a1 <- rbeta(rand,p.a,p.b)
        for (ii in 1:nn) {
          assign(paste("a", ii+1, sep = ""), rbeta(rand,
                                                   d.a[set2[ii]], d.b[set2[ii]]))
        }
        pmaxfunction <- paste0("pmax(", paste(paste0("a",
                                                     1:length(success_prob)), collapse = ","), ")")
        for (ii in 1:length(success_prob)) {
          assign(paste("c", ii, sep = ""), eval(parse(text = paste0("length(which(a",
                                                                    ii, "==", pmaxfunction, "))/rand"))))
        }
        if(p.cum.size < n_catchup){
          cc1 <- max(c1, 1/length(success_prob))
        }
        else{
          cc1 <- c1
        }
        for(ii in 1:nn){
          if(d.cum.size[set2[ii]] < n_catchup){
            max_val <- max(eval(parse(text = paste0("c", ii+1))), 1/length(success_prob))
            assign(paste0("cc", ii+1), max_val)
          }
          else{
            assign(paste0("cc", ii+1), eval(parse(text = paste0("c", ii+1))))
          }
        }
        for(ii in 1:length(success_prob)){
          assign(paste0("c", ii), eval(parse(text = paste0("cc", ii))))
        }
        sum_c <- eval(parse(text = paste0("sum(", paste("c", 1:length(success_prob), sep = "", collapse = ", "), ")")))
        for(ii in 1:length(success_prob)){
          assign(paste0("c", ii), eval(parse(text = paste0("c", ii, "/sum_c"))))
        }
        if (control_arm == "fixed") {
          previous_rem_sum <- 1 - c1
          c1 <- 1/length(success_prob)
          if(previous_rem_sum == 0){
            for (xx in 2:length(success_prob)) {
              assign(paste("c", xx, sep = ""),
                     (1-c1)/(length(success_prob)-1))
            }
          }
          else{
            after_rem_sum <- 1 - c1
            for (xx in 2:length(success_prob)) {
              assign(paste("c", xx, sep = ""),
                     eval(parse(text = paste0("c", xx, "/previous_rem_sum*after_rem_sum"))))
            }
          }

        }
        lambda <- power_c
        lambdafunction1 <- paste(paste0("c", 1:length(success_prob),
                                        "^lambda"), collapse = "+")
        lambdafunction2 <- paste(paste0("c", 2:length(success_prob),
                                        "^lambda"), collapse = "+")
        if (control_arm == "fixed") {
          rem_sum <- 1 - c1
          r1 <- c1
          for (z in 2:length(success_prob)) {
            assign(paste("r", z, sep = ""),
                   eval(parse(text = paste0("(c", z, "^lambda)/(",
                                            lambdafunction2, ")*rem_sum"))))
          }
        }
        else {
          for (z in 1:length(success_prob)) {
            assign(paste("r", z, sep = ""),
                   eval(parse(text = paste0("(c", z, "^lambda)/(",
                                            lambdafunction1, ")"))))
          }
        }
        assign("r", eval(parse(text = paste("c(",
                                            paste(paste0("r", 1:length(success_prob)),
                                                  collapse = ","), ")"))))

        all_i <- c(1:length(r))
        max_i <- which.max(r)
        min_i <- which.min(r)
        else_i <- all_i[!all_i %in% c(max_i, min_i)]
        min_threshold <- lower_bound
        max_threshold <- 1 - (length(r) - 1) * lower_bound
        if (max(r) > max_threshold & min(r) < min_threshold) {
          assign(paste0("r", max_i), eval(parse(text = paste0("min(",
                                                              max_threshold, ",r", max_i, ")"))))
          assign(paste0("r", min_i), eval(parse(text = paste0("max(",
                                                              min_threshold, ",r", min_i, ")"))))
          for (i in else_i) {
            assign(paste0("r", i), eval(parse(text = paste0("max(",
                                                            min_threshold, ",r", i, ")"))))
            if (eval(paste0("r", i)) > min_threshold) {
              assign(paste0("r", i), min_threshold)
            }
          }
        }
        if (max(r) <= max_threshold & min(r) < min_threshold) {
          assign(paste0("r", min_i), eval(parse(text = paste0("max(",
                                                              min_threshold, ",r", min_i, ")"))))
          for (i in else_i) {
            assign(paste0("r", i), eval(parse(text = paste0("max(",
                                                            min_threshold, ",r", i, ")"))))
          }
          assign(paste0("r", max_i), eval(parse(text = paste0("1-r",
                                                              min_i, "-", paste(paste0("r", else_i,
                                                                                       collapse = "-"))))))
        }
        assign("tempprob", eval(parse(text = paste("c(",
                                                   paste(paste0("r", 1:length(success_prob)),
                                                         collapse = ","), ")"))))
        Rho = tempprob
        #NNN <- c(p.cum.size, d.cum.size[set2])
        #if (control_arm == "fixed") {
        #  phi.jk <- Rho * ((Rho/(NNN/(sum(NNN))))^2)
        #  phi <- phi.jk
        #  ratio <- phi.jk/sum(phi.jk)
        #  previous_rem_sum <- 1 - ratio[1]
        #  phi[1] <- 1/length(success_prob)
        #  after_rem_sum <- 1 - phi[1]
        #  for (zzz in 2:length(success_prob)) {
        #    phi[zzz] <- after_rem_sum * ratio[zzz]/previous_rem_sum
        #  }
        #}
        #else {
        #  phi.jk <- Rho * ((Rho/(NNN/(sum(NNN))))^2)
        #  phi <- phi.jk/sum(phi.jk)
        #}
        phi <- Rho
        all_i <- c(1:length(phi))
        max_i <- which.max(phi)
        min_i <- which.min(phi)
        else_i <- all_i[!all_i %in% c(max_i, min_i)]
        min_threshold <- lower_bound
        max_threshold <- 1 - (length(phi) - 1) * lower_bound
        for (i in all_i) {
          assign(paste0("phi", i), eval(parse(text = paste0("phi[",
                                                            i, "]"))))
        }
        if (max(phi) > max_threshold & min(phi) < min_threshold) {
          assign(paste0("phi", max_i), eval(parse(text = paste0("min(",
                                                                max_threshold, ",phi", max_i, ")"))))
          assign(paste0("phi", min_i), eval(parse(text = paste0("max(",
                                                                min_threshold, ",phi", min_i, ")"))))
          for (i in else_i) {
            assign(paste0("phi", i), eval(parse(text = paste0("max(",
                                                              min_threshold, ",phi", i, ")"))))
            if (eval(paste0("phi", i)) > min_threshold) {
              assign(paste0("phi", i), min_threshold)
            }
          }
        }
        if (max(phi) <= max_threshold & min(phi) < min_threshold) {
          assign(paste0("phi", min_i), eval(parse(text = paste0("max(",
                                                                min_threshold, ",phi", min_i, ")"))))
          for (i in else_i) {
            assign(paste0("phi", i), eval(parse(text = paste0("max(",
                                                              min_threshold, ",phi", i, ")"))))
          }
          assign(paste0("phi", max_i), eval(parse(text = paste0("1-phi",
                                                                min_i, "-", paste(paste0("phi",
                                                                                         else_i, collapse = "-"))))))
        }
        assign("phi", eval(parse(text = paste("c(",
                                              paste(paste0("phi", 1:length(phi)), collapse = ","),
                                              ")"))))





        #aaaa=matrix(0,nrow=rand,ncol=nn)
        #pppp=c()


        #for(jj in 1:nn){
        #  aaaa[,jj] = rbeta(rand,d.a[set2[jj]],d.b[set2[jj]])
        #}
        #pppp = rbeta(rand,p.a,p.b)

        #if(nn>2){
        #  for(jj in 1:nn){
        #    d.ind[jj] = mean( aaaa[,jj]>apply(aaaa[,-jj],1,max) )
        #  }
        #  p.ind = mean( pppp>apply(aaaa,1,max) )
        #}
        #if(nn==2){
        #  for(jj in 1:nn){
        #    d.ind[jj] = mean( aaaa[,jj]>aaaa[,-jj] )
        #  }
        #  p.ind = mean( pppp>apply(aaaa,1,max)  )
        #}
        #if(nn==1){
        #  p.ind = mean( pppp>apply(aaaa,1,max) )
        #  d.ind = 1-p.ind
        #}




        p.ind <- phi[1]
        d.ind <- phi[-1]

        ## begin adaptive randomization process and etc......
        ddd = sample( c(p.ind, d.ind),1,prob=phi,replace=T )
        if( any(d.ind==ddd) ) {
          # find out which dose is selected out now.
          i = which(d.ind==ddd )
          i = set2[i]
          #$cat("i is ",i,"\n")
          # assign 3 patients to this i-th dose
          d.cum.size[i] = d.cum.size[i] + csize2
          # observe the response of this dose
          d.resp = rbinom(1,csize2,p.d[i])
          # update posterior beta paremeters of this dose
          d.a[i] = d.a[i] + d.resp; d.b[i]= csize2 - d.resp + d.b[i]
          d.resp <- d.a - eff_a
          ## after adaptive randomization, still conduct safety testing
          tox.resp[ i ] = rbinom(1,csize2,p.tox[ i ])
          ## update the posterior beta parameters for toxicity study
          t.a[ i ] = tox.resp[ i ] + t.a[ i ]; t.b[ i ] = csize2 - tox.resp[ i ] + t.b[ i ]
          tox.resp <- t.a - tox_a
          ## futility rule in phase II.
          exclude = NULL
          for(jj in 1:length(set2)){
            if( ((1-pbeta(p.p,d.a[set2[jj]],d.b[set2[jj]]))<f.star)  )
              exclude = c(exclude,jj)
          }
          if(length(exclude)!=0){
            ##cat(exclude, "\n")
            ##readline()
            set2 = set2[-exclude]
          }
          ##cat("set2 ",set2,'\n')


        } else {
          # assign 1 patients to the placebo
          p.cum.size = p.cum.size + csize2
          # observe the response of the placebo
          p.resp = rbinom(1,csize2,p.p)
          # updae the posterior beta paremeters of placebo
          p.a = p.a + p.resp; p.b = csize2 - p.resp + p.b }
          p.resp <- p.a - plac_a

      } ## Ends the if(length(set2)==0) loop.

      ## evaluate if we meet the stop criteria
      #       t.cum.size = t.cum.size + csize + ifelse(length(set2)!=0,csize2,0)
      t.cum.size = sum(d.cum.size) + p.cum.size

      #       cat("t.cum.size is ",t.cum.size,"\n")

      if( (p.cum.size >= p.cs)||(any(d.cum.size >= d.cs))||(t.cum.size >= t.cs) ) {stop.now=1;}
      ##         if( t.cum.size > 156-3 ) {    stop.now=1;}##cat("bug2", t.cum.size)




      #         cat(set1, "\n")
      #         readline()
      if( length(c(set1,set2))==0)  { stop.now=1 }

    } ## Ends while( stop.now==0 ) loop.

    #     ## Exclude uneffective arms in set2.
    #     exclude = NULL
    #     for(i in set2){
    #       if( !((pbeta(pi_t,t.a[i],t.b[i]))>p.star2)  )
    #         exclude = c(exclude,i)
    #     }
    #     if(length(exclude)!=0){
    #       set2 = set2[-exclude]
    #     }
    #
    ## more strict safety and efficacy here for recommending dose to the phase III
    bbb <- FALSE
    if(length(set2)!=0){
      for(i in 1:length(set2)){
        bbb = ((1-pbeta(pi_e,d.a[ set2[i] ],d.b[ set2[i] ])) >q.star2)&&(pbeta(pi_t,t.a[set2[i]],t.b[set2[i]])>p.star2)
        ##bbb = ( (1-pbeta(pi_e,d.a[ set2[i] ],d.b[ set2[i] ])) >q.star2 )
        if(bbb) {select[set2[i]]=1;
        ###### counting toxic numbers ######
        #                 if( (p.tox[i]==0.3)||(p.tox[i]==0.5) ) over_toxic = over_toxic + 1  ## date:12-24-2012
        ####################################
        }
      }
    }

    ##################################################################################################################
    #if( any(select!= 0) ) {
    ############### various sceniros for power computations and type I error computations ##################

    ## Sc 1-12.
    ## if sc=1 or sc=7, any dose selection is a type I error.
    #if( (sc==1)||(sc==7) ){
    #count.typeI[j] = 1
    #}
    ## if sc=6, 12, doses 1,2,3 are truely desirable doses
    #if((sc==3)||(sc==9)||(sc==6)||(sc==12)){
    #count.typeI[j] = ifelse(sum(select[c(1,2, 3)])!=0, 1, 0)
    #}
    ## if sc=2 or sc=8, sc=4 or sc=10, just dose 2 and dose 3 are 'good' doses for computing power.
    #if( (sc==2)||(sc==8)||(sc==4)||(sc==10) ){
    #count.typeI[j] = ifelse( sum(select[c(2,3)])!=0, 1, 0  )
    #}
    ## if sc=3 or sc=9, sc=5 or sc=11, just dose 1 and dose 2 are 'good' doses for computing power.
    #if( (sc==5)||(sc==11) ){
    #count.typeI[j] = ifelse( sum(select[c(1,2)])!=0, 1, 0  )
    #}

    ## Sc 13-36.
    ## for type I error of the following 13,19,25,31 sceniros.
    #if( (sc==13)||(sc==19)||(sc==25)||(sc==31) ){
    #count.typeI[j] = 1
    #}
    ## for sc=14,20,26,32, just dose 2,3,4,5 are 'good' doses for computing power.
    #if( (sc==14)||(sc==20)||(sc==26)||(sc==32) ){
    #count.typeI[j] = ifelse( sum(select[c(2,3,4,5)])!=0, 1, 0  )
    #}
    ## for sc=15,21,27,33, just dose 1,2,3,4 are 'good' doses for computing power.
    #if( (sc==15)||(sc==21)||(sc==27)||(sc==33) ){
    #count.typeI[j] = ifelse( sum(select[c(1,2,3,4)])!=0, 1, 0  )
    #}
    ## for sc=16,22,28,34, just dose 2,3,4 are 'good' doses for computing power.
    #if( (sc==16)||(sc==22)||(sc==28)||(sc==34) ){
    #count.typeI[j] = ifelse( sum(select[c(2,3,4)])!=0, 1, 0  )
    #}
    ## for sc=17,23,29,35, just dose 1,2,4,5 are 'good' doses for computing power.
    #if( (sc==17)||(sc==23)||(sc==29)||(sc==35) ){
    #count.typeI[j] = ifelse( sum(select[c(1,2,4,5)])!=0, 1, 0  )
    #}
    ## for sc=18,24,30,36, all the doses 1,2,3,4,5 are 'good' doses for computing power.
    #if( (sc==18)||(sc==24)||(sc==30)||(sc==36) ){
    #count.typeI[j] = 1 ##ifelse( sum(select[c(1,2,3,4,5)])!=0, 1, 0  )
    #}
    #if( (sc==37) || (sc==38) ){
    #count.typeI[j] = ifelse( sum(select[c(2,3)])!=0, 1, 0  )
    #}


    #} ## Ends if ( any(select!=0) )r.
    ###################################################################################################################
    if(sum(select) > 1){
      ## toxicity
      tox_prob_hat <- (t.a-tox_a)/d.cum.size
      eff_prob_hat <- (d.a-eff_a)/d.cum.size
      utility_vec <- eff_prob_hat - weight1 * tox_prob_hat - weight2 * tox_prob_hat * (tox_prob_hat > pT)
      max_utility <- max(utility_vec[select == 1])
      optimal_dose <- which(utility_vec == max_utility)
      select[select == 1] <- 0
      select[optimal_dose] <- 1
    }
    dose.select[j,] = select

    #t.cum.size = sum(d.cum.size) + p.cum.size
    #     over_toxic_vector[j] = over_toxic ## date:12-24-2012
    #     over_toxic_rate[j] = over_toxic_vector[j]/t.cum.size  ## date:12-24-2012

    final.tol.sampsize[j] = t.cum.size
    #      cat("t.cum.size is ",t.cum.size,"\n")
    final.dose.sampsize[j,] = d.cum.size
    final.plac.sampsize[j] = p.cum.size
    final.count[j,set2] = 1  ##count

    #   a = rowSums(final.count)
    #   b = final.count/a
    tox.count[j,] = t.a-tox_a # the toxic responses upon 5 dose levels on the j-th simulation 2-2-2013



  } ## Ends the for(j in 1:Nsim) loop.
  final.tol.sampsize.avg <- mean(final.tol.sampsize)
  final.dose.sampsize.avg <- apply(final.dose.sampsize, 2, mean)
  final.plac.sampsize.avg <- mean(final.plac.sampsize)
  tox.count.avg <- apply(tox.count, 2, mean)
  dose.select.avg <- apply(dose.select, 2, mean)
  #tox_safety <- tox.count/final.dose.sampsize <= pT
  #tox_safety[is.na(tox_safety)] <- FALSE
  tox_safety <- p.tox <= pT
  if(any((p.d > p.p) & (p.tox <= pT))){
    dose.select_f <- dose.select * (matrix(rep((p.d > p.p), Nsim), nrow = Nsim, byrow = TRUE) & matrix(rep(tox_safety, Nsim), nrow = Nsim, byrow = TRUE))
    type1_power <- sum(apply(dose.select_f, 1, function(x){if(any(x == 1)){return(TRUE)} else{return(FALSE)}}))/Nsim
    return(list(power = type1_power, tot_sampsize_avg = final.tol.sampsize.avg,
                dose_sampsize_avg = final.dose.sampsize.avg,
                placebo_sampsize_avg = final.plac.sampsize.avg,
                selection_percentage = dose.select.avg,
                tox_count_avg = tox.count.avg))
  }
  else{
    type1_power <- sum(apply(dose.select, 1, function(x){any(x == 1)}))/Nsim
    return(list(typeI = type1_power, tot_sampsize_avg = final.tol.sampsize.avg,
                dose_sampsize_avg = final.dose.sampsize.avg,
                placebo_sampsize_avg = final.plac.sampsize.avg,
                selection_percentage = dose.select.avg,
                tox_count_avg = tox.count.avg))
  }
}
