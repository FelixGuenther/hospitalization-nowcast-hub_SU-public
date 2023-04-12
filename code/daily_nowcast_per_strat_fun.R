#' Prepare data
prepare_data = function(dat, now, start_date, max_delay) {
  # Reportgin triangle funs
  prep_rep_tri = function(age_gr="00-04", loc="DE", start_date, dat, now, max_delay) {
    rep_tri = dat %>% filter(location==loc, 
                             age_group==age_gr,
                             date>=start_date,
                             date<=now) %>% 
      arrange(date) %>% select(starts_with("value")) %>% 
      as.matrix()
    rep_tri_rest = rep_tri[,1:(max_delay+1)]
    rep_tri_rest[,max_delay+1] = rowSums(rep_tri[,(max_delay+1):ncol(rep_tri)], na.rm = T)
    rep_tri_rest[is.na(rep_tri_rest)] = 0
    rep_tri_rest
  }
  prep_rep_tri_7d = function(age_gr="00-04", loc="DE", start_date, dat, now, max_delay) {
    rep_tri = dat %>% filter(location==loc, 
                             age_group==age_gr,
                             date>=start_date - 6,
                             date<=now) %>% 
      arrange(date) %>% select(starts_with("value")) %>% 
      as.matrix()
    per_len = as.numeric(now-start_date+1)
    rep_tri_rest = rep_tri[,1:(max_delay+1)]
    rep_tri_rest[,max_delay+1] = rowSums(rep_tri[,(max_delay+1):ncol(rep_tri)], na.rm = T)
    rep_tri_7d = matrix(0, nrow=per_len, max_delay + 1)
    
    for(d in 0:max_delay) {
      if (d==0) {
        for(t in 1:(per_len-d)) {
          for(k in t:(t+6)) {
            rep_tri_7d[t,d+1] = rep_tri_7d[t,d+1] + sum(rep_tri_rest[k,1:(t+6+1-k)])
          }
        }
      } else {
        for(t in 1:(per_len-d)) {
          for(k in t:(t+6)) {
            rep_tri_7d[t,d+1] = rep_tri_7d[t,d+1] + ifelse(t+6+1-k+d<=max_delay+1, rep_tri_rest[k,min(t+6+1-k+d, max_delay+1)], 0)
          }
        }
      }
    }
    rep_tri_7d[is.na(rep_tri_7d)] = 0
    rep_tri_7d
  }
  
  rep_tri_age_list = lapply(c("00-04", "05-14", "15-34", "35-59", "60-79", "80+"), 
                            function(age) {
                                prep_rep_tri(age_gr=age,
                                             loc="DE",
                                             start_date = start_date, 
                                             dat = dat, now=now, max_delay = max_delay)})
  names(rep_tri_age_list) = c("00-04", "05-14", "15-34", "35-59", "60-79", "80+")
  rep_tri_7d_age_list = lapply(c("00-04", "05-14", "15-34", "35-59", "60-79", "80+"), 
                            function(age) {
                              prep_rep_tri_7d(age_gr=age,
                                           loc="DE",
                                           start_date = start_date, 
                                           dat = dat, now=now, max_delay = max_delay)})
  rep_tri_fed_list = lapply(c("DE-BB", "DE-BE", "DE-BW", "DE-BY", "DE-HB", "DE-HE", 
                             "DE-HH", "DE-MV", "DE-NI", "DE-NW", "DE-RP", "DE-SH", 
                             "DE-SL", "DE-SN", "DE-ST", "DE-TH"),
                           function(region) {
                             prep_rep_tri(age_gr = "00+",
                                          loc = region,
                                          start_date = start_date, 
                                          dat = dat, now=now, 
                                          max_delay = max_delay)})
  names(rep_tri_fed_list) = c("DE-BB", "DE-BE", "DE-BW", "DE-BY", "DE-HB", "DE-HE", 
                              "DE-HH", "DE-MV", "DE-NI", "DE-NW", "DE-RP", "DE-SH", 
                              "DE-SL", "DE-SN", "DE-ST", "DE-TH")
  
  rep_tri_7d_fed_list = lapply(c("DE-BB", "DE-BE", "DE-BW", "DE-BY", "DE-HB", "DE-HE", 
                              "DE-HH", "DE-MV", "DE-NI", "DE-NW", "DE-RP", "DE-SH", 
                              "DE-SL", "DE-SN", "DE-ST", "DE-TH"),
                            function(region) {
                              prep_rep_tri_7d(age_gr = "00+",
                                           loc = region,
                                           start_date = start_date, 
                                           dat = dat, now=now, 
                                           max_delay = max_delay)})
  r_age = abind::abind(rep_tri_age_list, rev.along=0)
  r_7d_age = abind::abind(rep_tri_7d_age_list, rev.along=0)
  
  r_fed = abind::abind(rep_tri_fed_list, rev.along=0)
  r_7d_fed = abind::abind(rep_tri_7d_fed_list, rev.along=0)
  
  r_ov = array(prep_rep_tri(age_gr="00+",
                      loc="DE",
                      start_date = start_date, 
                      dat = dat, now=now, max_delay = max_delay),
               c(as.numeric(now-start_date)+1, 
                 max_delay+1, 1))
  r_7d_ov = array(prep_rep_tri_7d(age_gr="00+",
                            loc="DE",
                            start_date = start_date, 
                            dat = dat, now=now, max_delay = max_delay),
               c(as.numeric(now-start_date)+1, 
                 max_delay+1, 1))
  days = seq(start_date, now, "1 day")
  x_lambda = model.matrix(~-1+wd, 
                          data = tibble(wd=as.factor(
                            wday(days))))[,-1]
  # Design array hazard model
  wdays <- 2:7
  x_haz_wd <- array(NA, dim=c(max_delay, length(days), length(wdays)),
                 dimnames=list(paste("delay",0:(max_delay-1),sep=""),
                               as.character(days), 
                               wdays))
  # Loop over all times and lags
  for (t in seq_len(length(days))) {
    for (w in seq_len(length(wdays))) {
      x_haz_wd[,t, w] <- as.numeric(lubridate::wday(days[t] + 0:(max_delay-1), label=FALSE) == wdays[w])
    }
  }
  
  x_haz_week = model.matrix(~week, data=tibble(week=factor(rev(rep(1:(length(days)/7),each=7)))),
                            contrasts.arg = list(week = "contr.sum"))[,-1]


  list(r_ov=r_ov,
       r_age=r_age,
       r_age_list=rep_tri_age_list,
       r_fed=r_fed,
       r_fed_list=rep_tri_fed_list,
       x_lambda=x_lambda,
       x_haz_wd=x_haz_wd, 
       x_haz_week=x_haz_week,
       r_7d_ov=r_7d_ov,
       r_7d_age=r_7d_age,
       r_7d_fed=r_7d_fed)
}


summarise_fit = function(model_path, now, start_date, max_delay) {
  group=strsplit(model_path, "_")[[1]][4]
  age_group = ifelse(group %in% c("00+", "00-04", "05-14", "15-34", "35-59", "60-79", "80+"), group, "00+")
  location = ifelse(group %in% c("00+", "00-04", "05-14", "15-34", "35-59", "60-79", "80+"), "DE", group)
  
  print(paste0(group))
  fit = readRDS(model_path)
  print(paste0("# divergent transitions: ",fit$sampler_diagnostics(format = "draws_df") %>% summarise(sum(divergent__)) %>% unlist()))
  
  res_file = posterior::summarise_draws(fit$draws("N_7d"),
                                        ~quantile(.x, probs = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)),
                                        mean) %>%
    mutate(target_end_date=rep(seq(now-max_delay+1, now, "1 day"), times=1),
           age_group = age_group,
           location=location,
           forecast_date=now,
           diff=as.numeric(target_end_date-forecast_date)) %>%
    filter(diff>=-28) %>%
    mutate(target=paste0(diff, " day ahead inc hosp")) %>%
    select(-variable) %>%
    pivot_longer(cols=c(`2.5%`, `10%`, `25%`, `50%`, `75%`, `90%`, `97.5%`,  mean), names_to = "type") %>%
    mutate(quantile=factor(type, levels = c("2.5%", "10%", "25%", "50%", "75%", "90%", "97.5%"),
                           labels = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)),
           type = factor(type, levels = c("2.5%", "10%", "25%", "50%", "75%", "90%", "97.5%", "mean"),
                         labels = c(rep("quantile", 7), "mean")),
           pathogen="COVID-19") %>%
    select(location, age_group, forecast_date, target_end_date,target,type,quantile,value,pathogen)
  res_file
}
