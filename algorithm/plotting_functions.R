### Plotting functions ###

# Author: Malte Lüken
# Date: 26.04.2020


plot_samples_time <- function(object, from = 0, to = Inf, display = c("x", "y", "vel", "acc", "angle"), add.event = F) {
  
  require(ggplot2)
  require(dplyr)
  
  data_long <- object$samples %>% 
    gather(all_of(display), key = metric, value = value, factor_key = T) %>% 
    mutate_at(vars(label, state), as.factor) %>%
    dplyr::filter(t > from, t <= to)

  p <- ggplot(data_long , aes(x = t, y = value)) + facet_wrap(vars(metric), scales = "free") +
    geom_path() + xlab("time in s") + ylab("metric") +
    geom_point(aes(color = label, shape = state))
  
  if(add.event) p <- p + geom_text(aes(label = event, hjust = 1, vjust = 1, size = 3))
    
  return(p)
}


plot_samples_xy <- function(object, from = 0, to = Inf, add.event = F, add.fixpos = F) {
  
  require(ggplot2)
  require(dplyr)
  
  data <- object$samples %>% 
    dplyr::filter(t > from, t <= to, label != 0) %>%
    mutate_at(vars(label, state), as.factor) %>% 
    group_by(event) %>%
    mutate(fix.x = ifelse(label == 1, mean(x[!is.na(vel)], trim = 0.2), NA),
           fix.y = ifelse(label == 1, mean(y[!is.na(vel)], trim = 0.2), NA),)
  
  p <- ggplot(data, aes(x, y)) + 
    geom_point(aes(color = label, shape = state)) +
    geom_path()
  
  if(add.event) p <- p + geom_text(aes(label = event), hjust = 1, vjust = 1, size = 3)
  
  if(add.fixpos) p <- p + geom_point(aes(fix.x, fix.y), color = "green1")
  
  return(p)
}


plot_posterior_time <- function(object, from = 0, to = Inf, prob = T) {
  
  require(ggplot2)
  require(dplyr)
  
  data <- cbind(object$samples, object$model@posterior[,-1]) %>%
    gather(contains("S"), key = "state", value = "post.prob") %>%
    dplyr::filter(state != "state", t > from, t <= to)
  
  entropy <- data %>% 
    group_by(t) %>% 
    summarise(entropy = -sum(post.prob*ifelse(post.prob == 0, 0, log(post.prob))))
    
  
  p <- ggplot(data, aes(x = t)) + xlab("Time in s")
    
  if(!prob) { 
    p <- p + geom_line(data = entropy, aes(x = t, y = entropy)) + ylab("Entropy") 
  } else {
    p <- p + geom_line(aes(y = post.prob, color = state)) + ylab("Posterior\nprobability") 
  }
  
  return(p)
}
