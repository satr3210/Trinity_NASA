library(tidyverse)

import_data <- function(groups) {
    df <- read.table('cleaned.csv',header=TRUE, sep=",") %>%
     filter( group %in% groups) %>%
     transmute(
          group=as.factor(group),
          id=as.factor(id),
          session,
          trial,
          progress = ifelse(session==4, trial/10, (trial + 10 * (session-1))/30),
          ls_skill,
          ls_success = ifelse(ls_performance == -1, FALSE, TRUE),
          mc_skill,
          mc_success = ifelse(mc_performance == -1, FALSE, TRUE),
          de_skill,
          de_success = ifelse(de_performance == -1, FALSE, TRUE),
          integrated_skill = ls_skill + mc_skill + de_skill,
          mission_success = (
               (ls_performance != -1) &
               (mc_performance != -1) & 
               (de_performance != -1) &
               !(
                    ls_performance == mc_performance &
                    mc_performance == de_performance & 
                    de_performance == 0
               )
          )
     )
}