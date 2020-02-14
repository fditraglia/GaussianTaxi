## code to prepare `taxi` dataset goes here
library(tidyverse)
library(haven)

taxi <- read_dta('./data-raw/taxi034.dta')
cpi <- read_csv('./data-raw/CPI.NY.NJ.PA.Urban.All.Items.csv')

taxi <- taxi %>%
  mutate(year = year(date), month = month(date))

cpi <- cpi %>%
  rename(year = Year) %>%
  mutate(month = as.numeric(str_replace(Period, "M", "")))

# Base period for price index: Jun '99 = 1
base_value <- cpi %>%
  filter(year == 1999, month == 6) %>%
  pull(Value)

cpi <- cpi %>%
  mutate(price_index = Value / base_value) %>%
  select(year, month, price_index)

taxi <- left_join(taxi, cpi)
taxi <- taxi %>%
  rename(price = price_index) %>%
  mutate(earnings = round(income, 2), # Earnings don't quite make sense - not not expressed in $/cents
         driving = ttripmin / 60,
         breaks = tbrmin / 60,
         waiting = twaitmin / 60,
         working = driving + waiting,
         leisure = 24 - working,
         wage = earnings / driving,
         consumption = earnings / price) %>%
  select(id, date, shift, earnings, driving, breaks, waiting, working, leisure, wage, consumption)

usethis::use_data(taxi, overwrite = TRUE)
rm(list = ls())
