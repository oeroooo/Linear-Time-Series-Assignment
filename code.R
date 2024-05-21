# Linear Time Series Assignment
# Monthly series of the fabrication of french weapons and ammo from 1990 to 2024.
# Data is accessible on the INSEE website.

# Load required libraries
library(dplyr)
library(kableExtra)
library(astsa)
library(forecast)
library(portes)
library(fUnitRoots)
library(tseries)


# Read and preprocess the time series data
Time_series <- read.csv("data/valeurs_mensuelles.csv", sep = ";") %>%
  filter(!row_number() %in% c(1, 2, 3))
names(Time_series) <- c("date", "value", "codes")
Time_series$date <- as.Date(paste(Time_series$date, "-01", sep=""))
Time_series$value <- as.numeric(Time_series$value)


# Plotting the initial series
plotly::plot_ly(x = Time_series$date, y = Time_series$value, type = "scatter", mode = "lines") %>%
  plotly::layout(title = "Monthly series of the fabrication of french weapons and ammo from 1990 to 2024",
                 xaxis = list(title = "Date"),
                 yaxis = list(title = "Value"))

# Check for stationarity and differentiate 
diff <- diff(Time_series$value)
plotly::plot_ly(x = Time_series$date[-1], y = diff, type = "scatter", mode = "lines") %>%
  plotly::layout(title = "Log return of the series",
                 xaxis = list(title = "Date"),
                 yaxis = list(title = "Value"))

# ACF and PACF plots
png("images/acf.png")
acf(diff)
dev.off()
png("images/pacf.png")
pacf(diff)
dev.off()

# Stationarity tests
adf <- adfTest(diff)
pp <- pp.test(diff)
kpss <- kpss.test(diff)

# Display results of stationarity tests
data.frame(Test = c("Augmented Dickey-Fuller", "Phillips-Perron", "KPSS"), 
           Statistics = c(adf@test$statistic, pp$statistic, kpss$statistic), 
           lag_order = c(adf@test$parameter, pp$parameter, kpss$parameter), 
           P_value = c(paste0("<= ", adf@test$p.value), paste0("<= ", pp$p.value), paste0(">= ", kpss$p.value)))

# Selection and validation of the ARMA model
# Generate AIC and BIC values for ARMA models from order (0,0) to (6,6)
AIC_BIC <- data.frame(p = integer(), q = integer(), AIC = numeric(), BIC = numeric())
for (i in 0:6) {
  for (j in 0:6) {
    model <- arima(diff, order = c(i, 0, j))
    AIC_BIC <- rbind(AIC_BIC, data.frame(p = i, q = j, AIC = AIC(model), BIC = BIC(model)))
  }
}

# Determine validity of ARMA models based on significance of coefficients and Ljung-Box test
validity <- data.frame(p = integer(), q = integer(), Valid_AR_MA = character(), Valid_Ljung_Box = character())
for (i in 0:6) {
  for (j in 0:6) {
    model <- arima(diff, order = c(i, 0, j))
    coef_model <- coef(model)
    vcov_model <- sqrt(diag(vcov(model)))
    AR_valid <- ifelse(abs(coef_model[paste0("ar", i)] / vcov_model[paste0("ar", i)]) > 1.96, "YES", "NO")
    MA_valid <- ifelse(abs(coef_model[paste0("ma", j)] / vcov_model[paste0("ma", j)]) > 1.96, "YES", "NO")
    Valid_AR_MA <- if(AR_valid == "YES" && MA_valid == "YES") "YES" else "NO"
    Ljung_Box_p_value <- Box.test(model$residuals, lag = max(5, i, j), type = "Ljung-Box")$p.value
    Valid_Ljung_Box <- ifelse(Ljung_Box_p_value > 0.05, "YES", "NO")
    validity <- rbind(validity, data.frame(p = i, q = j, Valid_AR_MA = Valid_AR_MA, Valid_Ljung_Box = Valid_Ljung_Box))
  }
}

# Filter and display valid models
valid_models <- validity %>%
  filter(Valid_AR_MA == "YES" && Valid_Ljung_Box == "YES")
kable(valid_models) %>%
  kableExtra::kable_styling(bootstrap_options = "striped", full_width = F)

# Model significance testing
model <- arima(diff, order = c(4, 0, 5))
p_value <- (1 - pnorm(abs(coef(model)) / sqrt(diag(model$var.coef)))) * 2
p_value %>% as.data.frame() %>% setNames(c("p_value")) %>% kable() %>% kableExtra::kable_styling(bootstrap_options = "striped", full_width = F)

model <- arima(diff, order = c(2, 0, 3))
p_value <- (1 - pnorm(abs(coef(model)) / sqrt(diag(model$var.coef)))) * 2
p_value %>% as.data.frame() %>% setNames(c("p_value")) %>% kable() %>% kableExtra::kable_styling(bootstrap_options = "striped", full_width = F)
