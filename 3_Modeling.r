#필요한 패키지 설치 
install.packages("BiocManager")
install.packages("glmnet")
install.packages('sqldf')
install.packages("RPMM")
install.packages("tidyverse")
install.packages("caret")
install.packages("xgboost")
install.packages("lightgbm")
install.packages("ggplot2")
install.packages("data.table")
BiocManager::install("WGCNA")
BiocManager::install("GEOquery")

# 필요한 패키지 로드
library(WGCNA)  # Weighted gene co-expression network analysis
library(sqldf)  # SQL 데이터베이스 쿼리를 사용하여 데이터 처리
library(dplyr)
library(tidyverse)
library(glmnet)
library(caret)
library(xgboost)
library(data.table)
library(lightgbm)
library(ggplot2)

setwd("./")


#데이터 불러오기 
txt_file_path <- './GSE72680_beta_values.txt'
data <- fread(txt_file_path)

# 제거할 특정 단어를 지정
word_to_remove <- "Detection"

# 특정 단어가 포함된 열을 제거
data_cleaned <- data %>%
  select(-matches(word_to_remove))

#열이름 변경
data_cleaned <- data_cleaned %>%
  rename(ProbeID= 1)

#추가데이터 csv파일로 저장 
new_csv_file_path <- "./GSE72680_beta_values.csv"
write_csv(data_cleaned, new_csv_file_path)


# 결측치가 있는 행 제거
dat0 <- data_cleaned[complete.cases(data_cleaned), ]
print(dat0)

#결측값 확인
na_count <- colSums(is.na(dat0))
print(na_count)


#열이름 맞추기 
start_number <- 1868036
target_cols <- 2:ncol(dat0)
new_col_names <- paste0("GSM", seq(start_number, length.out = length(target_cols)))
colnames(dat0)[target_cols] <- new_col_names
print(dat0)




#-------------------------------------------------------------------------------
# ---***---***---***---***---***모델링 준비***---***---***---***---***---***---


datSample = read.csv("sample_data.csv")


datSample <- datSample %>%
  filter(!(datSample[,1] == "GSM1868428.soft"))
datSample  <- datSample %>%
  filter(!(datSample[,1] == "GSM1868429.soft"))


# 데이터 프레임 정의 
meth26 <- dat0 %>% as.data.frame()

rownames(meth26)[1]<- 'id'
colnames(datSample)[1] <- 'ProbeID'
datSample$ProbeID <- gsub("\\.soft", "", datSample$ProbeID)


#데이터 전치 
rownames(meth26) <- meth26$ProbeID  # 'ProbeID' 컬럼 값을 새로운 행 이름으로 사용
meth26_transposed <- t(meth26)
meth26_transposed <- as.data.frame(meth26_transposed)



# Convert rownames to a column
meth26_transposed <- tibble::rownames_to_column(meth26_transposed, var = "ProbeID")

# Now try the join operation again
datAll <- datSample %>% right_join(meth26_transposed, by = 'ProbeID')

#datAll에서 cg로 시작하는 열과, 그 열의 age를 선택해서 넣어준다.
datSe <- datAll %>% select(age, starts_with("cg"))

# "Age" 열에서 결측값이 있는 행의 인덱스 찾기
na_indices <- which(is.na(datSe$age))
print(na_indices)

#결측값 있는 행들 제거 
datSe_no_na <- datSe[-na_indices, ]

#숫자형으로 변환 
datSe_no_na <- lapply(datSe_no_na, as.numeric) %>% as.data.frame()



# x_matrix에 Age를 제외한 열로 구성된 행렬 만들기.
x_matrix <- as.matrix(datSe_no_na[, -which(names(datSe_no_na) == 'age')])

# y_vector 변수에는 "Age" 열의 데이터가 저장
y_vector <- datSe_no_na$age 


# 데이터를 훈련 데이터와 테스트 데이터로 나눔 (70% 훈련, 30% 테스트)
set.seed(123)  # 재현성을 위한 시드 설정

index <- createDataPartition(y_vector, p = 0.7, list = FALSE)
x_train <- x_matrix[index, ]
y_train <- y_vector[index]
x_test <- x_matrix[-index, ]
y_test <- y_vector[-index]





########################## M O D E L I N G #############################
#------------Ridge------------
lambda_values <- seq(from = 0.1, to = 1, by = 0.1)
mse_values <- numeric(length(lambda_values))

for (i in seq_along(lambda_values)) {
  l <- lambda_values[i]
  
  fit <- glmnet(x_train, y_train, alpha = 0, lambda = l)
  
  pred <- predict(fit, newx = x_test, s = l)
  
  mse <- mean((y_test - pred)^2)
  
  print(paste("Ridge Regression - When lambda =", l,", Mean Squared Error:", mse))
  
}



#------------Lasso------------


# Lasso Regression
lambda_values <- seq(from = 0.1, to = 1, by = 0.1)
mse_values_lasso <- numeric(length(lambda_values))
result_lasso <- list()

for (i in seq_along(lambda_values)) {
  l <- lambda_values[i]
  
  fit_lasso <- glmnet(x_train, y_train, alpha = 1, lambda = l)
  
  pred_lasso <- predict(fit_lasso, newx = x_test, s = l)
  
  # Lasso 회귀 모델의 계수 확인
  coef(fit_lasso, s = l) 
  
  mse_lasso <- mean((y_test - pred_lasso)^2)
  
  print(paste("Lasso Regression - When lambda =", l,", Mean Squared Error:", mse_lasso))
  
  mse_values_lasso[i] <- mse_lasso
  
  result_lasso[[as.character(l)]] <- data.frame(pred = as.vector(pred_lasso), y_test = y_test)
}


#------------Elastic net------------

# Elastic net
lambda_values <- seq(from = 0.1, to = 1, by = 0.1)
mse_values_elastic <- numeric(length(lambda_values))
result_elastic <- list()

for (i in seq_along(lambda_values)) {
  l <- lambda_values[i]
  
  fit_elastic <- glmnet(x_train, y_train, alpha = 0.5, lambda = l)
  
  pred_elastic <- predict(fit_elastic, newx = x_test, s = l)
  
  # Lasso 회귀 모델의 계수 확인
  coef(fit_elastic, s = l) 
  
  mse_elastic <- mean((y_test - pred_elastic)^2)
  
  print(paste("Elastic net - When lambda =", l,", Mean Squared Error:", mse_elastic))
  
  mse_values_elastic[i] <- mse_elastic
  
  result_elastic[[as.character(l)]] <- data.frame(pred = as.vector(pred_elastic), y_test = y_test)
}



#------------LightGBM------------
#  데이터셋 생성
train_data <- lgb.Dataset(data = as.matrix(x_train), label = y_train)
test_data <- lgb.Dataset(data = as.matrix(x_test), label = y_test, reference = train_data)

# 샘플링해서 로딩
sample_size <- 100  # 원하는 샘플 크기 지정 가능
sample_size <- min(sample_size, nrow(x_train))  # 최대값 조정

sample_indices <- sample(1:nrow(x_train), size = sample_size)

x_train_sampled <- x_train[sample_indices, ]
y_train_sampled <- y_train[sample_indices]

# 데이터셋 생성
train_data_sampled <- lgb.Dataset(data = as.matrix(x_train_sampled), label = y_train_sampled)
test_data <- lgb.Dataset(data = as.matrix(x_test), label = y_test, reference = train_data_sampled)

# LightGBM 파라미터 설정
params <- list(
  objective = "regression",
  learning_rate = 0.1,
  max_depth = -1,
  metric = "mse",
  lambda_l1=0.01,
  lambda_l2=0.01
)

# LightGBM 모델 학습 (샘플링된 데이터 활용)
lgb_model <- lgb.train(params,
                       data=train_data_sampled,
                       nrounds=100)
# 예측 수행 (전체 테스트 세트 활용)
pred <- predict(lgb_model, data = as.matrix(x_test))

mse <- mean((y_test - pred)^2)
print(paste("LightGBM - Mean Squared Error:", mse))



#------------Lasso with Best Lambda through Cross-Validation ------------


# 학습  
fit_lasso <- glmnet(x_train, y_train, alpha = 1)  

# Cross-validation for lambda selection
cv_fit <- cv.glmnet(x_train, y_train, alpha = 1)
plot(cv_fit, main = "Cross-validation for Lambda Selection(Lasso)")


# Optimal lambda where the cross-validation error is minimum
opt_lambda <- cv_fit$lambda.min
opt_lambda_rounded <- round(opt_lambda, 2)  # 소수점 둘째 자리까지 반올림
print(paste("Optimal Lambda:", opt_lambda_rounded))

# Coefficients at optimal lambda
opt_coef <- coef(fit_lasso, s = opt_lambda)


# 예측 수행 with optimal lambda
pred_lasso <- predict(fit_lasso, newx = x_test, s = opt_lambda)
mse <- mean((y_test - pred_lasso)^2)
print(paste("Mean Squared Error:", mse))

# Make a data frame for plotting
plot_df <- data.frame(pred = as.vector(pred_lasso), actual = y_test)

ggplot(plot_df, aes(x = pred, y = actual)) +
  geom_point(color = 'red') +
  stat_smooth(method = "lm", col = "black") +
  labs(x = "Predicted Age",
       y = "Chronological Age",
       title = paste0("Lasso Regression (lambda=", opt_lambda_rounded,", MSE=", round(mse, 2),")")) +
  theme(plot.title = element_text(hjust = 0.5)) 





#-----------결과 도출 ------------


# 전체 데이터셋에 대한 예측 계산
pred_all <- predict(fit_lasso, newx = x_matrix, s = opt_lambda)

#저장할 데이터셋 전처리
num_rows <- nrow(datSample)
datSample_no_na <- datSample[-((num_rows-34):num_rows), ] # 결측값 대상이었던 행을 제외 

datSample_no_na$DNAmAge_lasso  <- pred_alls

# CSV 파일로 저장
write.csv(datSample_no_na, file = "predicted(lasso).csv", row.names = FALSE)







