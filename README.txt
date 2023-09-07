1. 제출파일 목록
  1) README.txt : 코드 실행 방법 설명
  2) 1_data_download.r : 샘플의 나이와 건강 데이터셋 다운받는 과정을 담은 파일
  3) 2_preprocessing.ipynb : 샘플의 나이와 건강 데이터셋 전처리 과정을 담은 파일
  4) 3_Modeling.r : DNA 메틸화 데이터 전처리 과정과 머신러닝 모델을 돌려보는 파일
  5) 4_BDI_and_AgeDiff_Correlation_Analysis.ipynb : 예측 데이터를 통한 분석 과정을 담은 파일
  6) GSE72680_beta_values.txt : DNA 메틸화 데이터에 대한 파일
  7) Sample_data.zip : 서버에 업로드 하기위해 크롤링한 raw데이터들을 압축한 파일
  8) sample_data.csv : raw 데이터들을 데이터프레임화한 샘플의 나이와 건강정보를 담은 파일

2. Data 획득 방법
우리는 GEO에서 공개하는 Grady Trauma 프로젝트를 통한 아프리카계 미국인의 DNA 메틸화 데이터를 사용한다.
  1)첫번째 data(GSE72680_beta_values.txt)
 - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72680에서 DNA 메틸화 데이터인 beta_values 데이터를 다운로드 받고 파일로 저장한다.
 - 위의 방법으로 데이터 다운로드하면 시간이 오래걸리기 때문에, 미리 다운받은 데이터 또한 제출파일에 포함한다.
 2)두번째 data(sample_data.csv)
 - R 코드(1_data_download.r)를 통해 샘플의 건강정보를 가져온다.
 - 파이썬 코드(2_preprocessing.ipynb)를 통해 전처리한 샘플의 건강정보를 가져온다.
 - 위의 방법으로 데이터 다운로드하면 시간이 오래걸리기 때문에, 미리 다운받은 데이터 또한 제출파일에 포함한다.

3. 실행 필요 프로그램과 라이브러리
  1) python
  2) R
  3) RStudio
  4) jupyter notebook 혹은 jupyter lab
  5) 코드 실행에 필요한 Python 라이브러리
     - csv
     - os
     - pandas
     - seaborn 
     - matplotlib.pyplot
     - collections
     - scipy.stats
     - scipy
    파이썬 라이브러리들은 코드를 실행하면 설치할 수 있음(별도 설치 불필요)
  6)코드 실행에 필요한 R 라이브러리
     - WGCNA 
     - sqldf
     - dplyr
     - tidyverse
     - glmnet
     - caret
     - data.table
     - BiocManager
     - GEOquery
     - ggplot2
    R 라이브러리들은 코드를 실행하면 설치할 수 있음(별도 설치 불필요)

4. 실행방법

*1,2번 파일의 경우 데이터의 크기가 매우 커서 실행에 많은 시간이 소요되기 때문에 처리가 완료된 파일을 같이 첨부한다. 
 -> sample_data.csv 는 1번을 통해 크롤링된 raw 데이터들을 2번을 통해 분석을 위해 데이터프레임화한 파일이다. 

  제출한 데이터 파일(GSE72680_beta_values.txt, sample_data.csv)을 사용
  1) 코드 파일(*.R)과 데이터 파일(GSE72680_beta_values.txt,sample_data.csv)들이 동일한 작업 폴더 안에 위치한지 확인한다.
  2) RStudio를 실행하고,
      - 3_Modeling.r
      을 실행한다.
  3) 4_BDI_and_AgeDiff_Correlation_Analysis.ipynb 의 경우 이전 단계에서 완성된 predicted(lasso).csv 파일을 이용해 분석을 수행한다. 


5. 참조 논문 https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-10-r115#Abs1
