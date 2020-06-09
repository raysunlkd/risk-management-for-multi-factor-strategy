Our risk analysis result is in the file 'Risk_output_Antivirus.xlsx'!!!

1. At first, I use jupyter notebook 'Get_data.ipynb' get our basic datas for CSI300 and its constituent stocks.

2. I get data by tushare and it is hard to get all the great data we need, hence my original model can not be accomplished 
  and I deal the data with a simpler model (We just do the PCA analysis on the original 21 alpha factors and simplely selected 6 factors from them).

3. I process those basic datas in jupyter notebook 'Data_preprocess.ipynb' to get the datas for R codes and Risk_output.xlsx.

4. Also, I select 6 alpha factors from 21 alpha factor, get the composited factor, 
    caculate the factor scores and so on in jupyter notebook 'Data_preprocess.ipynb'

5. I use the R file 'Get_factor_returns.R' (which is just the same with the 'testStudent.R' shared by TA) 
   to caculate the factor return that is used in risk analysis.

6. I accomplish our risk analysis by the almost same way lecturer did in Rmarkdown file 'strategy.rmd'.
   Also, I accomplish the Risk_ouput.xlsx by the same way lecturer did to finish our risk analysis.

7. Maybe 'strategy.html' and those jupyter notebook files are easier to read.




