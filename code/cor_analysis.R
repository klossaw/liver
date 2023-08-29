library(corrplot)

cor(x,y=NULL,use="everything",method= c("pearson","kendall","spearman"))


cor_matrix <- t(matrix_C1_neg)
tmp <- cor(cor_matrix,method = 'pearson') %>% t() 
tmp
corrplot(tmp[1:100,1:100],type="upper")

# Update of correlation analysis
cor_res <- psych::corr.test(t(avg_samp), method = 'pearson', adjust = BH)
cor_mtx <- cor_res@r
corrplot::corrplot(cor_mtx, method = 'circle', type = 'upper', )
