library(argparse)
library(lmerTest)
library(tidyverse)

parser <- ArgumentParser(description="reQTL-mapping")

parser$add_argument("-i", "--input", help="each snp gene pair file", required=TRUE)
parser$add_argument("-o", "--output", help="output file", required=TRUE)

args <- parser$parse_args()

# 读取数据
data_need <- read.table(args$input, header=TRUE, sep="\t",stringsAsFactors = T)

gene <- colnames(data_need)[3]
snp <- colnames(data_need)[2]

# 修改列名
colnames(data_need) <- c('sample','ciseqtl','cisegene','cov1','cov2','cov3','cov4','cov5','cov6','cov7','cov8',
                    'cov9','cov10','cov11','cov12','cov13','cov14','cov15','env')


# 建立混合线性模型
model <- lmer(data = data_need,
              formula = cisegene ~ ciseqtl+env+cov1+cov2+cov3+cov4+cov5+cov6+cov7+cov8+cov8+cov10+cov11+cov12+cov13+cov14+cov15+
                cov1*env+cov2*env+cov3*env+cov4*env+cov5*env+cov6*env+cov7*env+cov8*env+cov9*env+cov10*env+cov11*env+cov12*env+
                cov13*env+cov14*env+cov15*env+ciseqtl*env+(1|sample))

p <- anova(model)["ciseqtl:env", "Pr(>F)"]

# 输出结果
write.table(data.frame(gene, snp, p), file=args$output, sep="\t", quote=FALSE, row.names=FALSE,col.names = FALSE)
