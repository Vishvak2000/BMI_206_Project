# Load necessary libraries
library(dplyr)
library(SCENT)
library(glmnet)
library(Matrix)

# Read in your data
fibroblast_SCENT <- readRDS("SCENT_2_covariate.rds")
rna_matrix <- fibroblast_SCENT@rna
atac_matrix <- fibroblast_SCENT@atac
atac_matrix <- t(atac_matrix)

# Set seed for reproducibility
set.seed(123)
merged_df <- do.call(rbind, fibroblast_SCENT@peak.info.list)
print(dim(merged_df))

peak_gene_counts <- table(merged_df$gene)
peak_gene_counts_df <- as.data.frame(peak_gene_counts)

colnames(peak_gene_counts_df) <- c("gene", "count")
filtered_genes <- peak_gene_counts_df %>% filter(count >= 10) %>% pull(gene)

merged_df <- merged_df %>% filter(gene %in% filtered_genes)
print(dim(merged_df))

unique_genes <- unique(merged_df$gene)
sampled_genes <- sample(unique_genes, 1000)

# Define the process_gene function
process_gene <- function(gene) {
  
  
  test_rna <- rna_matrix[which(row.names(rna_matrix) == gene),]
  test_atac <- atac_matrix[,which(colnames(atac_matrix) %in% merged_df[merged_df$gene == gene,]$peak)]
  
  print(gene)
  test_rna <- as.vector(test_rna)
  dense_test_atac <- as.matrix(test_atac)
  print(dim(dense_test_atac))


  if (ncol(dense_test_atac) < 2) {
    message("Skipping gene ", gene, " due to insufficient peaks.")
    return(NULL)
  }
  
  data_frame <- data.frame(y = test_rna, (dense_test_atac))
  glm_fit <- glm(y ~ ., data = data_frame, family = "poisson")
  
  cvfit <- cv.glmnet(y = test_rna, x = (dense_test_atac), family = "poisson")
  # plot(cvfit)  # Uncomment if you want to plot the CV curve
  
  glm_coefs <- coef(glm_fit)
  lasso_coefs <- coef(cvfit, s = "lambda.min")
  
  glm_predictions <- predict(glm_fit, newdata = data_frame, type = "response")
  lasso_predictions <- predict(cvfit, newx = (dense_test_atac), s = "lambda.min", type = "response")
  
  mse_glm <- mean((test_rna - glm_predictions)^2)
  mse_lasso <- mean((test_rna - lasso_predictions)^2)
  
  num_nonzero_coefs <- sum(lasso_coefs != 0)
  
  result <- list(
    gene = gene,
    glm_coefs = glm_coefs,
    lasso_coefs = lasso_coefs,
    mse_glm = mse_glm,
    mse_lasso = mse_lasso,
    num_nonzero_coefs = num_nonzero_coefs
  )
  
  return(result)
}

# Set up parallel backend
results_list <- list()

# Get the list of unique genes from sampled_df

# Loop through each gene and process
counter <- 0 
for (gene in sampled_genes) {

result <- process_gene(
gene = gene
)

counter = counter + 1

if (counter %% 100 == 0) {
    message("Processed ", counter, " genes...")
  }

results_list[[gene]] <- result
}

saveRDS(results_list,"glm_result.rds")

