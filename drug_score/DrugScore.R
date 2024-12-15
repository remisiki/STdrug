source("/nfs/dcmb-lgarmire/yangiwen/workspace/common/Utils.R")

calcDeg <- function(
  normal_samples,
  tumor_samples,
  use_col = "condition",
  min_cells = 3,
  min_domain_size = 3,
  test_use = "limma"
) {
  degs <- NULL
  if ((ncol(normal_samples) > min_domain_size) &&
      (ncol(tumor_samples) > min_domain_size)) {
    # Merge normal and tumor samples, subset by common genes
    domain_samples <- mergeSeuratList(c(normal_samples, tumor_samples))
    # DE analysis
    if (test_use == "limma") {
      X <- Seurat::GetAssayData(domain_samples)
      # Filter genes expressed with min_cells threshold
      X <- X[Matrix::rowSums(X > 0) >= min_cells, ]
      # Find differentially expressed genes with limma
      model_frame <- data.frame(
        batches = colnames(domain_samples),
        condition = domain_samples@meta.data[[use_col]]
      )
      model_matrix <- model.matrix(~ 0 + condition, data = model_frame)
      fit <- limma::lmFit(X, model_matrix)
      contrasts <- limma::makeContrasts(conditionT - conditionN, levels = colnames(coef(fit)))
      fit <- limma::contrasts.fit(fit, contrasts = contrasts)
      fit <- limma::eBayes(fit)
      de_out <- limma::topTable(fit, sort.by = "P", number = nrow(fit))
      degs <- data.frame(
        row.names = rownames(de_out),
        score = de_out$t,
        pval = de_out$P.Value,
        adj_pval = de_out$adj.P.Val
      )
    } else {
      # FIXME This should be run before FindMarkers on SCT assay, however a
      # bug prevent this from running. Find markers on Spatial assay instead.
      # See <https://github.com/satijalab/seurat/issues/9112>
      # domain_samples <- Seurat::PrepSCTFindMarkers(domain_samples)
      markers <- Seurat::FindMarkers(
        domain_samples,
        ident.1 = "T",
        ident.2 = "N",
        test.use = test_use,
        group.by = use_col,
        assay = "Spatial"
      )
      degs <- data.frame(
        row.names = rownames(markers),
        score = markers$avg_log2FC,
        pval = markers$p_val,
        adj_pval = markers$p_val_adj
      )
    }
  } else {
    message("Domain size smaller than threshold, ignored.")
  }
  return(degs)
}

rankDrugs <- function(drug_scores, domains = c("all"), fdr_threshold = NULL, score_threshold = NULL) {
  drug_scores <- drug_scores[drug_scores$domain %in% domains, ]
  drug_scores <- drug_scores[order(drug_scores$score, decreasing = T), ]
  n <- nrow(drug_scores)
  drug_scores$stand <- 1 - seq_len(n) / n
  if (!is.null(fdr_threshold)) {
    drug_scores <- drug_scores[drug_scores$fdr < fdr_threshold, ]
  }
  if (!is.null(score_threshold)) {
    drug_scores <- drug_scores[drug_scores$score > 0, ]
  }
  rownames(drug_scores) <- NULL
  return(drug_scores)
}

calcCciRatio <- function(
  object,
  batch_key = "batch",
  domain_key = "stads_domain"
) {
  data_input <- Seurat::GetAssayData(object)
  meta <- setNames(object@meta.data[c(batch_key, domain_key)], c("batches", "labels"))
  meta[] <- lapply(meta, as.factor)
  batches <- levels(meta$batches)
  spatial_locs <- Reduce(rbind, lapply(batches, function(batch) {
    return(object@images[[batch]]@coordinates[, c("imagerow", "imagecol")])
  }))
  spot_size <- 65
  spatial_factors <- Reduce(rbind, lapply(batches, function(batch) {
    spot_diameter_fullres <- object@images[[batch]]@scale.factors$spot
    return(data.frame(ratio = spot_size / spot_diameter_fullres, tol = spot_size / 2))
  }))
  rownames(spatial_factors) <- batches
  
  # Create cellchat object
  cellchat <- CellChat::createCellChat(
    object = data_input,
    meta = meta,
    group.by = "labels",
    datatype = "spatial",
    coordinates = spatial_locs,
    spatial.factors = spatial_factors
  )
  
  # Set the ligand-receptor interaction database
  cellchat@DB <- CellChat::CellChatDB.human
  
  # Preprocessing the expression data for cell-cell communication analysis
  # This step is necessary even if using the whole database
  cellchat <- CellChat::subsetData(cellchat)
  # future::plan("multisession", workers = 4)
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  
  # Compute the communication probability and infer cellular communication network
  cellchat <- CellChat::computeCommunProb(
    cellchat,
    type = "truncatedMean",
    trim = 0.1,
    distance.use = FALSE,
    interaction.range = 250,
    scale.distance = NULL,
    contact.dependent = TRUE,
    contact.range = 100
  )
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = 10)
  
  # Infer the cell-cell communication at a signaling pathway level
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  
  # Calculate the aggregated cell-cell communication network
  cellchat <- CellChat::aggregateNet(cellchat)
  
  return(cellchat@net$weight)
}

calcCegPval <- function(degs, l1000) {
  # Match common gene
  common_genes <- intersect(rownames(degs), unique(l1000$gene_info$pr_gene_symbol))
  common_gene_ids <- as.character(l1000$gene_info[match(common_genes, l1000$gene_info$pr_gene_symbol), "pr_gene_id"])
  response <- as.data.frame(l1000$response)[common_gene_ids, ]
  degs <- degs[common_genes, ]
  # Rank L1000 response and DEGs
  response <- apply(response, 2, function(x) {
    rank(-x, ties.method = "first")
  })
  degs$rank <- rank(degs$score, ties.method = "first")
  deg_rank <- as.matrix(replicate(ncol(response), degs$rank))
  # Calculate pval of type 1 CEG
  order_stat <- pmin(response, deg_rank)
  p_min <- 1 - pbeta((order_stat - 1) / nrow(order_stat), 1, 2, lower.tail = T)
  # Calculate pval of type 2 CEG
  order_stat <- pmax(response, deg_rank)
  p_max <- pbeta(order_stat / nrow(order_stat), 2, 1, lower.tail = T)
  # Calculate sum of z scores
  p_score <- pmin(p_min, p_max)
  return(p_score)
}

calcGeneImportance <- function(ceg_pval, l1000, drug_annotation) {
  train <- t(ceg_pval)
  sig_info <- l1000$sig_info
  gene_info <- l1000$gene_info
  drug_names <- sig_info[match(rownames(train), sig_info$sig_id), "pert_iname"]
  labels <- as.numeric(drug_annotation[match(drug_names, drug_annotation$drug_name), "response_gpt"] == "True")
  xgb <- xgboost::xgboost(
    train,
    labels,
    nrounds = 500,
    objective = "binary:logistic",
    verbose = 1
  )
  gene_importance <- data.frame(gene_id = colnames(train), score = 0)
  xgb_importance <- xgboost::xgb.importance(model = xgb)[, c("Feature", "Gain")]
  gene_importance[match(xgb_importance$Feature, gene_importance$gene_id), "score"] <- xgb_importance$Gain
  gene_importance$gene_symbol <- gene_info[match(gene_importance$gene_id, gene_info$pr_gene_id), "pr_gene_symbol"]
  gene_importance <- gene_importance[match(rownames(ceg_pval), gene_importance$gene_id), "score"]
  return(gene_importance)
}

calcDrugFdr <- function(drugs, ceg_pval, gene_importance, sig_info) {
  ceg_mask <- ceg_pval <= 0.05
  ceg_mask <- gene_importance * ceg_mask
  z_score <- qnorm(ceg_pval, lower.tail = F)
  z_score <- z_score * ceg_mask
  os_score <- colSums(z_score)
  # Calculate pval using K-S test
  df <- data.frame(matrix(ncol = 2, nrow = length(drugs)), row.names = drugs)
  colnames(df) <- c("pval", "fdr")
  for (drug in drugs) {
    sig_ids <- sig_info[sig_info$pert_iname == drug, "sig_id"]
    drug_os <- os_score[sig_ids]
    rest_os <- os_score[!names(os_score) %in% sig_ids]
    pval <- (ks.test(drug_os, rest_os, alternative = "less"))$p.value
    df[drug, "pval"] <- pval
  }
  df$fdr <- p.adjust(df$pval, method = "fdr")
  return(df)
}

calcDrugEfficacyScore <- function(gdsc, tissue, drugs) {
  scores <- gdsc[gdsc$tissue == tissue, ]
  scores <- 1 - scores[match(drugs, scores$drug), "normalized_ic50"]
  replace(scores, is.na(scores), 0)
}

calcDrugSafetyScore <- function(sider, drugs) {
  scores <- 1 - sider[match(drugs, sider$drug), "normalized_freq"]
  replace(scores, is.na(scores), 0)
}

calcDrugScore <- function(
  tissue,
  tumor_samples,
  normal_samples,
  domains,
  l1000,
  drug_annotation,
  gdsc,
  sider,
  domain_key = "stads_domain",
  checkpoint_path = NULL
) {
  # Calculate cluster proportions in diseased tissue
  domain_props <- prop.table(table(tumor_samples@meta.data[[domain_key]]))
  
  # Calculate CCI score
  cci_ratio_tumor <- calcCciRatio(tumor_samples)
  cci_ratio_normal <- calcCciRatio(normal_samples)
  cci_score <- exp(as.matrix(cci_ratio_tumor - cci_ratio_normal))
  cci_score <- rowSums(cci_score) - diag(cci_score)
  if (!is.null(checkpoint_path)) {
    write.csv(cci_ratio_tumor, file.path(checkpoint_path, paste0("cci_ratio_", patient, "T.csv")))
    write.csv(cci_ratio_normal, file.path(checkpoint_path, paste0("cci_ratio_", patient, "N.csv")))
  }
  
  # FDR
  drugs <- unique(l1000$sig_info[match(colnames(l1000$response), l1000$sig_info$sig_id), "pert_iname"])
  fdr_map <- list()
  for (domain in domains) {
    message(paste0("STADS domain ", domain))
    # Find DEGs
    message("Find DEG")
    # Get normal tissues in target domain
    domain_normal_samples <- subsetSeurat(normal_samples, domain_key, domain)
    # Get tumor tissues in target domain
    domain_tumor_samples <- subsetSeurat(tumor_samples, domain_key, domain)
    domain_degs <- calcDeg(domain_normal_samples, domain_tumor_samples)
    # Find CEGs
    message("Find CEG")
    ceg_pval <- calcCegPval(domain_degs, l1000)
    # Calculate gene importance
    message("Calculate gene importance")
    gene_importance <- calcGeneImportance(ceg_pval, l1000, drug_annotation)
    # Calculate OS score and FDR
    message("Calculate drug fdr")
    fdr_map[[domain]] <- calcDrugFdr(drugs, ceg_pval, gene_importance, l1000$sig_info)
  }
  
  # Calculate patient therapeutic score
  drug_scores <- data.frame()
  for (drug in drugs) {
    drug_score <- 0
    
    for (domain in domains) {
      fdr <- fdr_map[[domain]][drug, "fdr"]
      domain_score <- domain_props[[domain]] * (-log(fdr + 1e-10))
      if (!is.null(cci_score)) {
        domain_score <- domain_score * cci_score[[domain]]
      }
      drug_score <- drug_score + domain_score
      drug_scores <- rbind(drug_scores,
                           data.frame(
                             drug = drug,
                             domain = domain,
                             score = domain_score
                           ))
    }
    
    drug_scores <- rbind(drug_scores,
                         data.frame(
                           drug = drug,
                           domain = "all",
                           score = drug_score
                         ))
  }
  
  # Save drug results with each domain score
  if (!is.null(checkpoint_path)) {
    write.csv(drug_scores, file.path(checkpoint_path, paste0("drug_scores_", patient, ".csv")))
  }
  
  # Combine domain drugs, calculate ic50 and side effects
  drug_scores <- rankDrugs(drug_scores)
  drug_scores$score <- drug_scores$score + calcDrugEfficacyScore(gdsc, tissue, drug_scores$drug) + calcDrugSafetyScore(sider, drug_scores$drug)
  drug_scores <- rankDrugs(drug_scores)
  
  # Save final drug results
  write.csv(drug_scores, file.path(output_dir, paste0("drugs_", patient, ".csv")))
  
  drug_scores
}
