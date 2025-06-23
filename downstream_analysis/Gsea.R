source("/nfs/dcmb-lgarmire/yangiwen/workspace/stads/stable/drug_score/DrugScore.R")
source("/nfs/dcmb-lgarmire/yangiwen/workspace/stads/stable/drug_score/DataLoader.R")
source("/nfs/dcmb-lgarmire/yangiwen/workspace/stads/stable/drug_score/L1000Loader.R")

calcReversedGeneEnrichedPathways <- function(
  tumor_samples,
  normal_samples,
  domains,
  drug,
  l1000,
  use_cluster = "stads_domain",
  use_condition = "condition",
  ...
) {
  # Subset drug signatures
  drug_sigs <- rownames(l1000$sig_info[l1000$sig_info$pert_iname == drug, ])
  # Combine signatures using median perturbation value
  drug_perturbation <- matrixStats::rowMedians(as.matrix(l1000$response[, drug_sigs]))
  # Determine drug's perturbation direction in L1000 (up/down regulated)
  l1000_direction <- drug_perturbation > 0
  l1000_gene_ids <- names(l1000_direction)
  l1000_gene_names <- l1000$gene_info[match(l1000_gene_ids, l1000$gene_info$pr_gene_id), "pr_gene_symbol"]
  names(l1000_gene_names) <- l1000_gene_ids
  
  # Calculate pathways using enrichR
  library(enrichR)
  enrichR::setEnrichrSite("Enrichr")
  reversed_genes_list <- list()
  pathway_table <- data.frame()
  for (domain in domains) {
    message(paste0("STADS domain ", domain))
    # Get normal tissues in target domain
    domain_normal_samples <- subsetSeurat(normal_samples, use_cluster, domain)
    # Get tumor tissues in target domain
    domain_tumor_samples <- subsetSeurat(tumor_samples, use_cluster, domain)
    # Find DEGs
    domain_degs <- calcDeg(domain_normal_samples, domain_tumor_samples, use_col = use_condition, ...)
    # Filter degs by pval and logfc/limma score
    domain_degs <- domain_degs[(domain_degs$pval < 0.05) & abs(domain_degs$score) > 0, ]
    # Filter top 1000 differentially expressed genes in tissue
    domain_degs <- domain_degs[order(-abs(domain_degs$score)), ]
    domain_degs <- head(domain_degs, 1000)
    # Determine deg's direction in tissue (up/down regulated)
    deg_direction <- domain_degs$score > 0
    # Reorder drug's gene perturbation direction same as domain degs
    matched_l1000_direction <- l1000_direction[match(rownames(domain_degs), l1000_gene_names)]
    # Xor deg and drug gene direction to find reversed genes
    reversed_genes <- names(which(xor(deg_direction, matched_l1000_direction)))
    reversed_gene_names <- l1000_gene_names[reversed_genes]
    # Store reversed gene to node information in Cytoscape
    reversed_genes_list[[domain]] <- data.frame(
      Node2 = reversed_genes,
      Expression = domain_degs[match(reversed_gene_names, rownames(domain_degs)), "score"]
    )
    # Enrichr
    gse <- enrichR::enrichr(rownames(domain_degs), c(
      "KEGG_2021_Human",
      "MSigDB_Hallmark_2020"
    ))
    gse <- do.call(rbind, lapply(names(gse), function(n) cbind(gse[[n]], db = n)))
    gse <- gse[gse$Adjusted.P.value < 0.05, ]
    gse <- gse[grep("signaling pathway|Apoptosis|Cell cycle", gse$Term, ignore.case = T), ]
    # Store pathway information
    if (length(gse > 0)) {
      pathway_table <- rbind(pathway_table, data.frame(
        description = gse$Term,
        pval = gse$Adjusted.P.value,
        score = gse$Combined.Score,
        cluster = domain,
        member_genes = gse$Genes,
        db = gse$db
      ))
    }
  }
  list(
    l1000_gene_names = l1000_gene_names,
    l1000_direction = l1000_direction,
    drug_perturbation = drug_perturbation,
    reversed_genes_list = reversed_genes_list,
    pathway_table = pathway_table
  )
}

buildCytoscapeInput <- function(drug, reversedGeneEnrichedPathways) {
  l1000_gene_names <- reversedGeneEnrichedPathways$l1000_gene_names
  l1000_gene_ids <- names(l1000_gene_names)
  l1000_direction <- reversedGeneEnrichedPathways$l1000_direction
  drug_perturbation <- reversedGeneEnrichedPathways$drug_perturbation
  reversed_genes_list <- reversedGeneEnrichedPathways$reversed_genes_list
  pathway_table <- reversedGeneEnrichedPathways$pathway_table
  # Combine logfc/limma score for genes in all domains using median
  reversed_genes <- do.call(rbind, reversed_genes_list)
  reversed_genes <- aggregate(Expression ~ Node2, data = reversed_genes, FUN = median)
  # Value of the edge is the amount that drug can reverse
  reversed_genes$Value <- drug_perturbation[match(reversed_genes$Node2, l1000_gene_ids)]
  # Direction of edge is 1 if gene up-regulated in L1000
  reversed_genes$Direction <- ifelse(l1000_direction[match(reversed_genes$Node2, l1000_gene_ids)], 1, -1)
  # Pattern of node is -1 if gene up-regulated in L1000, opposite of edge direction
  reversed_genes$Pattern <- ifelse(l1000_direction[match(reversed_genes$Node2, l1000_gene_ids)], -1, 1)
  # Finalize other information for reversed genes
  reversed_genes$Node1 <- drug
  reversed_genes$Type <- "Gene"
  reversed_genes$Node2 <- l1000_gene_names[match(reversed_genes$Node2, l1000_gene_ids)]
  reversed_genes$Node <- reversed_genes$Node2
  # Split reversed genes information into Cytoscape network node and edges
  network_data <- reversed_genes[, c("Node1", "Node2", "Value", "Direction")]
  network_node <- reversed_genes[, c("Node", "Type", "Expression", "Pattern")]
  # Append pathway-cluster edges
  pathway_cluster_edges <- pathway_table[, c("description", "cluster")]
  colnames(pathway_cluster_edges) <- c("Node1", "Node2")
  pathway_cluster_edges$Value <- 1
  pathway_cluster_edges$Direction <- 2
  network_data <- rbind(pathway_cluster_edges, network_data)
  # Append pathway-gene edges
  pathway_gene_edges <- as.data.frame(tidyr::separate_rows(pathway_table, "member_genes", sep = ";"))
  pathway_gene_edges <- pathway_gene_edges[, c("description", "member_genes")]
  # Filter pathway member genes to include only reversed genes
  pathway_gene_edges <- pathway_gene_edges[pathway_gene_edges$member_genes %in% unique(network_node$Node), ]
  # Assign other attributes
  colnames(pathway_gene_edges) <- c("Node1", "Node2")
  pathway_gene_edges$Value <- 1
  pathway_gene_edges$Direction <- 0
  network_data <- rbind(pathway_gene_edges, network_data)
  # Append pathway nodes
  pathway_nodes <- data.frame(Node = unique(pathway_table$description), Type = "Pathway", Expression = 200, Pattern = 0)
  cluster_nodes <- data.frame(Node = domains, Type = "Cluster", Expression = 200, Pattern = 0)
  drug_nodes <- data.frame(Node = drug, Type = "Drug", Expression = 200, Pattern = 0)
  # Append other nodes
  network_node <- rbind(pathway_nodes, network_node)
  network_node <- rbind(cluster_nodes, network_node)
  network_node <- rbind(drug_nodes, network_node)
  list(
    network_node = network_node,
    network_data = network_data
  )
}

plotPathwayBubble <- function(pathway_table, output_dir = NULL) {
  if (!is.null(output_dir)) {
    png(
      file = file.path(output_dir, "gsea_bubble.jpg"),
      width = 1024,
      height = 512,
      res = 120
    )
  }
  print(
    ggplot2::ggplot(
      pathway_table,
      ggplot2::aes(
        x = cluster,
        y = description,
        size = -log10(pval),
        colour = "red"
      )
    ) +
      ggplot2::geom_point(alpha = 1) +
      ggplot2::guides(color = FALSE) +
      ggplot2::theme(legend.position = "top") +
      ggplot2::scale_size(range = c(1, 5),name = "-log10(pval)") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  )
  if (!is.null(output_dir)) {
    dev.off()
  }
}

runGsea <- function(
  drug,
  output_dir,
  ...
) {
  reversedGeneEnrichedPathways <- calcReversedGeneEnrichedPathways(drug = drug, ...)
  pathway_table <- reversedGeneEnrichedPathways$pathway_table
  cytoscapeInput <- buildCytoscapeInput(drug, reversedGeneEnrichedPathways)
  network_node <- cytoscapeInput$network_node
  network_data <- cytoscapeInput$network_data
  # Save results
  mkdir(output_dir)
  write.csv(network_data, file.path(output_dir, "network_data.csv"))
  write.csv(network_node, file.path(output_dir, "network_node.csv"))
  write.csv(pathway_table, file.path(output_dir, "pathway_table.csv"))
  # Plot pathway bubbles
  plotPathwayBubble(pathway_table, output_dir)
}

library(Seurat)

data_name <- "prostate"
tissue <- "prostate"
lincs_drug_response_phase1_path <- "/nfs/dcmb-lgarmire/yangiwen/workspace/stads/data/l1000/GSE92742"
lincs_drug_response_phase2_path <- "/nfs/dcmb-lgarmire/yangiwen/workspace/stads/data/l1000/GSE70138"
cluster_output_path <- file.path("/nfs/dcmb-lgarmire/yangiwen/workspace/stads/output", data_name, "stads/stable")
drug <- "mitoxantrone"
output_dir <- file.path("/nfs/dcmb-lgarmire/yangiwen/workspace/stads/output", data_name, "gsea", drug)

l1000 <- loadLincsTwoPhasesDrugResponse(lincs_drug_response_phase1_path, lincs_drug_response_phase2_path, tissue)

data_factory <- loadSampleData(cluster_output_path = cluster_output_path, data_name = data_name)
patients <- data_factory$patients
data_list <- data_factory$data
domains <- data_factory$domains

for (patient in patients) {
  message(patient)
  tumor_samples <- data_list[[patient]]$tumor
  normal_samples <- data_list[[patient]]$normal
  runGsea(
    drug,
    tumor_samples = tumor_samples,
    normal_samples = normal_samples,
    domains = domains,
    l1000 = l1000,
    output_dir = file.path(output_dir, patient)
  )
}
