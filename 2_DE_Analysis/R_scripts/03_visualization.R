# 발현 행렬과 메타데이터를 이용해 PCA plot과 loading matrix를 저장하는 함수
plot_pca_with_condition_and_batch <- function(
  expr_matrix,
  col_data,
  result_path,
  file_name = "PCA_plot.png",
  plot_title = "PCA: TMM + Length normalized log expression",
  condition_col_1 = "sample",      # 색상으로 사용할 컬럼
  condition_col_2 = "Condition",   # 모양으로 사용할 컬럼
  gene_annotation = NULL,          # 매트릭스에 gene_name이 없을 때 fallback용
  point_alpha = 0.55,               # 포인트 투명도(0~1)
  sample_legend_ncol = 1,
  sample_legend_cex = 0.85,
  legend_width_ratio = 1.2,
  pc_num_1 = 1,
  pc_num_2 = 2,
  show_sample_labels = TRUE,
  sample_label_cex = 0.7,
  sample_label_pos = 3,
  sample_label_offset = 0.5
  ) {

  # 0) expr_matrix에서 annotation 컬럼(gene_id, gene_name) 분리
  annot_col_names <- intersect(c("gene_id", "gene_name"), colnames(expr_matrix))
  if (length(annot_col_names) > 0) {
    gene_annot_from_matrix <- expr_matrix[, annot_col_names, drop = FALSE]
    expr_numeric <- expr_matrix[, setdiff(colnames(expr_matrix), annot_col_names), drop = FALSE]
  } else {
    gene_annot_from_matrix <- NULL
    expr_numeric <- expr_matrix
  }

  # 1) 필수 컬럼 확인
  req_cols <- c(condition_col_1, condition_col_2, "sample")
  miss <- setdiff(req_cols, colnames(col_data))
  if (length(miss) > 0) {
    stop("col_data에 다음 컬럼이 필요합니다: ", paste(miss, collapse = ", "))
  }

  # 2) 샘플 일치성 및 정렬 (수치 컬럼만 기준)
  samples_expr <- colnames(expr_numeric)
  samples_meta <- as.character(col_data$sample)
  if (!all(samples_expr %in% samples_meta)) {
    stop("expr_matrix의 모든 샘플 컬럼명이 col_data$sample에 존재해야 합니다. ",
         "누락: ", paste(setdiff(samples_expr, samples_meta), collapse = ", "))
  }
  col_data <- col_data[match(samples_expr, samples_meta), , drop = FALSE]

  # 3) 색상: condition_col_1
  tol15 <- c(
    "#332288", "#88CCEE", "#44AA99", "#117733", "#999933",
    "#DDCC77", "#CC6677", "#882255", "#AA4499", "#817066",
    "#E69F00", "#D55E00", "#53377A", "#7F180D", "#00538A"
  )

  condition_1_vals <- as.character(col_data[[condition_col_1]])
  condition_1_levels <- unique(condition_1_vals)
  n_cond_1 <- length(condition_1_levels)
  palette_1 <- if (n_cond_1 <= length(tol15)) tol15[seq_len(n_cond_1)] else grDevices::colorRampPalette(tol15)(n_cond_1)
  condition_1_colors <- setNames(palette_1, condition_1_levels)
  condition_1_colors_alpha <- sapply(condition_1_colors, function(col) grDevices::adjustcolor(col, alpha.f = point_alpha))

  # 4) 모양: condition_col_2
  condition_2_vals <- as.character(col_data[[condition_col_2]])
  condition_2_levels <- unique(condition_2_vals)
  filled_pchs <- c(21, 22, 23, 24, 25, 0, 2, 5, 6, 13)
  condition_2_pch_vals <- rep(filled_pchs, length.out = length(condition_2_levels))
  condition_2_pch <- setNames(condition_2_pch_vals, condition_2_levels)

  # 5) PCA (수치 컬럼만 사용)
  pca_res <- prcomp(t(expr_numeric), center = TRUE, scale. = TRUE)
  imp <- summary(pca_res)$importance[2, c(pc_num_1, pc_num_2)]

  # 6) 저장 시작
  png(filename = file.path(result_path, file_name), width = 12, height = 8, units = "in", res = 300)
  on.exit(dev.off(), add = TRUE)

  # 7) 레이아웃: [플롯 | 범례]
  opar <- par(no.readonly = TRUE); on.exit(par(opar), add = TRUE)
  layout(matrix(c(1, 2), nrow = 1), widths = c(3, legend_width_ratio))

  ## 좌측: PCA 플롯
  par(mar = c(5, 5, 4, 1))
  point_pch <- condition_2_pch[condition_2_vals]
  point_bg  <- condition_1_colors_alpha[condition_1_vals]
  point_text_col <- condition_1_colors[condition_1_vals]

  plot(
    x = pca_res$x[, pc_num_1], y = pca_res$x[, pc_num_2],
    pch = point_pch,
    bg  = point_bg,
    col = "black", cex = 1.6,
    xlab = paste0("PC", pc_num_1, " (", round(100 * imp[1], 1), "%)"),
    ylab = paste0("PC", pc_num_2, " (", round(100 * imp[2], 1), "%)"),
    main = plot_title
  )

  # Sample 이름 텍스트
  if (show_sample_labels) {
    text(
      pca_res$x[, pc_num_1], pca_res$x[, pc_num_2],
      labels = as.character(col_data$sample),
      col = point_text_col,
      cex = sample_label_cex,
      pos = sample_label_pos,
      offset = sample_label_offset
    )
  }

  ## 우측: 범례 전용 패널
  par(mar = c(2, 1, 2, 1))
  plot.new()

  legend(
    "topleft",
    title  = paste0(condition_col_1, " (color)"),
    legend = condition_1_levels,
    pch    = 21,
    pt.bg  = condition_1_colors_alpha[condition_1_levels],
    col    = "black",
    pt.cex = 1.5,
    cex    = sample_legend_cex,
    bty    = "n",
    ncol   = sample_legend_ncol,
    x.intersp = 0.6, y.intersp = 0.9
  )

  legend(
    "bottomleft",
    title  = paste0(condition_col_2, " (shape)"),
    legend = condition_2_levels,
    pch    = condition_2_pch[condition_2_levels],
    pt.bg  = "grey85",
    col    = "black",
    pt.cex = 1.5,
    cex    = sample_legend_cex,
    bty    = "n",
    x.intersp = 0.8, y.intersp = 0.9
  )

  # 8) PCA loading matrix 저장
  loading_df <- as.data.frame(pca_res$rotation)

  # gene_id: 매트릭스에 포함된 경우 사용, 없으면 rownames에서 추출
  if (!is.null(gene_annot_from_matrix) && "gene_id" %in% colnames(gene_annot_from_matrix)) {
    loading_df$gene_id <- gene_annot_from_matrix$gene_id
    loading_df <- loading_df %>% relocate(gene_id, .before = 1)

    if ("gene_name" %in% colnames(gene_annot_from_matrix)) {
      loading_df$gene_name <- gene_annot_from_matrix$gene_name
      loading_df <- loading_df %>% relocate(gene_name, .after = gene_id)
    }
  } else {
    # fallback: rownames를 gene_id로 사용
    loading_df <- rownames_to_column(loading_df, var = "gene_id")
  }

  loading_file_name <- paste0(tools::file_path_sans_ext(file_name), "_loading.csv")
  write.csv(loading_df, file = file.path(result_path, loading_file_name), row.names = FALSE)
}

# 기본 volcano plot을 PDF로 저장하는 함수
draw_volcanoplot = function(path,data,method,FC,pval,Comparison_name){
    make_folder(path,"1.volcano/")
    mycolors <- c("blue", "red", "black")
    names(mycolors) <- c("DOWN", "UP", "NO")
    Log10Pvalue = -log10(data[[method]])

    pdf(file= paste0(path,"1.volcano/",Comparison_name,".pdf"), width=7, height=5)

    plot = ggplot(data=data, aes(x=log2FoldChange, y=Log10Pvalue, col=diffexpressed)) + 
                geom_point() + 
                theme_minimal() +
                #geom_text_repel() +
                scale_color_manual(values=mycolors) +
                geom_vline(xintercept=c(-FC, FC), col="red") +
                geom_hline(yintercept=-log10(pval), col="red") + annotate("rect", xmin=Inf, xmax=Inf, ymin=Inf, ymax=Inf, fill="white", alpha=0.5) +
                annotate("text", x=Inf, y=Inf, label=paste("UP:", sum(data$diffexpressed == "UP")), vjust=3, hjust=1, size=3) + # sum이 원래 count 였음 24/11/30
                annotate("text", x=Inf, y=Inf, label=paste("DOWN:", sum(data$diffexpressed == "DOWN")), vjust=5, hjust=1, size=3) +
                ggtitle(Comparison_name) + 
                theme(plot.title = element_text(hjust = 0.5)) + 
                ylab("-Log10 adj P-value") + 
                xlab("Log2 Fold Change")
    print(plot)
    dev.off()
}

# DEG 라벨을 포함한 volcano plot을 PDF로 저장하는 함수
draw_volcanoplot_annot = function(path,data,method,FC,pval,Comparison_name){
    make_folder(path,"1.volcano_annot/")
    mycolors <- c("blue", "red", "black")
    names(mycolors) <- c("DOWN", "UP", "NO")
    Log10Pvalue = -log10(data[[method]])
    pdf(file= paste0(path,"1.volcano_annot/",Comparison_name,".pdf"), width=7, height=5)

    plot = ggplot(data=data, aes(x=log2FoldChange, y=Log10Pvalue, col=diffexpressed, label=delabel)) + ## , label=delabel 이거 붙으면 text 들어가는거야.
                geom_point() + 
                theme_minimal() +
                geom_text_repel() +
                scale_color_manual(values=mycolors) +
                geom_vline(xintercept=c(-FC, FC), col="red") +
                geom_hline(yintercept=-log10(pval), col="red") + annotate("rect", xmin=Inf, xmax=Inf, ymin=Inf, ymax=Inf, fill="white", alpha=0.5) +
                annotate("text", x=Inf, y=Inf, label=paste("UP:", sum(data$diffexpressed == "UP")), vjust=3, hjust=1, size=3) + # sum이 원래 count 였음 24/11/30
                annotate("text", x=Inf, y=Inf, label=paste("DOWN:", sum(data$diffexpressed == "DOWN")), vjust=5, hjust=1, size=3) +
                ggtitle(paste0(Comparison_name)) + 
                theme(plot.title = element_text(hjust = 0.5)) + 
                ylab("-Log10 adj P-value") + 
                xlab("Log2 Fold Change")
    print(plot)
    dev.off()
}

# 특정 유전자만 라벨링한 volcano plot을 PDF로 저장하는 함수
draw_volcanoplot_genes = function(path, data, method, FC, pval, Comparison_name, specific_genes = NULL){
    make_folder(path,"1.volcano_specific_genes/")
    mycolors <- c("blue", "red", "black")
    names(mycolors) <- c("DOWN", "UP", "NO")
    Log10Pvalue = -log10(data[[method]])
    data$specific_gene <- NA
    if (!is.null(specific_genes)) {
        data$specific_gene[data$gene_id %in% specific_genes] <- data$gene_id[data$gene_id %in% specific_genes]
    }
    Saving_file_name = paste0(Comparison_name)
    pdf(file= paste0(path,"1.volcano_specific_genes/",Saving_file_name,".pdf"), width=7, height=5)

    plot = ggplot(data=data, aes(x=log2FoldChange, y=Log10Pvalue, col=diffexpressed, label=specific_gene)) +
                geom_point() + 
                theme_minimal() +
                geom_text_repel() +
                scale_color_manual(values=mycolors) +
                geom_vline(xintercept=c(-FC, FC), col="red") +
                geom_hline(yintercept=-log10(pval), col="red") + 
                annotate("rect", xmin=Inf, xmax=Inf, ymin=Inf, ymax=Inf, fill="white", alpha=0.5) +
                annotate("text", x=Inf, y=Inf, label=paste("UP:", sum(data$diffexpressed == "UP")), vjust=3, hjust=1, size=3) +
                annotate("text", x=Inf, y=Inf, label=paste("DOWN:", sum(data$diffexpressed == "DOWN")), vjust=5, hjust=1, size=3) +
                ggtitle(paste0(Comparison_name)) + 
                theme(plot.title = element_text(hjust = 0.5)) + 
                ylab("-Log10 P-value") + 
                xlab("Log2 Fold Change")

    print(plot)
    dev.off()
}
