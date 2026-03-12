### DEG enrichment test

# SYMBOL을 ENTREZID로 변환하여 KEGG/Reactome 분석 입력에 맞추는 함수
convert_symbol_to_entrez <- function(symbols, orgdb,keytype) {
  conversion <- bitr(
    symbols,
    fromType = keytype,
    toType = "ENTREZID",
    OrgDb = orgdb
  )
  return(conversion)
}

# enrichment 결과가 없을 때 반환할 빈 데이터프레임을 생성하는 함수
empty_enrichment_df <- function() {
  data.frame(
    ID = character(),
    Description = character(),
    GeneRatio = character(),
    BgRatio = character(),
    pvalue = numeric(),
    p.adjust = numeric(),
    qvalue = numeric(),
    geneID = character(),
    Count = integer(),
    stringsAsFactors = FALSE
  )
}

# enrichment 결과 각 컬럼의 의미를 설명하는 데이터프레임을 생성하는 함수
description_enrichment_df <- function() {
  data.frame(
    Column = c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count"),
    Meaning = c(
        "Term ID",
        "Term name",
        "Ratio of DEGs annotated to this term (Intersected DEGs / QuerySize)",
        "Ratio of background genes annotated to this term (TermSize / Background gene set size)",
        "Raw enrichment p-value with Hypergeometric test",
        "Adjusted p-value (Benjamini-Hochberg method)",
        "Adjusted p-value (Estimated false discovery rate (FDR))",
        "List of DEGs in this term",
        "Number of matched DEGs with term"
    )
  )
}

# GO BP/CC/MF enrichment를 수행하고 ontology별 결과를 리스트로 반환하는 함수
perform_go_enrichment <- function(genes, keytype, orgdb) {
  ontologies <- c("BP", "CC", "MF")
  go_results <- list()

  for (ont in ontologies) {
    ego <- enrichGO(
      gene          = genes,
      OrgDb         = orgdb,
      ont           = ont,
      keyType       = keytype,
      pAdjustMethod = "BH",
      readable      = TRUE
    )

    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      result <- as.data.frame(ego)
      result_sorted <- result[order(result$p.adjust), ]  # p.adjust 기준 오름차순 정렬
    } else {
      result_sorted <- empty_enrichment_df()
    }
    go_results[[ont]] <- result_sorted
  }
  return(go_results)
}

# SYMBOL 입력을 ENTREZID로 변환한 뒤 KEGG enrichment를 수행하는 함수
perform_kegg_enrichment <- function(genes, keytype, orgdb, Org_code) {
  conv <- convert_symbol_to_entrez(genes, orgdb, keytype)
  if (nrow(conv) == 0) {
    return(empty_enrichment_df())
  }

  ekegg <- enrichKEGG(
    gene = unique(conv$ENTREZID),
    organism = Org_code,
    keyType = "ncbi-geneid",
    pAdjustMethod = "BH"
  )

  if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
    ekegg_readable <- setReadable(ekegg, OrgDb = orgdb, keyType = "ENTREZID")
    result <- as.data.frame(ekegg_readable)
    result <- result[order(result$p.adjust), ]  # p.adjust 기준 정렬
  } else {
    result <- empty_enrichment_df()
  }
  return(result)
}

# SYMBOL 입력을 ENTREZID로 변환한 뒤 Reactome enrichment를 수행하는 함수
perform_reactome_enrichment <- function(genes, keytype, orgdb, organism) {
  conv <- convert_symbol_to_entrez(genes, orgdb, keytype)
  if (nrow(conv) == 0) {
    return(empty_enrichment_df())
  }

  erct <- enrichPathway(
    gene          = unique(conv$ENTREZID),
    organism      = organism, #human, mouse
    pAdjustMethod = "BH",
    readable      = TRUE
  )
  if (!is.null(erct) && nrow(as.data.frame(erct)) > 0) {
    result <- as.data.frame(erct)
    result_sorted <- result[order(result$p.adjust), ]  # p.adjust 기준 오름차순 정렬
  } else {
    result_sorted <- empty_enrichment_df()
  }

  return(result_sorted)
}

# 사용자 지정 GMT gene set에 대해 enrichment를 수행하는 함수
perform_custom_geneset_enrichment <- function(genes, keytype, Custom_geneset_gmt) {

  if (length(genes) == 0) return(empty_enrichment_df())

  # Custom Geneset gmt 불러오기기
  Custom_geneset <- read.gmt(Custom_geneset_gmt)

  # enrichment 수행
  enrich_res <- enricher(
    gene = genes,
    TERM2GENE = Custom_geneset,
    pAdjustMethod = "BH"
  )

  # 정리
  if (!is.null(enrich_res) && nrow(as.data.frame(enrich_res)) > 0) {
    enrich_df <- as.data.frame(enrich_res)
    return(enrich_df[order(enrich_df$p.adjust), ])
  } else {
    return(empty_enrichment_df())
  }
}

# DEG enrichment 결과를 bubble plot 형태로 저장하는 함수
plot_DEG_enrichment_analysis <- function(df, Top_n, path, prefix, direction, genes, prefix_suffix = "GO_BP") {
  try({
    # 데이터 전처리
    df$Intersection <- as.numeric(str_split(df$GeneRatio, "/", simplify = TRUE)[, 1])
    df$QuerySize <- as.numeric(str_split(df$GeneRatio, "/", simplify = TRUE)[, 2])
    df$TermSize <- as.numeric(str_split(df$BgRatio, "/", simplify = TRUE)[, 1])
    df$generatio <- df$Intersection / df$QuerySize

    # 상위 GO term 선택
    max <- min(Top_n, nrow(df))
    topn_df <- df[order(df$p.adjust), ][1:max, ]

    # 플롯 크기 조정
    plot_size <- adjust_plot_size(topn_df)

    # PDF 저장
    pdf(file = paste0(path, "2.DEG_enrichment/1.figure/", prefix, "_", prefix_suffix, ".pdf"),
        width = plot_size$width + 5, height = plot_size$height)

    # ggplot 생성
    p <- ggplot(topn_df, aes(x = generatio, 
                              y = reorder(paste0(Description, " (size: ", TermSize, ") "), -log10(p.adjust)), 
                              size = Count, color = -log10(p.adjust))) +
      geom_point(alpha = 0.6) +
      scale_color_gradient(low = "blue", high = "red") +
      theme_minimal() +
      labs(title = paste0(prefix, "_", prefix_suffix), 
           x = "Intersection Size / Term Size", 
           y = "Term", 
           size = "Intersection Size", 
           color = "-Log p-value") +
      annotate("text", x = Inf, y = Inf, label = paste(direction, " regulated DEGs: ", length(genes)), 
               vjust = 2, hjust = 1, size = 3) +
      theme(axis.text = element_text(size = 16),
            plot.title = element_text(hjust = 0.5))

    # 플롯 출력
    print(p)
    dev.off()
  })
}

# DEG 결과에 대해 GO/KEGG/Reactome/Hallmark enrichment를 일괄 수행하는 함수
GO_DEG_analysis = function(path,data,Organism,Comparison_name,Keytype,Top_n = 40, Custom_geneset_gmt = NULL) {
    make_folder(path,"2.DEG_enrichment/")
    make_folder(path,"2.DEG_enrichment/0.excel/")
    make_folder(path,"2.DEG_enrichment/1.figure/")
    if (Organism == "human"){
        OrgDb = org.Hs.eg.db
        Org_code = "hsa"
        Hallmark = "/octopus/yeongjun.kim/Reference/Hallmark/h.all.v2025.1.Hs.symbols.gmt"
    } else if (Organism == "mouse") {
        OrgDb = org.Mm.eg.db
        Org_code = "mmu"
        Hallmark = "/octopus/yeongjun.kim/Reference/Hallmark/mh.all.v2025.1.Mm.symbols.gmt"
    }


    sets = c("UP", "DOWN")
    for (direction in sets) {
        all_results <- list()
        genes <- data$gene_id[data$diffexpressed == direction]
        prefix <- paste(Comparison_name, direction, sep = "__")

        # 각각 enrichment 수행
        go_res <- perform_go_enrichment(genes, Keytype, OrgDb)
        kegg_res <- perform_kegg_enrichment(genes, Keytype, OrgDb, Org_code)
        reactome_res <- perform_reactome_enrichment(genes, Keytype, OrgDb, Organism)
        hallmark_res <- perform_custom_geneset_enrichment(genes, Keytype, Hallmark)
        for (ont in names(go_res)) {
            GO_name = paste0("GO_", ont)
            all_results[[GO_name]] <- go_res[[ont]]
        }
        all_results[["KEGG"]] <- kegg_res
        all_results[["Reactome"]] <- reactome_res
        all_results[["Hallmark"]] <- hallmark_res
        if (!is.null(Custom_geneset_gmt)) {
            custom_res <- perform_custom_geneset_enrichment(genes, Keytype, Custom_geneset_gmt)
            all_results[["Custom"]] <- custom_res
        }
        all_results[["Column_Description"]] = description_enrichment_df()
        write_xlsx(all_results, paste0(path,"2.DEG_enrichment/0.excel/",prefix,".xlsx"))
        # 시각화 예시
        # GO_BP 시각화
        prefix_suffix = "GO_BP"
        plot_DEG_enrichment_analysis(all_results[[prefix_suffix]], Top_n, path, prefix, direction, genes,  prefix_suffix)
    }
}
