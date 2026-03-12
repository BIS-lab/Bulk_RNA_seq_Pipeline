
script_dir = "./R_scripts"

# 필요한 섹션별 스크립트를 순서대로 로드합니다.
source(file.path(script_dir, "00_libraries.R"))
source(file.path(script_dir, "01_deg_based_enrichment_analysis.R"))
source(file.path(script_dir, "02_gsea.R"))
source(file.path(script_dir, "03_visualization.R"))
source(file.path(script_dir, "04_pipeline.R"))


root_path = "./"

# 전체 설정 리스트화
project_settings <- list(
    Project_name = "test_run5/",
    group_col = "Condition",  # 비교군 기준 컬럼 (예: Condition, Group, Treatment)
    comb_set = list(
                    c("ConditionA", "ConditionB")
    ),
    tx2gene_path = "../1_Preprocessing/0_result/star_salmon/tx2gene.tsv",
    metadata_path = "./example_salmon_metadata.csv"
)

# 전체 파이프라인을 실행하고 주요 결과 객체를 꺼내는 코드입니다.
pipeline_result <- run_deg_pipeline(
  project_settings = project_settings,
  root_path = root_path
)


result_path <- pipeline_result$result_path
DGE_result <- pipeline_result$DGE_result

# 각 comparison에 대해 downstream analysis를 순회 실행하는 코드입니다.
for (comparison_name in names(DGE_result)) {
  analysis_data(
    path = result_path,
    Data = DGE_result[[comparison_name]],
    method = "padj",
    FC = 1,
    pval = 0.05,
    Comparison_name = comparison_name,
    Sample = "human",
    specific_geneset = "../3_References/Test_Specific_geneset.v2025.1.Hs.gmt",
    specific_genes = c("KRAS","ERBB2"),
    Keytype = "SYMBOL",
    norm_method = "edgeR"
  )
}
