# Bulk RNA-seq Pipeline

이 폴더는 Bulk RNA-seq 전처리부터 DEG, enrichment, GSEA, 시각화까지 한 번에 관리하기 위한 파이프라인입니다.

현재 구성은 크게 두 단계입니다.

1. '1_Preprocessing/': 'nf-core/rnaseq' 기반 전처리
2. '2_DE_Analysis/': 'tximport' + 'edgeR' 기반 DEG 및 downstream 분석

---

## 전체 폴더 구조

- '0_env/'
- '1_Preprocessing/'
- '2_DE_Analysis/'
- '3_References/'

---

## 전체 분석 흐름

이 파이프라인은 아래 순서로 사용합니다.

1. FASTQ 메타데이터 준비
2. 'nf-core/rnaseq'로 alignment 및 quantification 수행
3. 'star_salmon' 결과에서 'tx2gene.tsv'와 Salmon output 사용
4. '2_DE_Analysis'에서 DEG 계산
5. Volcano plot, PCA, DEG enrichment, GSEA 생성

---

## 1. Preprocessing 단계

전처리 실행 스크립트:

- '1_Preprocessing/run_rnaseq_preprocessing.sh'

이 스크립트는 'nf-core/rnaseq'를 사용해 RNA-seq preprocessing을 수행합니다.

현재 설정 요약:

- Nextflow 버전: '25.04.0'
- Pipeline: 'nf-core/rnaseq'
- Revision: '3.22.2'
- Profile: 'docker'
- Aligner: 'star_salmon'
- Output directory: '1_Preprocessing/0_result'
- Reference FASTA: '3_References/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa'
- Reference GTF: '3_References/Homo_sapiens.GRCh38.112.gtf'

실행 파일 내용 기준 핵심 명령은 다음과 같습니다.

'''bash
NXF_VER=25.04.0 nextflow run nf-core/rnaseq -r 3.22.2 \
--input example_fastq_metadata.csv \
-profile docker \
--genome null \
-c my.config \
--igenomes_ignore true \
--save_reference \
--outdir './0_result' \
--aligner 'star_salmon' \
--fasta '../3_References/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa' \
--gtf '../3_References/Homo_sapiens.GRCh38.112.gtf' \
-resume
'''

### Preprocessing 입력 파일

- '1_Preprocessing/example_fastq_metadata.csv'
- '1_Preprocessing/my.config'
- '3_References/' 아래 reference FASTA/GTF

### Preprocessing 주요 결과

주요 산출물은 '1_Preprocessing/0_result/star_salmon/' 아래 생성되며, DE 분석에서 특히 다음 파일을 사용합니다.

- 'tx2gene.tsv'
- sample별 Salmon 결과 파일

---

## 2. DE Analysis 단계

DE 분석 실행 예시 파일:

- '2_DE_Analysis/run_pipeline_code_example.R'

이 파일은 'R_scripts' 안의 함수 파일들을 순서대로 불러온 뒤, 설정값을 지정하고 전체 분석을 실행하는 예시입니다.

### 실행 예시 코드 개요

1. 'R_scripts' 내부 함수 로드
2. 'project_settings' 정의
3. 'run_deg_pipeline()' 실행
4. 각 comparison에 대해 'analysis_data()' 실행

예시 코드:

'''r
script_dir = "./R_scripts"

source(file.path(script_dir, "00_libraries.R"))
source(file.path(script_dir, "01_deg_based_enrichment_analysis.R"))
source(file.path(script_dir, "02_gsea.R"))
source(file.path(script_dir, "03_visualization.R"))
source(file.path(script_dir, "04_pipeline.R"))

root_path = "./"

project_settings <- list(
    Project_name = "test_run5/",
    group_col = "Condition",
    comb_set = list(
                    c("ConditionA", "ConditionB")
    ),
    tx2gene_path = "../1_Preprocessing/0_result/star_salmon/tx2gene.tsv",
    metadata_path = "./example_salmon_metadata.csv"
)

pipeline_result <- run_deg_pipeline(
  project_settings = project_settings,
  root_path = root_path
)

result_path <- pipeline_result$result_path
DGE_result <- pipeline_result$DGE_result

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
'''

---

## '2_DE_Analysis/R_scripts' 구성

### '00_libraries.R'
분석에 필요한 라이브러리를 로드합니다.

주요 패키지:

- 'DESeq2'
- 'tximport'
- 'RUVSeq'
- 'sva'
- 'clusterProfiler'
- 'fgsea'
- 'ReactomePA'
- 'org.Hs.eg.db'
- 'org.Mm.eg.db'
- 'ggplot2'
- 'ggrepel'
- 'dplyr'
- 'tibble'
- 'writexl'
- 'openxlsx'
- 'stringr'

### '01_deg_based_enrichment_analysis.R'
DEG 기반 enrichment 분석 함수 모음입니다.

주요 함수:

- 'convert_symbol_to_entrez()'
- 'empty_enrichment_df()'
- 'description_enrichment_df()'
- 'perform_go_enrichment()'
- 'perform_kegg_enrichment()'
- 'perform_reactome_enrichment()'
- 'perform_custom_geneset_enrichment()'
- 'plot_DEG_enrichment_analysis()'
- 'GO_DEG_analysis()'

### '02_gsea.R'
GSEA 분석 및 GSEA figure 생성 함수 모음입니다.

주요 함수:

- 'convert_vector_to_entrez()'
- 'empty_gsea_df()'
- 'description_gsea_df()'
- 'filter_gsea_result()'
- 'adjust_plot_size()'
- 'select_top_GSEA_abs()'
- 'draw_gsea_summary_plot()'
- 'perform_go_gsea()'
- 'perform_kegg_gsea()'
- 'perform_reactome_gsea()'
- 'perform_custom_geneset_gsea()'
- 'run_specific_gsea_analysis()'
- 'GSEA_analysis()'
- 'draw_Specific_GSEA_plot()'

### '03_visualization.R'
PCA 및 volcano plot 생성 함수 모음입니다.

주요 함수:

- 'plot_pca_with_condition_and_batch()'
- 'draw_volcanoplot()'
- 'draw_volcanoplot_annot()'
- 'draw_volcanoplot_genes()'

### '04_pipeline.R'
DE 분석 핵심 파이프라인 함수가 들어 있습니다.

주요 함수:

- 'make_folder()'
- 'TMM_to_DESeq2_format()'
- 'annotate_data()'
- 'analysis_data()'
- 'run_deg_pipeline()'

'run_deg_pipeline()' 내부 주요 단계:

- 입력 준비
- 'tximport' 결과 로드
- TPM 저장
- edgeR용 object 생성
- TMM normalized count 저장
- filtering
- PCA 생성
- contrast별 DEG 계산

---

## DE 분석 입력 파일

### '2_DE_Analysis/example_salmon_metadata.csv'
Salmon 결과 파일 경로와 샘플 메타데이터를 포함하는 CSV입니다.

필수 컬럼:

- 'sample'
- 'file'
- 비교 그룹 컬럼 예: 'Condition'

### 'project_settings'
예시 실행 파일에서 아래 항목을 설정합니다.

- 'Project_name'
- 'group_col'
- 'comb_set'
- 'tx2gene_path'
- 'metadata_path'

### 주요 인자 설명

- 'Sample': 현재 'human' 또는 'mouse'
- 'Keytype': 예: 'SYMBOL', 'ENSEMBL', 'ENTREZID'
- 'specific_geneset': 사용자 정의 GMT 파일 경로
- 'specific_genes': volcano plot에서 강조할 유전자 목록

---

## 주요 출력 결과

DE 분석 결과 폴더('Project_name') 아래 일반적으로 다음 결과가 생성됩니다.

### 루트 결과 파일

- 'TPM_gene.csv'
- 'TPM_gene_log2p1.csv'
- 'TMM_length_normalized_count.csv'
- 'log_TMM_length_normalized_count.csv'
- 'TMM_length_normalized_count_filtered.csv'
- 'PCA_plot.png'
- PCA loading CSV

### 하위 결과 폴더

- '0.DGE/'
- '1.volcano/'
- '1.volcano_annot/'
- '1.volcano_specific_genes/'
- '2.DEG_enrichment/'
- '3.GSEA/'

---

## 빠른 사용 순서

### A. 전처리 실행

'1_Preprocessing'에서 Nextflow 스크립트를 실행합니다.

### B. DE 분석 실행

'2_DE_Analysis'에서 예시 파일을 참고하여 R 스크립트를 실행합니다.

예:

'''r
source("run_pipeline_code_example.R")
'''

또는 파일 내용을 그대로 복사하여 단계별 실행해도 됩니다.

---

## 사용 시 주의사항

1. 'group_col'은 메타데이터에 반드시 존재해야 합니다.
2. 'metadata_path'의 'file' 경로는 실제 Salmon 결과 위치와 정확히 일치해야 합니다.
3. 'tx2gene_path'는 preprocessing 결과의 'star_salmon/tx2gene.tsv'를 사용해야 합니다.
4. 'Sample' 값은 organism DB 선택에 직접 영향을 줍니다.
5. 'Keytype'은 입력 유전자 ID 형식과 일치해야 합니다.
6. 사용자 정의 GMT 파일은 존재하는 경로여야 합니다.
7. 일부 분석은 organism annotation 패키지와 'clusterProfiler', 'ReactomePA'에 의존합니다.

---

## 관련 주요 파일

- '1_Preprocessing/run_rnaseq_preprocessing.sh'
- '2_DE_Analysis/run_pipeline_code_example.R'
- '2_DE_Analysis/DEG_analysis_separate_sample_real copy.ipynb'
- '2_DE_Analysis/R_scripts/'

---

## 문서 목적

이 README는 preprocessing 스크립트, DE 분석 예시 코드, R 함수 스크립트 구조를 한 문서에서 확인할 수 있도록 정리한 전체 파이프라인 안내서입니다.
