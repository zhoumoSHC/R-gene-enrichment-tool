# R-gene-enrichment-tool

这个R脚本用于进行基因集富集分析（GSEA），可以将您的基因列表与预定义的生物学通路进行比较，并进行富集分析。主要功能包括：

1. 基因集富集分析：
   - 比较您的基因列表与已知的基因集（如代表特定生物学功能的基因组），识别统计学上显著富集的基因集。
   - 基因集可在GSEA网站（https://www.gsea-msigdb.org/gsea/downloads.jsp ）下载。
   - 注意：您的输入基因列表中的基因应以gene symbol形式提供。

2. 通路注释：
   - 为您列表中的每个基因标注它们所属的生物学通路。

3. 多重检验校正：
   - 使用多种方法（如Bonferroni, Benjamini-Hochberg, FDR）进行校正。

输出：
脚本生成两个主要输出文件：
1. 富集分析的详细结果
2. 每个基因的通路注释信息

结果文件中的列名与DAVID工具保持一致，具体解释如下：
 - count：富集到此通路的基因数
 - list total：输入的基因列表与基因集的交集数
 - Population total（pop total）：基因集中的总基因数
 - Pop Hits（pop hits）：基因集中属于当前通路的基因数
 - Fold Enrichment：观察到的富集比例与期望富集比例的比值，计算公式为：(count/list total)/(pop hits/pop total)

注：本代码由克劳德老师协助完成。


This R script performs Gene Set Enrichment Analysis (GSEA), comparing your gene list with predefined biological pathways and conducting enrichment analysis. The main features include:

1. Gene Set Enrichment Analysis:
   - Compares your gene list with known gene sets (such as gene sets representing specific biological functions) to identify statistically significant enriched gene sets.
   - Gene sets can be downloaded from the GSEA website (https://www.gsea-msigdb.org/gsea/downloads.jsp).
   - Note: The genes in your input list should be provided as gene symbols.

2. Pathway Annotation:
   - Annotates each gene in your list with the biological pathways it belongs to.

3. Multiple Testing Correction:
   - Applies various methods (such as Bonferroni, Benjamini-Hochberg, FDR) for correction.

Output:
The script generates two main output files:
1. Detailed results of the enrichment analysis
2. Pathway annotation information for each gene

The column names in the result files are consistent with the DAVID tool. Here's what they mean:
 - count: Number of genes enriched in this pathway
 - list total: Number of genes in the intersection of your input gene list and the gene set
 - Population total (pop total): Total number of genes in the gene set
 - Pop Hits (pop hits): Number of genes in the gene set that belong to the current pathway
 - Fold Enrichment: The ratio of observed enrichment proportion to expected enrichment proportion, calculated as: (count/list total)/(pop hits/pop total)

Note: This code was developed with assistance from Claude AI.
