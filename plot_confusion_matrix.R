#!/usr/bin/env Rscript

# 加载所需的库
library(argparse)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(gtable)
library(grid)

# 加载patchwork包用于布局对齐
library(patchwork)

# 确保get_legend函数可用
if(!exists("get_legend")) {
    get_legend <- function(myggplot){
        tmp <- ggplot_gtable(ggplot_build(myggplot))
        leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
        if(length(leg) > 0) {
            legend <- tmp$grobs[[leg]]
            return(legend)
        } else {
            return(NULL)
        }
    }
}

# 尝试加载aricode包
if (!require(aricode, quietly = TRUE)) {
    cat("警告: aricode包未安装，将使用备用函数\n")
    # 定义备用的sortPairs函数
    sortPairs <- function(true, pred) {
        # 简化的sortPairs实现，基于列联表计算
        cont_table <- table(true, pred)
        
        # 获取行和列的名称用于配对
        true_levels <- rownames(cont_table)
        pred_levels <- colnames(cont_table)
        
        # 为每个真实类别和预测类别创建配对
        nij <- as.vector(cont_table)
        ni. <- rowSums(cont_table)
        n.j <- colSums(cont_table)
        
        # 创建配对索引
        pair_grid <- expand.grid(1:length(true_levels), 1:length(pred_levels))
        pair_c1 <- pair_grid[,1] - 1  # 转为0索引
        pair_c2 <- pair_grid[,2] - 1  # 转为0索引
        
        list(
            pair_c1 = pair_c1,
            pair_c2 = pair_c2,
            nij = nij,
            ni. = ni.,
            n.j = n.j
        )
    }
} else {
    cat("aricode包已加载\n")
}

# =============================================================================
# 自定义辅助函数占位符
# =============================================================================
# 用户需要在此处粘贴 adjusted_wallance_indices 函数的定义
# 这是一个关键依赖，必须在使用 cross_table_plot 之前定义
#
# 计算ARI和每个类别/簇的同质性/完整性指标
adjusted_wallance_indices <- function(true=NULL, pred=NULL, contigency_res=NULL){  
    # 创建列联表
    cont_table <- table(true, pred)
    
    # 总样本数
    n <- sum(cont_table)
    
    # 行列和
    row_sums <- rowSums(cont_table)
    col_sums <- colSums(cont_table)
    
    # 计算ARI
    # ARI的标准公式
    term1 <- sum(choose(cont_table, 2), na.rm=TRUE)
    term2 <- sum(choose(row_sums, 2), na.rm=TRUE) * sum(choose(col_sums, 2), na.rm=TRUE) / choose(n, 2)
    term3 <- 0.5 * (sum(choose(row_sums, 2), na.rm=TRUE) + sum(choose(col_sums, 2), na.rm=TRUE))
    
    if (term3 - term2 == 0) {
        ari <- 0
    } else {
        ari <- (term1 - term2) / (term3 - term2)
    }
    
    # 计算全局同质性 (Homogeneity)
    h_total <- -sum((col_sums/n) * log2(col_sums/n + 1e-10))
    h_conditional <- 0
    for(i in 1:nrow(cont_table)) {
        if(row_sums[i] > 0) {
            h_conditional <- h_conditional + (row_sums[i]/n) * 
                           (-sum((cont_table[i,]/row_sums[i]) * log2(cont_table[i,]/row_sums[i] + 1e-10)))
        }
    }
    global_homogeneity <- ifelse(h_total == 0, 1, 1 - h_conditional/h_total)
    
    # 计算全局完整性 (Completeness)
    h_classes <- -sum((row_sums/n) * log2(row_sums/n + 1e-10))
    h_clusters_given_classes <- 0
    for(j in 1:ncol(cont_table)) {
        if(col_sums[j] > 0) {
            h_clusters_given_classes <- h_clusters_given_classes + (col_sums[j]/n) * 
                                      (-sum((cont_table[,j]/col_sums[j]) * log2(cont_table[,j]/col_sums[j] + 1e-10)))
        }
    }
    global_completeness <- ifelse(h_classes == 0, 1, 1 - h_clusters_given_classes/h_classes)
    
    # 为每个预测簇计算单独的同质性分数
    avj <- numeric(ncol(cont_table))
    for(j in 1:ncol(cont_table)) {
        if(col_sums[j] > 0) {
            # 该簇内各类别的纯度
            cluster_entropy <- -sum((cont_table[,j]/col_sums[j]) * log2(cont_table[,j]/col_sums[j] + 1e-10))
            # 同质性 = 1 - 归一化熵
            max_entropy <- log2(sum(cont_table[,j] > 0))  # 最大可能熵
            if(max_entropy > 0) {
                avj[j] <- 1 - cluster_entropy/max_entropy
            } else {
                avj[j] <- 1  # 完全纯净
            }
        } else {
            avj[j] <- 0
        }
    }
    
    # 为每个真实类别计算单独的完整性分数
    awi <- numeric(nrow(cont_table))
    for(i in 1:nrow(cont_table)) {
        if(row_sums[i] > 0) {
            # 该类别被分散到各簇的程度
            class_entropy <- -sum((cont_table[i,]/row_sums[i]) * log2(cont_table[i,]/row_sums[i] + 1e-10))
            # 完整性 = 1 - 归一化熵
            max_entropy <- log2(sum(cont_table[i,] > 0))  # 最大可能熵
            if(max_entropy > 0) {
                awi[i] <- 1 - class_entropy/max_entropy
            } else {
                awi[i] <- 1  # 完全集中
            }
        } else {
            awi[i] <- 0
        }
    }
    
    return(list("AW"=global_homogeneity, "AV"=global_completeness, "ARI"=ari, 
               "AW2"=global_homogeneity, "AV2"=global_completeness, "ARI2"=ari,
               "Awi"=awi, "Avj"=avj))
}

# =============================================================================
# 主要绘图函数
# =============================================================================
cross_table_plot <- function(ground_truth, clusterings, a=1.3, b=5.7, c=2, m=0, n=0.2){
    # Define color palettes used inside the function
    my_color <- c("white", "#feebe2", "#fcc5c0", "#fa9fb5", "#f768a1", "#dd3497", "#ae017e", "#7a0177", "#49006a")
    my_color_4 <- c("white", "#ffffd4", "#fed98e", "#fe9929", "#d95f0e", "#993404")
    
    x <- unique(ground_truth)
    y <- as.factor(unique(clusterings))
    data <- expand.grid(X=x, Y=y)
    cross_count <- table(ground_truth, clusterings) 

    # cell count in the cross_count table
    data$Z1 <- as.numeric(apply(data, 1, function(x){cross_count[as.character(x[["X"]]),as.character(x[["Y"]])]}))
    # log transform Z1
    data$Z2 <- apply(data, 1, function(x){log(cross_count[as.character(x[["X"]]),as.character(x[["Y"]])]+1)})
    # row normalize of Z1
    data <- data %>% group_by(X) %>% mutate(Z3 = 100*Z1/sum(Z1))

    top_cluster <- data %>% group_by(X) %>% top_n(1, Z1)

    unselected_Y <- setdiff(unique(data$Y), unique(top_cluster$Y))
    top_cluster <- top_cluster[order(top_cluster$X),]
    new_levels <- as.numeric(c(as.character(unique(top_cluster$Y)),unselected_Y))
    data$Y <- factor(data$Y, levels=new_levels)

    res <- adjusted_wallance_indices(ground_truth, clusterings)

    # 获取列联表用于数据对齐
    cont_table <- table(ground_truth, clusterings)
    
    # 创建完整性(Completeness)数据框 - 对应真实类别
    # 确保与热图行顺序完全一致
    true_class_names <- rownames(cont_table)  # 从列联表获取真实顺序
    df_awi <- data.frame(
        cell_type = true_class_names,
        awi = res$Awi,  # 完整性分数作为条形图长度
        stringsAsFactors = FALSE
    )
    # 计算每个真实类别的细胞比例（用于颜色映射）
    true_counts <- table(ground_truth)
    df_awi$frac <- as.numeric(true_counts[df_awi$cell_type]) / length(ground_truth)
    # 确保顺序与热图y轴一致（从上到下，不需要反转）
    df_awi$cell_type <- factor(df_awi$cell_type, levels = true_class_names)

    # 创建同质性(Homogeneity)数据框 - 对应预测簇
    # 确保与热图列顺序完全一致
    cluster_names <- colnames(cont_table)  # 从列联表获取真实顺序
    df_avj <- data.frame(
        cell_type = cluster_names,
        avj = res$Avj,  # 同质性分数作为条形图长度
        stringsAsFactors = FALSE
    )
    # 计算每个预测簇的细胞比例（用于颜色映射）
    cluster_counts <- table(clusterings)
    df_avj$frac <- as.numeric(cluster_counts[df_avj$cell_type]) / length(clusterings)
    
    # 确保预测簇按照热图的顺序排列（使用new_levels保持一致性）
    df_avj$cell_type <- factor(df_avj$cell_type, levels = new_levels)

    main <- ggplot(data, aes(Y, X, fill= Z3)) + 
    geom_tile(colour="gray80", linewidth=0.3) + 
    scale_fill_gradientn(colours = my_color, breaks=seq(0,100,10), 
                        guide = guide_colourbar(title.position = "top", title = "Cells per true class %")) +
    scale_x_discrete(position = "bottom") +
    scale_y_discrete(position = "left") +
    labs(x="ATAC cluster", y="True class") + 
    theme_minimal() +
    theme(legend.direction = "horizontal", 
          legend.position = "bottom", 
          legend.key.width= unit(1.5, 'cm'),
          legend.text=element_text(size=12, face="bold", color="black"),
          legend.title=element_text(size=14, face="bold", color="black"),
          panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          panel.border = element_rect(color = "gray30", fill = NA, linewidth = 0.5),
          axis.title.x = element_text(size = 14, margin = margin(10,0,0,0), face="bold", color="black"), 
          axis.text.x = element_text(size = 12, face="bold", angle = 45, hjust = 1, color="black"),
          axis.title.y = element_text(size = 14, margin = margin(0,10,0,0), face="bold", color="black"), 
          axis.text.y = element_text(size = 12, face="bold", color="black"),
          axis.ticks = element_line(color="gray30"),
          plot.margin = unit(c(0, 2, 2, 2), "mm"))

    # 计算颜色范围 - 基于实际的细胞比例
    frac_range <- range(c(df_avj$frac, df_awi$frac))
    
    # 同质性条形图(顶部) - 每个柱子代表一个预测簇
    # 长度 = 同质性分数(avj), 颜色 = 簇大小比例(frac)
    bp.x <- ggplot(data = df_avj, aes(x = cell_type, y = avj, fill = frac)) + 
    geom_bar(stat = "identity", colour = "white", linewidth = 0.2) + 
    scale_x_discrete(labels = levels(df_avj$cell_type)) +  # 确保x轴标签显示
    scale_fill_gradientn(colours = my_color_4, 
                        limits = frac_range,
                        guide = "none") +  # 移除图例，避免重复
    ylim(0, 1) +  # 确保y轴范围一致
    theme_minimal() +
    theme(
          axis.title.x = element_blank(), 
          axis.text.x = element_blank(),  # 隐藏x轴标签，因为主热图有
          axis.ticks.x = element_blank(), 
          axis.text.y = element_text(size = 12, face = "bold", color = "black"), 
          axis.title.y = element_text(size = 14, margin = margin(0,5,0,0), face = "bold", color = "black"),
          legend.position = "none",
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          panel.border = element_rect(color = "gray30", fill = NA, linewidth = 0.5),
          plot.margin = unit(c(2, 2, 0, 2), "mm")) + 
    labs(y = "Homogeneity")

    # 完整性条形图(右侧) - 每个柱子代表一个真实类别
    # 长度 = 完整性分数(awi), 颜色 = 类别大小比例(frac)
    bp.y <- ggplot(data = df_awi, aes(x = cell_type, y = awi, fill = frac)) +
    geom_bar(stat = "identity", colour = "white", linewidth = 0.2) + 
    coord_flip() +
    scale_y_continuous(limits = c(0, 1)) +  # 确保x轴范围一致
    scale_fill_gradientn(colours = my_color_4, 
                        limits = frac_range,
                        guide = "none") +  # 隐藏此图的图例，使用共享图例
    theme_minimal() +
    theme(axis.title.x = element_text(size = 14, margin = margin(5,0,0,0), face = "bold", color = "black"), 
            axis.text.x = element_text(size = 12, face = "bold", color = "black"),
            axis.text.y = element_blank(),  # 隐藏y轴标签，因为主热图已经有
            axis.title.y = element_blank(), 
            axis.ticks.y = element_blank(), 
            axis.ticks.x = element_line(color = "gray30"),
            legend.position = "none",
            panel.grid = element_blank(),
            panel.background = element_rect(fill = "white", color = NA),
            plot.background = element_rect(fill = "white", color = NA),
            panel.border = element_rect(color = "gray30", fill = NA, linewidth = 0.5),
            plot.margin = unit(c(2, 2, 2, 0), "mm")) + 
    labs(x = "Completeness")

    df_hm <- data.frame(cols = numeric(0), value = numeric(0))

    # 创建ARI显示面板
    ari_text <- paste0("ARI = ", round(res$ARI, 3), "\n", "ARI2 = ", round(res$ARI2, 3))
    
    gg_empty <- ggplot() +
    geom_blank() +
    xlim(0, 1) + ylim(0, 1) +
    annotate("text", x = 0.5, y = 0.5, label = ari_text, 
             size = 8, fontface = "bold", color = "black", 
             hjust = 0.5, vjust = 0.5) +
    theme_void() +
    theme(panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          panel.border = element_rect(color = "gray30", fill = NA, linewidth = 0.5),
          plot.margin = unit(c(2, 2, 2, 2), "mm"))

    # 移除顶部条形图的图例，保留主热图的图例在底部
    bp.x_no_legend <- bp.x + theme(legend.position = "none")
    
    # 使用patchwork进行精确对齐
    # 使用patchwork实现完美对齐 - 2x2布局
    # 顶部行: 同质性条形图 + ARI文本
    # 底部行: 主热图 + 完整性条形图
    top_row <- bp.x_no_legend + gg_empty
    bottom_row <- main + bp.y
    
    plot <- top_row / bottom_row + 
            plot_layout(
                heights = c(0.8, 3),  # 顶部较短，底部较高
                widths = c(4, 1)      # 左列更宽，右列更窄
            )
    return(plot)
}

# =============================================================================
# 命令行参数解析
# =============================================================================
parser <- ArgumentParser(description = "生成混淆矩阵热图，比较真实标签与聚类结果")
parser$add_argument("--clusters", required = TRUE, 
                    help = "输入聚类结果TSV文件路径")
parser$add_argument("--truth", required = TRUE, 
                    help = "输入真实标签TSV文件路径")
parser$add_argument("--output", required = TRUE, 
                    help = "输出PNG图片文件路径")

args <- parser$parse_args()

# =============================================================================
# 主执行逻辑
# =============================================================================
cat("正在读取输入文件...\n")

# 读取输入TSV文件
clusters_df <- read.csv(args$clusters, sep = "\t", stringsAsFactors = FALSE)
truth_df <- read.csv(args$truth, sep = "\t", stringsAsFactors = FALSE)

cat("聚类文件包含", nrow(clusters_df), "行，", ncol(clusters_df), "列\n")
cat("真实标签文件包含", nrow(truth_df), "行，", ncol(truth_df), "列\n")

# 合并数据框（假设都有"cell"列）
merged_df <- merge(clusters_df, truth_df, by = "cell")
cat("合并后数据包含", nrow(merged_df), "行\n")

# 提取真实标签和聚类结果作为因子
ground_truth <- as.factor(merged_df$truth_label)
clusterings <- as.factor(merged_df$cluster)

cat("真实标签类别数:", length(levels(ground_truth)), "\n")
cat("聚类结果类别数:", length(levels(clusterings)), "\n")

# 调用绘图函数
cat("正在生成混淆矩阵图...\n")
plot_obj <- cross_table_plot(ground_truth, clusterings)

# 保存图片
cat("正在保存图片到:", args$output, "\n")
ggsave(filename = args$output, plot = plot_obj, 
       width = 10, height = 12, dpi = 300, units = "in")

cat("图片成功保存到:", args$output, "\n")
cat("脚本执行完成！\n")