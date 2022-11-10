
LXvolcano <- function (gene_data,Gene_FC,meta_data,Meta_FC,group1,group2){

  all_packages <- data.frame(installed.packages())

  pack <- data.frame(c("ggplot2","openxlsx","stringr","plyr","psych","dplyr","ggrepel",
                       "patchwork","raster","png","janitor","conflicted"))

  pack$type <- pack[,1] %in% all_packages$Package

  for (i in 1:nrow(pack)){if (!requireNamespace(pack[i,1], quietly=TRUE))
    install.packages(pack[i,1],update = F,ask = F)
  }
  rm(i)

  packages <- as.character(pack[,1])

  for(i in packages){
    library(i, character.only = T)
  }
  rm(i)


  conflict_scout()

  conflict_prefer("%+%", "psych")

  conflict_prefer("filter", "dplyr")
  conflict_prefer("lag", "dplyr")
  conflict_prefer("arrange", "dplyr")
  conflict_prefer("select", "dplyr")
  conflict_prefer("summarise", "dplyr")
  conflict_prefer("summarize", "dplyr")
  conflict_prefer("rename", "dplyr")
  conflict_prefer("mutate", "dplyr")
  conflict_prefer("count", "dplyr")
  conflict_prefer("failwith", "dplyr")
  conflict_prefer("id", "dplyr")
  conflict_prefer("combine", "dplyr")
  conflict_prefer("desc", "dplyr")
  conflict_prefer("collapse", "dplyr")
  conflict_prefer("slice", "dplyr")

  conflict_prefer("Position", "ggplot2")


#----- Gene volcano function----------------

  gene_volcano <- function(data_file) {

    df <- read.xlsx(data_file)

    df <- na.omit(df)

    group_gene <- colnames(df)[c(5:ncol(df))]
    group_df <- gsub("\\d+$","", group_gene)

    if(tolower(group1) %in% tolower(group_df)==FALSE)
    { g_gene <- paste(group_gene,collapse = ",")
      group_g <- paste("The groups in the gene data file are",g_gene)
      group_j <- paste("'",group1,"'","is not found in the gene data file. Please check it")
      print(group_g)
      print(group_j)
      print("-----------------------------------------------------------")} else
      group_df <- group_df

    if(tolower(group2) %in% tolower(group_df)==FALSE)
    { g_gene <- paste(group_gene,collapse = ",")
      group_g <- paste("The groups in the gene data file are",g_gene)
      group_j <- paste("'",group2,"'","is not found in the gene data file. Please check it")
      print(group_g)
      print(group_j)
      print("-----------------------------------------------------------")} else
      rm(group_df)

    ####################################################
    df$sum <-rowSums(df[,c(5:ncol(df))])
    df <- df[df[,ncol(df)]>0,]
    df <- df[,-ncol(df)]

    df$sum57 <- rowMeans(df[,c(5:7)])
    df$sum810 <- rowMeans(df[,c(8:10)])

    group_m <- colnames(df)[c(5:10)]
    group_m <- gsub("\\d+$","", group_m)

    colnames(df)[11] <- group_m[1]
    colnames(df)[12] <- group_m[length(group_m)]

    df_max <- df[df[,2]==max(df[,2]),]

    n1 <- ncol(df_max)-1
    n2 <- ncol(df_max) %>% as.numeric()

    if(df_max[1,n1]<df_max[1,n2])
      df_max$high <- colnames(df_max)[n2] else
        df_max$high <- colnames(df_max)[n1]

    if(tolower(group1) %in% tolower(df_max[1,ncol(df_max)]))
      title_gene <- c(paste("Gene volcano graphics","(",group1,"VS",group2,")")) else
        title_gene <- c(paste("Gene volcano graphics","(",group2,"VS",group1,")"))

    df <- df[,c(1:4)]
    #########################################

    for(i in 1:nrow(df)){
      if(substring(tolower(df[i,2]),1,3)=="inf")
        df[i,2]=NA}
    rm(i)

    df <- na.omit(df)

    df[,2] <- as.numeric(df[,2])
    df[,3] <- as.numeric(df[,3])
    df[,4] <- as.numeric(df[,4])

    colnames(df) <- c("ID","log2FC","pvalue","qvalue")

    df <- df[order(-abs(df[,2])),]

    data <- data.frame(distinct(df, ID, .keep_all = TRUE))


    if(is.null(Gene_FC)==TRUE)
    {data$type <- ifelse(data$log2FC>0 & data$pvalue<0.05,"Up",
                         ifelse(data$log2FC<0 & data$pvalue<0.05,"Down","Not sig")) } else
                         {data$type <- ifelse(data$log2FC>=log2(Gene_FC) & data$pvalue<0.05,"Up",
                                              ifelse(data$log2FC<=-log2(Gene_FC) & data$pvalue<0.05,"Down","Not sig"))}

    change_n <- table(data$type)
    change_nm <-data.frame(change_n)

    Down <- c(paste("Down",c(change_nm[1,2])))#Down
    No <- c(paste("Not sig",c(change_nm[2,2])))#Not sig
    UP <- c(paste("Up",c(change_nm[3,2])))#Up

    data$changes <- ifelse(data$type=="Up",UP,
                           ifelse(data$type=="Down",Down,No))

    mytheme<-theme_bw()+
      theme(text=element_text(family = "sans",colour ="black",face="bold",size =14),
            #panel.border = element_rect (linewidth = 0.8,color = "gray30"),
            axis.line = element_blank(),
            axis.ticks = element_line(linewidth = 0.6,colour = "gray30"),
            axis.ticks.length = unit(1.5,units = "mm"))+
      theme(plot.title = element_text(hjust = 0.5))+
      theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5),"cm"))

    xytheme <-theme(axis.text.x = element_text(face="bold",color="black",size=12,angle =0,hjust=1))+
      theme(axis.text.y = element_text(face="bold",color="black",size=12))+
      theme(legend.text=element_text(face="bold",color="black",size=12))

    p1 <- ggplot(data,aes(x=log2FC,y=-log10(pvalue)))+
      geom_point(aes(color=changes))+
      scale_color_manual(values = c("#00bfff","#c0c0c0","#ff4500"))+
      geom_hline(yintercept = -log10(0.05),linetype="dashed",color="#808080")+
      geom_vline(xintercept = c(-1,1),linetype="dashed",color="#808080")+
      labs(x="log2 (Fold Change)",y="-log10(p value)",title = title_gene) +
      mytheme+xytheme

    p1

    data01 <- data

    data01$p_log10 <- (-log10(data$pvalue))

    df_order <- data01[data01[,7]>1.301029996,]

    df_order <- df_order[order(-df_order[,7]),]

    i <- case_when(nrow(df_order)>1000 ~0.01,
                   nrow(df_order)>500 ~0.03,
                   nrow(df_order)>200 ~0.04,
                   nrow(df_order)>100 ~0.05,
                   TRUE ~0.1)

    show_n <- round(nrow(df_order)*i,0)

    df_n <- case_when(show_n >25 ~25,
                      TRUE ~show_n)

    df_new <- df_order[c(1:df_n),]

    df_show <- df_new[abs(df_new$log2FC)>=1,]

    p2<- p1+
      geom_label_repel(data=df_show,aes(x=log2FC,y=-log10(pvalue), label=ID),
                       label.size =0.1,size=2, box.padding = 0.4, max.overlaps =show_n)


    data$log2FoldChange <- abs(data$log2FC)


    label01=paste("       *",Down)
    label02=paste("          *",No)
    label03=paste("*",UP)

    label_all <- paste(label01,"\n", label02,"\n",label03)


    mytheme3<-theme_bw()+
      theme(text=element_text(family = "sans",colour ="black",face="bold",size =14),
            axis.line = element_blank(),
            axis.ticks = element_line(size = 0.6,colour = "gray30"),
            axis.ticks.length = unit(1.5,units = "mm"))+
      theme(plot.title = element_text(hjust = 0.5))+
      theme(plot.margin = unit(c(0.3, 5.5, 0.5, 1.0), "cm"))

    xytheme3 <-theme(axis.text.x = element_text(face="bold",color="black",size=12,angle =0,hjust=1))+
      theme(axis.text.y = element_text(face="bold",color="black",size=12))+
      theme(legend.text=element_text(face="bold",color="black",size=12),
            legend.position = c(1.2,0.40))

    # ------------------------------
    tag_theme3 <- theme(plot.tag = element_text(size =15,colour = "black"),
                        plot.tag.position = c(1.10,0.83))

    p3 <- ggplot(data,aes(x=log2FC,y=-log10(pvalue)))+
      geom_point(aes(color=-log10(pvalue),size=log2FoldChange))+
      scale_color_gradientn(values=seq(0,1,0.2),
                            colors=c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
      scale_size_continuous(rang=c(1,3))+
      geom_hline(yintercept = -log10(0.05),linetype="dashed",color="#808080")+
      geom_vline(xintercept = c(-1.2,1.2),linetype="dashed",color="#808080")+
      labs(x="log2 (Fold Change)",y="-log10(p value)",title = title_gene,tag =label_all)+
      mytheme3+xytheme3+tag_theme3
    #annotate("text",x=-Inf,y=Inf,hjust=-0.1,vjust=2,label=paste("???",Down),family="serif",size=4,colour="#006400")+
    #annotate("text",x=-Inf,y=Inf,hjust=-0.1,vjust=4,label=paste("???",No),family="serif",size=4,colour="#006400")+
    #annotate("text",x=-Inf,y=Inf,hjust=-0.15,vjust=6,label=paste("???",UP),family="serif",size=4,colour="#006400")

    p3

    p4 <- p3+
      geom_label_repel(data=df_show,aes(x=log2FC,y=-log10(pvalue), label=ID),#只显示绝对值log2FC>=1的ID
                       label.size =0.1,size=2, box.padding = 0.4, max.overlaps =show_n)

    p4

    if(dir.exists("analysis result")==FALSE)
      dir.create("analysis result")

    if(is.null(Gene_FC)==TRUE)
    {gene_df <- data[,c(1:5)]
    gene_df <- gene_df[order(gene_df[,5]),]
    write.xlsx(gene_df,"analysis result/Gene volcano data (P0.05, unlimited FC)_all.xlsx")

    gene_dd <- gene_df
    for(i in 1:nrow(gene_dd)){
      if(substring(tolower(gene_dd[i,5]),1,3)=="not")
        gene_dd[i,5] <- NA}
    rm(i)
    gene_dd <- na.omit(gene_dd)
    write.xlsx(gene_dd,"analysis result/Gene volcano data (P0.05, unlimited FC)_diferent.xlsx")

    ggsave("analysis result/Gene Volcano graphics 01 (P0.05, unlimited FC).png",p1, width=1200, height =1000, dpi=180,units = "px")
    ggsave("analysis result/Gene Volcano graphics 02 (P0.05, unlimited FC).png",p2, width=1200, height =1000, dpi=180,units = "px")
    ggsave("analysis result/Gene Volcano graphics 03 (P0.05, unlimited FC).png",p3, width=1400, height =1000, dpi=180,units = "px")
    ggsave("analysis result/Gene Volcano graphics 04 (P0.05, unlimited FC).png",p4, width=1400, height =1000, dpi=180,units = "px")} else
    {gene_df <- data[,c(1:5)]
    gene_df <- gene_df[order(gene_df[,5]),]
    write.xlsx(gene_df,"analysis result/Gene volcano data (P0.05, FC2)_all.xlsx")

    gene_dd <- gene_df
    for(i in 1:nrow(gene_dd)){
      if(substring(tolower(gene_dd[i,5]),1,3)=="not")
        gene_dd[i,5] <- NA}
    rm(i)
    gene_dd <- na.omit(gene_dd)
    write.xlsx(gene_dd,"analysis result/Gene volcano data (P0.05, FC2)_diferent.xlsx")

    ggsave("analysis result/Gene Volcano graphics 01 (P0.05, FC2).png",p1, width=1200, height =1000, dpi=180,units = "px")
    ggsave("analysis result/Gene Volcano graphics 02 (P0.05, FC2).png",p2, width=1200, height =1000, dpi=180,units = "px")
    ggsave("analysis result/Gene Volcano graphics 03 (P0.05, FC2).png",p3, width=1400, height =1000, dpi=180,units = "px")
    ggsave("analysis result/Gene Volcano graphics 04 (P0.05, FC2).png",p4, width=1400, height =1000, dpi=180,units = "px")}

    #p_gene <- p1+p2+p3+p4

    #dev.new(width=14,height=9, noRStudioGD = T)

    #p_gene

  }
  #---------------------------------------------------------

  #----- Metabolite volcano function----------------

  meta_volcano <- function(data_file) {

    df <- data_file

    colnames(df) <- c("ID","log2FC","pvalue","qvalue")

    data <- data.frame(distinct(df, ID, .keep_all = TRUE))

    if(is.null(Meta_FC)==TRUE)
    {data$type <- ifelse(data$log2FC>0 & data$pvalue<0.05,"Up",
     ifelse(data$log2FC<0 & data$pvalue<0.05,"Down","Not sig")) } else
                         {data$type <- ifelse(data$log2FC>=log(Meta_FC) & data$pvalue<0.05,"Up",
                        ifelse(data$log2FC<=-log(Meta_FC) & data$pvalue<0.05,"Down","Not sig"))}

    change_n <- table(data$type)
    change_nm <-data.frame(change_n)

    if("Not sig" %in% data$type){
      Down <- c(paste("Down",c(change_nm[1,2])))
      No <- c(paste("Not sig",c(change_nm[2,2])))
      UP <- c(paste("Up",c(change_nm[3,2])))} else
      {Down <- c(paste("Down",c(change_nm[1,2])))
      UP <- c(paste("Up",c(change_nm[2,2])))}


    if("Not sig" %in% data$type){
      data$changes <- ifelse(data$type=="Up",UP,
                             ifelse(data$type=="Down",Down,No))} else
                               data$changes <- ifelse(data$type=="Up",UP,Down)

    title_meta <- c(paste("Metabolite volcano graphics","(",group1,"VS",group2,")"))

    mytheme<-theme_bw()+
      theme(text=element_text(family = "sans",colour ="black",face="bold",size =14),
            axis.line = element_blank(),
            axis.ticks = element_line(size = 0.6,colour = "gray30"),
            axis.ticks.length = unit(1.5,units = "mm"))+
      theme(plot.title = element_text(hjust = 0.5))+
      theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

    xytheme <-theme(axis.text.x = element_text(face="bold",color="black",size=12,angle =0,hjust=1))+
      theme(axis.text.y = element_text(face="bold",color="black",size=12))+
      theme(legend.text=element_text(face="bold",color="black",size=12))

    # 指定图外，文字大小，颜色，以及位???
    tag_theme <- theme(plot.tag = element_text(size =14,colour = "blue"),
                       plot.tag.position = c(0.80,0.9)) #图形大小：长和高均为1

    p01 <- ggplot(data,aes(x=log2FC,y=-log10(pvalue)))+
      geom_point(aes(color=changes))+
      scale_color_manual(values = c("#00bfff","#c0c0c0","#ff4500"))+
      geom_hline(yintercept = -log10(0.05),linetype="dashed",color="#808080")+
      geom_vline(xintercept = c(-1,1),linetype="dashed",color="#808080")+
      labs(x="log2 (Fold Change)",y="-log10(p value)",title = title_meta) +
      mytheme+xytheme

    data01 <- data

    data01$p_log10 <- (-log10(data$pvalue))

    df_order <- data01[data01[,7]>1.301029996,]

    df_order <- df_order[order(-df_order[,7]),]

    i <- case_when(nrow(df_order)>1000 ~0.01,
                   nrow(df_order)>500 ~0.03,
                   nrow(df_order)>200 ~0.04,
                   nrow(df_order)>100 ~0.05,
                   TRUE ~0.1)

    show_n <- round(nrow(df_order)*i,0)

    df_n <- case_when(show_n >25 ~25,
                      TRUE ~show_n)

    df_new <- df_order[c(1:df_n),]

    df_show <- df_new[abs(df_new$log2FC)>=1,]

    p02<- p01+
      geom_label_repel(data=df_show,aes(x=log2FC,y=-log10(pvalue), label=ID),
                       label.size =0.1,size=2, box.padding = 0.4, max.overlaps =show_n)

    data$log2FoldChange <- abs(data$log2FC)

    # 设定文字标签（图形之外）
    label01=paste("     *",Down)
    label02=paste("      *",No)
    label03=paste("*",UP)

    label_all <- paste(label01,"\n", label02,"\n",label03)

    mytheme03<-theme_bw()+
      theme(text=element_text(family = "sans",colour ="black",face="bold",size =14),
            axis.line = element_blank(),
            axis.ticks = element_line(size = 0.6,colour = "gray30"),
            axis.ticks.length = unit(1.5,units = "mm"))+
      theme(plot.title = element_text(hjust = 0.5))+
      theme(plot.margin = unit(c(0.3, 5.5, 0.5, 1.0), "cm"))#对应着 上、右、下、左 的顺???

    xytheme03 <-theme(axis.text.x = element_text(face="bold",color="black",size=12,angle =0,hjust=1))+
      theme(axis.text.y = element_text(face="bold",color="black",size=12))+
      theme(legend.text=element_text(face="bold",color="black",size=12),
            legend.position = c(1.2,0.40))

    # 指定图外，文字大小，颜色，以及位???
    tag_theme03 <- theme(plot.tag = element_text(size =15,colour = "black"),
                         plot.tag.position = c(1.10,0.83)) #图形大小：长和高均为1

    p03 <- ggplot(data,aes(x=log2FC,y=-log10(pvalue)))+
      geom_point(aes(color=-log10(pvalue),size=log2FoldChange))+
      scale_color_gradientn(values=seq(0,1,0.2),
                            colors=c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
      scale_size_continuous(rang=c(1,3))+
      geom_hline(yintercept = -log10(0.05),linetype="dashed",color="#808080")+
      geom_vline(xintercept = c(-1.2,1.2),linetype="dashed",color="#808080")+
      labs(x="log2 (Fold Change)",y="-log10(p value)",title = title_meta,tag =label_all)+
      mytheme03+xytheme03+tag_theme03

    p04 <- p03+
      geom_label_repel(data=df_show,aes(x=log2FC,y=-log10(pvalue), label=ID),#只显示绝对值log2FC>=1的ID
                       label.size =0.1,size=2, box.padding = 0.4, max.overlaps =show_n)

    if(dir.exists("analysis result")==FALSE)
      dir.create("analysis result")

    if(is.null(Meta_FC)==TRUE)
    {meta_df <- data[,c(1:5)]
    meta_df <- meta_df[order(meta_df[,5]),]
    write.xlsx(meta_df,"analysis result/Metabolite volcano data (P0.05, unlimited FC)_all.xlsx")

    meta_dd <- meta_df
    for(i in 1:nrow(meta_dd)){
      if(substring(tolower(meta_dd[i,5]),1,3)=="not")
        meta_dd[i,5] <- NA}
    rm(i)
    meta_dd <- na.omit(meta_dd)
    write.xlsx(meta_dd,"analysis result/Metabolite volcano data (P0.05, unlimited FC)_diferent.xlsx")

    ggsave("analysis result/Metabolite Volcano graphics 01 (P0.05, unlimited FC).png",p01, width=1200, height =1000, dpi=180,units = "px")
    ggsave("analysis result/Metabolite Volcano graphics 02 (P0.05, unlimited FC).png",p02, width=1200, height =1000, dpi=180,units = "px")
    ggsave("analysis result/Metabolite Volcano graphics 03 (P0.05, unlimited FC).png",p03, width=1400, height =1000, dpi=180,units = "px")
    ggsave("analysis result/Metabolite Volcano graphics 04 (P0.05, unlimited FC).png",p04, width=1400, height =1000, dpi=180,units = "px")} else
    {meta_df <- data[,c(1:5)]
    meta_df <- meta_df[order(meta_df[,5]),]
    write.xlsx(meta_df,"analysis result/Metabolite volcano data (P0.05, FC2)_all.xlsx")

    meta_dd <- meta_df
    for(i in 1:nrow(meta_dd)){
      if(substring(tolower(meta_dd[i,5]),1,3)=="not")
        meta_dd[i,5] <- NA}
    rm(i)
    meta_dd <- na.omit(meta_dd)
    write.xlsx(meta_dd,"analysis result/Metabolite volcano data (P0.05, FC2)_diferent.xlsx")

    ggsave("analysis result/Metabolite Volcano graphics 01 (P0.05, FC2).png",p01, width=1200, height =1000, dpi=180,units = "px")
    ggsave("analysis result/Metabolite Volcano graphics 02 (P0.05, FC2).png",p02, width=1200, height =1000, dpi=180,units = "px")
    ggsave("analysis result/Metabolite Volcano graphics 03 (P0.05, FC2).png",p03, width=1400, height =1000, dpi=180,units = "px")
    ggsave("analysis result/Metabolite Volcano graphics 04 (P0.05, FC2).png",p04, width=1400, height =1000, dpi=180,units = "px")}

    #p_meta <- p01+p02+p03+p04

    # dev.new(width=14,height=9, noRStudioGD = T)

    #p_meta

  }

  #---------gene volcano -----------------------------------

  if(is.null(gene_data)==TRUE)
    print("The gene_data is not found!") else

      gene_volcano(data_file=gene_data)


  #------metabolite volcano-------------------
  if(is.null(meta_data)==TRUE)
    print("The meta_data is not found!") else

    {meta <- read.xlsx(meta_data)
    meta <- na.omit(meta)
    colnames(meta)[1]<- c("ID")

    group_df <- meta[,3]
    group_meta <- data.frame(group_df)
    colnames(group_meta) <- "names"
    group_meta <- distinct(.data = group_meta,.keep_all = T)
    group_meta[,1] <- gsub("[[:punct:]]", "", group_meta[,1])
    group_meta <- dplyr::filter(group_meta,nchar(group_meta[,1])>1)
    group_meta <- t(group_meta)  %>% as.character()

    if(tolower(group1) %in% tolower(group_df)==FALSE)
    {g_meta <- paste(group_meta,collapse = ",")
     group_m <- paste("The groups in the metabolite data file are",g_meta)
     group_j <- paste("'",group1,"'","is not found in the metabolite data file. Please check it")
     print(group_m)
     print(group_j)
     print("-----------------------------------------------------------")} else
     group_df <- meta[,3]

    if(tolower(group2) %in% tolower(group_df)==FALSE)
    {  g_meta <- paste(group_meta,collapse = ",")
       group_m <- paste("The groups in the metabolite data file are",g_meta)
       group_j <- paste("'",group2,"'","is not found in the metabolite data file. Please check it")
       print(group_m)
       print(group_j)
       print("-----------------------------------------------------------")} else
      rm(group_df)

    colnames(meta)[c(4:5)] <- c("pvalue","pajd")

    meta$log2FC <- NA

    for(i in 1:nrow(meta)){
      if(substring(tolower(meta[i,2]),1,3)=="inf" & tolower(meta[i,3])==tolower(group1) )
        meta[i,2] <- 100 }
    rm(i)

    for(i in 1:nrow(meta)){
      if(substring(tolower(meta[i,2]),1,3)=="inf" & tolower(meta[i,3])==tolower(group2))
        meta[i,2] <- 100}
    rm(i)

    meta[,2] <- as.numeric(meta[,2])


    i=1


    for(i in 1:nrow(meta)){
      if(tolower(meta[i,3])==tolower(group1))
        meta[i,6] <- log2(meta[i,2]) else
          meta[i,6] <- log2(1/meta[i,2])}
    rm(i)

    meta <- meta[order(-abs(meta[,6])),]

    meta <- data.frame(distinct(meta, ID, .keep_all = TRUE))

    meta_data <- dplyr::select(.data = meta,ID,log2FC,pvalue,pajd)

    meta_volcano(data_file=meta_data)
    }

  print("===========================================================")
  print("Please see the results in the folder of <analysis result> ")


  if(is.null(Meta_FC)==TRUE)
  {if(file.exists("analysis result/Gene Volcano graphics 01 (P0.05, unlimited FC).png")==TRUE)
  {plot_png <- readPNG("analysis result/Gene Volcano graphics 01 (P0.05, unlimited FC).png")
  r <- nrow(plot_png)/ncol(plot_png)
  plot(c(0,1),c(0,r),type="n",xlab="",ylab="",asp=1)
  rasterImage(plot_png,0,0,1,r)}} else
  {if(file.exists("analysis result/Gene Volcano graphics 01 (P0.05, FC2).png")==TRUE)
  {plot_png <- readPNG("analysis result/Gene Volcano graphics 01 (P0.05, FC2).png")
  r <- nrow(plot_png)/ncol(plot_png)
  plot(c(0,1),c(0,r),type="n",xlab="",ylab="",asp=1)
  rasterImage(plot_png,0,0,1,r)}}



  if(is.null(Meta_FC)==TRUE)
  {if(file.exists("analysis result/Metabolite Volcano graphics 01 (P0.05, unlimited FC).png")==TRUE)
  {plot_png <- readPNG("analysis result/Metabolite Volcano graphics 01 (P0.05, unlimited FC).png")
  r <- nrow(plot_png)/ncol(plot_png)
  plot(c(0,1),c(0,r),type="n",xlab="",ylab="",asp=1)
  rasterImage(plot_png,0,0,1,r)}} else
  {if(file.exists("analysis result/Metabolite Volcano graphics 01 (P0.05, FC2).png")==TRUE)
  {plot_png <- readPNG("analysis result/Metabolite Volcano graphics 01 (P0.05, FC2).png")
  r <- nrow(plot_png)/ncol(plot_png)
  plot(c(0,1),c(0,r),type="n",xlab="",ylab="",asp=1)
  rasterImage(plot_png,0,0,1,r)}}


}

