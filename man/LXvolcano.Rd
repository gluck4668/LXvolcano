\name{LXvolcano}
\alias{LXvolcano}
\title{LXvolcano:for the analysis of volcano}
\usage{
LXvolcano(gene_data,Gene_FC,meta_data,Meta_FC,group1,group2)
}
\examples{

install.packages("devtools")

library(devtools)

install_github("gluck4668/LXvolcano")

library(LXvolcano)

#--------------------------------

data(gene_data_example)
data(meta_data_example)

#--------------------------------

rm(list=ls())

setwd("D:/R-lin study/R packages/LXvolcano")

gene_data = "gene_df.xlsx" #(if there is no gene data to analyze,  here it should be NULL)

Gene_FC = 2 #(The value should be 2 or NULL )

meta_data = "meta_df.xlsx" #(if there is no gene data to analyze,  here it should be NULL)

Meta_FC = 2 #(The Fold Change should be 2 or NULL )

group1= "Model"  # (Treatment is in front and control is in rear)

group2 = "Normal"

LXvolcano(gene_data,Gene_FC,meta_data,Meta_FC,group1,group2)

}
