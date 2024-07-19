# PTSD and Genetic Associations Analysis

## Overview

This project aims to investigate the relationship between genetic variations and PTSD prevalence, as well as explore the impact of healthcare costs and ethnic fractionalization on PTSD. The analysis utilizes multiple datasets and employs data cleaning, transformation, and visualization techniques to derive insights. The ultimate aim of this project is to create a descriptive and predictive model for the development of PTSD. For a full project presentation please visit: https://docs.google.com/presentation/d/e/2PACX-1vQ2gcAm9AQUy83cOErSa0kt4HMHL2HEyXR7kCngY5GE3dee-vh5tTXXAwP0KqZ5GeTLPGLBymkkk_Yu/pub?start=false&loop=false&delayms=3000

## Project Proposal / Research Question

The primary research question of this project is: How do genetic variations, healthcare costs, and ethnic fractionalization influence the prevalence of PTSD in different populations?

## Data Sources

The usVeteransFKBP5 dataset was manually extracted from Lei Zhang et al.'s paper, "Genetic association of FKBP5 with PTSD in US service members deployed to Iraq and Afghanistan," published in the Journal of Psychiatric Research in 2020. This dataset includes genotype frequencies for specific SNPs associated with PTSD among US service members. The Healthcare Price Index dataset was sourced from the World Bank's Healthcare Price Index from the International Comparison Program, providing healthcare cost indices for various countries. The HIEF dataset was taken from Lenka Drazanova's "Historical Index of Ethnic Fractionalization Dataset (HIEF)," available on the Harvard Dataverse. The global PTSD data was obtained from K. C. Koenen et al.'s study, "Posttraumatic stress disorder in the World Mental Health Surveys," published in Psychol Med in 2017. Lastly, the trauma exposure data was collected from C. Benjet et al.'s research, "The epidemiology of traumatic event exposure worldwide: results from the World Mental Health Survey Consortium."


## Data Cleaning and Transformation

The data cleaning and transformation process involved several steps to ensure consistency and prepare the datasets for analysis. Initially, country names were standardized using the function standardizeCountryNames(df, country_col, country_scope), which mapped various representations of country names to standardized names. Genotype frequencies in the usVeteransFKBP5 dataset were extracted from parentheses and converted to numeric format using regex. The datasets were then merged based on common columns such as country names to create a comprehensive dataset for analysis. Data was categorized into bins for easier analysis and visualization, such as categorizing trauma exposure and PTSD prevalence. The data was reshaped to a long format suitable for plotting using techniques like melting and concatenation. These steps ensure that other researchers can replicate the analysis using the provided code and detailed explanations in the script.

## Visualizations

To contextualize the findings, heatmaps were generated to understand the relationships between PTSD prevalence, trauma exposure, healthcare costs, and ethnic fractionalization. For example, correlation matrices were visualized using seaborn's heatmap function to display the correlations clearly. Scatter plots were created to visualize genotype frequencies across control and PTSD groups, aiding in the comparison of genetic variations. Context-specific visualizations, such as choropleth maps, were used to display the geographical distribution of PTSD prevalence and trauma exposure. Combined category maps highlighted regions with significant findings by displaying combined categories of ethnic fractionalization and PTSD prevalence.

## Key Insights

The analysis confirmed the association between specific SNPs in the FKBP5 gene and PTSD prevalence in US service members. The correlation between healthcare costs and PTSD prevalence is 0.41. Given that higher values in the healthcare price index indicate more expensive healthcare, this suggests a slight negative correlation between the accessibility of healthcare and PTSD prevalence. More accessible (less expensive) healthcare tends to correlate with lower PTSD rates. The inclusion of the HIEF data was to investigate if countries with high PTSD rates but low trauma exposure rates, and high homogeneity, might indicate a genetic factor at play. The correlation between ethnic fractionalization and PTSD prevalence is -0.17. Lower values for the HIEF indicate more homogeneous countries, suggesting that more ethnically diverse countries tend to have lower PTSD prevalence. Analysis of genotype frequencies between control and PTSD groups generally showed a mutation frequency difference of +/- 1-7%. PTSD prevalence varies significantly across different countries, with certain regions exhibiting notably higher rates of trauma exposure and PTSD.


## Specific Country Insights

The United Kingdom (GBR) has the highest PTSD rate with a lifetime PTSD prevalence of 8.8%, trauma exposure of 60.6%, and a healthcare price index of 136.42. China (CHN) has the lowest PTSD rate with a lifetime PTSD prevalence of 0.3%, trauma exposure of 52.5%, and a healthcare price index of 61.64. Ukraine (UKR) has the highest trauma rate with a lifetime PTSD prevalence of 4.8%, trauma exposure of 84.6%, and a healthcare price index of 14.37. Bulgaria (BGR) has the lowest trauma rate with a lifetime PTSD prevalence of 1.9%, trauma exposure of 28.6%, and a healthcare price index of 32.03. Israel (ISR) has the highest healthcare cost with a lifetime PTSD prevalence of 1.6%, trauma exposure of 74.8%, and a healthcare price index of 186.3. Ukraine (UKR) again appears with the lowest healthcare cost, with a lifetime PTSD prevalence of 4.8%, trauma exposure of 84.6%, and a healthcare price index of 14.37.
