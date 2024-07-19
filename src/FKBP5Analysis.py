import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.express as px
import os

# Get the absolute paths to your datasets.
scriptDir = os.path.dirname(os.path.abspath(__file__))
hiefDataPath = os.path.join(scriptDir, '..', 'data', 'hiefData.csv')
usVeteransFkbp5Path = os.path.join(scriptDir, '..', 'data', 'usVeteransFKBP5.csv')
ptsdWorldMentalHealthSurveyPath = os.path.join(scriptDir, '..', 'data', 'ptsdWorldMentalHealthSurvey.csv')
healthCarePriceIndexPath = os.path.join(scriptDir, '..', 'data', 'healthCarePriceIndex.csv')
# Load datasets.
hiefData = pd.read_csv(hiefDataPath, comment='#')
usVeteransFkbp5 = pd.read_csv(usVeteransFkbp5Path, comment='#')
ptsdWorldMentalHealthSurvey = pd.read_csv(ptsdWorldMentalHealthSurveyPath, comment='#')
healthCarePriceIndex = pd.read_csv(healthCarePriceIndexPath, comment='#')
# Define country scope.
countryScope = {
    "COL": ["Colombia", "Colombia (Medellin)"],
    "IRQ": ["Iraq"],
    "PER": ["Peru"],
    "CHN": ["PRC China", "China"],
    "UKR": ["Ukraine"],
    "BRA": ["Brazil"],
    "BGR": ["Bulgaria"],
    "MEX": ["Mexico"],
    "ROU": ["Romania"],
    "ZAF": ["South Africa"],
    "AUS": ["Australia"],
    "BEL": ["Belgium"],
    "DEU": ["Germany", "German Federal Republic"],
    "ISR": ["Israel"],
    "ITA": ["Italy"],
    "JPN": ["Japan"],
    "NZL": ["New Zealand"],
    "GBR": ["United Kingdom", "Northern Ireland"],
    "PRT": ["Portugal"],
    "ESP": ["Spain", "Spain (Murcia)"],
    "NLD": ["The Netherlands", "Netherlands"],
    "USA": ["The USA", "USA", "United States of America"]
}
# Define color map for combined categories.
traumaToPtsdColorMap = {
    'Low Trauma-Low PTSD': '#1b8a5a',
    'Low Trauma-Moderate PTSD': '#9ACD32',
    'Low Trauma-High PTSD': '#FFD700',
    'Moderate Trauma-Low PTSD': '#1d4877',
    'Moderate Trauma-Moderate PTSD': '#ffd73e',
    'Moderate Trauma-High PTSD': '#FF8C00',
    'High Trauma-Low PTSD': '#f68838',
    'High Trauma-Moderate PTSD': '#ee3e32',
    'High Trauma-High PTSD': '#b40426'
}
# Define color map for EFindex and PTSD combined categories.
efToPtsdColorMap = {
    'High Homogeneity-Low PTSD': '#1b8a5a',
    'High Homogeneity-Moderate PTSD': '#9ACD32',
    'High Homogeneity-High PTSD': '#FFD700',
    'Moderate Homogeneity-Low PTSD': '#1d4877',
    'Moderate Homogeneity-Moderate PTSD': '#ffd73e',
    'Moderate Homogeneity-High PTSD': '#FF8C00',
    'Low Homogeneity-Low PTSD': '#f68838',
    'Low Homogeneity-Moderate PTSD': '#ee3e32',
    'Low Homogeneity-High PTSD': '#FF0000'
}
# Function to standardize country names to their ISO codes.
def standardizeCountryNames(df, countryCol, countryScope):
    # Creating a dictionary to map country names to ISO codes.
    standardizedCountries = {}
    # Looping through the countryScope dictionary to fill the standardizedCountries dictionary.
    for isoCode, names in countryScope.items():
        for name in names:
            standardizedCountries[name] = isoCode
    # Mapping the country names in the DataFrame to their corresponding ISO codes.
    # Using .fillna() to ensure any country name not found in the dictionary remains unchanged.
    df[countryCol] = df[countryCol].map(
        standardizedCountries).fillna(df[countryCol])
    return df
# Function to create a correlation heatmap.
def createCorrelationHeatmap(df, columns, title):
    # Creating a copy of the DataFrame and renaming columns for better readability in the heatmap.
    dfCopy = df[columns].copy()
    dfCopy.rename(columns={columns[0]: 'Lifetime PTSD Prevalence', columns[1]                  : 'Trauma Exposure', columns[2]: 'Healthcare Price Index'}, inplace=True)
    # Creating the correlation matrix from the DataFrame.
    corrMatrix = dfCopy.corr()
    # Setting the figure size for the heatmap.
    plt.figure(figsize=(10, 8))
    # Generating the heatmap using seaborn.
    sns.heatmap(corrMatrix, annot=True, cmap='coolwarm', vmin=-1, vmax=1)
    plt.title(title)
    plt.xticks(rotation=0)
    plt.yticks(rotation=90)
    plt.show()
# Function to categorize trauma and PTSD data.
def categorizeData(df, traumaBins, ptsdBins):
    # Categorizing trauma exposure using pd.cut() into bins labeled 'Low Trauma', 'Moderate Trauma', and 'High Trauma'.
    # The float('inf') is used to represent infinity, capturing all values above the specified bin thresholds.
    df['traumaCategory'] = pd.cut(df['Trauma exposure (%)'], bins=traumaBins, labels=[
                                  'Low Trauma', 'Moderate Trauma', 'High Trauma']).astype(str)
    # Categorizing PTSD prevalence similarly.
    df['ptsdCategory'] = pd.cut(df['Lifetime prevalence of PTSD in total sample (%)'], bins=ptsdBins, labels=[
                                'Low PTSD', 'Moderate PTSD', 'High PTSD']).astype(str)
    # Creating a combined category column by concatenating the trauma and PTSD categories.
    # .fillna('') is used to replace NaN values with empty strings to avoid issues during concatenation.
    df['combinedCategory'] = df['traumaCategory'].fillna(
        '') + '-' + df['ptsdCategory'].fillna('')
    return df
# Function to create a choropleth map.
def createChoroplethMap(df, colorMap, title):
    # Creating a column for country codes, assuming they are the same as the 'Country' column.
    df['countryCode'] = df['Country']
    # Mapping the combined category to the color map.
    df['Color'] = df['combinedCategory'].map(colorMap)
    # Creating the choropleth map using Plotly.
    fig = px.choropleth(df, locations="countryCode", color="combinedCategory", hover_name="Country",
                        hover_data=[
                            "Lifetime prevalence of PTSD in total sample (%)", "Trauma exposure (%)"],
                        color_discrete_map=colorMap, title=title)
    # Updating geographic layout for better visualization.
    fig.update_geos(showframe=True, showcoastlines=True,
                    projection_type='orthographic')
    fig.show()
# Function to calculate average EFindex and create correlation heatmap.
def calculateAverageEfindexAndPlot(df1, df2, countryCol, efindexCol, ptsdCol, title):
    # Calculating the average EFindex for each country.
    avgEfindex = df2.groupby(countryCol)[efindexCol].mean().reset_index()
    # Merging the PTSD data with the average EFindex data.
    mergedDf = pd.merge(df1[[countryCol, ptsdCol]], avgEfindex, on=countryCol)
    # Renaming columns for better readability in the heatmap.
    mergedDf.rename(columns={ptsdCol: 'Lifetime PTSD Prevalence',
                    efindexCol: 'Ethnic Homogeneity'}, inplace=True)
    # Creating the correlation matrix.
    corrMatrix = mergedDf[['Lifetime PTSD Prevalence',
                           'Ethnic Homogeneity']].corr()
    # Setting the figure size for the heatmap.
    plt.figure(figsize=(10, 8))
    # Generating the heatmap using seaborn.
    sns.heatmap(corrMatrix, annot=True, cmap='coolwarm', vmin=-1, vmax=1)
    plt.title(title)
    plt.xticks(rotation=0)
    plt.yticks(rotation=90)
    plt.show()
# Function to extract SNP frequencies and reshape data.
def extractAndReshapeSnpData(df):
    # Extracting numeric values from the genotype frequency columns within parentheses.
    df['Frequency 1:1'] = df['1:1'].str.extract(
        r'\((\d+\.\d+)\)').astype(float)
    df['Frequency 1:2'] = df['1:2'].str.extract(
        r'\((\d+\.\d+)\)').astype(float)
    df['Frequency 2:2'] = df['2:2'].str.extract(
        r'\((\d+\.\d+)\)').astype(float)
    # Adding detailed labels for each genotype.
    df['Snp_1:1'] = df['SNP'] + ' 1:1'
    df['Snp_1:2'] = df['SNP'] + ' 1:2'
    df['Snp_2:2'] = df['SNP'] + ' 2:2'
    # Reshaping the data for plotting.
    reshapedData11 = df[['Snp_1:1', 'Frequency 1:1', 'Group']].rename(
        columns={'Snp_1:1': 'SnpGenotype', 'Frequency 1:1': 'Frequency'})
    reshapedData12 = df[['Snp_1:2', 'Frequency 1:2', 'Group']].rename(
        columns={'Snp_1:2': 'SnpGenotype', 'Frequency 1:2': 'Frequency'})
    reshapedData22 = df[['Snp_2:2', 'Frequency 2:2', 'Group']].rename(
        columns={'Snp_2:2': 'SnpGenotype', 'Frequency 2:2': 'Frequency'})
    # Combining the reshaped data into a single DataFrame.
    return pd.concat([reshapedData11, reshapedData12, reshapedData22])
# Function to create scatter plot for SNP data.
def createScatterPlot(df, title):
    # Generating scatter plot using Plotly.
    fig = px.scatter(df, x='SnpGenotype', y='Frequency', color='Group', title=title,
                     labels={'SnpGenotype': 'SNP and Genotype',
                             'Frequency': 'Frequency'},
                     symbol='Group', size_max=10)
    # Updating layout for better aesthetics.
    fig.update_layout(title_font_size=24, xaxis_title_font_size=18, yaxis_title_font_size=18,
                      legend_title_font_size=16, legend_font_size=14, xaxis_tickangle=-90,
                      xaxis_tickfont_size=12, yaxis_tickfont_size=14, width=1200, height=800,
                      margin=dict(l=40, r=40, t=80, b=40))
    fig.show()
# Standardize country names.
ptsdWorldMentalHealthSurvey = standardizeCountryNames(ptsdWorldMentalHealthSurvey, 'Country', countryScope)
healthCarePriceIndex = standardizeCountryNames(healthCarePriceIndex, 'Country', countryScope)
# Merge healthcare price index with PTSD and trauma exposure data.
ptsdAndHealthcare = ptsdWorldMentalHealthSurvey.merge(healthCarePriceIndex[['Country','Healthcare price index world average = 100']], on='Country', how='inner')
# Verify that the merged DataFrame is not empty and create correlation heatmap.
if not ptsdAndHealthcare.empty:
    createCorrelationHeatmap(ptsdAndHealthcare,
                             ['Lifetime prevalence of PTSD in total sample (%)', 'Trauma exposure (%)',
                              'Healthcare price index world average = 100'],
                             "Lifetime Prevalence Of PTSD Compared To Trauma Exposure And The Cost Of Healthcare")
# Categorize and create combined categories.
ptsdWorldMentalHealthSurveyCopy = categorizeData(ptsdWorldMentalHealthSurvey.copy(),traumaBins=[-1, 55, 70, float('inf')], ptsdBins=[-1, 2.7, 5.3, float('inf')])
# Create choropleth map for PTSD and trauma exposure.
createChoroplethMap(ptsdWorldMentalHealthSurveyCopy, traumaToPtsdColorMap,
                    "Lifetime PTSD Prevalence & Trauma Exposure by Country")
# Standardize country names in hiefData.
hiefData = standardizeCountryNames(hiefData, 'Country', countryScope)
# Calculate average EFindex and create correlation heatmap.
calculateAverageEfindexAndPlot(ptsdWorldMentalHealthSurvey, hiefData,
                               'Country', 'EFindex', 'Lifetime prevalence of PTSD in total sample (%)',
                               'Correlation Between Lifetime PTSD Prevalence and EFindex')
# Merge PTSD data with EFindex data.
mentalHealthAndHiefComparisonData = pd.merge(ptsdWorldMentalHealthSurveyCopy, hiefData, on='Country')
# Categorize EFindex and PTSD data.
mentalHealthAndHiefComparisonData['efindexCategory'] = pd.cut(mentalHealthAndHiefComparisonData['EFindex'], bins=[-1, 0.3, 0.6, float('inf')], labels=['High Homogeneity', 'Moderate Homogeneity', 'Low Homogeneity']).astype(str)
mentalHealthAndHiefComparisonData['ptsdCategory'] = pd.cut(mentalHealthAndHiefComparisonData['Lifetime prevalence of PTSD in total sample (%)'], bins=[-1, 2.7, 5.3, float('inf')], labels=['Low PTSD', 'Moderate PTSD', 'High PTSD']).astype(str)
mentalHealthAndHiefComparisonData['combinedCategory'] = mentalHealthAndHiefComparisonData['efindexCategory'].fillna('') + '-' + mentalHealthAndHiefComparisonData['ptsdCategory'].fillna('')
# Create choropleth map for EFindex and PTSD combined categories.
createChoroplethMap(mentalHealthAndHiefComparisonData, efToPtsdColorMap,"Ethnic Fractionalization Index & PTSD Prevalence by Country")
# Extract and reshape SNP data for plotting.
reshapedSnpData = extractAndReshapeSnpData(usVeteransFkbp5)
# Create scatter plot for SNP data.
createScatterPlot(reshapedSnpData, "Genotype Frequencies in Control and PTSD Groups")
