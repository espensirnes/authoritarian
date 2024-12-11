import pandas as pd
from fuzzywuzzy import process
from matplotlib import pyplot as plt 
import statsmodels.api as sm
from scipy.stats import ttest_ind
import numpy as np
from statsmodels.stats.diagnostic import het_breuschpagan, het_white


MINVAL = 0.01
def main():
	#reformat()
	df = pd.read_csv('LEGITIMACY_GILEY_format.txt',
				  names = ['Country', 'Unweighted','Weighted', 'Attitudinal'] , sep='\t') 
	df_polity = pd.read_csv("data_new.csv")
	df_polity = df_polity.rename(columns={'ifhpol':'Polity2'})
	df_recent = df_polity[df_polity['year'] >= 2015]
	df_polity_avg = df_recent.groupby('country')['Polity2'].mean().reset_index()
	df = fuzzymatch(df, df_polity_avg)
	df = df.dropna()

	df['lnPolity2'] = np.log(df['Polity2']+MINVAL)
	df['lnPolity2sq'] = df['lnPolity2'] **2
	df['lnPolity2cube'] = df['lnPolity2'] **3

	fig, ax = plt.subplots()

	ax.scatter(df['Polity2'] , df['Attitudinal'] )
	ax.set_xlabel('Polity2')
	ax.set_ylabel('Attitudinal')
	
	
	

	b = reg(df).params

	x = np.exp(np.linspace(min(df['lnPolity2']), max(df['lnPolity2']),100))
	y = b[0] + b[1]*np.log(x+MINVAL)+ b[2]*np.log(x+MINVAL)**2 
	ax.plot(x,y, color='darkred')

	fig.savefig('results/Attitudinal.jpg')
	fig.show()

	pass



def reg(df):
	# Ensure 'Polity2' is in df, if not, merge it as done previously

	# Set up the variables
	X = df[['lnPolity2', 'lnPolity2sq']]  # Independent variables
	y = df['Attitudinal']  # Dependent variable

	# Add a constant (intercept) to the model
	X = sm.add_constant(X)

	# Fit the OLS model
	model = sm.OLS(y, X).fit()

	# Print the regression summary
	print("Standard OLS Regression Summary:")
	print(model.summary())

	# Perform Breusch-Pagan test for heteroskedasticity
	bp_test = het_breuschpagan(model.resid, model.model.exog)
	print("\nBreusch-Pagan Test:")
	print(f"LM Statistic: {bp_test[0]}, p-value: {bp_test[1]}")

	# Perform White's test for heteroskedasticity
	white_test = het_white(model.resid, model.model.exog)
	print("\nWhite's Test:")
	print(f"LM Statistic: {white_test[0]}, p-value: {white_test[1]}")

	# Re-run the regression with robust standard errors (Huber-White)
	robust_model = model.get_robustcov_results(cov_type='HC0')

	# Print the regression summary with robust standard errors
	print("\nRegression with Robust Standard Errors:")
	print(robust_model.summary())
	return model

def ttest(df):
	

	# Group 1: where Polity2 = 0
	group_1 = df[df['Polity2'] == 0]['Attitudinal']

	# Group 2: where 0 < Polity2 < 7
	group_2 = df[(df['Polity2'] > 0) & (df['Polity2'] < 7)]['Attitudinal']

	# Perform an independent t-test
	t_stat, p_value = ttest_ind(group_1, group_2, equal_var=False)  # Use Welch's t-test (unequal variance)
	print("Group 1: df['Polity2'] == 0")
	print("Group 2: df['Polity2'] > 0) & (df['Polity2'] < 7")
	print("t-statistic Group 1 -  Group 2:", t_stat)
	print("p-value:", p_value)

	# Interpretation
	if p_value < 0.05:
		print("Group 1 has significantly higher 'Attitudinal'")
	else:
		print("Group 1 does not have significantly higher 'Attitudinal'")

def reformat():
	a = ''
	heading = ''
	with open('LEGITIMACY_GILEY.txt', 'r') as f:
		s = f.read()
	for l in s.split('\n'):
		r =l.split(' ')
		if isfloat(r[-1]):
			a += ' '.join(r[isfloat(r[0]):-3])+'\t' + '\t'.join(r[-3:]) + heading + '\n'
		else:
			heading = l
	
	with open('LEGITIMACY_GILEY_format.txt', 'w') as f:
		f.write(a)

def fuzzymatch(df, df_polity_avg):

	# Step 2: Use fuzzy matching to map country names in df to those in df_polity_avg
	def get_closest_match(country_name, country_list):
		match, score = process.extractOne(country_name, country_list)
		return match if score > 80 else None  # Adjust score threshold as needed

	# Create a mapping of closest matches
	country_mapping = {country: get_closest_match(country, df_polity_avg['country'].tolist()) for country in df['Country'].unique()}

	# Map the country names in df to the names in df_polity_avg
	df['MappedCountry'] = df['Country'].map(country_mapping)

	# Step 3: Merge on the mapped country names
	final_df = df.merge(df_polity_avg, left_on='MappedCountry', right_on='country', how='left').drop(columns=['MappedCountry', 'country'])

	return final_df

def isfloat(x):
	try:
		float(x)
		return True
	except:
		return False

main()
