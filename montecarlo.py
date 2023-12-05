
import markov
import pandas as pd
import numpy as np

def monte_carlo(df, iterations, slices, name):
	df_ = df.copy(deep = True)
	n = len(df)
	m = []
	for i in range(iterations):
		df_['ifhpol'] = np.array(df['ifhpol'])[np.argsort(np.random.rand(n))]
		markov_matrix, tr, t = markov.markov(df_, slices, name)
		m.append(markov_matrix)
	
	mean = np.mean(m, axis=0)
	std = np.std(m, axis=2)

	return mean, std

