import pandas as pd

from scipy.optimize import minimize
import numpy as np
from matplotlib import pyplot as plt
import paneltime as pt
import statsmodels.api as sm
import os
os.chdir(os.path.dirname(__file__))
import networkdrawing as nw
import markov
import montecarlo as mt


def main():
    df = handle_data()
    #res = OLS(df, X)
    pres = panel_reg(df)
    plotting(df, pres)
    t0, t1  = mesh(df, 'ifhpol')
    markov_matrix, tr, t = markov.markov(df, [t0,t1], 'ifhpol')
    mean, std = mt.monte_carlo(df, 1000, [1/3,2/3], 'ifhpol')
    print(mean)
    print(std)
    print(pd.DataFrame(np.round(100*markov_matrix,1)))
    nw.draw_network(markov_matrix)

    markov.save_transitions(df)
    #res = OLS(df, X)



def plotting(df, res):
    df = df[['normalized_rank', 'change_pred']]
    df = df.sort_values('normalized_rank')
    df.plot('normalized_rank', 'change_pred')
    plt.show()
    x = np.linspace(0,1,100)
    a = res.ll.args.args_d
    plt.plot(x,a['beta'][1,0]+ a['beta'][2,0]*np.abs(x-0.5)**2)
    plt.title('prediction of normalized_rank and normalized_rank_sq only')
    plt.show()
    a=0
    
def handle_data():
    df = pd.read_csv("datasimplified.csv")
    df = get_annual_rank(df)

    df = pd.DataFrame(df)
    df['normalized_rank_sq'] = np.abs(df['normalized_rank']-0.5)**2
    df['demdummy'] = (df['normalized_rank']==1.0)*1.0
    df['semidemdummy'] = (df['normalized_rank']>0.5)*1.0
    return df
    
def panel_reg(df=None):
    pt.options.pqdkm.set([1,1,0,1,1])
    pt.options.fixed_random_time_eff.set(2)
    pt.options.fixed_random_group_eff.set(2)
    pt.options.tolerance.set(0.000001)
    
    res = pt.execute('abs_one_year_change~normalized_rank+normalized_rank_sq+semidemdummy+demdummy', df, 
               ID = 'country', T = 'year'
               )
    print(res)
    df['change_pred']  = res.predict()['in-sample predicted Y']
    return res
    
    
def OLS(df, X):
    X = sm.add_constant(df[['normalized_rank', 'normalized_rank_sq', 'semidemdummy', 'demdummy']])
    model = sm.OLS(df['abs_one_year_change'], X)
    results = model.fit()
    df['change_pred'] = results.predict()
    print(results.summary()) 
    return results

def get_annual_rank(df):
    df = pd.DataFrame(df)
    df = df.dropna()
    df = df.sort_values(['year', 'ifhpol'])
    # Group the data by 'year'
    grouped = df.groupby('year')
    # Add a new column 'rank' with the ranking for each identifier within each year
    df['rank'] = grouped['ifhpol'].rank(ascending=True)

    # Calculate the maximum rank for each year
    max_rank = grouped['rank'].transform('max')

    # Calculate the normalized rank for each identifier within each year
    df['normalized_rank'] = df['rank'] / max_rank
    df['max_rank'] = max_rank
    

    df = df.sort_values(['country', 'year'])

    grouped = df.groupby('country')

    # Calculate the one-year change in the 'variable' column
    df['rank_from'] = df['rank'].shift(1)
    df['ifhpol_from'] = df['ifhpol'].shift(1)
    df['totdurny2_from'] = df['totdurny2'].shift(1)
    df['abs_one_year_change'] = np.abs(grouped['normalized_rank'].diff())
    df['normalized_rank_from'] = df['rank_from'] / df['max_rank']

    # Drop the first row for each group as it will have a NaN value for the change
    df = df.dropna()
    # Print the resulting DataFrame
    return df


def mesh(df, name):
    d = {}
    x = 0.01
    mint = np.min(df[name])
    maxt1 =0.3
    mint2 = 0.9
    for t1 in np.arange(mint,maxt1,x):
        for t2 in np.arange(mint2,1.0,x):
            markov_matrix, tr, t = markov.markov(df, [t1,t2], name)
            obj = np.std(np.diag(markov_matrix), ddof=1)
            if not np.isnan(obj): 
                print(f"{obj}: {np.diag(markov_matrix)}; {(t1,t2)}")
                d[obj] =[t1, t2]

    t1, t2 = d[min(d)]
    markov_matrix, tr, t = markov.markov(df, [t1,t2], name)
    return t1, t2 

main()