import pandas as pd


#df = pd.read_csv('/Users/jrobertson/Desktop/cgmlst-profiles.csv',header=0, index_col=0)
#df.to_hdf('/Users/jrobertson/Desktop/cgmlst-profiles.hdf',key='cgmlst')

df = pd.read_hdf('/Users/jrobertson/Desktop/references/update/data/cgmlst/cgmlst-profiles.hdf',key='cgmlst')
df.to_csv('/Users/jrobertson/Desktop/cgmlst-profiles.csv')