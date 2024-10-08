from os import RWF_SYNC
import pickle
import pandas as pd

print('test')

with open('ames_qsar_model.pkl', 'rb') as f:
  rfc = pickle.load(f)

print(rfc)

# for key in rfc:
#   print('key:' + key + ' => ' + rfc[key])