import pickle

import pandas as pd
import numpy as np

from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
#from rdkit.Chem import Draw

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, roc_auc_score

print('test')
# Get dataset
ames_dataset = pd.read_csv(
'smiles_cas.smi',
  names=['SMILES', 'CASRN', 'Ames'],
  delimiter='\t'
)

print(ames_dataset)

# Convert SMILES to RDKit Mols
mols = [ Chem.MolFromSmiles(smi) for smi in ames_dataset['SMILES'] ]

print('mols')

# Convert RDKit Mols to fingerprints
mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
mfps = [ mfpgen.GetFingerprintAsNumPy(mol) for mol in mols ]

print('conversion')

# Split dataset into training and test sets
X_train, X_test, y_train, y_test = train_test_split(mfps, ames_dataset['Ames'], test_size=0.2, random_state=42, stratify=ames_dataset['Ames'])

print('split')

# Train random forest
rfc = RandomForestClassifier(n_estimators=100, random_state=42)
rfc.fit(X_train, y_train)
y_train_pred = rfc.predict_proba(X_train)[:,1]
acc_train = accuracy_score(y_train, np.round(y_train_pred))

print('Training accuracy:', acc_train * 100)

# Test random forest
y_test_pred = rfc.predict_proba(X_test)[:,1]
acc_test = accuracy_score(y_test, np.round(y_test_pred))
rocauc_test = roc_auc_score(y_test, y_test_pred)

print('Test accuracy:', acc_test * 100)
print('Test ROCAUC:', rocauc_test * 100)

# Save model
with open('ames_qsar_model.pkl', 'wb') as f:
  pickle.dump(rfc, f)