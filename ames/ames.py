import pickle
import time
from time import strftime

from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

#from rdkit.Chem import Draw   # Dependency problems within Replit re: libxrender, libexpat

with open('ames/ames_qsar_model.pkl', 'rb') as f:
  rfc = pickle.load(f)


def AmesPredictor(SMILES_string):
  output = {}

  # Populate metadata in output dictionary
  output['Input received'] = strftime('%d %b %Y %H:%M:%S (%Z%z)')
  output['Input SMILES'] = SMILES_string

  # Convert input SMILES string into RDKit Mol object
  mol = Chem.MolFromSmiles(SMILES_string)

  # Render the RDKit Mol as an image
  #img = Draw.MolToImage(mol, size=(100, 100))

  # Compute Morgan fingerprint
  mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
  mfp = mfpgen.GetFingerprintAsNumPy(mol).reshape(-1, 1024)

  # Make prediction with our trained RandomForestClassifier
  prediction = rfc.predict(mfp)

  # Convert machine prediction to human-readable prediction
  if prediction == 1:
    output['Ames prediction'] = 'Positive'
  if prediction == 0:
    output['Ames prediction'] = 'Negative'

  # Return output
  return output