from rdkit import Chem
from rdkit.Chem import Draw
import base64
import io
import pandas as pd
import cirpy

def generate_image(smiles):
    m = Chem.MolFromSmiles(smiles)
    image = Draw.MolToImage(m)
    data = io.BytesIO()
    image.save(data, "JPEG")
    encoded_img_data = base64.b64encode(data.getvalue())
    return encoded_img_data


def csv_to_dict(file):
    data = pd.read_csv(file)
    data_dict = data.to_dict(orient='records')
    return data_dict

def cas_to_smiles(cas):
    smiles = cirpy.resolve(cas, 'smiles')
    return smiles

print(cas_to_smiles('108-95-2'))