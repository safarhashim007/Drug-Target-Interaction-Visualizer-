import os
import urllib.request
import urllib.parse
import json

from flask import Flask, request, jsonify, render_template
from rdkit import Chem  # type: ignore
from rdkit.Chem import Descriptors  # type: ignore
from rdkit.Chem import AllChem  # type: ignore
from rdkit.Chem import QED  # type: ignore
from rdkit.Chem.Draw import rdMolDraw2D  # type: ignore
from rdkit.Chem.rdMolDescriptors import CalcMolFormula, CalcNumRings  # type: ignore

app = Flask(__name__)

def get_compound_name(smiles, mol):
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{urllib.parse.quote(smiles)}/property/Title/JSON"
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req, timeout=2) as response:
            data = json.loads(response.read().decode())
            return data['PropertyTable']['Properties'][0]['Title']
    except Exception:
        pass
        
    try:
        formula = CalcMolFormula(mol)
        return f"Compound ({formula})"
    except Exception:
        return "Custom Compound"

def detect_interactions(mol):
    # Mapping the previous detect_interactions to the new target format
    targets = []
    
    # Check for Oxygen (8) or Nitrogen (7) atoms
    has_h_bond_acceptor_donor = any(atom.GetAtomicNum() in (7, 8) for atom in mol.GetAtoms())
    if has_h_bond_acceptor_donor:
        targets.append({ "name": "H-Bond Dependent Targets", "cls": "Kinases / Receptors", "aff": 0.65, "conf": "med", "col": "#378ADD" })
        
    # Check for Carbon chains or non-polar groups
    has_c = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if has_c >= 4:
        targets.append({ "name": "Lipophilic Pockets", "cls": "GPCR", "aff": 0.55, "conf": "med", "col": "#1D9E75" })
        
    # Check for halogens for halogen bonding
    has_halogen = any(atom.GetAtomicNum() in (9, 17, 35, 53) for atom in mol.GetAtoms())
    if has_halogen:
        targets.append({ "name": "Halogen-binding Sites", "cls": "Enzyme", "aff": 0.75, "conf": "high", "col": "#D85A30" })
        
    # Check for aromatic rings for Pi-Pi stacking
    has_aromatic = any(atom.GetIsAromatic() for atom in mol.GetAtoms())
    if has_aromatic:
        targets.append({ "name": "Aromatic Stacking Targets", "cls": "Transcription Factor", "aff": 0.45, "conf": "low", "col": "#7F77DD" })
        
    if not targets:
        targets.append({ "name": "Low Specificity", "cls": "Various", "aff": 0.20, "conf": "low", "col": "#9c9a92" })

    return targets

def get_2d_svg(mol):
    d2d = rdMolDraw2D.MolDraw2DSVG(320, 190)
    d2d.drawOptions().addStereoAnnotation = True
    d2d.drawOptions().clearBackground = False
    d2d.DrawMolecule(mol)
    d2d.FinishDrawing()
    return d2d.GetDrawingText()

def get_3d_mblock(mol):
    try:
        mol_3d = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol_3d, randomSeed=42)
        try:
            AllChem.MMFFOptimizeMolecule(mol_3d)
        except Exception:
            pass
        return Chem.MolToMolBlock(mol_3d)
    except Exception:
        return None

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/api/analyze', methods=['POST'])
def analyze():
    data = request.get_json()
    smiles = data.get('smiles', '').strip()
    if not smiles:
        return jsonify({'error': 'No SMILES provided'}), 400
        
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return jsonify({'error': 'Invalid SMILES string'}), 400
        
    try:
        mw = round(Descriptors.MolWt(mol), 2)
        logp = round(Descriptors.MolLogP(mol), 2)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        tpsa = round(Descriptors.TPSA(mol), 2)
        rot_bonds = Descriptors.NumRotatableBonds(mol)
        qed = round(QED.qed(mol), 3)
        rings = CalcNumRings(mol)
    except Exception as e:
        return jsonify({'error': f'Error calculating properties: {str(e)}'}), 500
    
    targets = detect_interactions(mol)
    svg2d = get_2d_svg(mol)
    mblock = get_3d_mblock(mol)
    compound_name = get_compound_name(smiles, mol)
    
    result = {
        'name': compound_name,
        'mw': mw,
        'logp': logp,
        'hbd': hbd,
        'hba': hba,
        'tpsa': tpsa,
        'rot_bonds': rot_bonds,
        'qed': qed,
        'rings': rings,
        'targets': targets,
        'svg2d': svg2d,
        'mblock': mblock
    }
    return jsonify(result)

if __name__ == '__main__':
    app.run(debug=True, port=5000)
