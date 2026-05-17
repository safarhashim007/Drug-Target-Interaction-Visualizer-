import os
import urllib.request
import urllib.parse
import json
import ssl

# Create an unverified SSL context to bypass corporate proxy/SSL errors
ctx = ssl._create_unverified_context()

from flask import Flask, request, jsonify, render_template
from rdkit import Chem  # type: ignore
from rdkit.Chem import Descriptors  # type: ignore
from rdkit.Chem import AllChem  # type: ignore
from rdkit.Chem import QED  # type: ignore
from rdkit.Chem.rdMolDescriptors import CalcMolFormula, CalcNumRings, CalcNumAromaticRings, CalcFractionCSP3  # type: ignore
from rdkit.Chem.GraphDescriptors import BertzCT  # type: ignore
from rdkit.Chem.Draw import rdMolDraw2D  # type: ignore

app = Flask(__name__)

def get_compound_name(smiles, mol):
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{urllib.parse.quote(smiles)}/property/Title/JSON"
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req, timeout=2, context=ctx) as response:
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
    # Mapping the previous detect_interactions to the new target format with detailed descriptions
    targets = []
    
    # Check for Oxygen (8) or Nitrogen (7) atoms
    has_h_bond_acceptor_donor = any(atom.GetAtomicNum() in (7, 8) for atom in mol.GetAtoms())
    if has_h_bond_acceptor_donor:
        targets.append({
            "name": "H-Bond Dependent Targets",
            "cls": "Kinases / Receptors",
            "aff": 0.65,
            "conf": "med",
            "col": "#378ADD",
            "desc": "Targets containing polar amino acid residues (like Ser, Thr, Asn) in their active sites that strongly coordinate with the molecule's Nitrogen or Oxygen atoms via hydrogen bonds.",
            "conf_desc": "Presence of generic H-bond donors/acceptors indicates binding potential, but exact 3D spatial alignment determines true selectivity."
        })
        
    # Check for Carbon chains or non-polar groups
    has_c = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if has_c >= 4:
        targets.append({
            "name": "Lipophilic Pockets",
            "cls": "GPCR",
            "aff": 0.55,
            "conf": "med",
            "col": "#1D9E75",
            "desc": "Deep hydrophobic cavities typically found in G-protein coupled receptors (GPCRs) that accommodate the carbon framework via favorable van der Waals interactions.",
            "conf_desc": "Sufficient carbon skeleton size confirms fit into lipophilic regions, though flexible chains may lead to multi-target promiscuity."
        })
        
    # Check for halogens for halogen bonding
    has_halogen = any(atom.GetAtomicNum() in (9, 17, 35, 53) for atom in mol.GetAtoms())
    if has_halogen:
        targets.append({
            "name": "Halogen-binding Sites",
            "cls": "Enzyme",
            "aff": 0.75,
            "conf": "high",
            "col": "#D85A30",
            "desc": "Specific enzymatic pockets where highly directional halogen bonds interact strongly with Lewis bases (e.g., protein backbone carbonyl oxygens).",
            "conf_desc": "Halogens form exceptionally specific and highly directional interactions, strongly anchoring the drug to target enzyme active sites."
        })
        
    # Check for aromatic rings for Pi-Pi stacking
    has_aromatic = any(atom.GetIsAromatic() for atom in mol.GetAtoms())
    if has_aromatic:
        targets.append({
            "name": "Aromatic Stacking Targets",
            "cls": "Transcription Factor",
            "aff": 0.45,
            "conf": "low",
            "col": "#7F77DD",
            "desc": "Planar binding interfaces facilitating parallel or T-shaped π-π stacking interactions with aromatic protein residues (Phe, Tyr, Trp).",
            "conf_desc": "Aromatic rings are extremely common structural motifs; true binding requires specific flanking residues to prevent steric clashes."
        })
        
    if not targets:
        targets.append({
            "name": "Low Specificity",
            "cls": "Various",
            "aff": 0.20,
            "conf": "low",
            "col": "#9c9a92",
            "desc": "Lacks highly distinct functional binding handles, leading to transient or non-specific surface interactions.",
            "conf_desc": "Absence of strong electrostatic anchors or defined lipophilic cores reduces selective binding probability."
        })

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

@app.route('/api/compound')
def get_compound():
    """Proxy PubChem compound lookup by name or SMILES — avoids browser CORS issues."""
    query = request.args.get('query', '').strip()
    
    if not query:
        return jsonify({'error': 'No query provided'}), 400

    try:
        encoded = urllib.parse.quote(query)
        
        # Heuristic: if it contains typical SMILES structural characters, try 'smiles' first.
        # Otherwise, try 'name' first.
        if any(c in "=#[]()@" for c in query):
            namespaces = ["smiles", "name"]
        else:
            namespaces = ["name", "smiles"]
            
        data = None
        used_namespace = ""
        
        for ns in namespaces:
            props_url = (
                f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{ns}/{encoded}"
                f"/property/Title,CanonicalSMILES,IsomericSMILES,MolecularWeight,XLogP,TPSA,"
                f"HBondDonorCount,HBondAcceptorCount,RotatableBondCount/JSON"
            )
            req = urllib.request.Request(props_url, headers={'User-Agent': 'DrugTargetVisualizer/2.0'})
            try:
                with urllib.request.urlopen(req, timeout=10, context=ctx) as resp:
                    data = json.loads(resp.read().decode())
                    used_namespace = ns
                    break
            except urllib.error.HTTPError as e:
                if e.code == 404:
                    continue # Try the next namespace
                raise # Propagate other errors (like 502, 500)
                
        if not data:
            return jsonify({'error': f'Compound "{query}" not found in PubChem'}), 404

        props = data['PropertyTable']['Properties'][0]
        smiles = props.get('CanonicalSMILES') or props.get('IsomericSMILES') or props.get('ConnectivitySMILES') or ''
        
        if not smiles and used_namespace == "smiles":
            smiles = query # Fallback if PubChem didn't return one explicitly
            
        if not smiles:
            return jsonify({'error': f'No SMILES found for "{query}"'}), 404

        # Also fetch CID using the successful namespace
        cid_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{used_namespace}/{encoded}/cids/JSON"
        req2 = urllib.request.Request(cid_url, headers={'User-Agent': 'DrugTargetVisualizer/2.0'})
        cid = None
        try:
            with urllib.request.urlopen(req2, timeout=10, context=ctx) as resp2:
                cid_data = json.loads(resp2.read().decode())
                cid = cid_data['IdentifierList']['CID'][0]
        except Exception:
            pass

        return jsonify({
            'name': props.get('Title', query),
            'smiles': smiles,
            'cid': cid,
            'mw': props.get('MolecularWeight'),
            'logp': props.get('XLogP'),
            'tpsa': props.get('TPSA'),
            'hbd': props.get('HBondDonorCount'),
            'hba': props.get('HBondAcceptorCount'),
            'rot': props.get('RotatableBondCount'),
        })

    except Exception as e:
        return jsonify({'error': f'PubChem lookup failed: {str(e)}'}), 502


@app.route('/api/sdf')
def get_sdf():
    """Proxy PubChem 3D SDF fetch by CID — avoids browser CORS issues."""
    cid = request.args.get('cid', '').strip()
    if not cid:
        return jsonify({'error': 'No CID provided'}), 400
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/record/SDF/?record_type=3d"
        req = urllib.request.Request(url, headers={'User-Agent': 'DrugTargetVisualizer/2.0'})
        with urllib.request.urlopen(req, timeout=15, context=ctx) as resp:
            sdf_text = resp.read().decode()
        from flask import Response
        return Response(sdf_text, mimetype='text/plain')
    except urllib.error.HTTPError as e:
        return jsonify({'error': f'3D SDF not available (HTTP {e.code})'}), 404
    except Exception as e:
        return jsonify({'error': f'SDF fetch failed: {str(e)}'}), 502


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
        aromatic_rings = CalcNumAromaticRings(mol)
        fsp3 = round(CalcFractionCSP3(mol), 3)
        complexity = round(BertzCT(mol), 1)
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
        'aromatic_rings': aromatic_rings,
        'fsp3': fsp3,
        'complexity': complexity,
        'targets': targets,
        'svg2d': svg2d,
        'mblock': mblock
    }
    return jsonify(result)

if __name__ == '__main__':
    app.run(debug=True, port=5001)
