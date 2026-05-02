import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import py3Dmol
from stmol import showmol
import base64
from io import BytesIO

def detect_interactions(mol):
    interactions = []
    
    # Check for Oxygen (8) or Nitrogen (7) atoms
    has_h_bond_acceptor_donor = any(atom.GetAtomicNum() in (7, 8) for atom in mol.GetAtoms())
    if has_h_bond_acceptor_donor:
        interactions.append("💧 Hydrogen Bond Possible")
        
    # Check for Carbon chains or non-polar groups (heuristic: has at least one Carbon atom)
    has_c = any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    if has_c:
        interactions.append("🌿 Hydrophobic Interaction Possible")
        
    # Optional enhancement: check for halogens for halogen bonding
    has_halogen = any(atom.GetAtomicNum() in (9, 17, 35, 53) for atom in mol.GetAtoms())
    if has_halogen:
        interactions.append("⚡ Halogen Bonding Possible")
        
    # Check for aromatic rings for Pi-Pi stacking
    has_aromatic = any(atom.GetIsAromatic() for atom in mol.GetAtoms())
    if has_aromatic:
        interactions.append("💍 Pi-Pi Stacking Possible")

    return interactions

def generate_3d_view(mol):
    mol_3d = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol_3d, randomSeed=42)
    try:
        AllChem.MMFFOptimizeMolecule(mol_3d)
    except Exception:
        pass
        
    mblock = Chem.MolToMolBlock(mol_3d)
    
    view = py3Dmol.view(width=400, height=400)
    view.addModel(mblock, 'mol')
    view.setStyle({'stick': {}, 'sphere': {'radius': 0.3}})
    # Clean light off-white inner background for high contrast
    view.setBackgroundColor('#FFFFE3')
    view.zoomTo()
    return view

st.set_page_config(page_title="Drug-Target Interaction Visualizer", layout="wide", page_icon="🧬")

def inject_custom_css():
    st.markdown("""
    <style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap');

    :root {
        --bg-main: #4A4A4A;
        --glass-bg: rgba(255, 255, 255, 0.05);
        --glass-border: rgba(255, 255, 255, 0.15);
        --glass-shadow: 0 8px 32px 0 rgba(0, 0, 0, 0.37);
        --text-main: #FFFFE3;
        --accent-1: #6D8196;
        --accent-2: #CBCBCB;
    }

    html, body, [class*="css"], .stApp {
        font-family: 'Inter', sans-serif !important;
        background-color: transparent !important;
        color: var(--text-main) !important;
    }
    
    .stApp > header {
        background-color: transparent !important;
    }
    
    .stApp {
        background-color: var(--bg-main) !important;
    }

    /* Ambient Background Shapes for Glassmorphism */
    .bg-shape {
        position: fixed;
        border-radius: 50%;
        filter: blur(100px);
        z-index: 0;
        pointer-events: none;
    }
    .shape1 {
        width: 40vw;
        height: 40vw;
        background: rgba(109, 129, 150, 0.8);
        top: -10%;
        left: -10%;
    }
    .shape2 {
        width: 45vw;
        height: 45vw;
        background: rgba(203, 203, 203, 0.5);
        bottom: -20%;
        right: -10%;
    }
    .shape3 {
        width: 30vw;
        height: 30vw;
        background: rgba(255, 255, 227, 0.4);
        top: 20%;
        left: 40%;
    }

    [data-testid="stHeader"] {
        background-color: transparent !important;
        display: none !important;
    }

    /* Main UI Container */
    .block-container {
        background: rgba(74, 74, 74, 0.4) !important;
        backdrop-filter: blur(30px) !important;
        -webkit-backdrop-filter: blur(30px) !important;
        border: 1px solid var(--glass-border) !important;
        border-radius: 24px !important;
        box-shadow: var(--glass-shadow) !important;
        max-width: 1100px !important;
        padding: 3rem 4rem !important;
        margin-top: 2rem !important;
        margin-bottom: 2rem !important;
        z-index: 1;
        position: relative;
    }

    /* Typography */
    h1 {
        color: var(--accent-1) !important;
        font-size: 38px !important;
        font-weight: 700 !important;
        margin-bottom: 0.5rem !important;
        display: flex;
        align-items: center;
        gap: 12px;
        text-shadow: 0 2px 10px rgba(109, 129, 150, 0.3);
    }
    h2 {
        color: var(--accent-1) !important;
        font-size: 26px !important;
        font-weight: 600 !important;
        margin-bottom: 1rem !important;
    }
    h3 {
        color: var(--accent-1) !important;
        font-size: 22px !important;
        font-weight: 600 !important;
        margin-bottom: 1rem !important;
    }
    p, label, span {
        color: var(--text-main) !important;
    }

    /* Inputs & Dropdowns */
    .stTextInput > div > div > input, 
    div[data-baseweb="select"] > div {
        background: rgba(255, 255, 255, 0.05) !important;
        backdrop-filter: blur(20px) !important;
        border: 1px solid var(--glass-border) !important;
        border-radius: 12px !important;
        color: var(--text-main) !important;
        box-shadow: inset 0 2px 4px rgba(0,0,0,0.2) !important;
        transition: all 0.3s ease !important;
    }
    
    .stTextInput > div > div > input {
        padding: 12px !important;
    }

    .stTextInput > div > div > input:focus, 
    div[data-baseweb="select"] > div:focus-within {
        border-color: var(--accent-2) !important;
        box-shadow: 0 0 0 2px rgba(203, 203, 203, 0.3) !important;
        background: rgba(255, 255, 255, 0.1) !important;
    }

    /* Dropdown text color specifically */
    div[data-baseweb="select"], 
    div[data-baseweb="select"] * {
        color: var(--text-main) !important;
    }

    div[data-baseweb="popover"] > div {
        background: rgba(74, 74, 74, 0.85) !important;
        backdrop-filter: blur(25px) !important;
        border: 1px solid var(--glass-border) !important;
        border-radius: 12px !important;
        box-shadow: var(--glass-shadow) !important;
    }
    li[role="option"] {
        background: transparent !important;
        color: var(--text-main) !important;
    }
    li[role="option"]:hover, li[aria-selected="true"] {
        background: rgba(109, 129, 150, 0.8) !important;
    }

    /* Primary Button */
    .stButton > button {
        background: rgba(255, 255, 255, 0.1) !important;
        backdrop-filter: blur(10px) !important;
        border: 1px solid rgba(255,255,255,0.3) !important;
        border-radius: 12px !important;
        color: var(--accent-1) !important;
        font-weight: 600 !important;
        padding: 12px 24px !important;
        transition: all 0.3s ease !important;
        width: 100% !important;
        box-shadow: 0 4px 12px rgba(0,0,0,0.2) !important;
    }
    .stButton > button * {
        color: var(--accent-1) !important;
    }
    .stButton > button:hover {
        background: rgba(255, 255, 255, 0.2) !important;
        border: 1px solid var(--accent-1) !important;
        box-shadow: 0 6px 16px rgba(0,0,0,0.3), 0 0 15px rgba(109, 129, 150, 0.4) !important;
        color: var(--bg-main) !important;
    }
    .stButton > button:hover * {
        color: var(--bg-main) !important;
    }
    .stButton > button:active {
        transform: translateY(2px) !important;
    }
    
    /* Columns as Glass Sub-Panels */
    [data-testid="column"] {
        background: rgba(255, 255, 255, 0.03) !important;
        backdrop-filter: blur(20px) !important;
        -webkit-backdrop-filter: blur(20px) !important;
        border: 1px solid var(--glass-border) !important;
        border-radius: 16px !important;
        box-shadow: var(--glass-shadow) !important;
        padding: 24px !important;
    }
    
    /* Structure Display Panels (2D and 3D) */
    .square-panel, iframe {
        background: #FFFFE3 !important;
        border-radius: 24px !important;
        padding: 16px !important;
        display: flex !important;
        justify-content: center !important;
        align-items: center !important;
        width: 100% !important;
        height: 400px !important;
        box-sizing: border-box !important;
        box-shadow: 0 8px 24px rgba(0, 0, 0, 0.2), inset 0 2px 8px rgba(255, 255, 255, 0.3) !important;
        border: 1px solid rgba(255, 255, 255, 0.3) !important;
        margin-bottom: 12px !important;
        transition: transform 0.3s ease, box-shadow 0.3s ease !important;
        overflow: hidden !important;
    }
    
    .square-panel:hover, iframe:hover {
        transform: scale(1.02) !important;
        box-shadow: 0 12px 32px rgba(0, 0, 0, 0.3), inset 0 2px 8px rgba(255, 255, 255, 0.4) !important;
        z-index: 10 !important;
        position: relative !important;
    }

    .square-panel img {
        max-width: 100%;
        max-height: 100%;
        object-fit: contain;
        border-radius: 8px;
        mix-blend-mode: multiply;
    }
    
    /* Analysis Cards */
    .glass-card {
        background: rgba(255, 255, 255, 0.05);
        backdrop-filter: blur(25px);
        -webkit-backdrop-filter: blur(25px);
        border: 1px solid var(--glass-border);
        border-left: 5px solid var(--accent-2);
        border-radius: 12px;
        padding: 20px 24px;
        margin-bottom: 16px;
        color: var(--text-main);
        font-size: 16px;
        font-weight: 500;
        box-shadow: var(--glass-shadow);
        transition: transform 0.3s ease, box-shadow 0.3s ease;
    }
    .glass-card:hover {
        transform: translateY(-2px);
        box-shadow: 0 12px 40px 0 rgba(0, 0, 0, 0.4);
        background: rgba(255, 255, 255, 0.08);
    }
    
    hr {
        border-top: 1px solid var(--glass-border) !important;
        margin: 3rem 0 !important;
    }

    .stMarkdown {
        margin-bottom: 0 !important;
    }
    </style>
    """, unsafe_allow_html=True)

    # Inject background shapes
    st.markdown("""
    <div class="bg-shape shape1"></div>
    <div class="bg-shape shape2"></div>
    <div class="bg-shape shape3"></div>
    """, unsafe_allow_html=True)

def create_interaction_card(text):
    return f"""
    <div class="glass-card">
        {text}
    </div>
    """

inject_custom_css()

# 2. Header Section
st.markdown("<h1><span style='font-size: 32px;'>🧬</span> Drug-Target Interaction Visualizer</h1>", unsafe_allow_html=True)
st.markdown("<p style='margin-bottom: 3rem; font-size: 16px; line-height: 1.6;'>This application visualizes drug molecules from their SMILES representation and predicts basic drug-target interactions. Enter a SMILES string below to view its 2D and 3D structures, along with potential binding interactions based on its chemical properties.</p>", unsafe_allow_html=True)

# 3. Input Controls
st.markdown("<label style='font-size: 14px; font-weight: 500; margin-bottom: 8px; display: block;'>Select an example molecule or choose 'Custom...' to enter your own</label>", unsafe_allow_html=True)

examples = {
    "Aspirin (Pain relief)": "CC(=O)Oc1ccccc1C(=O)O",
    "Caffeine (Stimulant)": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "Ibuprofen (Anti-inflammatory)": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    "Paracetamol (Pain relief)": "CC(=O)Nc1ccc(O)cc1",
    "Dopamine (Neurotransmitter)": "C1=CC(=C(C=C1CCN)O)O",
    "Custom...": ""
}

selected_example = st.selectbox("Example selection", list(examples.keys()), label_visibility="collapsed")

st.markdown("<label style='font-size: 14px; font-weight: 500; margin-top: 20px; margin-bottom: 8px; display: block;'>Enter SMILES string</label>", unsafe_allow_html=True)
if selected_example == "Custom...":
    smiles_input = st.text_input("SMILES input", placeholder="e.g., CCO", label_visibility="collapsed")
else:
    smiles_input = st.text_input("SMILES input", value=examples[selected_example], label_visibility="collapsed")

st.markdown("<div style='margin-top: 24px;'></div>", unsafe_allow_html=True)

# Primary Action Button
if st.button("Generate Results", type="primary"):
    if not smiles_input:
        st.warning("Please enter a SMILES string.")
    else:
        mol = Chem.MolFromSmiles(smiles_input)
        
        if mol is None:
            st.error("Invalid SMILES string. Please check your input and try again.")
        else:
            st.markdown("<hr>", unsafe_allow_html=True)
            
            # 4. Visualization Section
            col1, col2 = st.columns(2)
            
            with col1:
                st.markdown("<h2>2D Structure</h2>", unsafe_allow_html=True)
                
                img = Draw.MolToImage(mol, size=(400, 400))
                buffered = BytesIO()
                img.save(buffered, format="PNG")
                img_b64 = base64.b64encode(buffered.getvalue()).decode()
                
                st.markdown(f"<div class='square-panel'><img src='data:image/png;base64,{img_b64}' style='max-width: 100%; border-radius: 8px; mix-blend-mode: multiply;'></div>", unsafe_allow_html=True)
                
            with col2:
                st.markdown("<h2>3D Structure</h2>", unsafe_allow_html=True)
                
                view = generate_3d_view(mol)
                try:
                    showmol(view, height=400, width=400)
                except Exception as e:
                    st.error(f"Error generating 3D view: {e}")
                
                st.markdown("<p style='font-size: 14px; text-align: center; margin-top: 12px; font-weight: 500;'>Interactive 3D Molecule: Click and drag to rotate, scroll to zoom</p>", unsafe_allow_html=True)

            st.markdown("<hr>", unsafe_allow_html=True)
            
            # 5. Analysis Section
            st.markdown("<h2>Predicted Drug-Target Interactions</h2>", unsafe_allow_html=True)
            st.markdown("<p style='margin-bottom: 24px; font-size: 16px;'>Based on standard rule-based chemical properties:</p>", unsafe_allow_html=True)
            
            interactions = detect_interactions(mol)
            
            if interactions:
                for interaction in interactions:
                    st.markdown(create_interaction_card(interaction), unsafe_allow_html=True)
            else:
                st.markdown(create_interaction_card("ℹ️ No Interactions Found"), unsafe_allow_html=True)

