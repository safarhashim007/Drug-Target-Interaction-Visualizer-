import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import py3Dmol
from stmol import showmol

def detect_interactions(mol):
    interactions = []
    
    # Check for Oxygen (8) or Nitrogen (7) atoms
    has_h_bond_acceptor_donor = any(atom.GetAtomicNum() in (7, 8) for atom in mol.GetAtoms())
    if has_h_bond_acceptor_donor:
        interactions.append("Hydrogen Bond Possible")
        
    # Check for Carbon chains or non-polar groups (heuristic: has at least one Carbon atom)
    has_c = any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    if has_c:
        interactions.append("Hydrophobic Interaction Possible")
        
    # Optional enhancement: check for halogens for halogen bonding
    has_halogen = any(atom.GetAtomicNum() in (9, 17, 35, 53) for atom in mol.GetAtoms())
    if has_halogen:
        interactions.append("Halogen Bonding Possible")
        
    # Check for aromatic rings for Pi-Pi stacking
    has_aromatic = any(atom.GetIsAromatic() for atom in mol.GetAtoms())
    if has_aromatic:
        interactions.append("Pi-Pi Stacking Possible")

    return interactions

def generate_3d_view(mol):
    mol_3d = Chem.AddHs(mol)
    # Generate 3D coordinates
    AllChem.EmbedMolecule(mol_3d, randomSeed=42)
    # Give it a quick geometry optimization
    try:
        AllChem.MMFFOptimizeMolecule(mol_3d)
    except Exception:
        pass # Better safe than sorry if optimization fails
        
    mblock = Chem.MolToMolBlock(mol_3d)
    
    view = py3Dmol.view(width=400, height=400)
    view.addModel(mblock, 'mol')
    view.setStyle({'stick': {}, 'sphere': {'radius': 0.3}})
    view.setBackgroundColor('#ffffff')
    view.zoomTo()
    return view

# --- UI Layout ---
st.set_page_config(page_title="Drug-Target Interaction Visualizer", layout="wide", page_icon="🧬")

def inject_custom_css():
    st.markdown("""
    <style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap');

    html, body, [class*="css"] {
        font-family: 'Inter', -apple-system, BlinkMacSystemFont, "SF Pro Display", "Segoe UI", Roboto, sans-serif !important;
        color: rgba(29, 29, 31, 0.85) !important;
    }

    /* Global App Background */
    .stApp, [data-testid="stAppViewContainer"] {
        background: linear-gradient(135deg, #e6e9f0 0%, #eef1f5 50%, #e0eafc 100%) !important;
        background-attachment: fixed !important;
    }
    
    /* Header/Top nav removal */
    [data-testid="stHeader"] {
        background-color: transparent !important;
        display: none !important;
    }

    /* Main Container */
    .block-container {
        max-width: 900px !important;
        padding: 2rem 3rem !important;
    }
    
    /* Typography */
    h1 {
        font-size: 32px !important;
        font-weight: 700 !important;
        color: #1d1d1f !important;
        letter-spacing: -0.5px !important;
        margin-bottom: 1rem !important;
    }
    h2, h3 {
        font-size: 20px !important;
        font-weight: 600 !important;
        color: #1d1d1f !important;
        margin-bottom: 1rem !important;
    }
    p {
        font-size: 16px !important;
        line-height: 1.6 !important;
    }

    /* 2D & 3D Columns as Glass Cards */
    [data-testid="column"] {
        background: rgba(255, 255, 255, 0.08) !important;
        backdrop-filter: blur(30px) !important;
        -webkit-backdrop-filter: blur(30px) !important;
        border-radius: 24px !important;
        border: 1px solid rgba(255, 255, 255, 0.4) !important;
        box-shadow: 0 8px 32px 0 rgba(0, 0, 0, 0.15) !important;
        padding: 24px !important;
        display: flex !important;
        flex-direction: column !important;
        align-items: center !important;
        margin-bottom: 24px !important;
    }
    
    /* Centering images in columns */
    [data-testid="column"] > div {
        width: 100% !important;
        display: flex !important;
        flex-direction: column !important;
        align-items: center !important;
    }

    /* Inputs */
    .stTextInput > div > div > input, 
    div[data-baseweb="select"] > div {
        background: rgba(255, 255, 255, 0.12) !important;
        border: 1px solid rgba(255, 255, 255, 0.5) !important;
        border-radius: 14px !important;
        box-shadow: inset 0 2px 4px rgba(0,0,0,0.02) !important;
        color: #1d1d1f !important;
        backdrop-filter: blur(20px) !important;
        transition: all 0.3s ease !important;
    }
    .stTextInput > div > div > input:focus, 
    div[data-baseweb="select"] > div:focus-within {
        border-color: #a8c0ff !important;
        box-shadow: 0 0 0 4px rgba(168, 192, 255, 0.3) !important;
    }

    /* Primary Button */
    .stButton > button {
        background: linear-gradient(135deg, #8ec5fc 0%, #e0c3fc 100%) !important;
        border: 1px solid rgba(255, 255, 255, 0.6) !important;
        border-radius: 16px !important;
        color: #1d1d1f !important;
        font-weight: 600 !important;
        font-size: 16px !important;
        padding: 10px 24px !important;
        box-shadow: 0 4px 16px rgba(142, 197, 252, 0.4) !important;
        transition: all 0.3s ease !important;
    }
    .stButton > button:hover {
        transform: scale(1.02) !important;
        box-shadow: 0 6px 20px rgba(142, 197, 252, 0.6) !important;
    }
    .stButton > button:active {
        transform: scale(0.98) !important;
    }

    /* Hide standard divider */
    hr {
        border-top: 1px solid rgba(255, 255, 255, 0.3) !important;
        margin: 2rem 0 !important;
    }
    
    /* Global Warning Alert (Missing SMILES) */
    [data-testid="stAlert"] {
        background: rgba(245, 158, 11, 0.1) !important; /* Muted amber */
        backdrop-filter: blur(20px) !important;
        border-radius: 16px !important;
        border: 1px solid rgba(245, 158, 11, 0.3) !important;
        color: #92400e !important;
        padding: 16px 20px !important;
    }
    </style>
    """, unsafe_allow_html=True)

def create_glass_card(icon, title, desc, bg_color):
    return f"""
    <div style="
        background: {bg_color};
        backdrop-filter: blur(20px);
        -webkit-backdrop-filter: blur(20px);
        border: 1px solid rgba(255, 255, 255, 0.4);
        border-radius: 16px;
        padding: 20px 24px;
        margin-bottom: 16px;
        box-shadow: 0 4px 16px rgba(0,0,0,0.05);
        transition: all 0.3s ease;
    " onmouseover="this.style.transform='translateY(-2px)'; this.style.boxShadow='0 8px 24px rgba(0,0,0,0.1)';" onmouseout="this.style.transform='translateY(0)'; this.style.boxShadow='0 4px 16px rgba(0,0,0,0.05)';">
        <div style="display: flex; align-items: center; margin-bottom: 8px;">
            <span style="font-size: 24px; margin-right: 12px;">{icon}</span>
            <h3 style="margin: 0; font-size: 18px; font-weight: 600; color: #1d1d1f;">{title}</h3>
        </div>
        <p style="margin: 0; font-size: 15px; color: rgba(29, 29, 31, 0.85); margin-left: 36px;">{desc}</p>
    </div>
    """

inject_custom_css()

st.title("🧬 Drug-Target Interaction Visualizer")
st.markdown("""
This application visualizes drug molecules from their **SMILES** representation and predicts basic drug-target interactions. 
Enter a SMILES string below to view its 2D and 3D structures, along with potential binding interactions based on its chemical properties.
""")

# Example molecules
examples = {
    "Aspirin (Pain relief)": "CC(=O)Oc1ccccc1C(=O)O",
    "Caffeine (Stimulant)": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "Ibuprofen (Anti-inflammatory)": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    "Paracetamol (Pain relief)": "CC(=O)Nc1ccc(O)cc1",
    "Dopamine (Neurotransmitter)": "C1=CC(=C(C=C1CCN)O)O",
    "Custom...": ""
}

selected_example = st.selectbox("Select an example molecule or choose 'Custom...' to enter your own:", list(examples.keys()))

if selected_example == "Custom...":
    smiles_input = st.text_input("Enter SMILES string:", placeholder="e.g., CCO")
else:
    smiles_input = st.text_input("Enter SMILES string:", value=examples[selected_example])

if st.button("Generate Results", type="primary"):
    if not smiles_input:
        st.warning("Please enter a SMILES string.")
    else:
        # Process the molecule
        mol = Chem.MolFromSmiles(smiles_input)
        
        if mol is None:
            st.error("Invalid SMILES string. Please check your input and try again.")
        else:
            col1, col2 = st.columns(2)
            
            with col1:
                st.subheader("2D Structure")
                img = Draw.MolToImage(mol, size=(400, 400))
                st.image(img, use_column_width=False, caption="2D Molecular Structure")
                
            with col2:
                st.subheader("3D Structure")
                # Generate 3D viewer
                view = generate_3d_view(mol)
                try:
                    showmol(view, height=400, width=400)
                    st.caption("Interactive 3D Molecular Structure (drag to rotate, scroll to zoom)")
                except Exception as e:
                    st.error(f"Error generating 3D view: {e}")

            st.divider()
            st.subheader("Predicted Drug-Target Interactions")
            st.markdown("Based on standard rule-based chemical properties:")
            
            interactions = detect_interactions(mol)
            
            if interactions:
                for interaction in interactions:
                    if "Hydrogen Bond" in interaction:
                        st.markdown(create_glass_card("💧", interaction, "Molecule contains Oxygen or Nitrogen atoms", "rgba(142, 197, 252, 0.1)"), unsafe_allow_html=True)
                    elif "Hydrophobic" in interaction:
                        st.markdown(create_glass_card("🌿", interaction, "Molecule contains Carbon chains", "rgba(167, 243, 208, 0.1)"), unsafe_allow_html=True)
                    elif "Halogen" in interaction:
                        st.markdown(create_glass_card("⚡", interaction, "Molecule contains Halogens like F, Cl, Br, I", "rgba(253, 230, 138, 0.1)"), unsafe_allow_html=True)
                    elif "Pi-Pi" in interaction:
                        st.markdown(create_glass_card("💍", interaction, "Molecule contains aromatic rings", "rgba(24bc, 165, 165, 0.1)"), unsafe_allow_html=True)
                    else:
                        st.markdown(create_glass_card("✨", interaction, "Based on basic chemical rules", "rgba(255, 255, 255, 0.1)"), unsafe_allow_html=True)
            else:
                st.markdown(create_glass_card("ℹ️", "No Interactions Found", "No specific interactions predicted based on the basic rules.", "rgba(255, 255, 255, 0.1)"), unsafe_allow_html=True)
