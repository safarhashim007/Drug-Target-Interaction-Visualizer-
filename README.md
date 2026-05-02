# Drug-Target Interaction Visualizer 🧬💊

A sleek, modern web application built with **Streamlit** and **RDKit** for visualizing drug molecules and predicting potential drug-target interactions. The UI features a premium, modern "glassmorphism" aesthetic with a responsive 2D and 3D molecular viewer.

## Features
- **SMILES Input:** Enter any valid SMILES string to visualize its molecular structure.
- **2D & 3D Visualization:** Seamlessly view molecules in both flat 2D structures and interactive 3D models using `py3Dmol` and `stmol`.
- **Chemical Analysis:** Extracts critical molecular properties (Molecular Weight, LogP, Hydrogen Bond Donors/Acceptors) using RDKit.
- **Lipinski's Rule of 5:** Automatically evaluates the drug-likeness of the input molecule.
- **Premium UI:** Designed with custom CSS utilizing modern styling, frosted glass effects, and a dark scientific theme.

## Local Development (Without Docker)

1. **Clone the repository:**
   ```bash
   git clone https://github.com/safarhashim007/chem.git
   cd chem
   ```

2. **Create a virtual environment:**
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows use `venv\Scripts\activate`
   ```

3. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```
   *Note: RDKit requires some underlying C++ system dependencies (like `libxrender1`). On macOS/Linux, these are usually handled, but you may need to install them if you encounter an `ImportError`.*

4. **Run the application:**
   ```bash
   streamlit run app.py
   ```
   The app will be available at `http://localhost:8501`.

## Docker Deployment 🐳

The easiest way to run the application locally without worrying about system dependencies is via Docker.

**Using Docker Compose (Recommended):**
```bash
docker compose up --build -d
```

**Using pure Docker:**
```bash
docker build -t chem-visualizer .
docker run -d -p 8501:8501 --name chem-app chem-visualizer
```
Your app will be live at `http://localhost:8501`.

## Cloud Deployment ☁️

This application is configured for seamless deployment on **Streamlit Community Cloud**. 

- The `packages.txt` file is included in the root directory specifically for Streamlit Cloud to automatically install the required Debian system dependencies (`libxrender1`, `libxext6`, etc.) for RDKit.
- Simply link your GitHub repository to Streamlit Community Cloud and set `app.py` as your main file.

## Technologies Used
- [Streamlit](https://streamlit.io/) - Web framework
- [RDKit](https://www.rdkit.org/) - Cheminformatics and machine learning software
- [py3Dmol](https://pypi.org/project/py3Dmol/) - 3D molecular visualization
- Custom HTML/Vanilla CSS for UI styling
