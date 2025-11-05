# 🚀 Quick Start Guide - BioMed Research Suite

## What You Got

A unified computational biology platform combining:
- **Molecular Docking** - Protein-ligand binding simulations
- **Cell Dynamics** - Real-time cell culture modeling

## 📦 Files Included

1. `unified_biomedical_suite.html` - Main web interface
2. `unified_backend.py` - Python Flask backend
3. `requirements.txt` - Python dependencies
4. `render.yaml` - Render.com deployment config
5. `Procfile` - Heroku deployment config
6. `runtime.txt` - Python version specification
7. `README.md` - Comprehensive documentation
8. `validate_system.py` - System testing script

## ⚡ 3-Minute Local Setup

```bash
# 1. Install dependencies
pip install -r requirements.txt

# 2. Validate system (optional but recommended)
python validate_system.py

# 3. Start backend
python unified_backend.py

# 4. Open frontend
# Open unified_biomedical_suite.html in your browser
```

That's it! You're running a full computational biology lab.

## 🌐 Deploy to Cloud (5 minutes)

### Option A: Render.com (Recommended)

1. Create GitHub repo with all files
2. Go to [render.com](https://render.com)
3. Click "New +" → "Blueprint"
4. Connect your repo
5. Deploy! ✨

### Option B: Heroku

```bash
# 1. Create app
heroku create your-app-name

# 2. Deploy
git add .
git commit -m "Deploy BioMed Suite"
git push heroku main

# 3. Open
heroku open
```

## 🧪 Quick Test

### Test Molecular Docking
1. Open the web interface
2. Click "Molecular Docking" module
3. Select "SARS-CoV-2 Mpro" protein
4. Choose "Remdesivir" ligand
5. Click "Run Docking Simulation"
6. View binding affinity results!

### Test Cell Dynamics
1. Click "Cell Dynamics" module
2. Select "HeLa" cell line
3. Click "Start Cell Culture Simulation"
4. Watch cells grow in real-time!

## 📊 What Each Module Does

### Molecular Docking
- Simulates drug-protein binding
- Predicts binding affinity (kcal/mol)
- Shows interaction types (H-bonds, hydrophobic, etc.)
- Analyzes binding site druggability

### Cell Dynamics
- Simulates cell culture growth
- Real-time visualization (800x600 μm field)
- Tracks viability and health
- Models cell cycle (G1/S/G2/M)

## 🎯 API Endpoints

```
GET  /api/health                  # System status
GET  /api/docking/proteins        # Available proteins
GET  /api/docking/ligands         # Available ligands
POST /api/docking/run             # Run docking
GET  /api/cells/cell-lines        # Available cell lines
POST /api/cells/simulate          # Run cell sim
POST /api/predict/drug-efficacy   # Predict drug effects
```

## 💡 Pro Tips

1. **For Research**: Export results as JSON for analysis
2. **For Teaching**: Use as interactive demo in classes
3. **For Development**: Extend with your own proteins/ligands
4. **For Production**: Upgrade to paid hosting for better performance

## ❓ Troubleshooting

**Backend won't start?**
→ Run: `pip install -r requirements.txt`

**Frontend can't connect?**
→ Check backend is running on port 5000

**Slow simulations?**
→ Reduce initial cell count or duration

## 📚 Next Steps

1. Read full `README.md` for detailed docs
2. Explore API endpoints
3. Try different proteins and ligands
4. Customize cell line parameters
5. Deploy to production!

## 🎓 Learn More

- Molecular Docking: Study drug-target interactions
- Cell Dynamics: Model cancer vs. normal cells
- PK/PD: Predict drug efficacy
- ADMET: Drug-like properties

---

**Happy Researching! 🧬**

Questions? Check README.md for comprehensive documentation.
