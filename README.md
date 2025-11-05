# BioMed Research Suite v2.0

A comprehensive computational biology platform combining molecular docking and cell dynamics simulation in one unified interface.

## 🧬 Features

### Molecular Docking Module
- Protein-ligand binding affinity prediction
- Multiple binding mode analysis
- Interaction fingerprinting (H-bonds, hydrophobic, π-π stacking)
- Binding site analysis (volume, druggability, flexibility)
- Support for major drug targets (HIV-1 Protease, COX-2, SARS-CoV-2 Mpro, EGFR Kinase)
- Comprehensive ligand library (Aspirin, Ibuprofen, Remdesivir, etc.)

### Cell Dynamics Module
- Real-time cell culture simulation
- Multiple cell line support (HeLa, MCF-7, A549, HEK293)
- Live visualization with canvas rendering
- Growth curve analysis
- Viability tracking
- Cell cycle modeling (G1/S/G2/M phases)
- Drug response prediction
- Metabolic state monitoring (ATP, glucose, oxygen)

## 🚀 Quick Start

### Local Development

1. **Install Dependencies**
   ```bash
   pip install -r requirements.txt
   ```

2. **Start the Backend**
   ```bash
   python unified_backend.py
   ```
   Server will start at `http://127.0.0.1:5000`

3. **Open the Frontend**
   Open `unified_biomedical_suite.html` in your web browser

### Deploy to Render.com

1. **Create a GitHub Repository**
   - Upload all files to your repository
   - Files needed: `unified_backend.py`, `requirements.txt`, `runtime.txt`, `render.yaml`, `unified_biomedical_suite.html`

2. **Connect to Render**
   - Go to [render.com](https://render.com)
   - Click "New +" → "Blueprint"
   - Connect your GitHub repository
   - Render will automatically detect `render.yaml` and deploy

3. **Access Your App**
   - Backend API: `https://your-app-name.onrender.com/api/health`
   - Frontend: Open `unified_biomedical_suite.html` in a browser
   - Update the API URL in the HTML if needed

### Deploy to Heroku

1. **Create Heroku App**
   ```bash
   heroku create your-app-name
   ```

2. **Deploy**
   ```bash
   git add .
   git commit -m "Deploy BioMed Research Suite"
   git push heroku main
   ```

3. **Open App**
   ```bash
   heroku open
   ```

## 📚 API Documentation

### Health Check
```
GET /api/health
Response: {
  "status": "healthy",
  "version": "2.0",
  "modules": ["molecular_docking", "cell_dynamics"]
}
```

### Molecular Docking
```
POST /api/docking/run
Body: {
  "proteinId": "6LU7",
  "ligandId": "remdesivir",
  "numModes": 9
}
Response: {
  "success": true,
  "data": {
    "protein": {...},
    "ligand": {...},
    "modes": [...],
    "best_affinity": -8.5
  }
}
```

### Get Proteins
```
GET /api/docking/proteins
Response: {
  "1HVH": { "name": "HIV-1 Protease", ... },
  ...
}
```

### Get Ligands
```
GET /api/docking/ligands
Response: {
  "aspirin": { "molecular_weight": 180.16, ... },
  ...
}
```

### Cell Simulation
```
POST /api/cells/simulate
Body: {
  "cellLineName": "HeLa",
  "experimentParams": {
    "initialCells": 50,
    "duration": 72,
    "timeInterval": 0.5
  },
  "environment": {
    "temperature": 37,
    "co2": 5,
    "humidity": 95
  }
}
Response: {
  "success": true,
  "data": [
    {
      "time": 0,
      "total": 50,
      "viable": 50,
      "viability": 100
    },
    ...
  ]
}
```

### Get Cell Lines
```
GET /api/cells/cell-lines
Response: {
  "HeLa": {
    "type": "Cancer",
    "doubling_time": 24,
    ...
  },
  ...
}
```

### Drug Efficacy Prediction
```
POST /api/predict/drug-efficacy
Body: {
  "cellLineName": "HeLa",
  "drugClass": "taxol",
  "concentration": 10
}
Response: {
  "ic50": 8.5,
  "predicted_efficacy": 65.2,
  "predicted_viability": 34.8
}
```

## 🔬 Usage Examples

### Example 1: Run Molecular Docking

1. Select a protein target (e.g., SARS-CoV-2 Mpro)
2. Choose a ligand (e.g., Remdesivir)
3. Click "Run Docking Simulation"
4. View binding modes, interactions, and affinity scores

### Example 2: Simulate Cell Culture

1. Select a cell line (e.g., HeLa)
2. Click "Start Cell Culture Simulation"
3. Watch real-time visualization of cell growth
4. Analyze growth curves and viability metrics

### Example 3: Combine Both Modules

1. Run docking to identify best drug-target affinity
2. Switch to cell dynamics module
3. Simulate drug effects on target cell line
4. Compare predicted efficacy with simulation results

## 🛠️ Technology Stack

### Frontend
- React 18
- Recharts for data visualization
- HTML5 Canvas for cell rendering
- Responsive CSS design

### Backend
- Flask 3.0 (Python web framework)
- NumPy & SciPy (scientific computing)
- Gunicorn (production server)
- Flask-CORS (cross-origin support)

### Deployment
- Render.com (recommended)
- Heroku
- Any platform supporting Python + Flask

## 📊 Data Models

### Protein Profile
- PDB ID, name, organism
- Resolution, binding site volume
- Flexibility score, druggability

### Ligand Profile
- SMILES structure, molecular weight
- LogP, hydrogen bond donors/acceptors
- Rotatable bonds

### Cell Line Profile
- Cell type, origin, doubling time
- Cell cycle parameters (G1/S/G2/M durations)
- Metabolic rates (glucose, oxygen, lactate)
- Drug sensitivity (IC50 values)

## 🎯 Future Enhancements

- [ ] 3D molecular visualization
- [ ] Advanced ADMET predictions
- [ ] Molecular dynamics (MD) simulation
- [ ] Multi-drug combination analysis
- [ ] Gene expression heatmaps
- [ ] Export to PDF reports
- [ ] Database integration for result storage
- [ ] User authentication
- [ ] Batch processing capabilities
- [ ] Machine learning models for prediction

## 📝 File Structure

```
biomed-research-suite/
├── unified_biomedical_suite.html  # Main frontend interface
├── unified_backend.py             # Unified Flask backend
├── requirements.txt               # Python dependencies
├── runtime.txt                    # Python version
├── Procfile                       # Heroku deployment
├── render.yaml                    # Render.com deployment
└── README.md                      # This file
```

## ⚙️ Configuration

### Environment Variables
- `PORT`: Server port (default: 5000)
- `FLASK_ENV`: Environment (development/production)
- `RENDER`: Flag for Render.com deployment
- `DYNO`: Flag for Heroku deployment

### CORS Configuration
The backend accepts requests from:
- localhost:* (development)
- 127.0.0.1:* (development)
- *.onrender.com (Render deployment)
- *.herokuapp.com (Heroku deployment)
- file:// (local HTML files)

## 🐛 Troubleshooting

### Backend won't start
- Ensure Python 3.11 is installed
- Install dependencies: `pip install -r requirements.txt`
- Check port 5000 is available

### Frontend can't connect to backend
- Verify backend is running
- Check API_BASE_URL in HTML file
- Enable CORS if needed

### Simulation runs slowly
- Reduce initial cell count
- Decrease simulation duration
- Upgrade to paid hosting plan

### Deployment fails
- Check Python version matches runtime.txt
- Verify all dependencies in requirements.txt
- Review deployment logs

## 📄 License

This is educational software for research and learning purposes.

## 🤝 Contributing

Contributions welcome! Areas for improvement:
- Algorithm optimization
- UI/UX enhancements
- Additional cell lines
- More drug compounds
- Advanced visualizations

## 📧 Support

For issues or questions:
- Review API documentation
- Check troubleshooting section
- Consult deployment logs

---

**BioMed Research Suite v2.0** - Advancing Computational Biology Research

Built with ❤️ for the scientific community
