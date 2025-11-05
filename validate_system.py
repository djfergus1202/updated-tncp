"""
BioMed Research Suite - System Validation
==========================================
Comprehensive testing script to validate all components
before deployment.
"""

import sys
import json
import time

print("="*70)
print("BioMed Research Suite - System Validation")
print("="*70)

# Test 1: Import Dependencies
print("\n[1/6] Testing Python Dependencies...")
try:
    import flask
    import flask_cors
    import numpy as np
    import scipy
    from scipy.integrate import odeint
    from scipy.spatial import distance_matrix
    print("✓ All dependencies imported successfully")
    print(f"  - Flask: {flask.__version__}")
    print(f"  - NumPy: {np.__version__}")
    print(f"  - SciPy: {scipy.__version__}")
except ImportError as e:
    print(f"✗ Dependency import failed: {e}")
    print("  Run: pip install -r requirements.txt")
    sys.exit(1)

# Test 2: Import Backend Module
print("\n[2/6] Testing Backend Module...")
try:
    import unified_backend as backend
    print("✓ Backend module imported successfully")
    print(f"  - Modules: {', '.join(backend.app.blueprints.keys()) if backend.app.blueprints else 'Core Flask app'}")
except Exception as e:
    print(f"✗ Backend module import failed: {e}")
    sys.exit(1)

# Test 3: Validate Data Structures
print("\n[3/6] Validating Data Structures...")
try:
    # Check protein database
    assert len(backend.PROTEIN_DATABASE) > 0, "Protein database is empty"
    print(f"✓ Protein database: {len(backend.PROTEIN_DATABASE)} entries")
    
    # Check ligand database
    assert len(backend.LIGAND_DATABASE) > 0, "Ligand database is empty"
    print(f"✓ Ligand database: {len(backend.LIGAND_DATABASE)} entries")
    
    # Check cell line database
    assert len(backend.CELL_LINE_DATABASE) > 0, "Cell line database is empty"
    print(f"✓ Cell line database: {len(backend.CELL_LINE_DATABASE)} entries")
    
except AssertionError as e:
    print(f"✗ Data validation failed: {e}")
    sys.exit(1)

# Test 4: Test Molecular Docking Engine
print("\n[4/6] Testing Molecular Docking Engine...")
try:
    # Get first protein and ligand
    protein_id = list(backend.PROTEIN_DATABASE.keys())[0]
    ligand_id = list(backend.LIGAND_DATABASE.keys())[0]
    
    print(f"  - Testing: {protein_id} + {ligand_id}")
    
    # Run docking
    start_time = time.time()
    results = backend.MolecularDockingEngine.run_docking(protein_id, ligand_id, num_modes=3)
    elapsed = time.time() - start_time
    
    # Validate results
    assert 'modes' in results, "Results missing 'modes'"
    assert len(results['modes']) == 3, "Expected 3 modes"
    assert 'best_affinity' in results, "Results missing 'best_affinity'"
    
    print(f"✓ Docking engine works (completed in {elapsed:.2f}s)")
    print(f"  - Generated {len(results['modes'])} binding modes")
    print(f"  - Best affinity: {results['best_affinity']:.2f} kcal/mol")
    
except Exception as e:
    print(f"✗ Docking engine test failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Test 5: Test Cell Simulation
print("\n[5/6] Testing Cell Simulation Engine...")
try:
    # Create test parameters
    params = {
        'cellLineName': 'HeLa',
        'experimentParams': {
            'initialCells': 10,
            'duration': 5,
            'timeInterval': 0.5
        },
        'environment': {
            'temperature': 37,
            'co2': 5,
            'humidity': 95
        },
        'treatment': {
            'type': 'none',
            'concentration': 0
        }
    }
    
    print(f"  - Testing: {params['cellLineName']} culture")
    
    # Run simulation
    start_time = time.time()
    sim = backend.CellCultureSimulation(params)
    results = sim.run(params['experimentParams']['duration'], 
                      params['experimentParams']['timeInterval'])
    elapsed = time.time() - start_time
    
    # Validate results
    assert len(results) > 0, "No results generated"
    assert 'time' in results[0], "Results missing 'time'"
    assert 'viable' in results[0], "Results missing 'viable'"
    
    print(f"✓ Cell simulation works (completed in {elapsed:.2f}s)")
    print(f"  - Generated {len(results)} data points")
    print(f"  - Final cell count: {results[-1]['total']}")
    print(f"  - Final viability: {results[-1]['viability']:.1f}%")
    
except Exception as e:
    print(f"✗ Cell simulation test failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Test 6: Test Flask Routes
print("\n[6/6] Testing Flask API Routes...")
try:
    app = backend.app
    client = app.test_client()
    
    # Test health endpoint
    response = client.get('/api/health')
    assert response.status_code == 200, f"Health check failed: {response.status_code}"
    data = json.loads(response.data)
    assert data['status'] == 'healthy', "Health status not 'healthy'"
    print("✓ /api/health")
    
    # Test proteins endpoint
    response = client.get('/api/docking/proteins')
    assert response.status_code == 200, f"Proteins endpoint failed: {response.status_code}"
    data = json.loads(response.data)
    assert len(data) > 0, "No proteins returned"
    print("✓ /api/docking/proteins")
    
    # Test ligands endpoint
    response = client.get('/api/docking/ligands')
    assert response.status_code == 200, f"Ligands endpoint failed: {response.status_code}"
    data = json.loads(response.data)
    assert len(data) > 0, "No ligands returned"
    print("✓ /api/docking/ligands")
    
    # Test cell lines endpoint
    response = client.get('/api/cells/cell-lines')
    assert response.status_code == 200, f"Cell lines endpoint failed: {response.status_code}"
    data = json.loads(response.data)
    assert len(data) > 0, "No cell lines returned"
    print("✓ /api/cells/cell-lines")
    
    # Test docking simulation endpoint
    docking_params = {
        'proteinId': list(backend.PROTEIN_DATABASE.keys())[0],
        'ligandId': list(backend.LIGAND_DATABASE.keys())[0],
        'numModes': 3
    }
    response = client.post('/api/docking/run',
                          data=json.dumps(docking_params),
                          content_type='application/json')
    assert response.status_code == 200, f"Docking endpoint failed: {response.status_code}"
    data = json.loads(response.data)
    assert data['success'] == True, "Docking did not succeed"
    print("✓ /api/docking/run")
    
    # Test cell simulation endpoint
    cell_params = {
        'cellLineName': 'HeLa',
        'experimentParams': {
            'initialCells': 5,
            'duration': 2,
            'timeInterval': 0.5
        },
        'environment': {
            'temperature': 37,
            'co2': 5,
            'humidity': 95
        },
        'treatment': {
            'type': 'none',
            'concentration': 0
        }
    }
    response = client.post('/api/cells/simulate',
                          data=json.dumps(cell_params),
                          content_type='application/json')
    assert response.status_code == 200, f"Cell simulation endpoint failed: {response.status_code}"
    data = json.loads(response.data)
    assert data['success'] == True, "Cell simulation did not succeed"
    print("✓ /api/cells/simulate")
    
    print("\n✓ All API routes functional")
    
except Exception as e:
    print(f"✗ API route test failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Summary
print("\n" + "="*70)
print("VALIDATION COMPLETE - ALL TESTS PASSED ✓")
print("="*70)
print("\nSystem is ready for deployment!")
print("\nNext steps:")
print("1. Start server: python unified_backend.py")
print("2. Open: unified_biomedical_suite.html")
print("3. Or deploy to Render/Heroku using provided configs")
print("="*70)
