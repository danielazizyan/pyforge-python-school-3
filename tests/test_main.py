import pytest
from fastapi.testclient import TestClient
from src.main import app, substructure_search, validate_smiles, smiles_db


# Fixture to reset smiles_db before each test
@pytest.fixture(autouse=True)
def reset_smiles_db():
    initial_state = smiles_db.copy()
    yield
    smiles_db.clear()
    smiles_db.extend(initial_state)


# Fixture for the FastAPI test client
@pytest.fixture
def client():
    return TestClient(app)


@pytest.mark.parametrize(
    "mols, substructure_smiles, expected",
    [
        (
            ["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"],
            "C",
            ["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"]
        ),
        (
            ["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"],
            "c1ccccc1",
            ["c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"]
        ),
        (
            ["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"],
            "O",
            ["CCO", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"]
        ),
        (["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"], "N", []),
        (["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"], "", ValueError),
        (["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"], "invalid", ValueError),
        ([], "C", [])
    ]
)
def test_substructure_search(mols, substructure_smiles, expected):
    """Test substructure search functionality."""
    if expected == ValueError:
        with pytest.raises(ValueError):
            substructure_search(mols, substructure_smiles)
    else:
        result = substructure_search(mols, substructure_smiles)
        assert result == expected


@pytest.mark.parametrize(
    "smiles, should_raise",
    [
        ("CCO", False),
        ("c1ccccc1", False),
        ("CC(=O)O", False),
        ("invalid_smiles", True),
        ("", True)
    ]
)
def test_validate_smiles(smiles, should_raise):
    """Test SMILES string validation."""
    if should_raise:
        with pytest.raises(ValueError):
            validate_smiles(smiles)
    else:
        validate_smiles(smiles)


def test_list_molecules(client):
    """Test listing all molecules."""
    response = client.get("/molecules")
    assert response.status_code == 200
    assert len(response.json()) == len(smiles_db)


def test_get_molecule(client):
    """Test retrieving a specific molecule by ID."""
    response = client.get("/molecules/1")
    assert response.status_code == 200
    assert response.json() == smiles_db[0]


def test_get_molecule_not_found(client):
    """Test retrieving a non-existent molecule by ID."""
    response = client.get("/molecules/999")
    assert response.status_code == 404
    assert response.json() == {"detail": "Molecule Not Found"}


def test_add_molecule(client):
    """Test adding a new molecule."""
    new_molecule = {"mol_id": 15, "name": "New Molecule", "smiles": "CCO"}
    response = client.post("/add", json=new_molecule)
    assert response.status_code == 201
    assert response.json() == new_molecule


def test_add_molecule_existing_id(client):
    """Test adding a molecule with an existing ID."""
    existing_molecule = {"mol_id": 1, "name": "Duplicate Ethanol", "smiles": "CCO"}
    response = client.post("/add", json=existing_molecule)
    assert response.status_code == 400
    assert response.json() == {"detail": "Molecule with this ID already exists."}


def test_add_molecule_invalid_smiles(client):
    """Test adding a molecule with invalid SMILES."""
    new_molecule = {"mol_id": 16, "name": "Invalid Molecule", "smiles": "invalid_smiles"}
    response = client.post("/add", json=new_molecule)
    assert response.status_code == 400
    assert response.json()["detail"] == "Invalid SMILES string: invalid_smiles"


def test_update_molecule(client):
    """Test updating an existing molecule."""
    updated_molecule = {"mol_id": 1, "name": "Updated Ethanol", "smiles": "CCO"}
    response = client.put("/molecules/1", json=updated_molecule)
    assert response.status_code == 200
    assert response.json() == updated_molecule


def test_update_molecule_not_found(client):
    """Test updating a non-existent molecule."""
    updated_molecule = {"mol_id": 999, "name": "Nonexistent Molecule", "smiles": "CCO"}
    response = client.put("/molecules/999", json=updated_molecule)
    assert response.status_code == 404
    assert response.json() == {"detail": "Molecule Not Found"}


def test_delete_molecule(client):
    """Test deleting an existing molecule."""
    response = client.delete("/molecules/1")
    assert response.status_code == 200
    assert response.json() == {"mol_id": 1, "name": "Ethanol", "smiles": "CCO"}
    response = client.get("/molecules/1")
    assert response.status_code == 404


def test_delete_molecule_not_found(client):
    """Test deleting a non-existent molecule."""
    response = client.delete("/molecules/999")
    assert response.status_code == 404
    assert response.json() == {"detail": "Molecule Not Found"}


def test_search_molecule(client):
    """Test searching molecules by substructure."""
    response = client.get("/search?substructure_smiles=c1ccccc1")
    assert response.status_code == 200
    assert len(response.json()) > 0
    expected_molecules = [
        {"mol_id": 2, "name": "Benzene", "smiles": "c1ccccc1"},
        {"mol_id": 4, "name": "Aspirin", "smiles": "CC(=O)Oc1ccccc1C(=O)O"},
        {"mol_id": 7, "name": "Toluene", "smiles": "Cc1ccccc1"},
        {"mol_id": 8, "name": "Phenol", "smiles": "c1ccc(cc1)O"},
        {"mol_id": 13, "name": "Styrene", "smiles": "c1ccccc1C=C"}
    ]
    assert response.json() == expected_molecules
