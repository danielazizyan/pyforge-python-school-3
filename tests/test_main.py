import pytest
from httpx import AsyncClient, ASGITransport
from src.main import app
from src.molecules.dao import MoleculeDAO
from sqlalchemy.exc import IntegrityError


@pytest.mark.asyncio(loop_scope="module")
async def test_get_server_id():
    async with AsyncClient(
        transport=ASGITransport(app=app),
        base_url="http://testserver"
    ) as ac:
        response = await ac.get("/")
    assert response.status_code == 200
    data = response.json()
    assert "server_id" in data


@pytest.mark.asyncio(loop_scope="module")
async def test_get_all_molecules():
    molecules_data = [
        {"name": "Molecule1", "smiles": "C=C"},
        {"name": "Molecule2", "smiles": "C#C"}
    ]

    # Clean up and add new molecules
    for molecule in molecules_data:
        await MoleculeDAO.delete(smiles=molecule["smiles"])
        await MoleculeDAO.add_molecule(**molecule)

    async with AsyncClient(
        transport=ASGITransport(app=app),
        base_url="http://testserver"
    ) as ac:
        response = await ac.get("/molecules/")
        assert response.status_code == 200
        data = response.json()
        assert isinstance(data, list)
        assert len(data) == 2, "There should be 2 molecules in the database."

    # Clean up after test
    for molecule in molecules_data:
        await MoleculeDAO.delete(smiles=molecule["smiles"])


@pytest.mark.asyncio(loop_scope="module")
async def test_get_all_molecules_with_limit():
    molecules_data = [
        {"name": "Molecule1", "smiles": "C=C"},
        {"name": "Molecule2", "smiles": "C#C"},
        {"name": "Molecule3", "smiles": "C#N"}
    ]

    # Clean up and add new molecules
    for molecule in molecules_data:
        await MoleculeDAO.delete(smiles=molecule["smiles"])
        await MoleculeDAO.add_molecule(**molecule)

    async with AsyncClient(
        transport=ASGITransport(app=app),
        base_url="http://testserver"
    ) as ac:
        response = await ac.get("/molecules/", params={"limit": 2})
        assert response.status_code == 200
        data = response.json()
        assert isinstance(data, list)
        assert len(data) == 2, "There should be 2 molecules returned due to the limit."

    # Clean up after test
    for molecule in molecules_data:
        await MoleculeDAO.delete(smiles=molecule["smiles"])


@pytest.mark.asyncio(loop_scope="module")
async def test_add_molecule():
    molecule_data = {"name": "TestMolecule", "smiles": "C=C"}

    # Delete any existing molecule with the same SMILES
    await MoleculeDAO.delete(smiles=molecule_data["smiles"])

    try:
        async with AsyncClient(
            transport=ASGITransport(app=app),
            base_url="http://testserver"
        ) as ac:
            response = await ac.post("/molecules/add/", json=molecule_data)
            assert response.status_code == 200
            data = response.json()
            assert data["molecule"]["name"] == molecule_data["name"]
            assert data["molecule"]["smiles"] == molecule_data["smiles"]

        # Verify that the molecule is added
        molecule = await MoleculeDAO.find_one_or_none(name="TestMolecule")
        assert molecule is not None
        assert molecule.name == "TestMolecule"
        assert molecule.smiles == "C=C"

    except IntegrityError as e:
        pytest.fail(f"Unique constraint violated: {e}")

    finally:
        # Clean up the molecule after the test
        await MoleculeDAO.delete(name="TestMolecule")


@pytest.mark.asyncio(loop_scope="module")
async def test_get_molecule_by_id():
    molecule_data = {"name": "TestMoleculeByID", "smiles": "C#C"}

    # Ensure molecule is not already in the DB
    await MoleculeDAO.delete(smiles=molecule_data["smiles"])

    async with AsyncClient(
        transport=ASGITransport(app=app),
        base_url="http://testserver"
    ) as ac:
        # Add the molecule
        response = await ac.post("/molecules/add/", json=molecule_data)
        assert response.status_code == 200

        # Retrieve the molecule by name using DAO
        molecule = await MoleculeDAO.find_one_or_none(name=molecule_data["name"])
        assert molecule is not None, "Molecule should be present in the database"
        mol_id = molecule.mol_id

        # Fetch the molecule by its ID
        response = await ac.get(f"/molecules/{mol_id}")
        assert response.status_code == 200
        retrieved_molecule = response.json()
        assert retrieved_molecule["mol_id"] == mol_id
        assert retrieved_molecule["name"] == molecule_data["name"]
        assert retrieved_molecule["smiles"] == molecule_data["smiles"]

    # Clean up the molecule after the test
    await MoleculeDAO.delete(name=molecule_data["name"])


@pytest.mark.asyncio(loop_scope="module")
async def test_update_molecule_by_id():
    molecule_data = {"name": "MoleculeToUpdate", "smiles": "C#C"}
    updated_data = {"name": "UpdatedMolecule", "smiles": "C=C"}

    # Ensure molecule is not already in the DB
    await MoleculeDAO.delete(smiles=molecule_data["smiles"])

    async with AsyncClient(
        transport=ASGITransport(app=app),
        base_url="http://testserver"
    ) as ac:
        # Add the molecule
        response = await ac.post("/molecules/add/", json=molecule_data)
        assert response.status_code == 200

        # Retrieve the molecule by name using DAO
        molecule = await MoleculeDAO.find_one_or_none(name=molecule_data["name"])
        assert molecule is not None, "Molecule should be present in the database"
        mol_id = molecule.mol_id

        # Update the molecule
        response = await ac.put(f"/molecules/{mol_id}", json=updated_data)
        assert response.status_code == 200
        updated_molecule_response = response.json()
        expected_message = f"The molecule with id {mol_id} has been updated!"
        assert updated_molecule_response["message"] == expected_message

        # Verify the molecule was updated using DAO
        updated_molecule = await MoleculeDAO.find_one_or_none(mol_id=mol_id)
        assert updated_molecule is not None
        assert updated_molecule.name == updated_data["name"]
        assert updated_molecule.smiles == updated_data["smiles"]

    # Clean up the molecule after the test
    await MoleculeDAO.delete(name=updated_data["name"])


@pytest.mark.asyncio(loop_scope="module")
async def test_delete_molecule_by_id():
    molecule_data = {"name": "MoleculeToDelete", "smiles": "C#N"}

    # Ensure molecule is not already in the DB
    await MoleculeDAO.delete(smiles=molecule_data["smiles"])

    async with AsyncClient(
        transport=ASGITransport(app=app),
        base_url="http://testserver"
    ) as ac:
        # Add the molecule
        response = await ac.post("/molecules/add/", json=molecule_data)
        assert response.status_code == 200

        # Fetch the molecule by name using DAO
        molecule = await MoleculeDAO.find_one_or_none(name=molecule_data["name"])
        assert molecule is not None, "Molecule should be present in the database"
        mol_id = molecule.mol_id

        # Delete the molecule by its ID
        response = await ac.delete(f"/molecules/delete/{mol_id}")
        assert response.status_code == 200
        delete_response = response.json()
        assert delete_response["message"] == f"The molecule with id {mol_id} is deleted!"

        # Verify the molecule was deleted using DAO
        molecule = await MoleculeDAO.find_one_or_none(mol_id=mol_id)
        assert molecule is None, "Molecule should be deleted from the database"


@pytest.mark.asyncio(loop_scope="module")
async def test_search_molecule():
    molecule_data = {"name": "MoleculeToSearch", "smiles": "C=C"}

    # Clean up and add the molecule
    await MoleculeDAO.delete(smiles=molecule_data["smiles"])
    await MoleculeDAO.add_molecule(**molecule_data)

    async with AsyncClient(
        transport=ASGITransport(app=app),
        base_url="http://testserver"
    ) as ac:
        search_params = {"substructure_smiles": "C"}
        response = await ac.get("/molecules/search", params=search_params)
        assert response.status_code == 200
        data = response.json()
        assert len(data) > 0, "Molecules should be found with the substructure 'C'."

    # Clean up after test
    await MoleculeDAO.delete(smiles=molecule_data["smiles"])
