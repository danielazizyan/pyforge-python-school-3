import asyncio
import pytest
from httpx import AsyncClient, ASGITransport
from src.main import app
from sqlalchemy import select, delete
from src.database import async_session_maker
from src.molecules.models import Molecule


@pytest.fixture(scope="module")
def event_loop():
    loop = asyncio.new_event_loop()
    yield loop
    loop.close()


async def add_and_verify_molecule(ac, molecule_data):
    """Helper function to add a molecule and verify its addition."""
    response = await ac.post("/molecules/add/", json=molecule_data)
    assert response.status_code == 200
    data = response.json()
    assert data["molecule"]["name"] == molecule_data["name"]
    assert data["molecule"]["smiles"] == molecule_data["smiles"]

    async with async_session_maker() as session:
        result = await session.execute(
            select(Molecule).filter_by(name=molecule_data["name"])
        )
        molecule = result.scalar_one_or_none()
        assert molecule is not None
        assert molecule.name == molecule_data["name"]
        assert molecule.smiles == molecule_data["smiles"]
        return molecule.mol_id


async def delete_molecule_by_id(mol_id):
    """Helper function to delete a molecule by its ID."""
    async with async_session_maker() as session:
        await session.execute(delete(Molecule).where(Molecule.mol_id == mol_id))
        await session.commit()


@pytest.mark.asyncio(loop_scope="module")
async def test_get_server_id():
    async with AsyncClient(
        transport=ASGITransport(app=app), base_url="http://testserver"
    ) as ac:
        response = await ac.get("/")
    assert response.status_code == 200
    data = response.json()
    assert "server_id" in data


@pytest.mark.asyncio(loop_scope="module")
async def test_get_all_molecules():
    async with AsyncClient(
        transport=ASGITransport(app=app),
        base_url="http://testserver"
    ) as ac:
        response = await ac.get("/molecules/")
        assert response.status_code == 200
        data = response.json()
        assert isinstance(data, list)
        assert len(data) > 0  # Check that the list is not empty
        assert all("mol_id" in molecule for molecule in data)
        assert all("name" in molecule for molecule in data)
        assert all("smiles" in molecule for molecule in data)


@pytest.mark.asyncio(loop_scope="module")
async def test_add_molecule():
    molecule_data = {"name": "TestMolecule", "smiles": "C=C"}

    async with AsyncClient(
        transport=ASGITransport(app=app),
        base_url="http://testserver"
    ) as ac:
        mol_id = await add_and_verify_molecule(ac, molecule_data)
    await delete_molecule_by_id(mol_id)


@pytest.mark.asyncio(loop_scope="module")
async def test_get_molecule_by_id():
    molecule_data = {"name": "TestMoleculeByID", "smiles": "C#C"}

    async with AsyncClient(
        transport=ASGITransport(app=app),
        base_url="http://testserver"
    ) as ac:
        mol_id = await add_and_verify_molecule(ac, molecule_data)

        response = await ac.get(f"/molecules/{mol_id}")
        assert response.status_code == 200
        retrieved_molecule = response.json()
        assert retrieved_molecule["mol_id"] == mol_id
        assert retrieved_molecule["name"] == molecule_data["name"]
        assert retrieved_molecule["smiles"] == molecule_data["smiles"]

    await delete_molecule_by_id(mol_id)


@pytest.mark.asyncio(loop_scope="module")
async def test_update_molecule_by_id():
    molecule_data = {"name": "MoleculeToUpdate", "smiles": "C#C"}
    updated_data = {"name": "UpdatedMolecule", "smiles": "C=C"}

    async with AsyncClient(
        transport=ASGITransport(app=app),
        base_url="http://testserver"
    ) as ac:
        mol_id = await add_and_verify_molecule(ac, molecule_data)

        response = await ac.put(f"/molecules/{mol_id}", json=updated_data)
        assert response.status_code == 200
        updated_molecule_response = response.json()
        expected_message = f"The molecule with id {mol_id} has been updated!"
        assert updated_molecule_response["message"] == expected_message

    async with async_session_maker() as session:
        result = await session.execute(select(Molecule).filter_by(mol_id=mol_id))
        updated_molecule = result.scalar_one_or_none()
        assert updated_molecule is not None
        assert updated_molecule.name == updated_data["name"]
        assert updated_molecule.smiles == updated_data["smiles"]

    await delete_molecule_by_id(mol_id)


@pytest.mark.asyncio(loop_scope="module")
async def test_delete_molecule_by_id():
    molecule_data = {"name": "MoleculeToDelete", "smiles": "C#N"}

    async with AsyncClient(
        transport=ASGITransport(app=app),
        base_url="http://testserver"
    ) as ac:
        mol_id = await add_and_verify_molecule(ac, molecule_data)

        response = await ac.delete(f"/molecules/delete/{mol_id}")
        assert response.status_code == 200
        delete_response = response.json()
        assert delete_response["message"] == f"The molecule with id {mol_id} is deleted!"

        async with async_session_maker() as session:
            result = await session.execute(select(Molecule).filter_by(mol_id=mol_id))
            deleted_molecule = result.scalar_one_or_none()
            assert deleted_molecule is None


@pytest.mark.asyncio(loop_scope="module")
async def test_search_molecule():
    molecule_data1 = {"name": "Molecule1", "smiles": "C=C"}
    molecule_data2 = {"name": "Molecule2", "smiles": "C#N"}

    async with AsyncClient(
        transport=ASGITransport(app=app),
        base_url="http://testserver"
    ) as ac:
        mol_id1 = await add_and_verify_molecule(ac, molecule_data1)
        mol_id2 = await add_and_verify_molecule(ac, molecule_data2)

        search_params = {"substructure_smiles": "C"}
        response = await ac.get("/molecules/search", params=search_params)
        assert response.status_code == 200

        data = response.json()
        assert isinstance(data, list)
        assert len(data) > 0  # Ensure that at least one molecule was found

    await delete_molecule_by_id(mol_id1)
    await delete_molecule_by_id(mol_id2)
