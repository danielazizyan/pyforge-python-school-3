import asyncio
from sqlalchemy import text
from src.database import async_session_maker
from src.molecules.models import Molecule

smiles_db = [
    {"name": "Ethanol", "smiles": "CCO"},
    {"name": "Benzene", "smiles": "c1ccccc1"},
    {"name": "Acetic acid", "smiles": "CC(=O)O"},
    {"name": "Aspirin", "smiles": "CC(=O)Oc1ccccc1C(=O)O"},
    {"name": "Methane", "smiles": "C"},
    {"name": "Propane", "smiles": "CCC"},
    {"name": "Toluene", "smiles": "Cc1ccccc1"},
    {"name": "Phenol", "smiles": "c1ccc(cc1)O"},
    {"name": "Formaldehyde", "smiles": "C=O"},
    {"name": "Acetone", "smiles": "CC(=O)C"},
    {"name": "Cyclohexane", "smiles": "C1CCCCC1"},
    {"name": "Ethylene", "smiles": "C=C"},
    {"name": "Styrene", "smiles": "c1ccccc1C=C"},
    {"name": "Hexane", "smiles": "CCCCCC"}
]


async def truncate_and_populate_db():
    async with async_session_maker() as session:
        # Truncate the table
        await session.execute(text("TRUNCATE TABLE molecules RESTART IDENTITY CASCADE"))
        await session.commit()  # Commit the truncation

        # Populate the table
        async with session.begin():
            for molecule_data in smiles_db:
                molecule = Molecule(**molecule_data)
                session.add(molecule)
                print(f"Inserted molecule: {molecule.name}")
            await session.commit()
            print("Commit complete")

if __name__ == "__main__":
    asyncio.run(truncate_and_populate_db())
