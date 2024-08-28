from sqlalchemy import delete, update
from sqlalchemy.future import select
from src.molecules.models import Molecule
from src.database import async_session_maker
from src.dao.base import BaseDAO
from rdkit import Chem


class MoleculeDAO(BaseDAO):
    model = Molecule

    @classmethod
    async def find_all_molecules(cls):
        async with async_session_maker() as session:
            query = select(cls.model)
            molecules = await session.execute(query)
            return molecules.scalars().all()

    @classmethod
    async def find_full_data(cls, mol_id: int):
        async with async_session_maker() as session:
            query = select(cls.model).filter_by(mol_id=mol_id)
            result = await session.execute(query)
            molecule_info = result.scalar_one_or_none()
            if not molecule_info:
                return None
            molecule_data = molecule_info.__dict__
            return molecule_data

    @classmethod
    async def add_molecule(cls, **molecule_data: dict):
        async with async_session_maker() as session:
            async with session.begin():
                new_molecule = cls.model(**molecule_data)
                session.add(new_molecule)
                await session.flush()
                new_molecule_id = new_molecule.mol_id
                await session.commit()
                return new_molecule_id

    @classmethod
    async def delete_molecule_by_id(cls, mol_id: int):
        async with async_session_maker() as session:
            async with session.begin():
                query = select(cls.model).filter_by(mol_id=mol_id)
                result = await session.execute(query)
                molecule_to_delete = result.scalar_one_or_none()
                if not molecule_to_delete:
                    return None
                await session.execute(delete(cls.model).filter_by(mol_id=mol_id))
                await session.commit()
                return mol_id

    @classmethod
    async def update_molecule(cls, mol_id: int, **updated_data: dict):
        async with async_session_maker() as session:
            async with session.begin():
                query = (
                    update(cls.model)
                    .where(cls.model.mol_id == mol_id)
                    .values(**updated_data)
                    .execution_options(synchronize_session="fetch")
                )
                result = await session.execute(query)
                await session.commit()
                return result.rowcount

    @classmethod
    async def search_molecule(cls, substructure_smiles: str):
        async with async_session_maker() as session:
            substructure = Chem.MolFromSmiles(substructure_smiles)

            query = select(cls.model)
            results = await session.execute(query)
            molecules = results.scalars().all()

            matches = []
            for molecule in molecules:
                if Chem.MolFromSmiles(molecule.smiles).HasSubstructMatch(substructure):
                    matches.append(molecule)

            return matches
