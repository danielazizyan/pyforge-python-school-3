from sqlalchemy.future import select
from src.molecules.models import Molecule
from src.database import async_session_maker
from src.dao.base import BaseDAO
from rdkit import Chem
import logging

logger = logging.getLogger(__name__)


class MoleculeDAO(BaseDAO):
    model = Molecule

    @classmethod
    async def find_all_molecules(cls, limit: int = None):
        async with async_session_maker() as session:
            query = select(cls.model)
            if limit:
                query = query.limit(limit)
            logger.debug(f"Executing query: {query}")
            result = await session.execute(query)

            for molecule in result.scalars():
                yield molecule

    @classmethod
    async def find_full_data(cls, mol_id: int):
        molecule_info = await cls.find_one_or_none(mol_id=mol_id)
        if not molecule_info:
            return None
        molecule_data = molecule_info.__dict__
        logger.debug(f"Found molecule: {molecule_data}")
        return molecule_data

    @classmethod
    async def add_molecule(cls, **molecule_data: dict):
        new_molecule = await cls.add(**molecule_data)
        logger.debug(f"Added molecule: {new_molecule}")
        return new_molecule.mol_id

    @classmethod
    async def delete_molecule_by_id(cls, mol_id: int):
        result = await cls.delete(mol_id=mol_id)
        if result == 0:
            logger.debug(f"No molecule found with ID {mol_id} to delete.")
            return None
        logger.debug(f"Deleted molecule with ID {mol_id}")
        return mol_id

    @classmethod
    async def update_molecule(cls, mol_id: int, **updated_data: dict):
        result = await cls.update({'mol_id': mol_id}, **updated_data)
        if result == 0:
            logger.debug(f"No molecule found with ID {mol_id} to update.")
            return None
        logger.debug(f"Updated molecule with ID {mol_id}")
        return result

    @classmethod
    async def search_molecule(cls, substructure_smiles: str):
        async with async_session_maker() as session:
            substructure = Chem.MolFromSmiles(substructure_smiles)

            query = select(cls.model)
            logger.debug(f"Executing query: {query}")
            results = await session.execute(query)
            molecules = results.scalars().all()

            matches = []
            for molecule in molecules:
                if Chem.MolFromSmiles(molecule.smiles).HasSubstructMatch(substructure):
                    matches.append(molecule)

            logger.debug(
                f"Found {len(matches)} molecules matching substructure "
                f"{substructure_smiles}."
            )
            return matches
