from src.celery_worker import celery
from src.molecules.dao import MoleculeDAO
from src.molecules.schema import MoleculeResponse
from rdkit import Chem
import json
from src.config import redis_client
import logging

logger = logging.getLogger(__name__)


@celery.task(bind=True)
def substructure_search_task(self, substructure_smiles: str):

    cache_key = f"search:{substructure_smiles}"
    cached_result = redis_client.get(cache_key)
    if cached_result:
        logger.info(f"Returning cached result for substructure: {substructure_smiles}")
        return cached_result

    try:
        logger.info(f"Starting substructure search for: {substructure_smiles}")
        if Chem.MolFromSmiles(substructure_smiles) is None:
            logger.error(f"Invalid SMILES structure: {substructure_smiles}")
            raise ValueError("This SMILES has an invalid structure")

        matches = MoleculeDAO.search_molecule(substructure_smiles=substructure_smiles)
        if matches:
            result = [
                MoleculeResponse.model_validate(mol).model_dump()
                for mol in matches
            ]
            logger.info(
                f"Found {len(result)} matches for substructure: {substructure_smiles}"
            )
        else:
            logger.info(f"No matches found for substructure: {substructure_smiles}")
            result = []

        result_json = json.dumps(result)

        redis_client.setex(cache_key, 600, result_json)
        return result_json
    except Exception as e:
        logger.exception(f"An error occurred during substructure search: {e}")
        raise e
