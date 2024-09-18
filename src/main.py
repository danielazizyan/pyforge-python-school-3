from fastapi import FastAPI
from os import getenv
from src.molecules.router import router as molecule_router
import logging

app = FastAPI()


@app.get("/")
def get_server():
    """
    Retrieve the server ID.
    """
    return {"server_id": getenv("SERVER_ID", "1")}


app.include_router(molecule_router)


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler("app.log")
    ]
)

logger = logging.getLogger(__name__)

logger.info("FastAPI application startup complete.")
