from fastapi import FastAPI
from os import getenv
from src.molecules.router import router as molecule_router

app = FastAPI()


@app.get("/")
def get_server():
    """
    Retrieve the server ID.
    """
    return {"server_id": getenv("SERVER_ID", "1")}


app.include_router(molecule_router)
