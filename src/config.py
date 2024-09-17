import os
import logging
import redis
from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    DB_HOST: str
    DB_PORT: int
    DB_NAME: str
    DB_USER: str
    DB_PASSWORD: str
    REDIS_HOST: str = "redis"
    REDIS_PORT: int = 6379

    CELERY_BROKER_URL: str = f"redis://{REDIS_HOST}:{REDIS_PORT}/0"
    CELERY_RESULT_BACKEND: str = f"redis://{REDIS_HOST}:{REDIS_PORT}/0"
    CELERY_TASK_ALWAYS_EAGER: bool = False

    # This ensures that .env variables are automatically loaded
    model_config = SettingsConfigDict(
        env_file=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".env")
    )


settings = Settings()


def get_db_url():
    """Constructs the database URL from the settings."""
    return (
        f"postgresql+asyncpg://{settings.DB_USER}:{settings.DB_PASSWORD}@"
        f"{settings.DB_HOST}:{settings.DB_PORT}/{settings.DB_NAME}"
    )


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler("app.log")
    ]
)

if os.getenv('PYTEST_CURRENT_TEST'):
    settings.REDIS_HOST = 'localhost'
    settings.DB_HOST = 'localhost'
    settings.DB_PORT = 5433
    settings.CELERY_TASK_ALWAYS_EAGER = True

settings.CELERY_BROKER_URL = f"redis://{settings.REDIS_HOST}:{settings.REDIS_PORT}/0"
settings.CELERY_RESULT_BACKEND = f"redis://{settings.REDIS_HOST}:{settings.REDIS_PORT}/0"

redis_client = redis.Redis(
    host=settings.REDIS_HOST,
    port=settings.REDIS_PORT,
    db=0,
    decode_responses=True
)
