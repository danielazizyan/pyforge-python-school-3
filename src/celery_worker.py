from celery import Celery
from src.config import settings

celery = Celery(
    'tasks',
    broker=settings.CELERY_BROKER_URL,
    backend=settings.CELERY_RESULT_BACKEND,
)

celery.conf.update(
    task_track_started=True,
    task_always_eager=settings.CELERY_TASK_ALWAYS_EAGER,
    task_store_eager_result=True,
)
