import uuid
import asyncio
import threading
from typing import List, Dict, Optional
from sqlalchemy.orm import Session
from app import models
from app.services import genome_service
from app.database import SessionLocal


class TaskManager:
    def __init__(self):
        self._progress_callbacks: Dict[str, list] = {}
        self._loop: Optional[asyncio.AbstractEventLoop] = None

    def set_loop(self, loop: asyncio.AbstractEventLoop):
        self._loop = loop

    def register_callback(self, task_id: str, callback):
        if task_id not in self._progress_callbacks:
            self._progress_callbacks[task_id] = []
        self._progress_callbacks[task_id].append(callback)

    def unregister_callback(self, task_id: str, callback):
        if task_id in self._progress_callbacks:
            if callback in self._progress_callbacks[task_id]:
                self._progress_callbacks[task_id].remove(callback)

    async def notify_progress(self, task_id: str, completed: int, total: int, status: str):
        if task_id in self._progress_callbacks:
            for callback in self._progress_callbacks[task_id]:
                try:
                    await callback(task_id, completed, total, status)
                except Exception:
                    pass

    def notify_progress_sync(self, task_id: str, completed: int, total: int, status: str):
        if self._loop and task_id in self._progress_callbacks:
            coro = self.notify_progress(task_id, completed, total, status)
            asyncio.run_coroutine_threadsafe(coro, self._loop)


task_manager = TaskManager()


def create_batch_task(db: Session, query_sequences: List[str]) -> models.BatchTask:
    """Create a new batch alignment task."""
    task_uuid = str(uuid.uuid4())

    task = models.BatchTask(
        task_id=task_uuid,
        status="pending",
        total=len(query_sequences),
        completed=0,
    )
    db.add(task)
    db.flush()

    for i, seq in enumerate(query_sequences):
        item = models.BatchTaskItem(
            task_id=task.id,
            order_index=i,
            query_sequence=seq,
            status="pending",
        )
        db.add(item)

    db.commit()
    db.refresh(task)
    return task


def get_batch_task(db: Session, task_id: str) -> Optional[models.BatchTask]:
    """Get batch task by task_id."""
    return db.query(models.BatchTask).filter(
        models.BatchTask.task_id == task_id
    ).first()


def run_batch_task_sync(task_id: str):
    """Run a batch alignment task in a background thread."""
    db = SessionLocal()
    try:
        task = get_batch_task(db, task_id)
        if not task:
            return

        task.status = "running"
        db.commit()

        task_manager.notify_progress_sync(task_id, 0, task.total, "running")

        items = db.query(models.BatchTaskItem).filter(
            models.BatchTaskItem.task_id == task.id
        ).order_by(models.BatchTaskItem.order_index).all()

        references = genome_service.get_all_references(db)

        for idx, item in enumerate(items):
            item.status = "processing"
            db.commit()

            try:
                results = genome_service.align_query_all_references(db, item.query_sequence)
                if results:
                    item.alignment_id = results[0].id
                item.status = "completed"
            except Exception as e:
                item.status = "failed"

            task.completed = idx + 1
            db.commit()

            task_manager.notify_progress_sync(task_id, task.completed, task.total, "running")

        task.status = "completed"
        from datetime import datetime
        task.completed_at = datetime.utcnow()
        db.commit()

        task_manager.notify_progress_sync(task_id, task.total, task.total, "completed")

    finally:
        db.close()


def start_batch_task(task_id: str):
    """Start batch task in background thread."""
    thread = threading.Thread(target=run_batch_task_sync, args=(task_id,), daemon=True)
    thread.start()
