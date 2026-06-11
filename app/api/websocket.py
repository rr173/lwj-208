import json
from fastapi import APIRouter, WebSocket, WebSocketDisconnect
from app.services.batch_service import task_manager

router = APIRouter(prefix="/ws", tags=["websocket"])


@router.websocket("/batch/{task_id}")
async def websocket_endpoint(websocket: WebSocket, task_id: str):
    await websocket.accept()

    async def progress_callback(task_id_str, completed, total, status):
        await websocket.send_text(json.dumps({
            "task_id": task_id_str,
            "completed": completed,
            "total": total,
            "status": status,
        }))

    task_manager.register_callback(task_id, progress_callback)

    try:
        while True:
            data = await websocket.receive_text()
    except WebSocketDisconnect:
        task_manager.unregister_callback(task_id, progress_callback)
