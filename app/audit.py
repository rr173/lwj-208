import json
import re
from typing import Optional, Tuple
from fastapi import Request, Response
from starlette.middleware.base import BaseHTTPMiddleware
from starlette.types import ASGIApp
from starlette.responses import StreamingResponse
from app.database import SessionLocal
from app import models


WRITE_METHODS = {"POST", "PUT", "DELETE", "PATCH"}


RESOURCE_PATTERNS = [
    (re.compile(r"^/api/reference(?:/.*)?$"), "reference"),
    (re.compile(r"^/api/alignment(?:/.*)?$"), "alignment"),
    (re.compile(r"^/api/samples(?:/[^/]+)?(?:/alignments(?:/.*)?)?$"), "sample"),
    (re.compile(r"^/api/samples/[^/]+/spectrum$"), "sample_spectrum"),
    (re.compile(r"^/api/samples/frequency/.*$"), "population_frequency"),
    (re.compile(r"^/api/samples/compare/.*$"), "sample_comparison"),
    (re.compile(r"^/api/samples/hotspots/.*$"), "hotspot_analysis"),
    (re.compile(r"^/api/batch(?:/.*)?$"), "batch_task"),
    (re.compile(r"^/api/phylogeny(?:/.*)?$"), "phylogeny"),
    (re.compile(r"^/api/scoring/rules(?:/.*)?$"), "rule"),
    (re.compile(r"^/api/scoring/samples/[^/]+/score$"), "scoring"),
    (re.compile(r"^/api/ld(?:/.*)?$"), "ld_analysis"),
    (re.compile(r"^/api/primer(?:/.*)?$"), "primer_design"),
    (re.compile(r"^/api/transmission(?:/.*)?$"), "transmission_analysis"),
    (re.compile(r"^/api/domain(?:/.*)?$"), "protein_domain"),
    (re.compile(r"^/api/synteny(?:/.*)?$"), "synteny_analysis"),
    (re.compile(r"^/api/stats(?:/.*)?$"), "stats"),
]


def method_to_operation_type(method: str) -> str:
    mapping = {
        "POST": "create",
        "PUT": "update",
        "DELETE": "delete",
        "PATCH": "update",
    }
    return mapping.get(method, "other")


def get_resource_type(path: str) -> str:
    for pattern, resource_type in RESOURCE_PATTERNS:
        if pattern.match(path):
            return resource_type
    return "other"


def extract_resource_id(path: str) -> Optional[str]:
    patterns = [
        re.compile(r"^/api/samples/(\d+)(?:/|$)"),
        re.compile(r"^/api/alignment/(\d+)(?:/|$)"),
        re.compile(r"^/api/scoring/rules/(\d+)(?:/|$)"),
        re.compile(r"^/api/reference/([^/]+)(?:/|$)"),
        re.compile(r"^/api/batch/([^/]+)(?:/|$)"),
        re.compile(r"^/api/phylogeny/([^/]+)(?:/|$)"),
    ]
    for pattern in patterns:
        match = pattern.match(path)
        if match:
            return match.group(1)
    return None


def generate_summary(
    method: str,
    path: str,
    resource_type: str,
    resource_id: Optional[str],
    request_body: Optional[str],
    status_code: int,
) -> str:
    operation = method_to_operation_type(method)

    if resource_type == "reference":
        if "/fasta" in path:
            return "上传参考序列FASTA"
        elif "/gff" in path:
            return "上传基因注释GFF"
        elif resource_id:
            if operation == "delete":
                return f"删除参考序列 {resource_id}"
            elif operation == "update":
                return f"更新参考序列 {resource_id}"
            else:
                return f"查询参考序列 {resource_id}"
        return "参考序列操作"

    if resource_type == "sample":
        if "/alignments" in path:
            if operation == "create":
                return f"为样本 {resource_id} 关联比对结果"
            elif operation == "delete":
                return f"移除样本 {resource_id} 的比对关联"
        if operation == "create":
            try:
                if request_body:
                    data = json.loads(request_body)
                    name = data.get("name", "")
                    if name:
                        return f"创建样本 {name}"
            except (json.JSONDecodeError, TypeError):
                pass
            return "创建样本"
        elif operation == "update":
            return f"更新样本 {resource_id}"
        elif operation == "delete":
            return f"删除样本 {resource_id}"
        return f"样本 {resource_id} 操作"

    if resource_type == "alignment":
        if operation == "create":
            return "执行序列比对"
        elif operation == "delete":
            return f"删除比对结果 {resource_id}"
        return f"比对结果 {resource_id} 操作"

    if resource_type == "rule":
        if operation == "create":
            try:
                if request_body:
                    data = json.loads(request_body)
                    name = data.get("name", "")
                    if name:
                        return f"创建评分规则 {name}"
            except (json.JSONDecodeError, TypeError):
                pass
            return "创建评分规则"
        elif operation == "update":
            return f"更新评分规则 {resource_id}"
        elif operation == "delete":
            return f"删除评分规则 {resource_id}"
        return f"评分规则 {resource_id} 操作"

    if resource_type == "batch_task":
        if operation == "create":
            return "创建批量比对任务"
        return f"批量任务 {resource_id} 操作"

    if resource_type == "phylogeny":
        if operation == "create":
            return "创建系统发育树任务"
        return f"系统发育树 {resource_id} 操作"

    if resource_type == "scoring":
        if operation == "create":
            return f"对样本 {resource_id} 进行致病性评分"
        return f"样本评分 {resource_id} 操作"

    if resource_type == "ld_analysis":
        if operation == "create":
            return "执行连锁不平衡分析"
        return "连锁不平衡分析"

    if resource_type == "primer_design":
        if operation == "create":
            return "设计引物"
        return "引物设计"

    if resource_type == "protein_domain":
        if operation == "create":
            return "添加蛋白质结构域"
        return "蛋白质结构域操作"

    if resource_type == "synteny_analysis":
        if operation == "create":
            return "执行共线性分析"
        return "共线性分析"

    if resource_type == "transmission_analysis":
        if operation == "create":
            return "执行传播链分析"
        return "传播分析"

    status_part = "成功" if 200 <= status_code < 300 else "失败"
    return f"{operation} {resource_type} - {status_part}"


class AuditLogMiddleware(BaseHTTPMiddleware):
    def __init__(self, app: ASGIApp):
        super().__init__(app)

    async def dispatch(self, request: Request, call_next):
        method = request.method
        path = request.url.path

        if method not in WRITE_METHODS:
            response = await call_next(request)
            return response

        if path.startswith("/health") or path.startswith("/docs") or path.startswith("/openapi.json") or path.startswith("/ws") or path.startswith("/api/audit"):
            response = await call_next(request)
            return response

        request_body_bytes = await request.body()
        request_body = request_body_bytes.decode("utf-8") if request_body_bytes else None

        async def receive_with_body():
            return {
                "type": "http.request",
                "body": request_body_bytes,
                "more_body": False,
            }

        original_receive = request._receive
        request._receive = receive_with_body

        try:
            response = await call_next(request)
        except Exception as e:
            self._save_audit_log(
                method=method,
                path=path,
                request_body=request_body,
                status_code=500,
                response_body=None,
            )
            raise e
        finally:
            request._receive = original_receive

        response_body = None
        try:
            if hasattr(response, "body") and response.body:
                response_body = response.body.decode("utf-8") if isinstance(response.body, bytes) else str(response.body)
            else:
                chunks = []
                async for chunk in response.body_iterator:
                    if isinstance(chunk, str):
                        chunk = chunk.encode("utf-8")
                    chunks.append(chunk)
                if chunks:
                    response_body = b"".join(chunks).decode("utf-8")
                    from starlette.responses import Response as StarletteResponse
                    response = StarletteResponse(
                        content=b"".join(chunks),
                        status_code=response.status_code,
                        headers=dict(response.headers),
                        media_type=response.media_type,
                    )
        except Exception:
            pass

        self._save_audit_log(
            method=method,
            path=path,
            request_body=request_body,
            status_code=response.status_code,
            response_body=response_body,
        )

        return response

    def _save_audit_log(
        self,
        method: str,
        path: str,
        request_body: Optional[str],
        status_code: int,
        response_body: Optional[str] = None,
    ):
        try:
            resource_type = get_resource_type(path)
            resource_id = extract_resource_id(path)

            if not resource_id and response_body and 200 <= status_code < 300:
                try:
                    data = json.loads(response_body)
                    if isinstance(data, dict) and "id" in data:
                        resource_id = str(data["id"])
                    elif isinstance(data, list) and len(data) > 0 and isinstance(data[0], dict) and "id" in data[0]:
                        resource_id = str(data[0]["id"])
                except (json.JSONDecodeError, TypeError):
                    pass

            operation_type = method_to_operation_type(method)
            summary = generate_summary(
                method, path, resource_type, resource_id, request_body, status_code
            )

            db = SessionLocal()
            try:
                audit_log = models.AuditLog(
                    operation_type=operation_type,
                    resource_type=resource_type,
                    resource_id=str(resource_id) if resource_id else None,
                    summary=summary,
                    method=method,
                    path=path,
                    status_code=status_code,
                    request_body=request_body,
                )
                db.add(audit_log)
                db.commit()
            finally:
                db.close()
        except Exception:
            pass


def get_audit_logs(
    db,
    start_time=None,
    end_time=None,
    resource_type=None,
    operation_type=None,
    page: int = 1,
    page_size: int = 20,
):
    query = db.query(models.AuditLog)

    if start_time:
        query = query.filter(models.AuditLog.timestamp >= start_time)
    if end_time:
        query = query.filter(models.AuditLog.timestamp <= end_time)
    if resource_type:
        query = query.filter(models.AuditLog.resource_type == resource_type)
    if operation_type:
        query = query.filter(models.AuditLog.operation_type == operation_type)

    total = query.count()
    query = query.order_by(models.AuditLog.timestamp.desc())
    query = query.offset((page - 1) * page_size).limit(page_size)
    items = query.all()

    return {
        "total": total,
        "page": page,
        "page_size": page_size,
        "items": items,
    }
