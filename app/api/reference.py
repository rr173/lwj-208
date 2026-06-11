from fastapi import APIRouter, Depends, HTTPException, UploadFile, File, Form
from sqlalchemy.orm import Session
from typing import List
from app import models, schemas
from app.database import get_db
from app.services import genome_service

router = APIRouter(prefix="/api/reference", tags=["reference"])


@router.get("", response_model=List[schemas.ReferenceSequenceOut])
def list_references(db: Session = Depends(get_db)):
    refs = genome_service.get_all_references(db)
    return refs


@router.get("/{name}", response_model=schemas.ReferenceSequenceDetail)
def get_reference(name: str, db: Session = Depends(get_db)):
    ref = genome_service.get_reference_by_name(db, name)
    if not ref:
        raise HTTPException(status_code=404, detail="Reference sequence not found")
    return ref


@router.post("/fasta", response_model=List[schemas.ReferenceSequenceOut])
def upload_fasta(fasta_text: str = Form(...), db: Session = Depends(get_db)):
    sequences = genome_service.parse_fasta(fasta_text)
    if not sequences:
        raise HTTPException(status_code=400, detail="No valid sequences found in FASTA")

    results = []
    for name, seq in sequences:
        if len(seq) > 100000:
            raise HTTPException(status_code=400, detail=f"Sequence {name} exceeds 100,000 base limit")
        ref = genome_service.create_reference_sequence(db, name, seq)
        results.append(ref)

    return results


@router.post("/fasta/upload", response_model=List[schemas.ReferenceSequenceOut])
async def upload_fasta_file(file: UploadFile = File(...), db: Session = Depends(get_db)):
    content = await file.read()
    fasta_text = content.decode("utf-8")
    sequences = genome_service.parse_fasta(fasta_text)
    if not sequences:
        raise HTTPException(status_code=400, detail="No valid sequences found in FASTA")

    results = []
    for name, seq in sequences:
        if len(seq) > 100000:
            raise HTTPException(status_code=400, detail=f"Sequence {name} exceeds 100,000 base limit")
        ref = genome_service.create_reference_sequence(db, name, seq)
        results.append(ref)

    return results


@router.post("/gff", response_model=List[schemas.GeneAnnotationOut])
def upload_gff(gff_text: str = Form(...), db: Session = Depends(get_db)):
    annotations = genome_service.parse_gff(gff_text)
    if not annotations:
        raise HTTPException(status_code=400, detail="No valid annotations found in GFF")

    seq_names = set(a["seq_name"] for a in annotations)
    results = []

    for seq_name in seq_names:
        ref = genome_service.get_reference_by_name(db, seq_name)
        if not ref:
            raise HTTPException(status_code=404, detail=f"Reference sequence '{seq_name}' not found")

        seq_annotations = [a for a in annotations if a["seq_name"] == seq_name]
        created = genome_service.add_gene_annotations(db, ref.id, seq_annotations)
        results.extend(created)

    return results


@router.post("/gff/upload", response_model=List[schemas.GeneAnnotationOut])
async def upload_gff_file(file: UploadFile = File(...), db: Session = Depends(get_db)):
    content = await file.read()
    gff_text = content.decode("utf-8")
    annotations = genome_service.parse_gff(gff_text)
    if not annotations:
        raise HTTPException(status_code=400, detail="No valid annotations found in GFF")

    seq_names = set(a["seq_name"] for a in annotations)
    results = []

    for seq_name in seq_names:
        ref = genome_service.get_reference_by_name(db, seq_name)
        if not ref:
            raise HTTPException(status_code=404, detail=f"Reference sequence '{seq_name}' not found")

        seq_annotations = [a for a in annotations if a["seq_name"] == seq_name]
        created = genome_service.add_gene_annotations(db, ref.id, seq_annotations)
        results.extend(created)

    return results


@router.get("/{name}/annotations", response_model=List[schemas.GeneAnnotationOut])
def get_annotations(name: str, db: Session = Depends(get_db)):
    ref = genome_service.get_reference_by_name(db, name)
    if not ref:
        raise HTTPException(status_code=404, detail="Reference sequence not found")

    return db.query(models.GeneAnnotation).filter(
        models.GeneAnnotation.reference_id == ref.id
    ).all()
