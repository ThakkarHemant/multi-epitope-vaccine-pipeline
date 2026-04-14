from datetime import datetime
from sqlalchemy import Column, String, DateTime, ForeignKey, Boolean, JSON
from sqlalchemy.orm import relationship

from .base import Base

class PipelineJob(Base):
    __tablename__ = "pipeline_jobs"

    id = Column(String, primary_key=True, index=True)
    uniprot_id = Column(String, nullable=False, index=True)
    status = Column(String, nullable=False, default="pending")  # pending, running, completed, failed
    created_at = Column(DateTime, default=datetime.utcnow)
    completed_at = Column(DateTime, nullable=True)
    error_message = Column(String, nullable=True)

    results = relationship("StageResult", back_populates="job", cascade="all, delete-orphan")


class StageResult(Base):
    __tablename__ = "stage_results"

    id = Column(String, primary_key=True, index=True)
    job_id = Column(String, ForeignKey("pipeline_jobs.id"), nullable=False)
    stage_name = Column(String, nullable=False)
    result_json = Column(JSON, nullable=True)
    passed = Column(Boolean, nullable=True)  # True if passed decision gate, False if failed/flagged
    error_message = Column(String, nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow)

    job = relationship("PipelineJob", back_populates="results")
