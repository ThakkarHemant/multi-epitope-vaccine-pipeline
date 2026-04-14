import logging
from datetime import datetime
from sqlalchemy.ext.asyncio import AsyncSession

from database.models import PipelineJob, StageResult
from report_generator.controller import run_report_generation

from services import sequence_service
from services import antigenicity_service
from services.antigenicity_service import AntigenicityGateError
from services import physicochemical_service
from services import epitope_service
from services import safety_service
from services import population_coverage_service
from services import ifn_gamma_service
from services import pepfold_service
from services import docking_service
from services import model_evaluation_service
from services import construct_builder_service
from services import codon_optimization_service
from services import cloning_service

logger = logging.getLogger(__name__)


async def _save(db: AsyncSession, job_id: str, stage_name: str,
                result: dict | None, passed: bool, error: str | None = None):
    stage_id = f"{job_id}_{stage_name.replace(' ', '_').replace('-','').lower()[:40]}"
    db.add(StageResult(
        id=stage_id,
        job_id=job_id,
        stage_name=stage_name,
        result_json=result or {},
        passed=passed,
        error_message=error,
    ))
    await db.commit()


async def run_pipeline(job_id: str, uniprot_id: str, db: AsyncSession):
    """Execute the 13-step MEV Design Pipeline sequentially."""
    try:
        job = await db.get(PipelineJob, job_id)
        if job:
            job.status = "running"
            await db.commit()

        # ------------------------------------------------------------------
        # Step 1 — Sequence Retrieval
        # ------------------------------------------------------------------
        r1 = await sequence_service.run(sequence="", uniprot_id=uniprot_id)
        sequence = r1["sequence"]
        await _save(db, job_id, "Step 01 - Sequence Retrieval", r1, True)
        logger.info("[%s] Step 1 done — %d AA", job_id, len(sequence))

        # ------------------------------------------------------------------
        # Step 2 — Antigenicity Screening (HARD GATE)
        # ------------------------------------------------------------------
        try:
            r2 = await antigenicity_service.run(sequence=sequence)
            await _save(db, job_id, "Step 02 - Antigenicity Screening", r2, True)
            logger.info("[%s] Step 2 done — score=%.4f", job_id, r2["antigenicity_score"])
        except AntigenicityGateError as gate_err:
            await _save(db, job_id, "Step 02 - Antigenicity Screening",
                        {"antigenicity_score": None, "classification": "Non-antigen"},
                        False, str(gate_err))
            job = await db.get(PipelineJob, job_id)
            if job:
                job.status = "failed"
                job.error_message = str(gate_err)
                job.completed_at = datetime.utcnow()
                await db.commit()
            logger.warning("[%s] Pipeline TERMINATED at Step 2 — antigenicity gate failed", job_id)
            return  # ← HARD STOP, no further steps

        # ------------------------------------------------------------------
        # Step 3 — Physicochemical Analysis
        # ------------------------------------------------------------------
        r3 = await physicochemical_service.run(sequence=sequence)
        await _save(db, job_id, "Step 03 - Physicochemical Analysis", r3, True)
        logger.info("[%s] Step 3 done", job_id)

        # ------------------------------------------------------------------
        # Step 4 — Epitope Prediction (4A B-cell, 4B MHC-I, 4C MHC-II)
        # ------------------------------------------------------------------
        r4 = await epitope_service.run(sequence=sequence, uniprot_id=uniprot_id)
        await _save(db, job_id, "Step 04 - Epitope Prediction", r4, True)
        logger.info("[%s] Step 4 done — bcell=%d mhc1=%d mhc2=%d",
                    job_id, r4["bcell_count"], r4["mhc1_count"], r4["mhc2_count"])

        # ------------------------------------------------------------------
        # Step 5 — Safety & Antigenicity Filtering of Epitopes
        # ------------------------------------------------------------------
        r5 = await safety_service.run(
            bcell_epitopes=r4["bcell_epitopes"],
            mhc1_epitopes=r4["mhc1_epitopes"],
            mhc2_epitopes=r4["mhc2_epitopes"],
        )
        await _save(db, job_id, "Step 05 - Safety Filtering", r5, True)
        logger.info("[%s] Step 5 done — safe B=%d I=%d II=%d",
                    job_id, r5["safe_bcell_count"], r5["safe_mhc1_count"], r5["safe_mhc2_count"])

        # ------------------------------------------------------------------
        # Step 6 — Population Coverage Analysis
        # ------------------------------------------------------------------
        r6 = await population_coverage_service.run(
            bcell_epitopes=r5["safe_bcell_epitopes"],
            mhc1_epitopes=r5["safe_mhc1_epitopes"],
            mhc2_epitopes=r5["safe_mhc2_epitopes"],
        )
        await _save(db, job_id, "Step 06 - Population Coverage", r6, True)
        logger.info("[%s] Step 6 done — world coverage class1=%.1f%% class2=%.1f%%",
                    job_id, r6["world_class1_coverage_pct"], r6["world_class2_coverage_pct"])

        # ------------------------------------------------------------------
        # Step 7 — IFN-Gamma Epitope Prediction
        # ------------------------------------------------------------------
        r7 = await ifn_gamma_service.run(safe_mhc2_epitopes=r5["safe_mhc2_epitopes"])
        await _save(db, job_id, "Step 07 - IFN-Gamma Prediction", r7, True)
        logger.info("[%s] Step 7 done — ifn_gamma=%d", job_id, r7["ifn_gamma_count"])

        # ------------------------------------------------------------------
        # Step 8 — PEP-FOLD 3D Structure Modeling
        # ------------------------------------------------------------------
        r8 = await pepfold_service.run(safe_mhc1_epitopes=r5["safe_mhc1_epitopes"])
        await _save(db, job_id, "Step 08 - PEP-FOLD Modeling", r8, True)
        logger.info("[%s] Step 8 done — %d peptides modeled", job_id, r8["successful_models"])

        # ------------------------------------------------------------------
        # Step 9 — Molecular Docking (HDOCK)
        # ------------------------------------------------------------------
        r9 = await docking_service.run(
            safe_mhc1_epitopes=r5["safe_mhc1_epitopes"],
            safe_mhc2_epitopes=r5["safe_mhc2_epitopes"],
        )
        await _save(db, job_id, "Step 09 - Molecular Docking", r9, True)
        logger.info("[%s] Step 9 done — %d complexes docked", job_id, r9["total_docked"])

        # ------------------------------------------------------------------
        # Step 10 — Docking Model Evaluation (PyMOL)
        # ------------------------------------------------------------------
        r10 = await model_evaluation_service.run(docked_complexes=r9["docked_complexes"])
        await _save(db, job_id, "Step 10 - Model Evaluation", r10, True)
        logger.info("[%s] Step 10 done — %d best models", job_id, r10["best_model_count"])

        # ------------------------------------------------------------------
        # Step 11 — Vaccine Construct Design
        # ------------------------------------------------------------------
        r11 = await construct_builder_service.run(
            safe_mhc1_epitopes=r5["safe_mhc1_epitopes"],
            safe_mhc2_epitopes=r5["safe_mhc2_epitopes"],
            bcell_epitopes=r5["safe_bcell_epitopes"],
        )
        await _save(db, job_id, "Step 11 - Construct Design", r11, True)
        logger.info("[%s] Step 11 done — construct=%d AA", job_id, r11["construct_length_aa"])

        # ------------------------------------------------------------------
        # Step 12 — Codon Optimization
        # ------------------------------------------------------------------
        r12 = await codon_optimization_service.run(
            construct_sequence=r11["final_construct_sequence"]
        )
        await _save(db, job_id, "Step 12 - Codon Optimization", r12, True)
        logger.info("[%s] Step 12 done — DNA=%d bp CAI=%.2f GC=%.1f%%",
                    job_id, r12["dna_length_bp"], r12["cai_score"], r12["gc_content_pct"])

        # ------------------------------------------------------------------
        # Step 13 — In-Silico Cloning
        # ------------------------------------------------------------------
        r13 = await cloning_service.run(optimized_dna=r12["optimized_dna_sequence"])
        await _save(db, job_id, "Step 13 - In-Silico Cloning", r13, True)
        logger.info("[%s] Step 13 done — status=%s", job_id, r13["cloning_status"])

        # ------------------------------------------------------------------
        # Report generation + job completion
        # ------------------------------------------------------------------
        job = await db.get(PipelineJob, job_id)
        if job:
            job.status = "completed"
            job.completed_at = datetime.utcnow()
            await db.commit()
        await run_report_generation(job_id, db)
        logger.info("[%s] Pipeline completed successfully", job_id)

    except Exception as exc:
        logger.exception("[%s] Pipeline failed: %s", job_id, exc)
        job = await db.get(PipelineJob, job_id)
        if job:
            job.status = "failed"
            job.error_message = str(exc)[:500]
            job.completed_at = datetime.utcnow()
            await db.commit()
