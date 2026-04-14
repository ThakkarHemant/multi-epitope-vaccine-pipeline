import json
import os
import re
from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession
from fpdf import FPDF

from database.models import PipelineJob, StageResult


# ─────────────────────────────────────────────────────────────────────────────
# DISCLAIMER TEXT
# ─────────────────────────────────────────────────────────────────────────────

EPITOPE_DISCLAIMER = (
    "Disclaimer: The epitopes presented in this report have been selected based on "
    "population coverage analysis and promiscuous allele binding profiles. "
    "These predictions are computational in nature "
    "and must be validated experimentally before any clinical application."
)

def _sanitize(text: str) -> str:
    """Replace Unicode characters Helvetica cannot render."""
    return (
        text
        .replace("\u2013", "-").replace("\u2014", "-")
        .replace("\u2018", "'").replace("\u2019", "'")
        .replace("\u201c", '"').replace("\u201d", '"')
        .replace("\u2026", "...")
        .encode("latin-1", errors="replace")
        .decode("latin-1")
    )


def _get_stage_num(stage_name: str) -> int:
    """Extract integer step number for sorting ('Step 10' -> 10)."""
    match = re.search(r'Step\s*(\d+)', stage_name)
    return int(match.group(1)) if match else 999


# --------------------------------------------------------------------------
# JSON Report  — errors field included as-is
# --------------------------------------------------------------------------

async def generate_json_report(job_id: str, db: AsyncSession) -> dict:
    job = await db.get(PipelineJob, job_id)
    if not job:
        raise ValueError(f"Job {job_id} not found.")

    results_q = await db.execute(
        select(StageResult).where(StageResult.job_id == job_id)
    )
    stage_results = sorted(results_q.scalars().all(), key=lambda s: _get_stage_num(s.stage_name))

    stages = {}
    flags = []
    for sr in stage_results:
        stages[sr.stage_name] = {
            "result": sr.result_json,
            "passed": sr.passed,
            "error": sr.error_message,
        }
        if sr.passed is False:
            flags.append(sr.stage_name)

    # Candidacy assessment
    antigenicity_gate = stages.get("Step 02 - Antigenicity Screening", {}).get("passed", False)
    safe_mhc1 = (stages.get("Step 05 - Safety Filtering", {}).get("result") or {}).get("safe_mhc1_count", 0)
    safe_mhc2 = (stages.get("Step 05 - Safety Filtering", {}).get("result") or {}).get("safe_mhc2_count", 0)
    cloning_ok = (stages.get("Step 13 - In-Silico Cloning", {}).get("result") or {}).get("cloning_compatible", False)

    if antigenicity_gate and safe_mhc1 >= 3 and safe_mhc2 >= 3 and cloning_ok:
        candidacy = "STRONG MEV CANDIDATE"
    elif antigenicity_gate and (safe_mhc1 + safe_mhc2) >= 3:
        candidacy = "MODERATE MEV CANDIDATE"
    elif antigenicity_gate:
        candidacy = "WEAK MEV CANDIDATE"
    else:
        candidacy = "REJECTED — ANTIGENICITY GATE FAILED"

    return {
        "job_id": job_id,
        "uniprot_id": job.uniprot_id,
        "status": job.status,
        "created_at": job.created_at.isoformat() if job.created_at else None,
        "completed_at": job.completed_at.isoformat() if job.completed_at else None,
        "stages": stages,
        "flags": flags,
        "overall_mev_candidacy": candidacy,
    }


# --------------------------------------------------------------------------
# PDF Report
# --------------------------------------------------------------------------

# Keys in result_json that are suppressed from the PDF (still present in JSON)
_PDF_SKIP_KEYS = {"errors"}


class MEVReportPDF(FPDF):
    def header(self):
        self.set_font("Helvetica", "B", 14)
        self.cell(0, 10, "Multi-Epitope Vaccine (MEV) Design Report",
                  new_x="LMARGIN", new_y="NEXT", align="C")
        self.ln(4)

    def footer(self):
        self.set_y(-15)
        self.set_font("Helvetica", "I", 8)
        self.cell(0, 10, f"Page {self.page_no()}/{{nb}}", align="C")

    def section_title(self, title: str):
        self.set_font("Helvetica", "B", 11)
        self.set_fill_color(210, 230, 255)
        self.cell(0, 8, _sanitize(title), new_x="LMARGIN", new_y="NEXT", fill=True)
        self.ln(2)

    def key_value(self, key: str, value: str):
        self.set_font("Helvetica", "B", 9)
        self.cell(65, 6, _sanitize(key) + ":", new_x="RIGHT")
        self.set_font("Helvetica", "", 9)
        self.multi_cell(0, 6, _sanitize(str(value)), new_x="LMARGIN", new_y="NEXT")


def _draw_pipeline_flowchart(pdf: "MEVReportPDF", report: dict):
    """Draw a branching pipeline flowchart page with routing logic."""
    stages = report.get("stages", {})

    STEPS = [
        (1,  "Step 1: Sequence Retrieval",   lambda r: f"{r.get('length', '?')} AA"),
        (2,  "Step 2: Antigenicity (GATE)",  lambda r: f"Score: {r.get('antigenicity_score', '?')}"),
        (3,  "Step 3: Physicochemical",      lambda r: f"MW: {r.get('molecular_weight', '?')} Da"),
        (4,  "Step 4: Epitope Prediction",   lambda r: f"B:{r.get('bcell_count',0)} I:{r.get('mhc1_count',0)} II:{r.get('mhc2_count',0)}"),
        (5,  "Step 5: Safety Filtering",     lambda r: f"Safe B:{r.get('safe_bcell_count',0)} I:{r.get('safe_mhc1_count',0)} II:{r.get('safe_mhc2_count',0)}"),
        (6,  "Step 6: Population Coverage",  lambda r: f"W_Class1: {r.get('world_class1_coverage_pct','?')}%"),
        (7,  "Step 7: IFN-Gamma",            lambda r: f"{r.get('ifn_gamma_count', 0)} inducers"),
        (8,  "Step 8: PEP-FOLD 3D",          lambda r: f"{r.get('successful_models', 0)} structures"),
        (9,  "Step 9: Molecular Docking",    lambda r: f"{r.get('total_docked', 0)} complexes"),
        (10, "Step 10: Model Eval",          lambda r: f"{r.get('best_model_count', 0)} best models"),
        (11, "Step 11: Construct Design",    lambda r: f"{r.get('construct_length_aa', '?')} AA"),
        (12, "Step 12: Codon Opt",           lambda r: f"CAI: {r.get('cai_score','?')} GC: {r.get('gc_content_pct','?')}%"),
        (13, "Step 13: In-Silico Cloning",   lambda r: f"{r.get('cloning_status', '?')}"),
    ]

    pdf.add_page()
    pdf.set_auto_page_break(auto=False)

    pdf.set_font("Helvetica", "B", 13)
    pdf.cell(0, 10, "Pipeline Execution Flowchart", new_x="LMARGIN", new_y="NEXT", align="C")
    pdf.ln(1)

    page_w  = pdf.w - pdf.l_margin - pdf.r_margin
    center_x = pdf.l_margin + page_w / 2

    box_w  = 72
    box_h  = 13.0
    gap_y  = 5.5
    start_y = pdf.get_y() + 2

    def get_coords(step_num):
        if step_num <= 6:
            row = step_num - 1
            x = center_x - box_w / 2
        elif step_num == 7:
            row = 6
            x = 14
        elif step_num == 8:
            row = 6
            x = 124
        else:
            row = step_num - 2
            x = center_x - box_w / 2
        y = start_y + row * (box_h + gap_y)
        return x, y

    def draw_box(step_num, label, extractor):
        x, y = get_coords(step_num)
        step_key = [k for k in stages.keys() if k.startswith(f"Step {step_num:02d}")]
        if not step_key:
            stage_data = {}
            passed = None
            result = {}
            error  = None
        else:
            stage_data = stages[step_key[0]]
            passed = stage_data.get("passed", None)
            result = stage_data.get("result") or {}
            error  = stage_data.get("error")

        if passed is True:
            pdf.set_fill_color(200, 240, 200)
            badge, br, bg, bb = "TRUE", 34, 139, 34
        elif passed is False:
            pdf.set_fill_color(255, 210, 210)
            badge, br, bg, bb = "FALSE", 200, 30, 30
        else:
            pdf.set_fill_color(230, 230, 230)
            badge, br, bg, bb = "N/A", 120, 120, 120

        pdf.set_draw_color(100, 100, 100)
        pdf.set_line_width(0.3)
        pdf.rect(x, y, box_w, box_h, "DF")

        pdf.set_xy(x + 1, y + 1.5)
        pdf.set_font("Helvetica", "B", 7.5)
        pdf.set_text_color(30, 30, 30)
        pdf.cell(box_w - 20, 4, _sanitize(label))

        bx = x + box_w - 15
        by = y + 1.5
        pdf.set_fill_color(br, bg, bb)
        pdf.rect(bx, by, 13, 4, "DF")
        pdf.set_xy(bx, by)
        pdf.set_font("Helvetica", "B", 6)
        pdf.set_text_color(255, 255, 255)
        pdf.cell(13, 4, badge, align="C")

        try:
            res_txt = extractor(result)
        except Exception:
            res_txt = "-"
        pdf.set_xy(x + 2, y + 6)
        pdf.set_font("Helvetica", "", 7)
        pdf.set_text_color(60, 60, 60)
        pdf.cell(box_w - 4, 4, _sanitize(res_txt))

        if error:
            pdf.set_xy(x + 2, y + 9.5)
            pdf.set_font("Helvetica", "I", 6)
            pdf.set_text_color(180, 0, 0)
            pdf.cell(box_w - 4, 3, _sanitize(error[:50]))

    for step_num, label, ext in STEPS:
        draw_box(step_num, label, ext)

    def draw_straight_arrow(step_from, step_to):
        x1, y1 = get_coords(step_from)
        x2, y2 = get_coords(step_to)
        x_mid   = x1 + box_w / 2
        y_start = y1 + box_h
        y_end   = y2
        pdf.set_draw_color(100, 100, 100)
        pdf.set_line_width(0.4)
        pdf.line(x_mid, y_start, x_mid, y_end)
        pdf.line(x_mid, y_end, x_mid - 1.5, y_end - 2)
        pdf.line(x_mid, y_end, x_mid + 1.5, y_end - 2)

    for i in [1, 2, 3, 4, 5]:
        draw_straight_arrow(i, i + 1)
    for i in [9, 10, 11, 12]:
        draw_straight_arrow(i, i + 1)

    x6, y6 = get_coords(6)
    x7, y7 = get_coords(7)
    x8, y8 = get_coords(8)

    mid_y_6 = y6 + box_h + gap_y / 2
    pdf.set_draw_color(100, 100, 100)
    pdf.line(center_x, y6 + box_h, center_x, mid_y_6)

    pdf.line(center_x, mid_y_6, x7 + box_w / 2, mid_y_6)
    pdf.line(x7 + box_w / 2, mid_y_6, x7 + box_w / 2, y7)
    pdf.line(x7 + box_w / 2, y7, x7 + box_w / 2 - 1.5, y7 - 2)
    pdf.line(x7 + box_w / 2, y7, x7 + box_w / 2 + 1.5, y7 - 2)
    pdf.set_xy(x7 + box_w / 2 + 2, mid_y_6 - 3.5)
    pdf.set_text_color(100, 100, 100)
    pdf.set_font("Helvetica", "B", 6)
    pdf.cell(10, 3, "MHC-II")

    pdf.line(center_x, mid_y_6, x8 + box_w / 2, mid_y_6)
    pdf.line(x8 + box_w / 2, mid_y_6, x8 + box_w / 2, y8)
    pdf.line(x8 + box_w / 2, y8, x8 + box_w / 2 - 1.5, y8 - 2)
    pdf.line(x8 + box_w / 2, y8, x8 + box_w / 2 + 1.5, y8 - 2)
    pdf.set_xy(x8 + box_w / 2 - 12, mid_y_6 - 3.5)
    pdf.cell(10, 3, "MHC-I")

    x9, y9   = get_coords(9)
    mid_y_78 = y7 + box_h + gap_y / 2

    pdf.line(x7 + box_w / 2, y7 + box_h, x7 + box_w / 2, mid_y_78)
    pdf.line(x8 + box_w / 2, y8 + box_h, x8 + box_w / 2, mid_y_78)
    pdf.line(x7 + box_w / 2, mid_y_78, x8 + box_w / 2, mid_y_78)
    pdf.line(center_x, mid_y_78, center_x, y9)
    pdf.line(center_x, y9, center_x - 1.5, y9 - 2)
    pdf.line(center_x, y9, center_x + 1.5, y9 - 2)

    pdf.set_draw_color(220, 100, 0)
    pdf.set_line_width(0.5)

    bypass_start_x = x6 + box_w
    bypass_start_y = y6 + box_h / 2
    x11, y11       = get_coords(11)
    bypass_end_x   = x11 + box_w
    bypass_end_y   = y11 + box_h / 2
    far_x          = pdf.w - pdf.r_margin - 3

    pdf.line(bypass_start_x, bypass_start_y, far_x, bypass_start_y)
    pdf.line(far_x, bypass_start_y, far_x, bypass_end_y)
    pdf.line(far_x, bypass_end_y, bypass_end_x, bypass_end_y)
    pdf.line(bypass_end_x, bypass_end_y, bypass_end_x + 2, bypass_end_y - 1.5)
    pdf.line(bypass_end_x, bypass_end_y, bypass_end_x + 2, bypass_end_y + 1.5)

    pdf.set_text_color(220, 100, 0)
    pdf.set_font("Helvetica", "B", 7)
    pdf.set_xy(far_x - 11, bypass_start_y - 4)
    pdf.cell(10, 4, "B-cell")

    candidacy = report.get("overall_mev_candidacy", "N/A")
    cand_w    = 110
    cand_x    = pdf.l_margin + (page_w - cand_w) / 2
    cand_y    = start_y + 12 * (box_h + gap_y) + 4
    pdf.set_fill_color(50, 50, 120)
    pdf.rect(cand_x, cand_y, cand_w, 12, "DF")
    pdf.set_xy(cand_x, cand_y + 2.5)
    pdf.set_font("Helvetica", "B", 9)
    pdf.set_text_color(255, 255, 255)
    pdf.cell(cand_w, 7, _sanitize(candidacy), align="C")

    pdf.set_text_color(0, 0, 0)
    pdf.set_draw_color(0, 0, 0)
    pdf.set_auto_page_break(auto=True, margin=15)
    pdf.add_page()


async def generate_pdf_report(job_id: str, db: AsyncSession) -> bytes:
    report = await generate_json_report(job_id, db)

    pdf = MEVReportPDF()
    pdf.alias_nb_pages()
    pdf.add_page()

    # Metadata
    pdf.section_title("Protein Metadata")
    pdf.key_value("Job ID",          report["job_id"])
    pdf.key_value("UniProt ID",      report["uniprot_id"])
    pdf.key_value("Status",          report["status"])
    pdf.key_value("Created",         report.get("created_at", ""))
    pdf.key_value("Completed",       report.get("completed_at", ""))
    seq_data = (report.get("stages", {})
                .get("Step 01 - Sequence Retrieval", {})
                .get("result", {})) or {}
    pdf.key_value("Protein Name",    seq_data.get("protein_name", "N/A"))
    pdf.key_value("Organism",        seq_data.get("organism", "N/A"))
    pdf.key_value("Sequence Length", str(seq_data.get("length", "N/A")) + " AA")
    pdf.ln(3)

    # ── Metadata ─────────────────────────────────────────────────────────────
    pdf.section_title("Protein Metadata")
    pdf.key_value("Job ID",          report["job_id"])
    pdf.key_value("UniProt ID",      report["uniprot_id"])
    pdf.key_value("Status",          report["status"])
    pdf.key_value("Created",         report.get("created_at", ""))
    pdf.key_value("Completed",       report.get("completed_at", ""))
    seq_data = (report.get("stages", {})
                .get("Step 01 - Sequence Retrieval", {})
                .get("result", {})) or {}
    pdf.key_value("Protein Name",    seq_data.get("protein_name", "N/A"))
    pdf.key_value("Organism",        seq_data.get("organism", "N/A"))
    pdf.key_value("Sequence Length", str(seq_data.get("length", "N/A")) + " AA")
    pdf.ln(3)

    # ── Disclaimer (added after metadata) ────────────────────────────────────
    pdf.set_fill_color(255, 243, 205)           # soft amber background
    pdf.set_draw_color(200, 150, 0)             # amber border
    pdf.set_line_width(0.5)
    disclaimer_x = pdf.l_margin
    disclaimer_w = pdf.w - pdf.l_margin - pdf.r_margin

    # measure height needed before drawing
    pdf.set_font("Helvetica", "B", 9)
    title_h = 6
    pdf.set_font("Helvetica", "", 8.5)
    # draw the box — estimate ~22mm height for wrapped text
    pdf.rect(disclaimer_x, pdf.get_y(), disclaimer_w, 26, "DF")

    # title line
    pdf.set_xy(disclaimer_x + 3, pdf.get_y() + 2)
    pdf.set_font("Helvetica", "B", 9)
    pdf.set_text_color(160, 100, 0)
    pdf.cell(disclaimer_w - 6, title_h,
             "Epitope Selection Criteria & Disclaimer",
             new_x="LMARGIN", new_y="NEXT")

    # body text — multi_cell wraps automatically
    pdf.set_x(disclaimer_x + 3)
    pdf.set_font("Helvetica", "", 8)
    pdf.set_text_color(80, 60, 0)
    pdf.multi_cell(disclaimer_w - 6, 4.5,
                   _sanitize(EPITOPE_DISCLAIMER),
                   new_x="LMARGIN", new_y="NEXT")

    # reset colours
    pdf.set_text_color(0, 0, 0)
    pdf.set_draw_color(0, 0, 0)
    pdf.set_line_width(0.2)
    pdf.ln(4)
    # ── end disclaimer ────────────────────────────────────────────────────────

    # Flowchart page
    _draw_pipeline_flowchart(pdf, report)

    # Per-step sections
    for stage_name in sorted(report.get("stages", {}).keys(), key=_get_stage_num):
        stage      = report["stages"][stage_name]
        passed_str = "PASSED" if stage.get("passed") else "FLAGGED"
        pdf.section_title(f"{stage_name}  [{passed_str}]")

        result = stage.get("result") or {}
        if isinstance(result, dict):
            for k, v in result.items():
                if k in _PDF_SKIP_KEYS:          # ← skip 'errors' in PDF only
                    continue
                if isinstance(v, str) and len(v) > 300:
                    pdf.key_value(k, f"({len(v)} chars - see JSON report)")
                elif isinstance(v, (dict, list)):
                    pretty = json.dumps(v, indent=2)
                    pdf.key_value(k, pretty if len(pretty) <= 400 else "(complex - see JSON report)")
                else:
                    pdf.key_value(k, str(v))

        if stage.get("error"):
            pdf.key_value("Error", stage["error"])
        pdf.ln(2)

    # Overall assessment
    pdf.section_title("Overall MEV Candidacy Assessment")
    pdf.set_font("Helvetica", "B", 12)
    pdf.cell(0, 10, report["overall_mev_candidacy"],
             new_x="LMARGIN", new_y="NEXT", align="C")
    pdf.ln(3)

    if report.get("flags"):
        pdf.section_title("Flagged Steps")
        for flag in report["flags"]:
            pdf.set_font("Helvetica", "", 9)
            pdf.cell(0, 6, f"  - {flag}", new_x="LMARGIN", new_y="NEXT")

    return bytes(pdf.output())


# --------------------------------------------------------------------------
# Entry-point
# --------------------------------------------------------------------------

async def run_report_generation(job_id: str, db: AsyncSession):
    import shutil
    from datetime import datetime

    reports_dir = os.path.join(os.path.dirname(__file__), "..", "data", "reports")
    os.makedirs(reports_dir, exist_ok=True)

    backend_dir = os.path.join(os.path.dirname(__file__), "..")
    output_dir  = os.path.abspath(os.path.join(backend_dir, "output"))
    os.makedirs(output_dir, exist_ok=True)

    report_json = await generate_json_report(job_id, db)
    uniprot_id  = report_json.get("uniprot_id", job_id)
    ts          = datetime.utcnow().strftime("%Y%m%d_%H%M%S")
    base_name   = f"{uniprot_id}_{ts}"

    # JSON — internal + output/
    internal_json = os.path.join(reports_dir, f"{job_id}.json")
    with open(internal_json, "w", encoding="utf-8") as f:
        json.dump(report_json, f, indent=2)
    output_json = os.path.join(output_dir, f"{base_name}_report.json")
    shutil.copy2(internal_json, output_json)

    # PDF — internal + output/
    pdf_bytes    = await generate_pdf_report(job_id, db)
    internal_pdf = os.path.join(reports_dir, f"{job_id}.pdf")
    with open(internal_pdf, "wb") as f:
        f.write(pdf_bytes)
    output_pdf = os.path.join(output_dir, f"{base_name}_report.pdf")
    shutil.copy2(internal_pdf, output_pdf)

    import logging as _logging
    _logging.getLogger(__name__).info(
        "Reports written to output/%s_report.[json|pdf]", base_name
    )
    return {"json_path": output_json, "pdf_path": output_pdf}