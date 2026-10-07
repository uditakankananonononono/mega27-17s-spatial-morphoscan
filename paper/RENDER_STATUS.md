# Render status (updated 2026-10-08)

paper/main.pdf is rendered from paper/main.tex by the render_paper GitHub Actions workflow (xelatex).

**STALE: the committed PDF (51pp) does NOT include the P5-PLS1 subsection added to main.tex in 7bc26bd.** CI runs for 7bc26bd (cancelled after about 25 min) and 20521bd (about 30 min, still running at last check) both stalled in the "Install TeX" apt step. Cause unknown; no log access. The PDF includes P3-SN1 (6.19) but not P5-PLS1. The source (main.tex) is current. Rerun the workflow (workflow_dispatch, or touch paper/main.tex) when the runner recovers.
