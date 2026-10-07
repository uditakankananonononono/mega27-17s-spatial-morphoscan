# Render status (updated 2026-10-08)

paper/main.pdf is rendered from paper/main.tex by the render_paper GitHub Actions workflow (xelatex, texlive container, separate commit job). The PDF is current: 52pp, includes P3-SN1 (6.19) and P5-PLS1 (6.20), rendered at 513adb8.

History: runs for 7bc26bd and 20521bd stalled in the apt "Install TeX" step (cause unknown, no log access), so the workflow was changed to render inside the texlive/texlive container with no apt step. The container run took about a minute.

Known open item: paper wording in Section 6.23 says G-B1 "passes on the corrected-geometry rerun". Under docs/GATE_LEDGER.md, G-B1 stays FAIL (v2 is post-hoc, descriptive only). The paper text has not been reconciled yet.
