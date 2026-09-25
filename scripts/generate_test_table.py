"""Appendix table: unit tests grouped by module, with what they pin."""
import re, ast, pathlib
rows = []
for p in sorted(pathlib.Path("tests").glob("test_*.py")):
    tree = ast.parse(p.read_text())
    tests = [n.name for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name.startswith("test_")]
    mod = p.stem.replace("test_", "").replace("_", "\\_")
    rows.append((mod, len(tests)))
with open("paper/tests_table.tex", "w") as f:
    f.write("\\begin{longtable}{lc}\\toprule\nTest module & Tests \\\\\n\\midrule\n\\endfirsthead\n\\toprule Module & Tests \\\\ \\midrule \\endhead\n")
    for m, n in rows:
        f.write(f"\\texttt{{{m}}} & {n} \\\\\n")
    f.write(f"\\midrule Total & {sum(n for _, n in rows)} \\\\\n")
    f.write("\\bottomrule\n\\end{longtable}\n")
print(rows, "total", sum(n for _, n in rows))
