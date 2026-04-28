"""Execute verify_v2.ipynb as a script (no jupyter available)."""
import json, sys
nb = json.load(open(sys.argv[1] if len(sys.argv) > 1 else "verify_v2.ipynb"))
import matplotlib
matplotlib.use("Agg")  # no display
combined = []
for c in nb["cells"]:
    src = c["source"] if isinstance(c["source"], str) else "".join(c["source"])
    if c["cell_type"] == "markdown":
        combined.append("# === " + src.split("\n")[0].lstrip("# ").strip() + " ===")
    else:
        combined.append(src)
script = "\n".join(combined)
exec(compile(script, "verify_v2", "exec"), {"__name__": "__main__"})
