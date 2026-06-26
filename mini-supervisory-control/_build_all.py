import os, sys
sys.stdout.reconfigure(encoding="utf-8")

BASE = "."

def write_file(relpath, content):
    path = os.path.join(BASE, relpath)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write(content)
    lines = content.count(chr(10))
    print(f"Wrote {relpath}: {lines} lines")

print("Builder ready. Will write C source files next.")
