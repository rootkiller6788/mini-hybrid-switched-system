#!/usr/bin/env python3

"II" -- generate all files for mini-impulsive-system."i"
import os

BASE = os.path.dirname(os.path.abspath(__file__))

def w(relpath, content):
    path = os.path.join(BASE, relpath)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w', encoding='utf-8') as f:
        f.write(content)
    nc = content.count('\n') + 1
    print(f"  OK: {relpath} ({nc} lines)")

print("Starting build ...")

==========
== LIBERIES ==
==========