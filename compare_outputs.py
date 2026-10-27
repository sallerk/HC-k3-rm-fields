"""
compare_outputs.py -- compare the result files of a rerun with the reference outputs.

Usage:  python -I compare_outputs.py REFERENCE_DIR RERUN_DIR
Every .txt file in REFERENCE_DIR is compared with the file of the same name in RERUN_DIR, line by line, with
trailing white space and the line-end convention ignored, and with running times (a number followed by "s" or
" s" just before a closing parenthesis or bracket, as in "(0.1 s)" or "(exact count 2.2 s)") replaced by a
placeholder; such times occur only in c20_ttv_exact_traces.txt and c22_ttv_ntt_traces.txt.  The random tests of
the scripts use fixed seeds, so nothing else is set aside.  Each file is reported as "equal" or "DIFFERENT" (with
the first line that differs).  Expected last line: "ALL EQUAL".
"""

import os
import re
import sys

TIME = re.compile(r"\d+(?:\.\d+)? ?s(?=[\])])")
KINDS = (".txt",)


def text_lines(path):
    with open(path, encoding="utf-8") as f:
        return [TIME.sub("<time>", line.rstrip()) for line in f]


def first_difference(A, B):
    for i, (a, b) in enumerate(zip(A, B)):
        if a != b:
            return "line %d differs" % (i + 1)
    return "lengths %d and %d" % (len(A), len(B))


def main():
    ref, new = sys.argv[1], sys.argv[2]
    allequal = True
    names = sorted(n for n in os.listdir(ref) if n.endswith(KINDS))
    for name in names:
        a = os.path.join(ref, name)
        b = os.path.join(new, name)
        if not os.path.exists(b):
            print("%-24s missing in rerun" % name)
            allequal = False
            continue
        A, B = text_lines(a), text_lines(b)
        eq = A == B
        note = "" if eq else "  (%s)" % first_difference(A, B)
        print("%-24s %s%s" % (name, "equal" if eq else "DIFFERENT", note))
        allequal &= eq
    extra = sorted(n for n in os.listdir(new) if n.endswith(KINDS) and n not in names)
    for name in extra:
        print("%-24s only in rerun" % name)
    print("ALL EQUAL" if allequal and not extra else "DIFFERENCES FOUND")


if __name__ == "__main__":
    main()
