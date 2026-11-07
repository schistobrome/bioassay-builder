"""
M2 setup check for the Assay Strategy Planner.

Confirms that the language model used for matching assay descriptions to
assay-approach names works on your computer. Nothing is sent anywhere except
the one-time model download from huggingface.co.

Run it with the Python inside the .venv folder:
    .\\.venv\\Scripts\\python.exe m2_setup_check.py

It saves its output as m2_setup_results.txt next to this file.
Paste that file's contents back into the chat.
"""

import os
import platform
import sys
import time

LOG = []


def say(text=""):
    print(text)
    LOG.append(text)


MODEL = "sentence-transformers/all-MiniLM-L6-v2"

# Real-looking descriptions (taken from the WRN data) and approach names.
DESCRIPTIONS = [
    "Inhibition of WRN ATPase activity measured by ADP-Glo kinase assay based luminescence",
    "Inhibition of WRN helicase using fluor/quencher-labeled forked DNA, fluorescence based assay",
    "Binding affinity to WRN assessed as dissociation constant by SPR analysis",
    "Inhibition of WRN helicase activity by polyacrylamide gel electrophoresis",
    "Binding affinity to WRN assessed as melting temperature by SYPRO Orange thermal shift assay",
]
APPROACHES = [
    "ATPase activity, luminescent ADP detection",
    "Helicase unwinding, fluorogenic quenched DNA substrate",
    "Surface plasmon resonance binding kinetics",
    "Gel-based DNA unwinding",
    "Thermal shift / differential scanning fluorimetry",
]


def main():
    say("M2 setup check")
    say("Python: %s on %s" % (sys.version.split()[0], platform.platform()))
    try:
        import numpy  # noqa: F401
        import torch
        from sentence_transformers import SentenceTransformer
    except Exception as error:
        say("FAIL: could not load the packages: %s: %s" % (type(error).__name__, error))
        say("Make sure you ran this with the .venv Python (see the instructions at the top).")
        save()
        return
    say("torch %s, sentence-transformers loaded." % torch.__version__)

    say("Loading model %s (first run downloads about 90 MB)..." % MODEL)
    started = time.time()
    try:
        model = SentenceTransformer(MODEL)
    except Exception as error:
        say("FAIL: could not load the model: %s: %s" % (type(error).__name__, error))
        save()
        return
    say("Model ready in %.1f s." % (time.time() - started))

    started = time.time()
    d = model.encode(DESCRIPTIONS, normalize_embeddings=True)
    a = model.encode(APPROACHES, normalize_embeddings=True)
    scores = d @ a.T
    say("Encoding took %.2f s." % (time.time() - started))

    say("")
    say("For each description, the best-matching approach (similarity 0 to 1):")
    right = 0
    for i, text in enumerate(DESCRIPTIONS):
        order = list(scores[i].argsort()[::-1])
        best = order[0]
        mark = "ok " if best == i else "?? "
        right += best == i
        say("  %s%.2f  %s" % (mark, scores[i][best], APPROACHES[best]))
        say("        for: %s" % text[:80])
        say("        runner-up: %.2f  %s" % (scores[i][order[1]], APPROACHES[order[1]]))
    say("")
    say("%d of %d matched the approach I expected." % (right, len(DESCRIPTIONS)))
    save()


def save():
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "m2_setup_results.txt")
    try:
        with open(path, "w", encoding="utf-8") as handle:
            handle.write("\n".join(LOG) + "\n")
        print("\nSaved to: %s\nPlease paste its contents back into the chat." % path)
    except OSError as error:
        print("\nCould not save the file (%s). Please copy the text above." % error)


if __name__ == "__main__":
    main()
