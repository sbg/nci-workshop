# dna2protein

This example consists of two CWL tools (Transcribe, Translate) and one workflow (dna2protein), as well as the Dockerfile you can use to build a container with the tools inside.

Transcribe takes a TXT file containing a DNA Sequence as an input and produces a TXT file with an mRNA Sequence. Translate takes an mRNA Sequence, identifies the first ORF, and produces a TXT file with a peptide sequence as an output. dna2protein consists of *input -> Transcribe > Translate --> output*.

To run locally with Python:

1. `python scripts/transcribe.py -d data/input.txt > data/rna.txt`
2. `python scripts/translate.py -r data/rna.txt > data/protein.txt`

OR

`python scripts/transcribe.py -d data/input.txt | python scripts/translate.py > data/protein.txt`

To run locally with rabix, pull this folder into your local rabix directory, then enter:

`./rabix.sh -a dna2protein.cwl.json -i input_transcribe.json`

All CWL apps in this example were written using the Seven Bridges platform.