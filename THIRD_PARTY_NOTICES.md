# Grantham distance table

`data/grantham.csv` is copied from
https://github.com/patterninstitute/grantham/blob/main/data-raw/grantham_distance_matrix.csv
and contains the original Grantham (1974) amino-acid distance table.

Source package copyright (2021): Ramiro Magno. Contributors include Isabel
Duarte, Ana-Teresa Maia, University of Algarve, CINTESIS and Pattern Institute.

MIT License

Copyright (c) 2021 Ramiro Magno

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

# Genetic codes

`data/genetic_codes.tsv` records the NCBI genetic code mappings exposed by
Biostrings `GENETIC_CODE_TABLE` in the local R reference environment, in
TCAG codon order. Source: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
These biological data are included as a static snapshot; no Biostrings code
is required at runtime.
