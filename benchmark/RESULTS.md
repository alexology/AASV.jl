# Misure locali del 14 ottobre 2026

Windows, Julia 1.12.2 e R 4.5.0, sullo stesso computer. DECIPHER con
`processors=1`; nessuna esecuzione contemporanea dei due benchmark.
Mediana di tre chiamate successive alla prima, in secondi:

| Operazione | R | Julia |
|---|---:|---:|
| Allineamento di 20 proteine × 150 aa | 0.280 | 0.089 |
| Allineamento di 80 proteine × 150 aa | 1.000 | 0.527 |
| Analisi funzionale + classificazione di 10 ASV su tre livelli | 11.600 | 0.284 |

La pipeline usa 20 riferimenti già allineati e ASV di 390 nt. Include
input/output FASTA ed Excel; esclude BOLD e VSEARCH. Su questi input,
FI e classi coincidono per tutte le 30 combinazioni ASV/livello.

Prima chiamata della pipeline, esclusi avvio del processo e importazione dei
pacchetti: Julia 13.838 s, R 12.660 s. La compilazione iniziale Julia pesa
quindi molto rispetto alle chiamate successive.

Su questi input, il motore Julia è circa 3,1 volte più veloce con 20 proteine
e 1,9 volte con 80; la pipeline funzionale e di classificazione è anch'essa
più veloce a regime. Questi dati non stabiliscono quale intero workflow su un
dataset reale sia più veloce: i motori di allineamento sono differenti e i
dati qui sono sintetici.

Il profiling della versione precedente attribuiva quasi tutto il tempo del
caso da 80 sequenze al raffinamento: 32,66 s su circa 33,9 s. Il punteggio ora
scorre direttamente le stringhe ASCII e il raffinamento valuta soltanto le
coppie attraversate da ciascuna bipartizione, mantenendo in cache le altre.
La mediana a 80 sequenze è scesa da 33,913 s a 0,527 s. La selezione UPGMA usa
inoltre una coda minima delle distanze, con costo O(N² log N) anziché O(N³).

Il codice R originale si interrompe su una query priva di gap iniziali;
la pipeline è misurata su frammenti più corti dei riferimenti per evitare
quel caso senza modificare il sorgente originale. Vedere README.md.

Script riproducibili: `compare.jl`, `compare.R`, `summarize.jl`.
Misure grezze: `../results/benchmark/julia_times.csv`, `r_times.csv` e
`summary.csv` nella stessa cartella risultati.
