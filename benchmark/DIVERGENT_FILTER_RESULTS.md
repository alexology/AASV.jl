# Filtro divergente: misure locali del 15 ottobre 2026

Windows, Julia 1.12.2. Confronto del calcolo delle distanze prima e dopo
l'ottimizzazione, su sequenze sintetiche di 220 posizioni. Una chiamata dopo
il riscaldamento, con uguaglianza esatta dei risultati verificata prima della misura.

| Sequenze | Distinte | Prima (s) | Dopo (s) | Accelerazione | Allocazioni prima | Allocazioni dopo |
|---|---:|---:|---:|---:|---:|---:|
| 400 | 400 | 0.1255123 | 0.0200584 | 6.26× | 452,541,550 byte | 58,148 byte |
| 1.000 | 50 | 0.6708134 | 0.0003065 | 2.189× | 2,831,735,566 byte | 36,467 byte |

Nel secondo scenario ogni sequenza distinta compare 20 volte: ogni record
ha quindi un duplicato a distanza zero e non servono scansioni delle coppie.
Questo e un caso particolarmente favorevole, non una previsione sui dati reali.
Le allocazioni sono il totale cumulativo, non il picco di memoria residente.
Tempi esclusi: avvio, caricamento pacchetti, lettura/scrittura FASTA, checkpoint
e report. Riproduzione: `julia --project=. benchmark/divergent_filter.jl`.
