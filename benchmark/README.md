# Confronto locale R / Julia

Per confrontare il solo calcolo delle distanze del filtro divergente con
l'implementazione precedente:

```powershell
julia --project=. benchmark/divergent_filter.jl
```

Il benchmark verifica prima l'uguaglianza delle distanze, poi misura una chiamata
dopo il riscaldamento, riportando tempo e byte allocati. Include 400 sequenze
tutte diverse e 1.000 sequenze con 50 sequenze distinte, lunghe 220 posizioni.
Sono dati sintetici: le misure escludono lettura/scrittura FASTA, checkpoint e
report e non rappresentano il tempo dell'intera pipeline.
Le misure locali sono in `benchmark/DIVERGENT_FILTER_RESULTS.md`.

Eseguire in sequenza, evitando altri carichi sulla CPU:

```powershell
julia --project=. benchmark/compare.jl
& 'C:\Program Files\R\R-4.5.0\bin\Rscript.exe' benchmark/compare.R results/benchmark reference/aasv
```

Il primo script crea fixture identiche per entrambi i linguaggi sotto
`results/benchmark` e misura Julia; il secondo legge le fixture e misura il
codice R originale al commit indicato nel README del progetto.

Scenari: allineamento multiplo di 20 e 80 proteine sintetiche di 150 residui,
e analisi funzionale + classificazione di 10 ASV (390 nt) con 20 riferimenti
su specie, genere e famiglia. I profili della pipeline sono già allineati e
identici per entrambe le versioni; il tempo della pipeline include FASTA/Excel.
Non include download BOLD o VSEARCH.

Le ASV sono frammenti delle posizioni 13–142 delle proteine di riferimento,
che hanno 150 residui. Nel sorgente R originale, una query con zero gap
iniziali provoca un errore nel sottoinsieme `audit_results[-seq_len(0)]`;
il confronto della pipeline usa frammenti per evitare quel caso senza
modificare il codice R. I risultati non rappresentano quindi tutte le
possibili query. `--pipeline-only` permette di ripetere la pipeline conservando
le misure degli allineamenti, sia nello script Julia che come terzo argomento
dello script R.

Ogni scenario contiene la prima chiamata e tre chiamate successive; riportare
la mediana delle tre chiamate successive. La prima chiamata Julia include
eventuale compilazione delle funzioni, ma **non** l'avvio del processo o il
caricamento dei pacchetti, che avvengono prima della misura.

Il comando sopra avvia Julia con un solo thread, quindi DECIPHER usa
`processors=1` e ProfileAlignments usa un solo thread. Per misurare anche la
modalità parallela Julia, eseguire separatamente
`julia -t auto --project=. benchmark/compare.jl`. Si
confrontano le impostazioni normali dei due motori, che usano algoritmi,
pesi e parametri differenti. I tempi misurano queste implementazioni su
questi input, non le prestazioni intrinseche di R contro Julia. Non sono
benchmark di accuratezza biologica.

Risultati grezzi: `results/benchmark/julia_times.csv` e `r_times.csv`.
Il riepilogo dell'ultima misura locale è in `benchmark/RESULTS.md`.
