# AASV.jl

Implementazione Julia della pipeline [aasv](https://github.com/alexology/aasv)
per l'analisi funzionale delle ASV COI. Le formule scientifiche sono compatibili
con la versione R 0.1.0, commit `9960347948da3142714deb656fcf38124b4e38d8`;
l'architettura usa tipi, dispatch e convenzioni Julia.

Richiede Julia 1.10 o successiva. Non richiede R, Biostrings o DECIPHER.

```powershell
julia --project=. setup.jl
julia -t auto --project=. test/runtests.jl
julia -t auto --project=. examples/offline_pipeline.jl
```

L'esempio locale crea FASTA, tabelle dei codoni, report Excel e
`results/demo/analysis/classification.csv`.

Per installare la versione pubblicata su GitHub:

```julia
using Pkg
Pkg.develop(url="https://github.com/alexology/ProfileAlignments.jl")
Pkg.develop(url="https://github.com/alexology/AASV.jl")
using AASV, BioSequences

record = Sequence("COI", dna"ATGCCAGCA")
protein = best_translation(record.sequence)
```

## Pipeline con dati reali

```julia
using AASV, DataFrames, CSV

taxonomy = CSV.read("taxonomy.csv", DataFrame)
# Colonne obbligatorie: family, genus, species.
align_species_seq!(taxonomy; alignment_dir="references",
                   min_length=500, max_length=1000, genetic_code="5")

for family in unique(skipmissing(taxonomy.family))
    problems = check_alignment(family; alignment_dir="references", tax_lev=SpeciesRank)
    improve_alignment!(problems; log_dir="logs")
    dir = joinpath("references", "alignments", family, "species")
    paths = filter(p -> endswith(p, ".fasta"), readdir(dir; join=true))
    filter_divergent_seq!(paths; log_dir="logs", max_distance=0.10)
end

align_genus_seq!(; alignment_dir="references")
align_family_seq!(; alignment_dir="references")

asvs = CSV.read("asvs.csv", DataFrame)
# Colonne obbligatorie: family, genus, species, ASV (sequenza DNA).
# Colonne numeriche aggiuntive: abbondanze nei campioni.
revised = asv_functional_structure!(asvs;
    alignment_dir="references", output_dir="analysis", genetic_code="5")
verdicts = classify_asv("analysis")
CSV.write("classification.csv", verdicts)
```

Le opzioni Julia si passano dopo `;`. I valori statistici mancanti usano
`missing`; l'assenza di un risultato o di un'opzione usa `nothing`. CSV, Excel
e `DataFrame` sono adattatori di ingresso/uscita, mentre il nucleo usa strutture
concrete. Le sequenze sono vettori di `Sequence(id, sequence)`, con nomi univoci.
Gli argomenti `alignment_dir` e `output_dir` accettano percorsi Windows, Linux
e macOS. Usare `overwrite=true` per ricostruire i riferimenti già presenti.

## API Julia

`Sequence` è lo stesso tipo esportato da `ProfileAlignments`: gli allineamenti
non copiano più ogni record tra due rappresentazioni del dominio. La traduzione
tipizzata restituisce un solo tipo di risultato:

```julia
translation = find_translation(dna"ATGCCAGCA")
if !isnothing(translation)
    translation.amino_acids
    translation.frame
end
```

`assess_structure` e `assess_asv` restituiscono rispettivamente
`StructureAssessment` e `ASVAssessment`; le differenze amminoacidiche sono
`MutationObservation`. Queste API compongono il calcolo in memoria senza creare
tabelle intermedie. `check_asv_structure` e `calculate_asv_distance` conservano
l'output tabellare per i report esistenti.

Le opzioni tassonomiche accettano `SpeciesRank`, `GenusRank` e `FamilyRank`.
Una politica FI riutilizzabile si rappresenta con un valore immutabile:

```julia
policy = ClassificationConfig(quality_threshold=0.995,
                              include_hydro=false)
verdicts = classify_asv("analysis", policy)
```

Le funzioni che modificano file o cartelle terminano con `!`:
`align_species_seq!`, `align_genus_seq!`, `align_family_seq!`, `re_align!`,
`improve_alignment!`, `filter_divergent_seq!` e
`asv_functional_structure!`. I nomi precedenti restano come adattatori per la
compatibilità del codice esistente.

`filter_divergent_seq!` mostra l'avanzamento dei file e, ogni 5 secondi,
dei confronti tra sequenze nel file corrente. Per visualizzarlo anche con
`print_messages=false`, passare `show_progress=true`; la frequenza si regola
con `progress_interval` (in secondi).
Il calcolo delle distanze usa direttamente i byte delle sequenze ASCII, senza
copie delle righe per ogni coppia, e precalcola gli estremi senza gap. Le sequenze
identiche condividono i confronti, mantenendo tutti gli identificativi originali.
Le distanze e le regole di eliminazione restano esatte, inclusi i casi senza
sovrapposizione o con soli gap. Il guadagno aumenta in presenza di duplicati;
con sequenze tutte diverse il numero massimo di confronti resta quadratico.

Per elaborare tutte le famiglie senza scrivere cicli nello script:

```julia
species_bad = check_alignments!(; alignment_dir=reference_dir,
    tax_lev=SpeciesRank)
isempty(species_bad) || improve_alignment!(species_bad)

# Dopo la correzione vengono ricontrollate solo le famiglie modificate.
species_still_bad = check_alignments!(; alignment_dir=reference_dir,
    tax_lev=SpeciesRank)

filter_divergent_references!(; alignment_dir=reference_dir,
    tax_lev=SpeciesRank, max_distance=0.10)
```

Entrambe le funzioni scoprono le famiglie e saltano le cartelle prive di FASTA.
`families=["FamilyA", "FamilyB"]` limita la selezione. `resume=true` (predefinito)
riusa i risultati salvati in `alignment_dir/checkpoints` solo se nomi e contenuti
SHA-256 dei file, parametri e versione dell'algoritmo coincidono. Anche nuovi
download, file aggiunti o rimossi e riallineamenti invalidano i risultati quando
cambia il contenuto. `resume=false` forza una nuova elaborazione.

Su un cluster, `align_species_seq!`, `align_genus_seq!`, `align_family_seq!`,
`check_alignments!` e `filter_divergent_references!` accettano
`shard_index` e `shard_count`. Le famiglie vengono ordinate e ripartite in modo
deterministico, senza sovrapposizioni. Per esempio, un job array di 8 elementi
usa `shard_index=task_id, shard_count=8`. I job possono condividere
`alignment_dir`, perché ogni shard scrive famiglie e checkpoint distinti;
download o elaborazioni non suddivise non devono essere avviati in parallelo
sugli stessi file.

`asv_functional_structure!` usa invece task dinamici entro il nodo e accetta
`threads=:auto` oppure un limite intero. I calcoli delle ASV procedono in
parallelo in batch limitati al numero di thread, mentre i report aggregati
vengono composti nello stesso ordine dell'input.
`check_alignments!` restituisce anche i percorsi problematici delle famiglie
riutilizzate; supporta inoltre `GenusRank`, `FamilyRank` e `type="raw_sequences"`.

Il filtro registra l'impronta dei file **dopo** la modifica e conserva nei
checkpoint il resoconto delle sequenze eliminate. I report per famiglia sono
scritti in `alignment_dir/logs/divergent_references/<livello>/<famiglia>.tsv`
(personalizzabile con `log_dir`, oppure disattivabile con `log_dir=nothing`).
Le riprese non duplicano le righe dei report. Gli output vengono preparati
prima della sostituzione dei FASTA: una sostituzione interrotta puo essere
completata al rilancio senza ripetere il filtro. Se l'interruzione avviene
durante il calcolo, la famiglia in corso viene ricalcolata; quelle completate
vengono riutilizzate. I file di preparazione restano disponibili in caso di
interruzione. Il filtro opera sui file correnti: una soglia meno restrittiva
non ripristina sequenze eliminate in precedenza; per questo serve rigenerare
i riferimenti. `show_progress=false` disattiva i messaggi di avanzamento.

`align_species_seq!` mostra per impostazione predefinita le famiglie elaborate
o riutilizzate, l'avanzamento delle pagine BOLD, il filtraggio e la specie
in corso di allineamento con il numero di sequenze. Passare `show_progress=false`
per disattivare questi messaggi; `progress_interval` regola gli aggiornamenti
durante il filtraggio. I messaggi di allineamento sono emessi prima di ogni specie.
Le specie indipendenti della stessa famiglia vengono allineate in batch con task
del pacchetto separato `ProfileAlignments`; ogni batch è limitato al numero di
thread disponibile per contenere la memoria ed evitare parallelismo annidato.
Gli algoritmi di allineamento, i preset e le strutture numeriche restano nel
pacchetto `ProfileAlignments`; AASV gestisce download, tassonomia e file.

`align_genus_seq!` e `align_family_seq!` mostrano famiglie e taxa elaborati o
riutilizzati, il numero di sequenze nei profili da unire, le fusioni completate
con la percentuale per il taxon corrente e il tempo finale. Gli aggiornamenti
avvengono prima e dopo ogni fusione; `show_progress=false` li disattiva.

## Funzioni convertite

| Funzione | Operazione |
|---|---|
| `aa_similarity_matrix` | 4096 coppie di codoni con similarità Hamming |
| `align_species_seq!` | download BOLD COI-5P, filtri, traduzione, allineamenti |
| `align_genus_seq!`, `align_family_seq!` | fusione progressiva dei profili |
| `check_alignment` | percorsi degli allineamenti con gap interni |
| `check_alignments!` | controllo di tutte le famiglie con risultati persistenti |
| `filter_divergent_references!` | filtro per famiglia con ripresa e report persistenti |
| `re_align!` | riallineamento, opzionalmente dai FASTA originali |
| `improve_alignment!` | rimozione iterativa di sequenze problematiche e log |
| `filter_divergent_seq!` | filtro della distanza dal vicino più prossimo |
| `asv_functional_structure!` | analisi delle ASV ai tre livelli tassonomici |
| `classify_asv` | indice FI, classe, evidenze e metriche intermedie |

Sono disponibili anche `read_fasta`, `write_fasta`, `best_translation`,
`aa_table`, `align_sequences`, `align_profiles`, `calculate_asv_distance`,
`get_hydro_metrics` e `check_asv_structure`. `?nome_funzione` visualizza
la documentazione in Julia.

## Compatibilità e differenze

- La formula FI, i pesi (0.35/0.30/0.35), la soglia 0.90, il filtro qualità,
  l'entropia normalizzata e lo sconto di conservazione `1-exp(-n_ref/20)`
  seguono il codice R. Le distanze Grantham usano la tabella originale, non
  un'approssimazione. L'idrofobicità usa la scala Kyte–Doolittle.
- L'allineamento usa il pacchetto indipendente
  [ProfileAlignments.jl](../ProfileAlignments.jl/README.md), installato come
  progetto fratello e interamente Julia:
  BLOSUM62, gap affini (`gap_open=-10`, `gap_extend=-1`), distanze k-mer,
  albero guida, due passaggi progressivi e raffinamento iterativo.
  I gap terminali sono gratuiti; le query ASV vengono inserite nel profilo
  con `mode=:query`. `anchors=:fft` abilita l'accelerazione approssimata FFT,
  mentre il default usa la programmazione dinamica completa.
  Il precedente motore a penalità lineari è stato sostituito.
  Non riproduce DECIPHER o MAFFT: coordinate, gap, mutazioni e FI possono
  differire. La validazione biologica su benchmark reali resta necessaria.
- `Sequence`, `best_translation` e le metriche proteiche accettano direttamente
  i tipi BioSequences (`LongDNA`, `LongAA` e i literal `dna"..."`, `aa"..."`).
  La pipeline converte una sola volta questi input nella codifica interna ASCII,
  mantenendo compatibilità con i file e con l'API a stringhe.
- Gli argomenti specifici di DECIPHER (per esempio `processors`,
  `terminalGap`) non si trasferiscono direttamente; usare i parametri Julia
  sopra elencati. Gli argomenti precedenti `match`, `mismatch` e `gap`
  sono ancora accettati: `gap` imposta una penalità lineare.
  La distanza di filtro è la frazione di residui differenti
  nell'intersezione delle regioni coperte, escludendo i gap terminali.
- Codici genetici incorporati: i 20 codici della tabella Biostrings usata
  nel confronto R (1–6, 9–16, 21–26), senza reinterpretare i codoni
  iniziali. Per altri codici passare un dizionario completo codone→aminoacido.
  Le basi IUPAC traducibili senza ambiguità producono il residuo corrispondente;
  altrimenti `X`. La traduzione cerca esclusivamente i tre frame forward,
  ignora il codone terminale incompleto e non cerca il reverse complement.
- I profili identici vengono deduplicati per l'analisi strutturale, come in R.
  Nel ripristino di un allineamento interamente problematico, senza riferimento
  esterno viene mantenuta deterministicamente la prima sequenza, mentre R
  ne sceglie una casualmente.
- BOLD viene letto con l'API paginata anche per famiglie piccole. I retry
  mantengono lo stato del backoff; errori HTTP 4xx permanenti interrompono
  subito il download. Un marker `.complete` evita che una cartella creata
  durante un download fallito venga considerata completa. Le pagine scaricate
  vengono salvate in `alignment_dir/bold_downloads/<famiglia>` solo dopo la
  ricezione completa. Con `overwrite=false`, un rilancio riusa le pagine locali
  e richiede solo quelle mancanti della stessa query; `overwrite=true` avvia
  una nuova query senza mescolare pagine di download diversi. La cache viene
  mantenuta anche dopo il completamento. Se BOLD non rende piu disponibile
  la query salvata, occorre ripartire con `overwrite=true`. Per dati locali,
  passare `records_provider=family -> records` a `align_species_seq`.
- L'analisi registra anche confronti privi di mutazioni o con idrofobicità
  non calcolabile, evitando di perdere ASV dai report. In assenza della
  metrica idrofobica, il contributo FI è zero, come in R; il report contiene
  `UNAVAILABLE`. Uno stop o un gap incompatibile produce il sentinel
  `aa_pos=-1` e FI zero.
- I report mantengono nomi e colonne principali della versione R.
  `classify_asv` legge Excel o TSV. Il riesame della pipeline sostituisce
  i report precedenti nella cartella `aa_structure_results`.

L'indice FI misura la compatibilità con l'evoluzione di COI funzionale;
non è una probabilità né una diagnosi dell'origine NUMT/mitocondriale.

## Qualità delle letture con VSEARCH

```julia
asv_functional_structure(asvs;
    alignment_dir="references", output_dir="analysis",
    vsearch_path=raw"C:\tools\vsearch.exe",
    trimmed_folder="trimmed", ee=40)
```

I campioni sono selezionati dalle colonne numeriche con abbondanza >0,
cercando il nome del campione nel nome del FASTQ. Sono supportati FASTQ
non compressi (`.fastq`, `.fq`) con qualità Phred+33. Le probabilità di errore
sono `10^(-Q/10)`; il filtro `ee` usa la loro somma e la qualità del codone
è `1-somma(errori)` limitata a [0,1]. Gli errori di VSEARCH interrompono
l'analisi. Il calcolo delle probabilità è testato localmente; l'integrazione
con un eseguibile VSEARCH reale richiede la sua installazione.

## Verifiche locali

I test verificano FASTA, traduzione, gap, pulizia degli allineamenti,
distanze, qualità Phred, FI, retry BOLD simulati e pipeline completa offline.
Il confronto opzionale con R verifica ogni colonna numerica e la classe
su otto scenari, oltre alla traduzione, all'idrofobicità e ai codici genetici:

```powershell
julia --project=. test/parity.jl 'C:\Program Files\R\R-4.5.0\bin\Rscript.exe'
```

Questo controllo richiede R con `dplyr`, `readxl`, `Biostrings`, `Peptides`
e i sorgenti originali sotto `reference/aasv` (o un percorso alternativo
come secondo argomento). Il normale utilizzo e i test Julia non richiedono R.
La CI inclusa esegue i test su Windows e Linux, Julia 1.10 e versione stabile.
Il pacchetto di allineamento ha anche una propria suite, eseguibile da solo:

```powershell
julia --project=../ProfileAlignments.jl -e 'using Pkg; Pkg.test()'
```

## Provenienza

Il porting deriva da [aasv](https://github.com/alexology/aasv), di Alex Laini,
Samuele Voyron e Laura Gruppuso, distribuito con GPL >=2. Anche questo porting
è distribuito con GPL >=2. La tabella Grantham proviene dal pacchetto
[grantham](https://github.com/patterninstitute/grantham), licenza MIT;
vedere `THIRD_PARTY_NOTICES.md`. Le copie dei repository originali in
`reference/` sono materiale di confronto locale e sono escluse da Git.
