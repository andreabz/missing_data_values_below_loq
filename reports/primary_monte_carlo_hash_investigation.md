# Verifica forense della discrepanza SHA-256 di `results.rds`

Data dell'indagine: 22 novembre 2026. Tutti i controlli sui file sono stati eseguiti in sola lettura; gli oggetti R sono stati caricati in memoria e confrontati senza risalvarli.

## 1. Obiettivo e perimetro

L'indagine verifica la differenza tra l'hash attribuito a `results.rds` nella tabella del rapporto di audit e l'hash calcolato sul file corrente. Ho confrontato il contenuto corrente con i sei checkpoint, esaminato la produzione e i controlli previsti in `R/12_primary_monte_carlo.R`, `R/13_validate_primary_monte_carlo.R` e `R/14_analyze_primary_monte_carlo.R`, e ricostruito le evidenze disponibili nei rapporti.

Non ho eseguito runner, audit o script di analisi; non ho modificato né risalvato alcun artefatto primario. Non ho cercato copie fuori dal repository.

## 2. Stato iniziale del repository

- HEAD: `c29e92a3ea5b0116e46dbbdfb7465e60d96c88cf`, detached.
- Modifica non staged preesistente: `protocol/experimental_protocol.md`, con le due sostituzioni documentali `20261122` → `20261122L`.
- File locali non committati già presenti: `R/13_validate_primary_monte_carlo.R`, `R/14_analyze_primary_monte_carlo.R`, `reports/primary_monte_carlo_data_quality.md`, `reports/primary_monte_carlo_analysis.md` e `reports/figures/`.
- Era già presente anche `reports/primary_monte_carlo_scientific_review.md`.

Lo stato Git è stato consultato soltanto con comandi di lettura. La modifica al protocollo è rimasta intatta.

## 3. Hash registrati e correnti

L'hash “registrato nell'audit” è trascritto dalla tabella in `reports/primary_monte_carlo_data_quality.md`. Gli hash correnti sono stati calcolati direttamente sui percorsi indicati con `sha256sum`.

| File | Hash registrato nell'audit | Hash corrente | Corrispondenza |
|---|---|---|---|
| `n30_LOQ_10.rds` | `3ffcf53c0517e43554565872fc8b16507a04a2037643aba9555de770404eae7d` | `3ffcf53c0517e43554565872fc8b16507a04a2037643aba9555de770404eae7d` | Sì |
| `n30_LOQ_40.rds` | `a56a5adb05b80a1e440b7d45323aa50da880131cc67df99eaa34080fcde8f856` | `a56a5adb05b80a1e440b7d45323aa50da880131cc67df99eaa34080fcde8f856` | Sì |
| `n100_LOQ_10.rds` | `163410359b2d2ab34d20fe1ccc3a4b965a987cc0f4d478fa5e7262cfa46b143e` | `163410359b2d2ab34d20fe1ccc3a4b965a987cc0f4d478fa5e7262cfa46b143e` | Sì |
| `n100_LOQ_40.rds` | `397d32e6e1ceb89c71cb758f8a8fac802010888e51307e163ad687d9eee33ff0` | `397d32e6e1ceb89c71cb758f8a8fac802010888e51307e163ad687d9eee33ff0` | Sì |
| `n500_LOQ_10.rds` | `440dd278456def635563ab7bd8d5564c1ed1af259a96f7428efbf3e9afe4e5b0` | `440dd278456def635563ab7bd8d5564c1ed1af259a96f7428efbf3e9afe4e5b0` | Sì |
| `n500_LOQ_40.rds` | `529cd9fc4c48951e7b45ace8a07f10d986eddb9887f60b4accacd341adece905` | `529cd9fc4c48951e7b45ace8a07f10d986eddb9887f60b4accacd341adece905` | Sì |
| `results.rds` | `063c1719bf19c31b7bbfe8f8a8fac802010888e51307e163ad687d9eee33ff0` | `063c1719bf19c31b7bb16bc56483e644e6af2bd1d8a96e17140f06d71b2b8e7c` | **No** |
| `run_manifest.rds` | `5b05259fcc3c4ba530aab4d2e68999f347ee6eec52c1ef5c2908196a3be22161` | `5b05259fcc3c4ba530aab4d2e68999f347ee6eec52c1ef5c2908196a3be22161` | Sì |

La revisione scientifica riporta lo stesso hash discordante dell'audit come riferimento storico, e registra lo stesso hash corrente di questa indagine. Non è documentato altrove un secondo digest storico di `results.rds`.

## 4. Metadati dei file

Dimensioni e tempi correnti sono stati letti prima del caricamento R degli RDS in questa indagine. `atime` era disponibile. L'ora di modifica documentata nell'audit è riportata con precisione al secondo; per l'ora di accesso l'audit non riporta un valore storico. Tutti gli orari di filesystem sotto sono Europe/Rome (+0100).

| File | Byte | mtime nell'audit | mtime corrente | ctime corrente | atime osservato all'inizio dell'indagine | modo |
|---|---:|---|---|---|---|---:|
| `n30_LOQ_10.rds` | 1.301.726 | 2026-11-22 12:40:22 | 2026-11-22 12:40:22.659 | 2026-11-22 12:40:22.659 | 2026-11-22 12:50:02.414 | 664 |
| `n30_LOQ_40.rds` | 1.328.753 | 2026-11-22 12:46:12 | 2026-11-22 12:46:12.995 | 2026-11-22 12:46:12.997 | 2026-11-22 12:50:02.699 | 664 |
| `n100_LOQ_10.rds` | 1.364.688 | 2026-11-22 12:41:49 | 2026-11-22 12:41:49.254 | 2026-11-22 12:41:49.255 | 2026-11-22 12:50:01.502 | 664 |
| `n100_LOQ_40.rds` | 1.359.375 | 2026-11-22 12:47:24 | 2026-11-22 12:47:24.623 | 2026-11-22 12:47:24.624 | 2026-11-22 12:50:02.067 | 664 |
| `n500_LOQ_10.rds` | 1.361.700 | 2026-11-22 12:45:20 | 2026-11-22 12:45:20.751 | 2026-11-22 12:45:20.752 | 2026-11-22 12:50:02.976 | 664 |
| `n500_LOQ_40.rds` | 1.354.996 | 2026-11-22 12:49:37 | 2026-11-22 12:49:37.950 | 2026-11-22 12:49:37.951 | 2026-11-22 12:50:03.253 | 664 |
| `results.rds` | 8.779.726 | 2026-11-22 12:49:41 | 2026-11-22 12:49:41.110 | 2026-11-22 12:49:41.110 | 2026-11-22 12:50:03.522 | 664 |
| `run_manifest.rds` | 1.413 | 2026-11-22 12:49:38 | 2026-11-22 12:49:38.356 | 2026-11-22 12:49:38.356 | 2026-11-22 12:50:01.361 | 664 |

La dimensione corrente di `results.rds` (8.779.726 byte) coincide con quella documentata nell'audit. mtime e ctime correnti sono antecedenti all'indagine; non dimostrano da soli l'identità storica del contenuto. Nella directory dei risultati erano presenti soltanto `.gitkeep`, i sei checkpoint, `results.rds` e `run_manifest.rds`; la ricerca nel repository ha trovato un solo `results.rds` e nessuna copia/backup/temp con quel nome. Il manifest registra avvio `2026-11-22 11:39:18 UTC`, aggiornamento finale `2026-11-22 11:49:38 UTC`, `complete = TRUE` e 12.000 repliche completate. I rapporti documentano la data dell'audit e della revisione (22 novembre 2026), ma non l'ora di esecuzione; il rapporto di analisi non registra una data/ora di generazione.

## 5. Confronto strutturale

`results.rds` è un `data.table`/`data.frame` di 360.000 righe × 32 colonne. Ogni payload `checkpoint$results` è un `data.table`/`data.frame` di 60.000 × 32; i sei checkpoint coprono ciascuno 2.000 ID di replica, da 1 a 2.000. `results.rds` ha 60.000 record per scenario e 2.000 repliche per scenario. Le chiavi `(scenario_id, replicate_id, method, estimand)` non hanno duplicati.

Nomi, ordine e numero delle colonne, tipi e classi delle colonne coincidono tra il combined e la combinazione dei payload. Il combined ha 37.634 `NA` in `estimate` e nei due campi di errore; 12.000 `NA` in ciascun riferimento `sample_truth` e `population_target`; 300.000 `NA` in `probability`; e campi d'inferenza non disponibili (`interval_lower`, `interval_upper`, `ucl95` tutti `NA`). Sono valori mancanti attesi dagli estimand/status e dall'inferenza disabilitata. Non sono stati trovati `NaN` o `Inf` nei campi numerici. `warnings` e `error_message` sono `NA` in tutte le 360.000 righe; l'audit documenta zero fallimenti e zero warning/errori di esecuzione.

## 6. Confronto logico record per record

Ho caricato il combined e ciascun checkpoint separatamente, combinato in memoria i sei payload e confrontato i record usando la chiave indicata sopra. Ordinando copie dei due oggetti per chiave, `identical()` è risultato `TRUE`, inclusi valori, colonne, classi e attributi dell'oggetto confrontato. Anche il confronto a precisione zero, senza vincoli sugli attributi, è risultato esatto. Pertanto non sono emerse differenze per chiavi, stime, benchmark, errori, status, metadati o proporzioni empiriche di censura.

La sequenza originale delle righe non è identica alla concatenazione non ordinata dei checkpoint; il runner ordina il combined per scenario, replica, metodo ed estimand, mentre i checkpoint hanno il proprio ordine di accumulo. La differenza di sequenza scompare con l'ordinamento canonico per chiave, senza differenze record-valore. Questo è un confronto logico normalizzato, non un confronto binario dei file né la prova che il file storico avesse gli stessi byte.

Per gli status, risultano 48.000 record `success` per ciascuno dei sei metodi momentali (tutte le quattro statistiche × 12.000 repliche). Per Reverse KM sono registrati: media non identificabile 12.000; q2,5% e q10% teoricamente non identificabili 12.000 ciascuna; mediana non stimabile 415 e successi 11.585; q90% successi 12.000; q97,5% non restituita 4.000 e successi 8.000. Tali categorie coincidono con l'audit e l'analisi.

## 7. Analisi delle possibili cause

| Ipotesi | Evidenze favorevoli | Evidenze contrarie / verifica | Valutazione |
|---|---|---|---|
| Cambiamento sostanziale del contenuto corrente | Il digest riportato nell'audit non coincide con quello corrente. | Il combined corrente coincide esattamente, record per record e dopo ordinamento canonico, con i sei checkpoint i cui hash sono invariati; dimensione e mtime coincidono con i valori audit. | Nessuna evidenza di differenza logica corrente rispetto ai checkpoint. La storia dei byte precedenti non è ricostruibile. |
| Riserializzazione dello stesso oggetto | Una riserializzazione può cambiare i byte senza cambiare i valori. | Non esiste un hash storico verificabile direttamente sul vecchio file; nessuna evidenza documenta una riserializzazione successiva al run. | Possibile in astratto, ma non dimostrata e non necessaria per spiegare il dato. |
| Differenze di classe/attributi o precisione | Il confronto può essere sensibile a tipi e attributi. | Dopo normalizzazione dell'ordine, i due oggetti in memoria sono `identical()`; tipi/classi coincidono e confronto numerico a tolleranza zero passa. | Non supportata per il contenuto corrente. |
| Ordinamento o ricomposizione dai checkpoint | Il runner ricostruisce `results` da `scenario_results`, combina con `rbindlist`, ordina e salva il file completo. | Il confronto record per record conferma la stessa unione logica; l'ordinamento originale differisce ma la chiave normalizzata coincide esattamente. | Spiega il possibile ordine delle righe, non la discrepanza del digest testuale. |
| Errore di registrazione/trascrizione dell'hash nell'audit | `R/13_validate_primary_monte_carlo.R` calcola SHA-256 direttamente sui percorsi con `sha256sum` e stampa `file_info`; non genera/scrive il rapporto Markdown. Sette hash trascritti nella tabella corrispondono; uno no. L'analisi dichiara che gli otto hash prima/dopo sono rimasti invariati, e `R/14` verifica esplicitamente questa uguaglianza. | Non è disponibile l'output console o un log originale dell'audit con il digest di `results.rds`; non si può stabilire se l'errore sia avvenuto nella trascrizione o in un'altra fase. | Spiegazione plausibile e più coerente con le evidenze, ma non accertata. Confidenza moderata. |
| Altra causa | Nessun backup o file temporaneo corrispondente è presente nel repository. | Non sono disponibili copie esterne o log storici e non sono stati cercati fuori dal repository. | Non determinabile. |

Il flusso di `R/12` scrive i checkpoint e, a run completo, combina i payload, li ordina per scenario/replica/metodo/estimand e salva `results.rds` tramite salvataggio atomico. `R/13` calcola SHA-256 sul file fisico prima di caricarlo, poi compara record ordinati; produce output console e oggetto restituito in memoria, ma non scrive il rapporto Markdown presente. `R/14` calcola hash prima della lettura e dopo aver scritto report/figure, controlla che siano identici e verifica il payload dei checkpoint contro `results.rds` a tolleranza zero. Il rapporto di analisi contiene l'esito positivo del controllo prima/dopo, ma non i digest individuali né un log originale.

## 8. Conseguenze per l'analisi scientifica

Per il contenuto attualmente presente, il confronto record per record dimostra che `results.rds` è logicamente equivalente all'unione dei checkpoint conservati; non sono emerse differenze nei valori che alimentano le metriche o i grafici. Inoltre, l'analisi esistente documenta la riconciliazione del payload e l'invarianza degli hash durante la propria esecuzione. Perciò la discrepanza testuale dell'hash non fornisce evidenza di un cambiamento sostanziale nei dati usati dall'analisi corrente.

Questo non prova quale fosse il digest fisico del file nel momento dell'audit, né prova la provenienza storica di quei byte. La verifica di equivalenza con i checkpoint è forte per i record correnti, ma non risolve la riserva distinta sui valori latenti non conservati: `sample_truth` non può essere ricostruito indipendentemente dai file disponibili.

Non è necessario rigenerare l'analisi sulla base delle evidenze correnti. La discrepanza binaria/documentale deve però restare registrata finché l'hash originale non viene confrontato con un output console o log contemporaneo affidabile. Non ho corretto la tabella dell'audit o alcun altro rapporto.

## 9. Limiti dell'indagine

- Non esiste una copia storica di `results.rds` né un digest originale indipendente oltre al testo nella tabella audit.
- L'audit e la revisione riportano solo la data, senza ora esatta. L'orario di modifica del file e i timestamp del manifest attestano quando gli artefatti furono scritti/aggiornati, non quando fu osservato o trascritto il digest.
- L'uguaglianza binaria tra due versioni storiche del combined non è verificabile; il confronto svolto riguarda l'oggetto corrente e i checkpoint correnti.
- `sample_truth` è coerente tra le righe e con i riferimenti memorizzati, ma non è indipendentemente ricostruibile senza i campioni latenti.
- La lettura degli RDS può aggiornare `atime` secondo le regole del filesystem. Non sono state alterate dimensioni, contenuti, mtime o permessi intenzionalmente.

## 10. Giudizio conclusivo

**Discrepanza non spiegata, contenuto logico equivalente.** L'hash SHA-256 riportato per `results.rds` nell'audit è diverso dall'hash corrente e la provenienza della stringa storica non è recuperabile dai file disponibili. L'ipotesi di errore isolato di trascrizione/registrazione è plausibile, ma non provata. Al contempo, la corrispondenza record per record del combined corrente con i sei checkpoint è stata dimostrata a precisione esatta dopo normalizzazione dell'ordine; non sono state individuate differenze nei dati logici.

### Verifiche finali

- Hash di `results.rds` prima e dopo l'indagine: `063c1719bf19c31b7bb16bc56483e644e6af2bd1d8a96e17140f06d71b2b8e7c` (invariato durante questa indagine).
- Hash dei sei checkpoint e del manifest: identici ai valori registrati nell'audit e invariati tra inizio e fine dell'indagine.
- HEAD: `c29e92a3ea5b0116e46dbbdfb7465e60d96c88cf`, detached.
- `protocol/experimental_protocol.md`: modifica non staged `20261122` → `20261122L` conservata.
- Nessuna operazione Git modificativa, nessun commit e nessun push.
- L'unico file creato in questa indagine è `reports/primary_monte_carlo_hash_investigation.md`.
