# Revisione scientifica dell'analisi Monte Carlo principale

Data della revisione: 22 novembre 2026 (Europe/Rome)  
Commit esaminato: `c29e92a3ea5b0116e46dbbdfb7465e60d96c88cf` (detached HEAD)

## A. Sintesi esecutiva

La revisione indipendente di protocollo, runner, metodi, audit, script di analisi, rapporto, manifest/checkpoint e 19 grafici non ha rilevato errori nelle formule aggregate principali né scambi tra i due benchmark. Le metriche sono stratificate per scenario, metodo ed estimand; i confronti con il campione completo e con la popolazione sono distinti. Il codice esclude correttamente dalle metriche Reverse KM gli estimand non identificabili anche quando esiste un valore numerico.

Il rapporto è scientificamente prudente e non formula una graduatoria generale. Prima di considerare conclusa la fase, raccomando però modifiche interpretative minori: dichiarare in modo più diretto le assunzioni del log-ROS e della MLE lognormale, chiarire che le metriche Reverse KM sono condizionate alla restituzione/identificabilità, e ricordare che i boxplot nascondono i punti outlier. La mancanza dei valori latenti salvati resta una riserva sostanziale già correttamente dichiarata dall'audit e dal rapporto: non è possibile ricostruire indipendentemente le statistiche `sample_truth`.

**Giudizio complessivo: idoneo con modifiche minori.** Non ho trovato un problema bloccante nelle conclusioni numeriche riportate; la valutazione vale per il disegno lognormale simulato e non convalida prestazioni fuori da tale scenario.

## B. Aspetti metodologici verificati

- **Disegno:** il manifest e il runner corrispondono al protocollo: generazione lognormale con `meanlog = 0`, `sdlog = 1`; n = 30, 100, 500; LOQ fissati ai quantili teorici 10% e 40%; 2.000 repliche per ciascuna delle sei condizioni; seed intero `20261122L`; sette metodi; inferenza disabilitata.
- **Appaiamento e benchmark:** tutti i metodi applicabili ricevono i dati censurati derivati dal medesimo campione latente della replica. `error_vs_sample_truth` ed `error_vs_population` usano riferimenti separati. Il rapporto distingue l'errore rispetto alla statistica del campione completo dall'errore rispetto al parametro teorico.
- **Metriche:** il codice calcola media e mediana degli errori, deviazione standard campionaria degli errori, RMSE come `sqrt(mean(error^2))`, quantili di tipo 7 e, per il benchmark di popolazione, SD delle stime. Il bias è la media dell'errore di popolazione. Il bias relativo, essendo calcolato dividendo ogni errore per il parametro teorico fisso prima della media, coincide con bias/parametro nei gruppi considerati. Non ho rilevato denominatori teorici nulli nel disegno lognormale; i riferimenti nulli o prossimi a zero sono comunque esclusi dalle metriche relative e conteggiati.
- **Variabilità:** le distribuzioni delle statistiche del campione completo sono riportate separatamente; il testo non attribuisce tale variabilità ai metodi di censura. Il rapporto distingue bias, RMSE e dispersione delle stime e non sostiene un miglioramento uniforme con n.
- **Metodi e estimand:** complete case è presentato come analisi dei rilevati; le sostituzioni come statistiche di campioni costruiti; MLE come stima parametrica con momenti derivati, non come imputazione; Reverse KM è separato dalle statistiche momentali ricostruite. Le tabelle stratificano le metriche per estimand.
- **Reverse KM:** la media è esclusa perché strutturalmente non identificabile con LOQ comune. Gli stati teoricamente non identificabile, non stimabile nel campione e non restituito sono distinti. Valori finiti associati a non-identificabilità teorica non entrano nelle metriche. Il conteggio riporta 12.000 medie non identificabili, 12.000 q2,5% e 12.000 q10% teoricamente non identificabili, 415 mediane non stimabili e 4.000 q97,5% non restituite.
- **Censura empirica:** le proporzioni osservate derivano dai campioni e variano tra repliche. Medie e dispersioni sono coerenti con le probabilità teoriche; le correlazioni sono calcolate entro scenario, metodo, estimand e benchmark solo su coppie finite ed eleggibili, con n riportato. Il testo le qualifica come descrittive e non causali.
- **Integrità dei risultati:** l'audit documenta sei checkpoint da 2.000 repliche e 60.000 righe ciascuno, manifest completo, 360.000 righe aggregate riconciliate record per record, nessun fallimento computazionale dei sei metodi momentali e nessun errore/warning registrato. L'analisi riporta i conteggi di disponibilità e le categorie di esclusione.
- **Grafici:** sono presenti i 19 file previsti: 12 pannelli 2×2 per gli errori dei sei scenari rispetto ai due benchmark, sei pannelli per censura empirica/errore di popolazione e un pannello Reverse KM 2×3. Le routine condivise riportano scenario, benchmark o asse di errore, statistica, metodi e conteggi disponibili; Reverse KM resta separato. Le etichette indicano proporzione empirica e errore nelle unità della statistica, non percentuali di errore. Le figure non costruiscono una graduatoria aggregata.
- **Riserva di provenienza:** il rapporto ripete correttamente che i campioni latenti non sono stati conservati. La coerenza aritmetica degli errori e l'uguaglianza dei riferimenti tra metodi non sono una ricostruzione indipendente di `sample_truth`.

## C. Problemi individuati

### C1. Assunzioni del log-ROS non esplicitate con sufficiente precisione

- **Gravità:** minore.
- **File/sezione:** `reports/primary_monte_carlo_analysis.md`, “Definitions and eligibility” e conclusioni; metodo in `R/04_method_ros.R`.
- **Descrizione:** il rapporto identifica ROS e lo distingue dalle sostituzioni, ma non spiega che l'implementazione usa `NADA2::ros` con trasformazione log/exp. La documentazione primaria del pacchetto descrive la regressione dei quantili normali associati ai valori non censurati e il back-transform; quindi non è un metodo privo di assunzioni. In questo esperimento la generazione lognormale rende l'assunzione di normalità sulla scala log coerente con il modello generatore, ma non dimostra robustezza a misspecificazione.
- **Conseguenza:** il lettore potrebbe interpretare l'esito ROS come proprietà generale del metodo, anziché come risultato sotto un disegno che coincide con la specificazione lognormale adottata.
- **Correzione proposta:** aggiungere una frase che specifichi la trasformazione log-ROS, la regressione sui plotting positions degli uncensored e la dipendenza dei risultati dall'adeguatezza della forma lognormale. Riservare una valutazione di robustezza allo scenario Gamma già previsto separatamente. Riferimento tecnico: [documentazione CRAN di NADA2, sezione ROS](https://cran.r-project.org/web/packages/NADA2/refman/NADA2.html).

### C2. Metriche Reverse KM condizionate all'eleggibilità

- **Gravità:** minore.
- **File/sezione:** `reports/primary_monte_carlo_analysis.md`, disponibilità/metriche Reverse KM e figura `reverse_km_quantile_errors_vs_population.png`.
- **Descrizione:** il rapporto documenta che le metriche usano solo righe con stato `success`, e mostra numerosità e stati di disponibilità. Tuttavia, RMSE, bias e quantili di errore delle quantili restituite sono riassunti condizionatamente alle repliche in cui la quantile è identificabile nel campione e restituita dal pacchetto. Questo è particolarmente rilevante per q97,5%, con 8.000 successi su 12.000; la distribuzione dell'errore osservata non descrive le 4.000 mancate restituzioni.
- **Conseguenza:** un confronto basato solo sull'errore dei successi può apparire più favorevole o più stabile di una valutazione che includesse anche la disponibilità. Le statistiche sono valide come accuratezza condizionale, ma non come prestazione complessiva incondizionata.
- **Correzione proposta:** aggiungere esplicitamente “metriche condizionate all'eleggibilità; leggere insieme alla disponibilità; non rappresentano una valutazione incondizionata delle repliche richieste”. Mantenere separati mancata identificabilità e mancata restituzione, come già fa l'analisi.

### C3. Outlier non disegnati nei boxplot

- **Gravità:** minore.
- **File/sezione:** codice grafico in `R/14_analyze_primary_monte_carlo.R` (`outline = FALSE`) e i 13 grafici a boxplot.
- **Descrizione:** `outline = FALSE` sopprime la visualizzazione dei punti fuori dai baffi. Le tabelle includono quantili 2,5% e 97,5%, quindi le code sono parzialmente documentate, ma non sono necessariamente visibili nella figura. Il campo “n” sulle figure indica le osservazioni incluse, non il numero di punti mostrati.
- **Conseguenza:** lettori che si basano soltanto sulla figura possono non percepire stime/errori estremi o asimmetrie di coda.
- **Correzione proposta:** dichiarare questa scelta nella didascalia e richiamare i quantili tabellari; in una revisione grafica successiva si può valutare di mostrare outlier o aggiungere un pannello/riassunto di coda, senza cambiare ora i grafici.

### C4. Associazione tra censura ed errore: interpretazione da mantenere strettamente descrittiva

- **Gravità:** osservazione.
- **File/sezione:** `reports/primary_monte_carlo_analysis.md`, “Empirical censoring and error association”.
- **Descrizione:** le correlazioni Pearson/Spearman sono implementate correttamente per coppie utilizzabili entro scenario e riportano n. Il rapporto già specifica che la proporzione empirica è realizzata dal campione, non assegnata, e che le associazioni non sono causali. Resta utile esplicitare che campione, errore e tasso empirico di censura condividono la stessa estrazione; la correlazione non isola un effetto della sola censura e può riflettere composizione campionaria e dipendenze del disegno.
- **Conseguenza:** una lettura causale o un confronto della forza delle correlazioni tra n diversi sarebbe improprio.
- **Correzione proposta:** accompagnare la sezione con la precisazione che si tratta di associazioni entro scenario/dimensione campionaria e non di effetti causali né di una decomposizione dell'errore.

Non sono stati individuati problemi critici o maggiori nelle formule implementate e nella loro applicazione ai risultati descritti.

### C5. Discrepanza nell'hash registrato per `results.rds`

- **Gravità:** osservazione da chiarire prima di attestare l'invarianza degli artefatti.
- **File/sezione:** `reports/primary_monte_carlo_data_quality.md`, tabella degli hash; directory dei risultati primari.
- **Descrizione:** alla verifica finale, `results.rds` ha SHA-256 `063c1719bf19c31b7bb16bc56483e644e6af2bd1d8a96e17140f06d71b2b8e7c`, mentre l'hash riportato nell'audit e annotato prima di questa revisione è `063c1719bf19c31b7bbfe8f8a8fac802010888e51307e163ad687d9eee33ff0`. Una seconda lettura ha confermato l'hash corrente. Gli altri sette RDS (sei checkpoint e manifest) corrispondono agli hash annotati. Il timestamp di `results.rds` è 22 novembre 2026 alle 12:49:41 +0100, antecedente a questa revisione.
- **Conseguenza:** non è possibile attestare la corrispondenza di `results.rds` con l'hash riportato nell'audit; la differenza potrebbe derivare da un hash trascritto erroneamente nell'audit o da una diversa versione del file, ma questa revisione non determina quale spiegazione sia corretta. Non ho eseguito scritture sui risultati primari.
- **Correzione proposta:** confrontare l'hash con un record di provenienza indipendente/affidabile e correggere la documentazione soltanto dopo tale verifica. Non sovrascrivere né rigenerare `results.rds` sulla base della sola discrepanza.

## D. Risultati scientifici da approfondire

- **Censura empirica e n:** per LOQ teorico 10%, la SD della proporzione empirica scende da 5,408 punti percentuali a n=30 a 1,358 a n=500; per LOQ 40%, da 9,254 a 2,149 punti. Le frequenze realizzate restano variabili e discrete, soprattutto ai campioni piccoli. Ciò sostiene l'interpretazione che il tasso teorico non è imposto replica per replica.
- **Errore campionario e errore del metodo:** nel caso n=30/LOQ 10%, la SD della media dei campioni completi è 0,4044; il RMSE rispetto alla popolazione per i sei metodi momentali va da 0,4034 a 0,4727. Questo rende importante non leggere il RMSE totale come puro errore introdotto dalla censura: contiene anche la variabilità di estrazione.
- **Differenza tra bias e RMSE:** a n=500/LOQ 40%, il bias della media è 0,8096 per complete case, ma vicino a zero per ROS (-0,0028) e MLE (-0,0003). I relativi RMSE sono 0,8219, 0,0921 e 0,0897. È un contrasto forte nel singolo scenario, ma non giustifica una graduatoria generale: occorre leggerlo accanto a variabilità, disponibilità, estimand e specificazione dei metodi.
- **Reverse KM alle code:** la disponibilità non equivale a identificabilità. q2,5% e q10% hanno rispettivamente 76 e 2.705 valori numerici finiti pur essendo non identificabili secondo lo status teorico; q97,5% ha 4.000 mancate restituzioni. La corretta gestione degli status è quindi sostanziale, non meramente tecnica.
- **Campione truth:** la diminuzione della variabilità delle statistiche complete all'aumentare di n è coerente con il disegno; la fonte dei valori per replica non è però riproducibile indipendentemente senza i dati latenti. Le conclusioni che dipendono da `sample_truth` restano quindi subordinate alla correttezza del runner e alla verifica interna disponibile.

## E. Valutazione delle conclusioni

### Pienamente supportate dai risultati disponibili

- I risultati sono completi e riconciliati rispetto a manifest, checkpoint e disegno dichiarato.
- Errori rispetto alla statistica del campione completo e rispetto al parametro di popolazione rispondono a domande diverse e sono riportati separatamente.
- La variabilità delle statistiche del campione completo è separata dagli errori dei metodi; non viene presentata come errore causato dalla censura.
- Il Reverse KM non va classificato insieme ai metodi che stimano momenti ricostruiti; gli stati non identificabili/non stimabili/non restituiti non sono conteggiati come errori computazionali.
- Non è giustificata una graduatoria generale dei metodi e il rapporto non la formula.

### Supportate con riserve

- Le conclusioni quantitative basate su `sample_truth` sono internamente coerenti, ma non costituiscono verifica indipendente dei valori campionari latenti.
- Le prestazioni ROS e MLE sono risultati sotto un generatore lognormale coerente con le rispettive assunzioni. Non descrivono la robustezza a distribuzioni diverse.
- Gli errori Reverse KM descrivono solo le stime eleggibili restituite; la disponibilità deve rimanere parte dell'interpretazione.
- Le correlazioni con la censura empirica sono descrittive e condividono la variabilità di campionamento con gli errori.

### Da riformulare

- Qualsiasi frase futura che sintetizzi l'accuratezza Reverse KM senza specificare il numero di repliche eleggibili e la natura condizionale delle metriche.
- Qualsiasi descrizione di ROS come non parametrico o privo di assunzioni, o di MLE come generalmente robusta senza scenario di misspecificazione.

### Non supportate

- Generalizzazioni a popolazioni non lognormali, a LOQ diversi o a scenari applicativi non simulati.
- Affermazioni causali basate sulle correlazioni tra proporzione empirica censurata ed errore.
- Confronti o ranking globali che ignorino estimand, benchmark, scenario e disponibilità.

## F. Raccomandazioni

1. Integrare nel rapporto una breve descrizione delle assunzioni del log-ROS e del carattere parametrico lognormale della MLE, chiarendo che lo scenario corrente non valuta misspecificazione.
2. Esplicitare che le metriche Reverse KM sono condizionate alle quantili identificabili/restituite e vanno lette con disponibilità e conteggi degli status.
3. Segnalare che i boxplot sopprimono gli outlier e indirizzare alle tabelle di quantili per la lettura delle code.
4. Mantenere l'interpretazione delle correlazioni come descrittiva, entro ciascuno scenario, senza attribuzione causale o decomposizione dell'errore.
5. Conservare la riserva sulla non ricostruibilità indipendente di `sample_truth`. Per rendere verificabili in futuro le statistiche del campione completo, considerare in un progetto successivo la persistenza di statistiche latenti essenziali o dei campioni latenti, con una decisione esplicita su dimensioni e riproducibilità degli artefatti.
6. Trattare lo scenario Gamma separato come evidenza di robustezza/misspecificazione, senza incorporarlo retroattivamente nelle conclusioni di questo esperimento.

Non è stata applicata alcuna correzione al rapporto analitico, ai grafici, agli script, ai risultati, al manifest, ai checkpoint o al protocollo.

## G. Giudizio finale

**Idoneo con modifiche minori.** L'analisi è internamente coerente con il disegno e le statistiche descritte sono calcolate con definizioni appropriate. Le principali riserve non indicano errori numerici individuati: riguardano chiarezza delle assunzioni ROS/MLE, interpretazione condizionale delle metriche Reverse KM e visibilità delle code nei grafici. La dipendenza dal solo scenario lognormale e la mancata disponibilità dei campioni latenti delimitano la portata delle conclusioni. Nessuna correzione proposta è stata applicata in questa revisione.

## Verifiche di integrità e arresto

La revisione non ha eseguito simulazioni, bootstrap, inferenza, aggregazioni aggiuntive o scritture ai dati primari. I sei checkpoint e il manifest conservano gli hash annotati all'avvio della revisione; per `results.rds` permane la discrepanza descritta in C5, dunque la sua invarianza rispetto al riferimento documentato non è attestabile. Non sono state eseguite operazioni Git modificative; non sono stati creati commit né effettuati push. Il solo file aggiunto in questa fase è questo rapporto; la modifica locale al seed in `protocol/experimental_protocol.md` è stata conservata.
