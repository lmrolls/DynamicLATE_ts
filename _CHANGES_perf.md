# Registro cumulativo delle modifiche — repository `dynamiclate_ts`
# Branch corrente di questa cartella: `perf`

File privato, non tracciato da git (escluso in `.git/info/exclude`): non entra in nessun commit e non
lascia questo PC. Il registro è cumulativo: ogni branch riporta tutte le modifiche dei branch
precedenti nella catena `main → branch3-tryMM → perf → fixes → dp`, più le proprie. Il file del
branch k è un sottoinsieme del file del branch k+1.

Convenzione: "invarianti" = ciò che il branch garantisce rispetto al precedente; "verifica" = come
è stato controllato; hash = commit pubblicato su GitHub.

---------------------------------------------------------------------------------------------------
## 0. Base: `main` (1b88302, 2026-07-08, lmrolls)
---------------------------------------------------------------------------------------------------
Stato di partenza (nessuna modifica nostra):
- `Dynamic_LATE_Analysis.m`, `test_Dynamic_LATE_all_DGPs.m`, `README.md`, `.gitattributes` (radice).
- `Section4/ftest/`: `piLATE_Ftest.m` (test F sup, euristica di selezione), `runSection4_Simul.m`,
  dati `dataForMatlabAltEstTick.csv`, cartella risultati pubblicati
  `FinalResults_Paper_Exhibits_Results_5000sims_FixedCV_ALL/` (log, .mat, .txt, figure 5–6 .eps e
  script `codeForPlots/`).
- `Section4/cvals/`: `Compute_CVs.m` (cv simulati con selezione `S_Dates_Norm` sul moto browniano),
  `RunMc_CVs_Fstat.m`, cartelle `FinalVersion_Results_CV_Export*/` (cv pubblicati).
- `Section5/`: `est_2sls_pilate.m` (OLS/2SLS su Ŝ), `est_fgls_pilate.m` (FGLS, selezione propria
  `S_Dates_multi_greedy_double_eq`), `runSection5_Simul.m`.
- `Section6/`: `robust_LATE_inference.m` (AR/LM/CLR su Ŝ, selezione per M2), `runSection6_Simul.m`.
- `m_break/` (67 .m: datazione Bai–Perron `Procedure_piL_Full_Sample`, `nldat`, `dating`,
  `S_Dates_m.m` a cinque blocchi m = 1…5, `OLS_Robust`, `NeweyWest*`, driver legacy `Run_pi_LATE*.m`,
  `S_Dates_m_Heter.m`, `S_Dates.m`, dati in `results/`, `yearscount/`), `Mbreak_Multivariate/` (57 .m,
  datazione multivariata Qu–Perron per il FGLS).

---------------------------------------------------------------------------------------------------
## 1. `branch3-tryMM` (c653b59, 2026-10-20, lmrolls) rispetto a `main`
---------------------------------------------------------------------------------------------------
Scopo: importare i competitor delle simulazioni (stimatore a campione intero per la Sez. 5; test di
Magnusson–Mavroeidis 2014 per la Sez. 6) senza toccare il codice π-LATE, con un'opzione negli script.
Due commit: 8db1e80 "reinstate MM" (2026-07-09: codice MM in `Section6/helpersMM/`, script separato
`runSection6_simulMM.m`, quattro cartelle di log legacy) e c653b59 "Step 1" (riorganizzazione:
sottocartelle `competitors_*`, opzione `include_competitors`, log legacy rimossi, script MM separato
eliminato). Effetto netto rispetto a `main`: 15 file, +1583/−58 righe.

File nuovi
- `Section5/competitors_FS/est_iv_fullsample.m`: IV/2SLS a campione intero, β̂_FS delle Tab. 7–8
  (porting letterale di `estimate_full_sample`/`OLS_Robust_full` del pacchetto originale; errori
  standard HC1 o NW, non usati nelle tabelle). `README.md` con provenienza e uso.
- `Section6/competitors_MM/` (11 file): `computeMMtests2.m` (GMM-S, split-S (GMM), qLL-S, exp-S,
  ave-S), `computeMM_IVstats2.m` (split-S (IV), split-KLM, split-JKLM, split-CLR, rk-stat),
  `getMagnussonAsymptoticCV.m` (cv asintotici simulati di qLL/exp/ave-S), `split_CLR_cv_seed.m`
  (cv condizionale di split-CLR, seme `rng_seed+3`), helper HAC `HACest.m`, `HACest_pw.m`,
  `HACest_tr.m`, `HAC_Sigma_hat.m`, `HACC_Sigma_hat.m`, `lag0.m`; `README.md`. Tutti i dieci .m sono
  byte per byte (a meno dei fine riga) quelli del pacchetto legacy (`Monte_Carlo/section_6`,
  `Magnusson_Mavroeidis/functions`): convenzioni legacy conservate (kernel NW con lag ⌊T^{1/3}⌋,
  trimming 0.05, ramo NW degli split, funzionale del cv di exp-S, Ω scalare nel ramo HOM) e
  dichiarate come da rivedere nel branch `fixes`.

File modificati
- `Section5/runSection5_Simul.m` (+112/−27): terzo argomento `include_competitors` (default true);
  stage 4 (`parfor` su `est_iv_fullsample`) sui dati già generati negli stage 1–3; campi
  `Bias_FS`, `MSE_FS`; righe FS nelle Tab. 7–8 dell'export nell'ordine pubblicato (FS, OLS, FGLS);
  token `_FS` nei nomi di cartella e file. Con `false` output identico a `main`.
- `Section6/runSection6_Simul.m` (+557/−31): terzo argomento `include_competitors` (default true);
  nella stessa iterazione, dopo `robust_LATE_inference`, calcolo delle 9 statistiche MM e delle
  diagnostiche grezze (M2_AMS, M2_ST, rk-stat); `size`/`power` a 15 colonne (6 π-LATE + 9 MM) più
  `raw`; cv MM: χ² per GMM-S, split-S, split-KLM, split-JKLM; asintotici simulati per qLL/exp/ave-S
  (per la Tab. 9 risimulati per scenario dopo `rng(d·1000 + π0·100 + X)` come nel legacy
  `simulFun_TAB_9.m`, interruttore `mm_cv_seed_per_scenario`); split-CLR per replica con seme
  `rng_seed+3`; Fig. 7–8 a potenza size-adjusted per tutte le 15 statistiche; export `_MM`.
  Con `false` output identico a `main`.

Invarianti
- Nessuna funzione numerica π-LATE toccata: `robust_LATE_inference.m`, `est_2sls_pilate.m`,
  `est_fgls_pilate.m`, `piLATE_Ftest.m`, `Compute_CVs.m`, `m_break/*`, `Mbreak_Multivariate/*`
  identici a `main`.
- Negli script le chiamate π-LATE (`generate_DGP(…, rng_seed)`, `robust_LATE_inference(…, rng_seed)`,
  regola di rigetto, `prctile` per colonna) sono byte per byte quelle di `main`; ogni iterazione si
  risemina (`rng(rng_seed,"twister")` nel DGP e nel cv del CLR), quindi i competitor non possono
  alterare le colonne π-LATE.

Verifica
- MATLAB (autore): `runSection5_Simul(2,'tab8',false/true)`: colonne OLS/FGLS identiche nei due modi,
  righe FS aggiunte; `runSection6_Simul(2,'tab10',false/true)` complete (gli zeri di rigetto con 2
  repliche sotto H0 sono attesi).
- Diff `main..branch3-tryMM` riesaminato riga per riga il 21/10: nessuna modifica aritmetica.

---------------------------------------------------------------------------------------------------
## 2. `perf` (721cf2c, 2026-10-20, lmrolls) rispetto a `branch3-tryMM`
---------------------------------------------------------------------------------------------------
Scopo: eliminare i calcoli ripetuti della griglia di ricerca (π_low × π × m) e i cinque blocchi
copia-incolla di `S_Dates_m`, con risultati identici bit per bit a seme fissato. Nessuna
correzione di metodo. Commit intermedio 078851a (file di prova `PERF_BRANCH_TEST.txt`, rimosso in
721cf2c). Effetto netto: 7 file, +285/−528 righe.

Principio: la selezione `S_Dates_m(…, π, m)` dipende da π solo attraverso `Tpi = floor(π·T)` (idem
`S_Dates_Norm`, `S_Dates_m_M2_Core`), e F/M2 dipendono da π solo attraverso S. Le celle (Tpi, m)
sono calcolate una volta e lette dal ciclo su π_low (che prima rifaceva le celle dei π_low maggiori:
125 celle per 45 coppie (π, m) distinte). Nessun `inv()`→`\`, nessuna vettorizzazione, nessuna somma
riordinata, nessuna estrazione casuale spostata; `max` sugli stessi vettori nello stesso ordine.

File
- `m_break/S_Dates_m.m` (riscritto, intestazione originale conservata): ora chiama
  `S_Dates_m_regimes` + `S_Dates_m_select`.
- `m_break/S_Dates_m_regimes.m` (nuovo): statistiche |t| (q = 1) o F (q > 1) degli m+1 regimi
  delimitati da `datevec(1:m, m)`, tabella `[|stat|, inizio, fine, lunghezza cumulata]` ordinata
  per |stat| decrescente; dipende da (y, z, x, datevec, m), non da π.
- `m_break/S_Dates_m_select.m` (nuovo): selezione greedy fino a Tpi date e finestra scorrevole
  nell'ultimo regime, con la convenzione di finestra di tryMM (`Start+i : Start+i+Missing`,
  `i = 1…Length−Missing`) per garantire l'identità.
- `m_break/S_Dates_m_grid.m` (nuovo): tutte le celle (Tpi, m) di una griglia, regimi una volta per
  m e ogni Tpi distinto una volta; `S_Dates_m(…, π, m) ≡ S_Dates_m_grid(…, floor(π·T), m)`.
- `Section4/ftest/piLATE_Ftest.m`: celle via `S_Dates_m_grid`, F per cella una volta, ciclo su
  π_low in sola lettura; `S_m_opt` preso dalla cella. `S_Dates_m` 125 → 0 chiamate dirette,
  `F_ST_Stat_Legacy` 125 → 45 per chiamata.
- `Section4/cvals/Compute_CVs.m`: cache (Tpi, m) della selezione (`S_Dates_Norm` + `sortrows` +
  `merge_Dates`) per replica; la statistica `norm(vec_sum)^2/(qw·π)` resta calcolata per ogni
  (π_low, π, m) perché dipende da π, non solo da Tpi. `S_Dates_Norm` e `merge_Dates` invariate.
- `Section6/robust_LATE_inference.m`: `S_Dates_m_M2_Core` spezzata in `S_Dates_m_M2_regimes`
  (una volta per m) e `S_Dates_m_M2_select` (per Tpi); `computeStats_Internal` a campione intero
  una volta per cella (45 invece di 125); `S_seg_final` dalla cella; `S_Dates_m_M2_Core` mantenuta
  come punto d'ingresso per una cella. `LR_cv_legacy_exact` invariata.

Non modificato deliberatamente: finestre scorrevoli con una regressione per finestra, datazione
Bai–Perron, FGLS, `S_Dates_Norm`, `merge_Dates` (vettorizzabili ma non identici bit per bit).

Verifica
- Octave, alberi tryMM e perf affiancati sugli stessi dati: `S_Dates_m` 3240 casi (T ∈ {100, 200,
  300}, q ∈ {1, 2}, X ∈ {[], 1, [1 x]}, 9 π × 5 m), `piLATE_Ftest` 64 casi, `robust_LATE_inference`
  32 casi (inclusi 2 errori identici), `Compute_CVs` 25 configurazioni: 0 differenze.
- MATLAB (autore, livello A): Sez. 5 72 celle e Sez. 6 54 celle identiche ai riferimenti tryMM.
- Da fare alla fine (livello B): `compare_runs(dirA, dirB)` sugli output completi di `main` e
  `perf` (Sez. 4 incluse `runSection4_Simul` e `Compute_CVs`).
