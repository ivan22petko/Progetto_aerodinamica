Per runnare i bash su xfoil **usare il comando bash** e non run, altrimenti non funziona.
Tutti i file `.ps` bisogna aprirli su LINUX, ogni file `.ps` è poi stato convertito anche in `.png` in modo da poterlo aprire anche su Windows.

# HESS SMITH

In `hess_smith.m` sono presenti i risultati dell'analisi del profilo NACA 0012 a velocità incidente di 1m/s e angolo di incidenza di 2°. Il codice restituisce: 
- plot del profilo
- plot del coefficiente di pressione e confronto con quello proveniente da XFOIL, i cui dati sono caricati tramite il file estratto da XFOIL chiamato `CPWR_NACA_0012_2.txt` già presente nella cartella ma che può essere ottenuto usando i file bash riportati in seguito
- valori dei coefficienti di portanza e di momento rispetto al C.A. all'incidenza assegnata, confronto con i risultati provenienti da xfoil che sono contenuti nei file `polar_cl.txt` e `polar_cm.txt` 

Il codice `hess_smith_vett.m` riporta passaggi analoghi alla versione `hess_smith.m` ma per un range di alpha variabili, in modo tale da potere restituire:
- la curva Cl/alpha ed il confronto di tale curva con quella ottenuta da XFOIL, i cui dati sono contenuti nel file  `polar_cl.txt`

In `run_xfoil.d.sh` c’è il salvataggio dei valori del cp del NACA 0012, la polare usata per ottenere i dati del cl al variare di angoli di incidenza alfa (ogni 0.5°) e la polare da cui abbiamo preso i dati del CM al variare di angoli di incidenza alfa (ogni 0.25°).

# SEPARAZIONE E TRANSIONIZIONE

In questa cartella sono riportati tutti i bash che servono per calcolare i grafici del CF per vedere separazione e transizione. Inoltre è presente anche un file `.dat` dove sono contenuti i CF per ogni punto della corda. 
- Per CF1 si intende re=250000, alfa=1 (N=9)
- Per CF2 si intende re=250000, alfa=2 (N=9)
- Per CF3 si intende re=2500000, alfa=1 (N=9)
- Per CF4 si intende re=2500000, alfa=2 (N=9)
- Per CF5 si intende re=25000000, alfa=1 (N=9)
- Per CF6 si intende re=25000000, alfa=2 (N=9)

# THEODORSEN
Ci sono due cartelle:

**theodorsen_matlab**: deve essere compilato il codice `theodorsen_BL207_spline.m` che permette di calcolare l’angolo di progetto, l'angolo di zero lift e l'angolo di momento nullo rispetto all’origine del profilo BL207. Il grafico del cp generato dal codice è quello presente nel report.

**theodorsen_bash**: contiene due bash (da eseguire sul terminale usando bash), il primo è `run_theod.sh`, che esegue la simulazione inviscida del profilo a incidenza pari all’angolo di theodorsen e salva un’immagine con il grafico del cp e l’elenco di coordinate cp-x associate a quel grafico. Il secondo è `run_zerolift.sh` che fa la stessa simulazione ma all’incidenza di zero lift calcolata in matlab, salva grafico e coordinate, e serve per mostrare che il Cl calcolato a quell’incidenza vale 0.0033 come detto nel report.

# WEISSINGER
I codici `cessna_1surf_def.m`, `hawk_1surf_def.m` applicano il metodo di Weissinger sulle ali dei due velivoli, considerate isolate. I codici `cessna_2surf_def.m`, `hawk_2surf_def.m` applicano Weissinger ai due velivoli considerando sia ala che coda.

Il codice `confronto_wing_wingtail.m` permette di generare il grafico presente nel report che confronta le polari generate nei quattro codici precedenti. **Importante:** per essere compilato correttamente, devono prima essere compilati i 4 codici precedenti.

Il codice `plot_circolazioni.m` permette di generare il grafico con i confronti delle distribuzioni di circolazione presente nel report. **Importante:** per essere compilato correttamente devono prima essere compilati i primi 4 codici qui citati.