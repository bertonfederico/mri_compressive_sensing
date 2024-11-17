# Compressione delle Immagini MRI

Nel campo dell’imaging medico, in particolare con la Risonanza Magnetica (MRI), la necessità di compressione delle immagini nasce dalla grande quantità di dati che esse generano. Le immagini MRI richiedono infatti molto spazio di archiviazione, il che può risultare costoso e complicato da gestire, specialmente quando si tratta di scansioni ad alta risoluzione o di più sezioni (slice). Le tecniche di compressione diventano quindi fondamentali per ridurre i costi di archiviazione e trasmissione dei dati, mantenendo al contempo una qualità sufficiente per le finalità diagnostiche.

Le immagini MRI sono spesso acquisite in scala di grigi e presentano una certa simmetria intrinseca. Questa simmetria può essere sfruttata per la compressione. In particolare, un'immagine MRI, che rappresenta il dominio spaziale, può esibire proprietà simmetriche che consentono di ridurre i dati. Ad esempio, quando un'immagine MRI è in scala di grigi, la metà sinistra dell'immagine è spesso una riflessione della metà destra. Questa simmetria suggerisce che, in alcuni casi, i dati potrebbero essere ridotti concentrandosi solo su un lato dell'immagine. Tuttavia, la simmetria nelle immagini MRI non è sempre perfetta, quindi vengono esplorate altre modalità di compressione per ottenere una riduzione ottimale delle dimensioni dei dati.

## Tecniche di compressione sul k-space dell'MRI

Il k-space è il dato grezzo acquisito durante la scansione MRI, che contiene tutti i componenti di frequenza necessari per ricostruire l'immagine. L'idea alla base della compressione nel k-space è che non tutte le frequenze contribuiscono in modo uguale all'immagine finale. Alcune frequenze, in particolare quelle al centro del k-space, rappresentano i componenti a bassa frequenza, che generalmente contengono le informazioni più significative sull'immagine. Le componenti ad alta frequenza, situate ai margini del k-space, rappresentano spesso dettagli più fini, che possono essere eliminati più facilmente senza compromettere significativamente la qualità dell'immagine.

Per esplorare questa idea, possiamo applicare diverse tecniche di campionamento nel k-space dopo aver eseguito una trasformata FFT 2D dell'immagine MRI. Tre diverse modalità di campionamento vengono considerate per la compressione:

1. **Campionamento gaussiano**: in questa modalità, i valori di k-space vengono campionati seguendo una distribuzione gaussiana, concentrando il campionamento maggiormente verso il centro del k-space.
   
2. **Campionamento casuale**: in questa modalità, i valori di k-space vengono scelti in modo casuale, senza alcuna priorità, per creare una selezione di campioni che non segue un pattern specifico.

3. **Campionamento a soglia**: qui, vengono selezionati solo i coefficienti nel k-space che superano una certa soglia di ampiezza, escludendo i valori più bassi, che generalmente contengono informazioni meno significative per l'immagine.

Per confrontare queste modalità, viene presa la stessa percentuale (10%) dei valori più alti e gli altri vengono posti a zero. Questo approccio permette di osservare l'effetto di ciascun tipo di campionamento sulla qualità dell'immagine finale e sulla compressione dei dati. I risultati mostrano che la modalità ottimale è quella che seleziona i coefficienti più alti, poiché essa conserva le informazioni più rilevanti per la ricostruzione dell’immagine. Tuttavia, la modalità gaussiana non si discosta molto da essa, dato che la FFT2D nel k-space tende a concentrarsi principalmente attorno al centro, il che rende simile il campionamento gaussiano al campionamento dei valori più alti.

## Modalità di Condivisione dei Dati dopo il Campionamento della FFT2D

Una volta eseguito il campionamento del k-space, è importante considerare come i dati possano essere condivisi efficacemente. La FFT2D del k-space di un’immagine MRI rappresenta di fatto la forma grezza dei dati immediatamente ottenuta dalla risonanza. Pertanto, è possibile condividere direttamente il k-space dell'MRI, ma dopo il campionamento, diventa necessario conoscere le modalità con cui i dati sono stati campionati per poterli ricostruire correttamente.

- **Modalità gaussiana e casuale**: per queste due modalità, si potrebbe utilizzare una maschera di campionamento fissa. In altre parole, la stessa maschera di campionamento potrebbe essere utilizzata ogni volta, il che eliminerebbe la necessità di trasmettere informazioni aggiuntive su come i dati siano stati campionati. Questa soluzione semplificherebbe la gestione dei dati dopo il campionamento.

- **Modalità di Selezione dei Valori Più Alti**: nel caso del campionamento basato sulla selezione dei valori più alti, è necessario in qualche modo condividere informazioni aggiuntive riguardo alla posizione dei valori campionati nel k-space. Poiché i valori vengono selezionati in modo non uniforme, è necessario sapere esattamente quali coefficienti sono stati scelti per poter ricostruire correttamente l'immagine.

Per affrontare questa necessità, sono state testate due modalità per rappresentare i dati campionati:

1. **COO (Coordinate Format)**: questa modalità rappresenta i dati come una matrice sparsa, memorizzando solo i valori non nulli e le loro posizioni. L'indice delle righe e delle colonne viene memorizzato separatamente, riducendo notevolmente la quantità di dati da trasmettere rispetto a una matrice completa.

2. **RLE (Run-Length Encoding)**: l'RLE è una tecnica di compressione che rappresenta sequenze consecutive di valori uguali (run) come una coppia di valore e frequenza (cioè il valore stesso e il numero di volte che si ripete). Nel contesto del k-space, questa modalità può essere utile per comprimere aree in cui molti coefficienti sono zero o hanno valori simili, riducendo il numero di dati da trasmettere.

<table align="center">
  <tr>
    <th style="width: 300px;">Metric</th>
    <th style="width: 200px;">Value</th>
  </tr>
  <tr>
    <td>Original size of FFT2D (in bytes)</td>
    <td>1,048,576</td>
  </tr>
  <tr>
    <td>Size of COO sparse matrix (in bytes)</td>
    <td>157,320</td>
  </tr>
  <tr>
    <td>Size of RLE compressed data (in bytes)</td>
    <td>263,580</td>
  </tr>
  <tr>
    <td>Memory saved percentage (COO)</td>
    <td>15.00%</td>
  </tr>
  <tr>
    <td>Memory saved percentage (RLE)</td>
    <td>25.14%</td>
  </tr>
</table>

