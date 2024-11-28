# Sparsità del segnale MRI

La **rappresentazione sparsa** è un concetto fondamentale nell'elaborazione dei segnali, particolarmente utile per comprimere e rappresentare segnali che contengono una gran parte di informazioni in un numero ridotto di componenti significative. In generale, un segnale $x$ può essere rappresentato in una base $\Psi$ come una combinazione lineare di pochi elementi:

$$
x = \Psi s
$$

Dove:
- $x$ è il segnale originale.
- $\Psi$ è una matrice che rappresenta una **base ortonormale** (oppure, se la rappresentazione quadrata non è sufficiente a garantire la sparsità desiderata, una matrice con un numero maggiore di colonne rispetto alle righe).
- $s$ è un vettore di **coefficienti sparsi**, cioè un insieme di pesi che combinano le colonne di $\Psi$ per ricostruire $x$.


Il problema si traduce nell'ottimizzare la seguente funzione:

$$
\mathbf{s}^* = \arg\min_{\mathbf{s}} \|\mathbf{x} - \Psi \mathbf{s}\|_2^2 \quad \text{s.t.} \quad \|\mathbf{s}\|_0 \leq k,
$$

dove $k$ rappresenta la quantità massima di coefficienti non nulli che si desidera utilizzare nella forma sparsa $s$.

La risoluzione esatta del problema è **NP-hard**; per superare la complessità computazionale, si utilizza una forma rilassata del problema, sostituendo la norma $\ell_0$ (discontinua e non convessa) con la norma convessa. 

Nel caso di segnali naturali, come MRI, la loro rappresentazione sparsa è spesso ben approssimata in basi come **Fourier**, **Wavelet** o **DCT**.

# Trasformata di Fourier 2D

Nel caso delle immagini, la teoria generale della sparsità può essere applicata utilizzando la **Trasformata Discreta di Fourier (DFT)**, che è una base ortonormale per rappresentare un'immagine nel dominio delle frequenze.

Nel caso delle immagini MRI, i dati acquisiti non sono direttamente immagini spaziali, ma sono misurazioni in quello che viene chiamato il **k-space**. Il k-space è una rappresentazione della frequenza spaziale dell'immagine, che descrive come le varie frequenze (o modelli spaziali) sono distribuite nell'immagine. Per ottenere l'immagine finale, i dati nel k-space devono essere trasformati nel dominio spaziale tramite l'inversa della trasformata di Fourier.

La **trasformata di Fourier bidimensionale** (**DFT2D**) è una tecnica matematica utilizzata per trasformare un'immagine dal dominio spaziale al dominio delle frequenze spaziali. Questa trasformazione consente di rappresentare un'immagine come una combinazione di frequenze (componenti ad alta e bassa frequenza) anziché come una distribuzione di intensità spaziale.

In termini matematici, la **DFT2D** di un'immagine $f[x, y]$ (con $x$ e $y$ che rappresentano le coordinate spaziali dell'immagine) è data da:

$$
F[u, v] = \sum_{x=0}^{M-1} \sum_{y=0}^{N-1} f[x, y] e^{-j2\pi \left(\frac{ux}{M} + \frac{vy}{N}\right)}
$$

dove:
- $F[u, v]$ sono i coefficienti nel dominio delle frequenze,
- $f[x, y]$ è il valore del pixel nell'immagine spaziale,
- $M$ e $N$ sono le dimensioni dell'immagine,
- $u$ e $v$ rappresentano le frequenze spaziali nelle direzioni $x$ e $y$,
- $j$ è l'unità immaginaria.

Quando un'immagine è trasformata nel dominio delle frequenze con la DFT2D, essa viene rappresentata come una combinazione di frequenze basse e frequenze alte:

- **Frequenze basse**: Si trovano vicino al centro del dominio delle frequenze e rappresentano le informazioni globali dell'immagine, come forme, contorni e altre strutture larghe e morbide.
- **Frequenze alte**: Si trovano più lontano dal centro e rappresentano i dettagli fini, come i bordi, le texture e il rumore.

La DFT 2D è separabile lungo le righe e le colonne, il che significa che può essere calcolata come il prodotto di due DFT 1D. Le matrici Fourier separabili lungo le righe e le colonne, quindi, sono calcolate come:

$$
\Psi_{\text{row}}[i, j] = e^{-2\pi i \frac{ij}{n}} \quad \text{e} \quad \Psi_{\text{col}}[i, j] = e^{-2\pi i \frac{ij}{m}}
$$

dove \(i, j\) sono gli indici di riga e colonna. L'immagine $X$ è quindi rappresentata come:

$$
X = \Psi_{\text{row}}\ s\  \Psi_{\text{col}}^H
$$

dove:
- $\Psi_{\text{row}}$ è la matrice della trasformata di Fourier applicata alle righe dell'immagine.
- $\Psi_{\text{col}}^H$ è la matrice Hermitiana della trasformata di Fourier applicata alle colonne.
- $s$ è il vettore dei coefficienti sparsi nel dominio delle frequenze.

L'operazione di trasformazione del segnale in frequenza avviene tramite il prodotto matrice-coefficiente:

$$
s = \Psi_{\text{row}}^H\ X\ \Psi_{\text{col}}^T
$$

La matrice $s$ contiene i coefficienti delle frequenze, cioè l'intensità delle sinusoidi che compongono l'immagine. Se l'immagine ha una struttura semplice o contiene solo poche frequenze rilevanti, i coefficienti $s$ saranno sparsi.

Di seguito viene riportato il codice Python per eseguire tale formulazione matematica:

```python
##################################################################
####################### Fourier matrix (Ψ) #######################
##################################################################
n, m = image.shape
omega_n = np.exp(-2j * np.pi / n)
omega_m = np.exp(-2j * np.pi / m)

# Separable Fourier matrix along rows and columns.
Psi_row = np.array([[omega_n ** (i * j) for j in range(n)] for i in range(n)]) / np.sqrt(n)
Psi_col = np.array([[omega_m ** (i * j) for j in range(m)] for i in range(m)]) / np.sqrt(m)

# Sparse matrix: s = Ψ_row^H * x * Ψ_col^T
s = np.dot(Psi_row.conj().T, np.dot(image, Psi_col.T))

# Image reconstruction
image_reconstructed = np.real(np.dot(Psi_row, np.dot(s, Psi_col.conj().T)))
```
![Figure_13](https://github.com/user-attachments/assets/c0033707-04f4-4b09-957d-ed848eff6760)


# Trasformata di Wavelet 2D

La Trasformata Wavelet Discreta (DWT) è una tecnica potente per analizzare segnali e immagini, che divide un segnale in componenti a diverse risoluzioni. Essa consente di rappresentare i segnali in modo sparso, concentrando la maggior parte delle informazioni nelle basse frequenze, e riducendo l'importanza delle alte frequenze. 

In un contesto unidimensionale, la DWT è utilizzata per separare un segnale $x[n]$ in due componenti principali:
- approssimazione (a bassa frequenza), tramite un filtro passa-basso
$$A[n]=∑_k​x[k]h[2n−k]$$
- dettagli (a alta frequenza), tramite un filtro passa-alto
$$D[n]=∑_k​x[k]g[2n−k]$$

La DWT 2D estende la DWT 1D alle immagini bidimensionali (matrici di dati). In questa versione 2D, si applica la DWT separatamente alle righe e alle colonne dell'immagine:
1. DWT sulle righe: ogni riga dell'immagine viene sottoposta alla DWT 1D, creando quindi una matrice di approssimazione a bassa frequenza e una matrice di dettaglio ad alta frequenza per ogni riga
2. DWT sulle colonne: si applica la DWT 1D separatamente alle colonne dei risultati ottenuti dalla fase precedente (sia le approssimazioni che i dettagli orizzontali). Si creano quindi quattro sottobande
    - LL: approssimazione verticale nella matrice di approssimazione orizzontale - bassa frequenza
    - LH: approssimazione verticale nella matrice di dettaglio orizzontale - alta frequenza
    - HL: dettaglio verticale nella matrice di approssimazione orizzontale - alta frequenza
    - HH: dettaglio verticale nella matrice di dettaglio orizzontale - alta frequenza

Imposta un livello iterativo maggiore di 1, ad ogni livello viene utilizzata la sottobanda LL del livello precedente e viene decomposta nuovamente in quattro sottobande, conservando la struttura principale del segnale in LL e ottenendo dettagli sempre più fini in LH, HL, HH.

In generale, la DWT porta a una rappresentazione sparsa perché separa efficacemente le frequenze basse dalle alte, consentendo di mantenere solo i coefficienti più significativi per la ricostruzione del segnale.

```python
##################################################################
####################### Wavelet matrix (Ψ) #######################
##################################################################
coeffs = pywt.wavedec2(image, wavelet='db1', level=4)

# Creating the figure to display the results
fig, axs = plt.subplots(4, 4, figsize=(12, 12))

for i in range(1, len(coeffs)):
    # Extract details for each layer
    c = coeffs[i]
    
    # Extract LH, HL, HH
    if isinstance(c, tuple):
        LH, HL, HH = c
    else:
        LH = HL = HH = c
        
    # LL
    axs[0, i-1].imshow(coeffs[0], cmap='gray')
    axs[0, i-1].set_title(f"LL - Level {i}")
    axs[0, i-1].axis('off')
    # LH
    axs[1, i-1].imshow(LH, cmap='gray')
    axs[1, i-1].set_title(f"LH - Level {i}")
    axs[1, i-1].axis('off')
    # HL
    axs[2, i-1].imshow(HL, cmap='gray')
    axs[2, i-1].set_title(f"HL - Level {i}")
    axs[2, i-1].axis('off')
    # HH
    axs[3, i-1].imshow(HH, cmap='gray')
    axs[3, i-1].set_title(f"HH - Level {i}")
    axs[3, i-1].axis('off')

plt.tight_layout()
plt.show()
```
![Figure_11](https://github.com/user-attachments/assets/bda5d0b1-90de-4178-8be5-62bb49e10077)


# Trasformata Discreta del Coseno 2D

La Trasformata Discreta del Coseno (DCT) è una tecnica di trasformazione ampiamente utilizzata in elaborazione dei segnali e immagini, specialmente per la compressione. La DCT è simile alla trasformata di Fourier, ma utilizza solo funzioni coseno, senza i termini complessi di Fourier. 

La DCT per una sequenza di segnali discreti $x_0, x_1,…, x_N−1$ ​ è definita come:

$$X_k = \sum_{n=0}^{N-1} x_n \cdot \cos\left(\frac{\pi (2n + 1) k}{2N}\right), \quad k = 0, 1, \dots, N-1
$$

dove:
- $x_n$ ​sono i campioni del segnale originale.
- $X_k$ ​ sono i coefficienti della DCT.
- $N$ è il numero di campioni.
- $k$ è l'indice della frequenza nella trasformata DCT.

Nel caso delle immagini, la Trasformata Discreta del Coseno 2D è applicata separatamente sulle righe e sulle colonne di un'immagine.

$$X_{k_1, k_2} = \sum_{n_1=0}^{N_1-1} \sum_{n_2=0}^{N_2-1} x_{n_1, n_2} \cdot \cos\left(\frac{\pi (2n_1 + 1) k_1}{2N_1}\right) \cdot \cos\left(\frac{\pi (2n_2 + 1) k_2}{2N_2}\right)$$

Nell'elaborazione di immagini la DCT (come anche FFT e Wavelet) viene applicata su blocchi locali (tipicamente di dimensione $8 \cdot 8$) piuttosto che sull'intera immagine. Questa poichè:
- riduce la complessità computazionale: l'applicazione della DCT su un'intera immagine $M \times N$ richiede un costo computazionale pari a $O(MN \log(MN))$, mentre invece, suddividendo l'immagine in blocchi di dimensione $B \times B$, la complessità complessiva diventa $O(MN \log B)$, con $B \ll \min(M, N)$
- maggiore sparsità: l'applicazione della DCT su blocchi locali consente una rappresentazione più accurata e sparsa in ogni blocco, favorendo quindi il compressed sensing.


```python

##################################################################
########################## DCT transform #########################
##################################################################
h, w = image.shape
dct_image = np.zeros_like(image)

plt.figure(figsize=(12, 6))
plt.suptitle("DCT sparsity")
plt.subplot(1, 4, 1)
plt.imshow(image, cmap='gray')
plt.title("Original image")
plt.axis('off')

# Block iterations
def dct_transform(block_size, image_wind):
    for i in range(0, h, block_size):
        for j in range(0, w, block_size):
            # Extract the block
            block = image[i:i+block_size, j:j+block_size]
            # Apply DCT to the block
            dct_image[i:i+block_size, j:j+block_size] = dct(dct(block.T, norm='ortho').T, norm='ortho')
    plt.subplot(1, 4, image_wind)
    plt.imshow(np.log(np.abs(dct_image) + 1), cmap='hot')
    plt.title(f"DCT image - block_size {block_size}")
    plt.axis('off')

dct_transform(block_size=4, image_wind=2)
dct_transform(block_size=150, image_wind=3)
dct_transform(block_size=256, image_wind=4)
    
plt.show()
```

![Figure_12](https://github.com/user-attachments/assets/652e4226-5169-408e-bf5b-4d1885490f2e)


## Confronto tra FFT2D, Wavelet e DCT: sparsità e complessità computazionale

Nel contesto del sensing sparso e della compressione delle immagini, è stato effettuato un confronto tra tre trasformate: FFT2D, e DCT. L'obiettivo del confronto è analizzare le seguenti caratteristiche:

1. **Sparsità delle trasformate**
2. **Complessità computazionale**

Le tre trasformate sono state applicate all'immagine senza suddividerla in blocchi, essendo di fatto la modalità di esecuzione in ambito clinico.

### 1. **Sparsità delle Trasformate**

La sparsità di una trasformata è misurata come la percentuale di valori normalizzati che sono inferiori a una soglia predefinita. La normalizzazione avviene dividendo i valori assoluti della trasformata per il valore massimo della trasformata stessa, ottenendo così valori compresi tra 0 e 1. Le soglie utilizzate per il calcolo della sparsità sono $t = 0$ e $t = 0.00001$, e la **sparsità** viene calcolata come:

$$Sparsità = \frac{\text{Numero di elementi sotto la soglia}}{\text{Numero totale di elementi}} \times 100$$

### 2. **Complessità computazionale**

La complessità computazionale di ciascuna trasformata è stata misurata in termini di tempo di esecuzione. La complessità di ciascun metodo, espressa in notazione O-grande, è la seguente:

- **FFT2D**:
  La **FFT bidimensionale** è una versione ottimizzata della **Trasformata di Fourier**. La sua complessità computazionale è:

  $$O(n^2 \log n)$$

  Questo è dovuto al fatto che la FFT riduce il numero di operazioni necessarie separando i calcoli sulle righe e sulle colonne dell'immagine.

- **Wavelet**:
  La **trasformata Wavelet** decomprime l'immagine in una serie di coefficienti a diverse scale. La sua complessità computazionale è:

  $$O(n^2)$$

  La Wavelet è un processo iterativo che esegue una trasformata multi-risoluzione, ma richiede meno operazioni rispetto a una trasformata completa.

- **DCT**:
  La **Trasformata Coseno Discreta** è utilizzata in applicazioni di compressione, come nel JPEG. La sua complessità computazionale è simile a quella della Wavelet:

  $$O(n^2)$$

  La DCT esegue la trasformazione nel dominio spaziale dell'immagine, ma è generalmente più costosa rispetto alla FFT.

### Risultati

Come si può notare, la Wavelet presenta una peculiarità: con un threshold pari a zero, molti dei suoi coefficienti sono esattamente nulli, il che comporta un valore di sparsità relativamente alto anche quando la soglia è nulla. 

Al contrario, nelle trasformate FFT e DCT, con threshold pari a zero nessun coefficiente rientra in tale range. Quando si aumenta il valore del threshold, c'è un picco nella sparsità: molti coefficienti vengono eliminati, portando a una riduzione marcata della quantità di dati significativi e aumentando la sparsità in modo considerevole. 


| Method    | Threshold   | Below Threshold (%)   | Total Size  
| --------- | ----------- | --------------------- | ------------
| FFT2D     | 0           | 0.00                  | 65536                
| FFT2D     | 0.001       | 67.42                 | 65536                
| Wavelet2D | 0           | 35.46                 | 65536               
| Wavelet2D | 0.001       | 61.69                 | 65536               
| DCT2D     | 0           | 0.00                  | 65536                
| DCT2D     | 0.001       | 71.95                 | 65536                


![Figure_1](https://github.com/user-attachments/assets/fd645493-db7f-4d98-b5e2-0c8b6a5cf3dc)







