# Sparsità del segnale MRI

La **rappresentazione sparsa** è un concetto fondamentale nell'elaborazione dei segnali, particolarmente utile per comprimere e rappresentare segnali che contengono una gran parte di informazioni in un numero ridotto di componenti significative. In generale, un segnale $x$ può essere rappresentato in una base (o dizionario) $\Psi$ come una combinazione lineare di pochi elementi:

$$
x = \Psi s
$$

Dove:
- $x$ è il segnale originale.
- $\Psi$ è una matrice che rappresenta una **base ortonormale**.
- $s$ è un vettore di **coefficienti sparsi**, cioè un insieme di pesi che combinano le colonne di $\Psi$ per ricostruire $x$.

Nel caso di segnali naturali, come MRI, la loro rappresentazione sparsa è spesso ben approssimata in basi come **Wavelet** o **Fourier**.

# Applicazione alla trasformata di Fourier 2D

Nel caso delle immagini, la teoria generale della sparsità può essere applicata utilizzando la **Trasformata Discreta di Fourier (DFT)**, che è una base ortonormale per rappresentare un'immagine nel dominio delle frequenze.

Nel caso delle immagini MRI, i dati acquisiti non sono direttamente immagini spaziali, ma sono misurazioni in quello che viene chiamato il **k-space**. Il k-space è una rappresentazione della frequenza spaziale dell'immagine, che descrive come le varie frequenze (o modelli spaziali) sono distribuite nell'immagine. Per ottenere l'immagine finale, i dati nel k-space devono essere trasformati nel dominio spaziale tramite l'inversa della trasformata di Fourier.

La **trasformata di Fourier bidimensionale** (**DFT2D**) è una tecnica matematica utilizzata per trasformare un'immagine dal dominio spaziale al dominio delle frequenze spaziali (k-space per le immagini MRI). Questa trasformazione consente di rappresentare un'immagine come una combinazione di frequenze (componenti ad alta e bassa frequenza) anziché come una distribuzione di intensità spaziale.

In termini matematici, la **DFT2D** di un'immagine $f[x, y]$ (con $x$ e $y$ che rappresentano le coordinate spaziali dell'immagine) è data da:

$$
F[u, v] = \sum_{x=0}^{M-1} \sum_{y=0}^{N-1} f[x, y] e^{-j2\pi \left(\frac{ux}{M} + \frac{vy}{N}\right)}
$$

Dove:
- $F[u, v]$ sono i coefficienti nel dominio delle frequenze,
- $f[x, y]$ è il valore del pixel nell'immagine spaziale,
- $M$ e $N$ sono le dimensioni dell'immagine,
- $u$ e $v$ rappresentano le frequenze spaziali nelle direzioni $x$ e $y$,
- $j$ è l'unità immaginaria.

Quando un'immagine è trasformata nel dominio delle frequenze con la DFT2D, essa viene rappresentata come una combinazione di frequenze basse e frequenze alte:

- **Frequenze basse**: Si trovano vicino al centro del dominio delle frequenze e rappresentano le informazioni globali dell'immagine, come forme, contorni e altre strutture larghe e morbide.
- **Frequenze alte**: Si trovano più lontano dal centro e rappresentano i dettagli fini, come i bordi, le texture e il rumore.

Supponiamo di avere un'immagine $X$ di dimensioni $N \times M$. La sua rappresentazione in frequenza $s$ può essere ottenuta tramite la Trasformata di Fourier 2D. L'immagine $X$ è quindi rappresentata come:

$$
X = \Psi_{\text{row}} s \Psi_{\text{col}}^H
$$

dove:
- $\Psi_{\text{row}}$ è la matrice della trasformata di Fourier applicata alle righe dell'immagine.
- $\Psi_{\text{col}}^H$ è la matrice Hermitiana della trasformata di Fourier applicata alle colonne.
- $s$ è il vettore dei coefficienti sparsi nel dominio delle frequenze.

La matrice della DFT2D è separabile, ovvero possiamo applicare la trasformata di Fourier separatamente sulle righe e sulle colonne dell'immagine:

- La matrice $\Psi_{\text{row}}$ contiene le basi della trasformata di Fourier per le righe.
- La matrice $\Psi_{\text{col}}$ contiene le basi della trasformata di Fourier per le colonne.
  
L'operazione di trasformazione del segnale in frequenza avviene tramite il prodotto matrice-coefficiente:

$$
s = \Psi_{\text{row}}^H \cdot X \cdot \Psi_{\text{col}}^T
$$

dove:
- $\Psi_{\text{row}}^H$ è la trasposta coniugata di $\Psi_{\text{row}}$.
- $\Psi_{\text{col}}^T$ è la trasposta di $\Psi_{\text{col}}$.

La matrice $s$ contiene i coefficienti delle frequenze, cioè l'intensità delle sinusoidi che compongono l'immagine. Se l'immagine ha una struttura semplice o contiene solo poche frequenze rilevanti, i coefficienti $s$ saranno sparsi (molti di essi saranno zero o vicini a zero).

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

Come nel caso monodimensionale, anche nel caso bidimensionale la complessità computazionale risulta migliore utilizzando FFT2D. Se la matrice di pixel su cui applicare la trasformata è di dimensione, si ha una complessità:
- DFT2D: $O(M^2*N^2)$
- FFT2D: $O(M*N*log(M*N))$

Per questo motivo, in python viene usualmente utilizzata la libreria *scipy.fft*.

La dimensione della matrice $s$ risultante rimane invariata rispetto alla matrice dell'immagine di partenza, mentre la sparsità (ovvero la percentuale di valori al di sotto di $ϵ = 0.1$) risulta essere intorno al $90.0 % $.

# Applicazione alla trasformata di Wavelet

La Trasformata Wavelet Discreta (DWT) è una tecnica potente per analizzare segnali e immagini, che divide un segnale in componenti a diverse risoluzioni. Essa consente di rappresentare i segnali in modo sparso, concentrando la maggior parte delle informazioni nelle basse frequenze, e riducendo l'importanza delle alte frequenze. 

In un contesto unidimensionale, la DWT è utilizzata per separare un segnale $x[n]$ in due componenti principali:
- approssimazione (a bassa frequenza).
- dettagli (a alta frequenza).

Nel caso di un segnale unidimensionale $x[n]$, la DWT separa il segnale utilizzando due filtri:
- filtro passa-basso $h[n]$: cattura le componenti a bassa frequenza.

$$A[n]=∑_k​x[k]h[2n−k]$$

- filtro passa-alto $g[n]$: cattura le componenti ad alta frequenza.

$$D[n]=∑_k​x[k]g[2n−k]$$

La DWT 2D estende la DWT 1D alle immagini bidimensionali (matrici di dati), separando il segnale in componenti orizzontali e verticali. In questa versione 2D, si applica la DWT separatamente alle righe e alle colonne dell'immagine.
1. DWT sulle righe: ogni riga dell'immagine viene sottoposta alla DWT 1D, creando quindi una matrice di approssimazione a bassa frequenza e una matrice di dettaglio ad alta frequenza per ogni riga
2. DWT sulle colonne: si applica la DWT 1D separatamente alle colonne dei risultati ottenuti dalla fase precedente (sia le approssimazioni che i dettagli orizzontali). Si creano quindi quattro sottobande
    - LL: approssimazione verticale della matrice di approssimazione orizzontale - bassa frequenza
    - LH: approssimazione verticale della matrice di dettaglio orizzontale - alta frequenza
    - HL: dettaglio verticale della matrice di approssimazione orizzontale - alta frequenza
    - HH: dettaglio verticale della matrice di dettaglio orizzontale - alta frequenza

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

In sintesi, i risultati di questo algoritmo sono i seguenti:
- la complessità computazionale, in vantaggio rispetto a FFT2D, risulta $O(N*M*)$.
- la dimensione delle matrici di risultato (LL, LH, HL, HH) dimezza ad ogni iterazione.
- la sparsità di segnali MRI, rispetto a $ϵ = 0.1$, risulta essere intorno al $0.45 %$.