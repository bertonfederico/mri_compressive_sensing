## Campionamento dei dati per la compressione delle immagini

In molte applicazioni, specialmente nella **ricostruzione delle immagini MRI**, è comune ridurre drasticamente la quantità di dati da acquisire, riducendo la risoluzione spaziale dell'immagine originale. Questo approccio viene comunemente chiamato **undersampling** o **campionamento incompleto**. La riduzione dei dati si traduce in una compressione dell'immagine, mantenendo solo una parte delle informazioni originali, ma preservando la qualità dell'immagine finale dopo un opportuno processo di ricostruzione.

### Campionamento nel dominio di Fourier (k-space e FFT2D)

Nel caso delle immagini MRI, i dati acquisiti non sono direttamente immagini spaziali, ma sono misurazioni in quello che viene chiamato il **k-space**. Il k-space è una rappresentazione della frequenza spaziale dell'immagine, che descrive come le varie frequenze (o modelli spaziali) sono distribuite nell'immagine. Per ottenere l'immagine finale, i dati nel k-space devono essere trasformati nel dominio spaziale tramite l'inversa della trasformata di Fourier.

La **trasformata dapida di Fourier bidimensionale** (**FFT2D**) è una tecnica matematica utilizzata per trasformare un'immagine dal dominio spaziale al dominio delle frequenze spaziali (k-space per le immagini MRI). Questa trasformazione consente di rappresentare un'immagine come una combinazione di frequenze (componenti ad alta e bassa frequenza) anziché come una distribuzione di intensità spaziale.

In termini matematici, la **FFT2D** di un'immagine $f(x, y)$ (con $x$ e $y$ che rappresentano le coordinate spaziali dell'immagine) è data da:

$$
F(u, v) = \sum_{x=0}^{M-1} \sum_{y=0}^{N-1} f(x, y) e^{-j2\pi \left(\frac{ux}{M} + \frac{vy}{N}\right)}
$$

Dove:
- $F(u, v)$ sono i coefficienti nel dominio delle frequenze,
- $f(x, y)$ è il valore del pixel nell'immagine spaziale,
- $M$ e $N$ sono le dimensioni dell'immagine,
- $u$ e $v$ rappresentano le frequenze spaziali nelle direzioni $x$ e $y$,
- $j$ è l'unità immaginaria.

Quando un'immagine è trasformata nel dominio delle frequenze con la FFT2D, essa viene rappresentata come una combinazione di frequenze basse e frequenze alte:

- **Frequenze basse**: Si trovano vicino al centro del dominio delle frequenze e rappresentano le informazioni globali dell'immagine, come forme, contorni e altre strutture larghe e morbide.
- **Frequenze alte**: Si trovano più lontano dal centro e rappresentano i dettagli fini, come i bordi, le texture e il rumore.

La frequenza si riferisce alla quantità di variazione spaziale nell'immagine: basse frequenze corrispondono a variazioni lente e ampie, mentre le alte frequenze rappresentano cambiamenti rapidi, come bordi o dettagli molto piccoli.

### Campionamento random o Gaussiano rispetto al centro

Una tecnica comune per ridurre la quantità di dati acquisiti è campionare selettivamente il k-space, ossia acquisire solo alcuni punti di questo dominio invece di raccogliere tutti i dati. Questo processo è realizzato tramite l'uso di una **maschera di campionamento** che mantiene solo una frazione dei punti originali nel k-space.
Il campionamento nel k-space avviene in questi esempi con **campionamento random** o **campionamento gaussiano**. Entrambi questi metodi preservano la struttura centrale del k-space e rimuovono i punti più periferici. Questo approccio ha un'importante motivazione teorica legata alla natura delle immagini MRI e alla trasformata di Fourier bidimensionale (FFT2D).
Le informazioni principali di un'immagine (cioè le componenti spaziali più significative) sono concentrate principalmente intorno al centro del k-space. Pertanto, campionare in modo gaussiano intorno al centro è molto efficace per preservare le caratteristiche principali dell'immagine.





## Wavelet per la decompressione

Le **Wavelet** sono una famiglia di funzioni matematiche che consentono di analizzare un segnale a diverse risoluzioni, e permettono di esaminare il contenuto di un segnale (o immagine) sia a livello globale che nei dettagli più fini. Questo le rende particolarmente utili per applicazioni come il compressione dei dati e la ricostruzione di immagini.

Una **trasformazione Wavelet** è una trasformazione che permette di rappresentare un segnale in termini di componenti a diverse frequenze e risoluzioni. Mentre la FFT divide il segnale in componenti sinusoidali a frequenza fissa, la trasformazione Wavelet permette una localizzazione sia nel dominio delle frequenze che in quello del tempo o spazio. A differenza della FFT, che non fornisce informazioni temporali o spaziali (le componenti sono globali), le Wavelet sono in grado di localizzare le frequenze sia in tempo che in spazio. Le Wavelet possono essere utilizzate per analizzare i dettagli fini di un'immagine, concentrandosi su piccole scale o su grandi scale.

Nel contesto delle immagini, la **Trasformazione Wavelet Discreta** è la versione discreta e computazionalmente efficiente della trasformazione Wavelet. La DWT applica la funzione Wavelet su una griglia discreta di dati (come un'immagine digitale) ed è tipicamente utilizzata per decomporre l'immagine in diverse risoluzioni.
Una **DWT** di un'immagine bidimensionale (come un'immagine MRI) divide l'immagine in coefficenti associati a differenti livelli di risoluzione e dettagli. La decomposizione in DWT produce tipicamente quattro componenti:
- **LL (Low-Low)**: contiene le informazioni a bassa frequenza (strutture globali dell'immagine).
- **LH (Low-High)**: contiene le informazioni orizzontali ad alta frequenza (dettagli orizzontali).
- **HL (High-Low)**: contiene le informazioni verticali ad alta frequenza (dettagli verticali).
- **HH (High-High)**: contiene le informazioni ad alta frequenza (dettagli fini).

Uno degli utilizzi principali delle Wavelet nelle immagini è la compressione dei dati e il denoising. Nel contesto della ricostruzione delle immagini, le Wavelet vengono utilizzate per ottenere una rappresentazione sparsa dell'immagine. Le immagini vengono scomposte in componenti Wavelet, e successivamente si applica una tecnica di soft-thresholding per promuovere la sparsità. Questo è un approccio efficace per la ricostruzione compressa, dove si cerca di ottenere una buona approssimazione dell'immagine originale con pochi dati.





## ISTA (Iterative Shrinkage-Thresholding algorithm)

L'algoritmo **ISTA** è una tecnica iterativa utilizzata per risolvere problemi di ottimizzazione in cui si cerca di minimizzare una funzione obiettivo composta da due termini principali:

1. **Errore di adattamento (termine di data fidelity)**: rappresenta la discrepanza tra la stima e i dati osservati.
2. **Regolarizzazione (penalizzazione della sparsità)**: penalizza soluzioni non sparse, cioè soluzioni in cui molte variabili non sono zero.

L'ISTA è comunemente utilizzato per risolvere il seguente problema di ottimizzazione:

$$
\min_x \frac{1}{2} \| A x - y \|^2_2 + \lambda \| x \|_1
$$

Dove:
- $ A $ è la matrice di trasformazione
- $ y $ è il vettore di dati osservati
- $ x $ è la variabile da ricostruire
- $ \lambda $ è il parametro di regularizzazione che bilancia il trade-off tra l'adattamento ai dati e la sparsità della soluzione.

### Estensione di ISTA con penalizzazione Wavelet e funzione $ \Phi(x) $

L'ISTA tradizionale usa la penalizzazione $ \| x \|_1 $, ovvero la norma $ l_1 $ dei coefficienti dell'immagine. Tuttavia, come già citato il dominio che presenta miglior sparsità è Wavelet.
L'uso delle Wavelet consente di rappresentare l'immagine in termini di una serie di coefficienti che possono essere sparsi. Una trasformata Wavelet scompone un'immagine in una serie di coefficienti a diverse risoluzioni, che descrivono la struttura dell'immagine a diverse scale.

La funzione obiettivo dell'ISTA può quindi essere modificata per applicare la penalizzazione $ l_1 $ sui coefficienti Wavelet, piuttosto che sui valori dell'immagine. Questo porta al seguente problema:

$$
\min_x \frac{1}{2} \| A x - y \|^2_2 + \lambda \| W(x) \|_1
$$

dove $ W(x) $ è la funzione Wavelet applicata all'immagine $ x $, e quindi la norma $ l_1 $ è ora applicata sui coefficienti Wavelet.

Nella versione tradizionale di ISTA, il modello di osservazione è descritto da un sistema lineare, in cui i dati osservati y sono ottenuti applicando una matrice A a un vettore $x$, ovvero $y=Ax+r$, dove $r$ è il rumore. In questo contesto, la matrice $A$ rappresenta un operatore che mappa il vettore $x$ nel dominio delle osservazioni $y$.

La funzione di perdita in questo caso è
$$​L(x)=\frac{1}{2} \|Φ(x)−y\|_2^2$$
dove $Φ(x)$ rappresenta il processo che mappa l'immagine $x$ dal dominio spaziale al dominio delle frequenze tramite la FFT2D, applicando successivamente una maschera.

Dovendo considerare nell'algoritmo ISTA il grdiente dell'errore $r$ da minimizzare, si ha che
$$∇_x​L(x)=Φ_T(Φ(x)−y)$$
dove $Φ_T$ consiste nell'applicazione della maschera e la successiva inversione della FFT2D.

### Passi iterativi dell'algoritmo

Nel codice che hai fornito, vediamo che l'algoritmo segue una struttura molto simile a quella di ISTA con una penalizzazione basata sulle Wavelet:

  1. **Passo di gradiente discendente**: viene calcolato il residuo $ r = y - \Phi(x) $, dove $ \Phi $ rappresenta la trasformazione legata alla FFT e alla maschera di undersampling. Successivamente, viene effettuato un passo di gradiente discendente:
$$
x^{k+1} = x^k - \alpha \nabla f(x^k) = x^k + \alpha \  \Phi_T(r)
$$

  2. **Passo di soft-thresholding sui coefficienti Wavelet**: dopo il passo di gradiente, il codice esegue una decomposizione Wavelet sull'immagine aggiornata $ x^{k+1} $; questa trasformazione scompone l'immagine in una serie di coefficienti a diversi livelli di risoluzione. Successivamente, il soft-thresholding è applicato a ciascun coefficiente: a questo punto, i coefficienti di piccole dimensioni sono ridotti a zero, portando a una soluzione più sparsa.
  $$
  z^{k+1} = wavelet(x^{k+1})
  $$
  $$
  S_\lambda(z^{k+1}) = \text{sign}(z^{k+1}) \cdot \max(|z^{k+1}| - \lambda, 0)
  $$
  3. **Passo di ricostruzione dell'immagine**: una volta che i coefficienti Wavelet sono soggetti a soft-thresholding, l'immagine viene ricostruita attraverso l'inversa della trasformazione Wavelet:
  $$
  x^{k+1} = wavelet_T(S_\lambda(z^{k+1}))
  $$
  4. **Verifica della convergenza**: infine, il codice verifica la convergenza confrontando la differenza tra la soluzione corrente e quella precedente. Se la differenza è inferiore a una tolleranza $ \text{tol} $, o se il numero di iterazioni raggiunge una soglia, l'algoritmo si ferma. 

### Codice dell'algoritmo

```python
# Soft-thresholding function for sparsity regularization
def soft_thresholding(x, threshold):
    """
    Applies soft thresholding to promote sparsity, commonly used in Wavelet-based denoising.

    Parameters:
    - x: numpy array or tuple of arrays, the input data to threshold. If a tuple, each element is processed recursively.
    - threshold: float, the threshold level. Values with absolute magnitude below this will be set to zero.

    Returns:
    - Thresholded result, either as a numpy array or a tuple of arrays, with values shrunk towards zero.
    """
    # If x is a tuple (e.g., Wavelet coefficient structure), apply soft thresholding to each element
    if isinstance(x, tuple):
        return tuple(soft_thresholding(c, threshold) for c in x)
    
    # Apply soft thresholding for individual values in the array
    return np.sign(x) * np.maximum(np.abs(x) - threshold, 0)


# Function for MRI reconstruction using the ISTA algorithm
def ISTA_MRI_reconstruction(y, mask, lam, max_iter=1000, tol=1e-5, step_size = 0.001, Wavelet='db1', level=4):
    """
    Perform MRI image reconstruction using Iterative Shrinkage-Thresholding Algorithm (ISTA).

    Parameters:
    - y: numpy array, observed undersampled k-space data.
    - mask: numpy array, binary mask used for undersampling.
    - lam: float, regularization parameter for controlling sparsity.
    - max_iter: int, maximum number of iterations to perform.
    - tol: float, tolerance for stopping criterion based on convergence.
    - Wavelet: str, type of Wavelet to use for sparsity in the transform domain.
    - level: int, number of decomposition levels for the Wavelet transform.

    Returns:
    - x: numpy array, reconstructed MRI image.
    """
    
    # Initialize the reconstructed image with the zero-filled back-projection of y
    x_or = Phi_T(y, mask)
    x = np.copy(x_or)
    

    for k in range(max_iter):
        # Compute the residual in k-space and perform a gradient step
        r = y - Phi(x, mask)
        x_new = x + step_size * Phi_T(r, mask)

        # Wavelet decomposition of x_new and soft-thresholding of coefficients
        coeffs = pywt.wavedec2(x_new, Wavelet, level=level)
        coeffs_thresh = [soft_thresholding(c, lam * step_size) for c in coeffs]
        
        # Reconstruct the image from the thresholded Wavelet coefficients
        x_new = pywt.waverec2(coeffs_thresh, Wavelet)


        # Check for convergence by comparing the norm of the update difference to the tolerance
        if np.linalg.norm(x_new - x) < tol:
            break
        
        # Update the estimate for the next iteration
        x = x_new

    print(f"ISTA convergence reached at iteration {k+1}")
    print(y * mask - fft2c(x) * mask)
    return x
```






## ADMM con TOTAL VARIATION
L'Alternating Direction Method of Multipliers (ADMM) è una tecnica di ottimizzazione particolarmente adatta per problemi che combinano più termini di costo, ciascuno dei quali richiede una diversa forma di regolarizzazione. Questo approccio si applica bene alla ricostruzione di immagini MRI sottocampionate, dove vogliamo contemporaneamente preservare i dati misurati e ottenere un'immagine "pulita" e priva di artefatti. La penalizzazione di Variazione Totale (TV) è un metodo comunemente usato per questo scopo, in quanto minimizza le variazioni locali senza sfocare eccessivamente i bordi, preservando i dettagli strutturali importanti dell'immagine.
Consideriamo il problema della ricostruzione di un'immagine x da dati incompleti del k-spazio k. L'obiettivo è trovare un'immagine che:
- rispetti i dati misurati nel k-spazio (ovvero le frequenze osservate)
- riduca il rumore e gli artefatti che possono derivare dall'undersampling tramite la regolarizzazione.

Questo può essere formulato come un problema di ottimizzazione che minimizza la seguente funzione obiettivo: 


In forma primale, il problema si pone come
$$
\min_x \frac{1}{2} \| mask \ (FFT2D(x) - k) \|^2_2 + \lambda*TV(x)
$$
dove:
- $mask$ è la maschera di campionamento, una matrice che conserva i valori di k-spazio osservati e ignora quelli non misurati.
- $FFT2D(x)$ è la trasformata di Fourier dell'immagine $x$, che permette di confrontarla con i dati $k$ nel dominio k-spazio.
- $TV(x)$ rappresenta la variazione totale di $x$. Non essendo calcolabile in forma chiusa, anche la minimizzazione della TV viene effettuato iterativamente:

  - calcolo del gradiente nelle direzioni orizzontali e verticali;
  - calcolo della norma totale del gradiente e normalizzazione tramite esso dei gradienti precedentemente calcolati;
  - calcolo della divergenza sia orizzontale che verticale, e combinazione di esse per il calcolo della divergenza totale per ogni pixel;
  - aggiornamento dell'immagine a partire da tale divergenza, utilizzando una variabile peso desiderata;


In forma primale-duale, il problema diventa:
$$
\min_{x, z} \frac{1}{2} \| mask \ (FFT2D(x) - k) \|^2_2 + \lambda*TV(x)\ \ \ \ \ \ \ s.t\ \ \ x = z
$$
che si risolve in modalità iterativa tramite ADMM:
- aggiornamento di $x$: minimizza il termine di fedeltà ai dati mantenendo fisso $z$ e $u$.
- aggiornamento di $z$: minimizza il termine di variazione totale, effettuando una regolarizzazione TV.
- aggiornamento della variabile duale $u$: aggiorna la variabile duale per spingere $x$ e $z$ a coincidere, soddisfacendo progressivamente il vincolo $x=z$.


Anche lo stesso metodo di Total Variation si applica in modalità iterativa in modalità *gradient dessent*; iterativamente:
- calcola i gradienti in entrambe le direzioni per misurare la variazione locale.
- utilizza una norma dei gradienti per calcolare la divergenza e aggiornare l'immagine riducendo le variazioni eccessive, preservando i bordi.

```python
import numpy as np
from domain_transforms.fft_transforms import fft2c, ifft2c

def _total_variation(image, weight, num_iters=100, tau=0.125):
    """
    Applies denoising using Total Variation (TV) regularization with Gradient Descent.

    Parameters:
        image (numpy.ndarray): Input image to be denoised.
        weight (float): Weight for the TV penalty (higher values increase smoothing).
        num_iters (int): Number of gradient descent iterations.
        tau (float): Gradient descent step size.

    Returns:
        numpy.ndarray: Denoised/reconstructed image.
    """
    # Initialize the working image copy
    img = np.copy(image)
    
    # Iteratively apply gradient descent for TV minimization
    for _ in range(num_iters):
        
        # Calculate gradients in x and y directions (differences between neighboring pixels)
        grad_x = np.roll(img, -1, axis=1) - img  # Gradient along x-axis
        grad_y = np.roll(img, -1, axis=0) - img  # Gradient along y-axis

        # Compute the gradient magnitude for isotropic TV norm
        # Adding a small constant (1e-8) for numerical stability
        grad_norm = np.sqrt(grad_x**2 + grad_y**2 + 1e-8)
        
        # Calculate divergence (TV regularization term)
        # Dividing by grad_norm and computing divergence to achieve smoothness
        div_x = (grad_x / grad_norm) - np.roll(grad_x / grad_norm, 1, axis=1)
        div_y = (grad_y / grad_norm) - np.roll(grad_y / grad_norm, 1, axis=0)
        divergence = div_x + div_y
        
        # Update the image by moving in the direction of the negative gradient
        # The term (tau * weight * divergence) controls the degree of smoothing
        img = img + tau * weight * divergence

    return img


def admm_mri_tv_reconstruction(kspace_sub, mask, lam=0.1, num_iters=50):
    """
    Reconstructs an MRI image from undersampled k-space data using ADMM with Total Variation (TV) regularization.

    Parameters:
        kspace_sub (numpy.ndarray): Undersampled k-space data.
        mask (numpy.ndarray): Sampling mask for k-space (1 for sampled, 0 for unsampled).
        lam (float): Regularization parameter for TV denoising.
        num_iters (int): Number of ADMM iterations.

    Returns:
        numpy.ndarray: Reconstructed image in the spatial domain.
    """
    
    # Initial estimate for the image: inverse Fourier transform of the undersampled k-space
    img = ifft2c(kspace_sub)
    
    # Initialize auxiliary variable z and dual variable u
    z = np.copy(img)
    u = np.zeros_like(img)

    # ADMM Iterative reconstruction loop
    for _ in range(num_iters):
        
        # 1. Data fidelity step:
        # Apply the mask in the Fourier domain to retain sampled values and replace unsampled values
        img = ifft2c(mask * kspace_sub + (1 - mask) * fft2c(z - u))
        
        # 2. Regularization step:
        # Apply TV denoising to the real part of (img + u) to encourage spatial smoothness
        z = _total_variation(np.real(img + u), weight=lam)
        
        # Convert z back to complex type for further calculations
        z = z.astype(np.complex128)
        
        # 3. Dual variable update:
        # Adjust the dual variable u to account for the difference between img and z
        u += img - z

    # Return the final reconstructed image in the spatial domain
    return img
```





## ADMM con Wavelet
Mentre la regularizzazione TV è efficace nel preservare i contorni attraverso penalizzazioni sui gradienti e garantire una certa regolarità nell'immagine, la regularizzazione basata su Wavelet si dimostra particolarmente utile in applicazioni come la ricostruzione delle immagini MRI, dove le immagini tendono ad avere caratteristiche a più scale che possono essere rappresentate in modo sparso nel dominio delle Wavelet.

Le Wavelet forniscono un modo potente di rappresentare i dati dell'immagine con coefficienti sparsi. Il principale vantaggio di utilizzare una regularizzazione Wavelet è che, rispetto alla regularizzazione TV, la Wavelet non solo preserva i bordi, ma è anche capace di adattarsi meglio a strutture che cambiano a diverse risoluzioni (come i dettagli finissimi nelle immagini MRI).

L'algoritmo ADMM con regularizzazione Wavelet differisce da quello con regularizzazione TV  nella fase di regolarizzazione. Invece di applicare una penalizzazione sui gradienti dell'immagine (come nel caso del TV), applica una sogliatura ai coefficienti di Wavelet dell'immagine.

In forma primale, il problema si pone come
$$
\min_x \frac{1}{2} \| mask*(FFT2D(x) - y) \|^2_2 + \lambda \| W(x) \|_1\ \ \
$$
mentre in forma primale-duale
$$
\min_{x, z} \frac{1}{2} \| mask*(FFT2D(x) - y) \|^2_2 + \lambda \| z \|_1\ \ \ \ \ \ \ s.t\ \ \ W(x) = z
$$
che si risolve in modalità iterativa tramite ADMM:
- passo di fedeltà, come esposto nel caso TV.
- passo di regolarizzazione: in questo caso si applica il denoising basato su Wavelet. Questo passaggio sfrutta la decomposizione Wavelet dell'immagine, applica la sogliatura ai coefficienti di dettaglio e ricostruisce l'immagine utilizzando i coefficienti sogliati.
- aggiornamento della variabile duale: come nel caso del TV, la variabile duale viene aggiornata per correggere.

La principale differenza tra la regularizzazione TV e quella Wavelet sta nella natura della penalizzazione. Mentre la TV penalizza la variazione dei valori di pixel (favorendo l'ombreggiatura uniforme), la Wavelet penalizza i coefficienti di dettaglio ad alta frequenza, favorendo una rappresentazione sparsa delle strutture a più scale.

```python
import numpy as np
import pywt
from domain_transforms.fft_transforms import ifft2c, fft2c

def Wavelet_denoise(image, Wavelet='db1', level=4, threshold=0.1):
    """
    Applies soft-thresholding denoising to Wavelet coefficients.

    Parameters:
        image (numpy.ndarray): Input image to be denoised.
        Wavelet (str): Type of Wavelet to use for decomposition, e.g., 'db1' (Daubechies).
        level (int): Number of decomposition levels.
        threshold (float): Threshold value for soft-thresholding.

    Returns:
        numpy.ndarray: Denoised image.
    """
    # Decompose the image into Wavelet coefficients
    coeffs = pywt.wavedec2(image, Wavelet=Wavelet, level=level)
    
    # Apply soft-thresholding to the detail coefficients
    coeffs_thresh = [coeffs[0]]  # Keep approximation coefficients unchanged
    for detail_level in coeffs[1:]:
        # Apply soft-thresholding to each tuple of detail coefficients
        coeffs_thresh.append(tuple(pywt.threshold(c, threshold, mode='soft') for c in detail_level))
    
    # Reconstruct the image from the thresholded coefficients
    return pywt.waverec2(coeffs_thresh, Wavelet)


def admm_mri_Wavelet_reconstruction(kspace_sub, mask, num_iters=100, Wavelet='db1', level=4, threshold=0.1):
    """
    Reconstructs an MRI image from undersampled k-space data using ADMM with Wavelet-based denoising.

    Parameters:
        kspace_sub (numpy.ndarray): Undersampled k-space data.
        mask (numpy.ndarray): Sampling mask for k-space (1 for sampled, 0 for unsampled).
        num_iters (int): Number of ADMM iterations.
        Wavelet (str): Type of Wavelet for Wavelet denoising.
        level (int): Number of decomposition levels for Wavelet transform.
        threshold (float): Threshold for Wavelet coefficient denoising.

    Returns:
        numpy.ndarray: Reconstructed image in the spatial domain.
    """
    
    # Initial estimate for the image: inverse Fourier transform of the undersampled k-space
    img = ifft2c(kspace_sub)
    
    # Initialize auxiliary variable z and dual variable u
    z = np.copy(img)
    u = np.zeros_like(img)

    # ADMM iterative reconstruction loop
    for i in range(num_iters):
        
        # 1. Data fidelity step:
        # Update img by enforcing data consistency in the Fourier domain
        # This step retains sampled k-space values and fills in unsampled areas
        img = ifft2c(mask * kspace_sub + (1 - mask) * fft2c(z - u))

        # 2. Regularization step:
        # Apply Wavelet denoising to the real part of (img + u) to encourage sparsity in the Wavelet domain
        z = Wavelet_denoise(np.real(img + u), Wavelet=Wavelet, level=level, threshold=threshold)

        # Convert z back to complex type for compatibility with further calculations
        z = z.astype(np.complex128)

        # 3. Dual variable update:
        # Adjust the dual variable u based on the difference between img and z
        # This step enforces consistency between the current estimate and regularized solution
        u += img - z

    # Return the final reconstructed image in the spatial domain
    return img
```