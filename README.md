# Compressive sensing for MRI

Il progetto si propone di studiare la compressione delle immagini di risonanza magnetica (MRI), focalizzandosi sulla modalità in cui queste immagini vengono inizialmente acquisite, ovvero nel formato k-space. Il k-space rappresenta una matrice complessa che contiene informazioni spaziali delle immagini acquisite, ma non è direttamente interpretabile come un'immagine visibile. Questo tipo di dato viene successivamente trasformato, tramite una trasformata di Fourier, in un'immagine spaziale. L'obiettivo principale del progetto è sviluppare tecniche per comprimere efficacemente le immagini MRI partendo dal k-space e cercando di mantenere una buona qualità visiva con una riduzione del volume di dati. Il progetto è stato suddiviso in tre sotto-progetti principali: lo studio iniziale della sparsità delle immagini mediche, il campionamento della FFT2 delle immagini MRI, e infine le tecniche di ottimizzazione per la ricostruzione approssimata delle immagini originali.

## Studio della sparsità nelle immagini mediche
Il primo sotto-progetto esplora le possibili modalità per generare la sparsità nelle immagini mediche. In particolare, sono state testate diverse trasformazioni, come la FFT bidimensionale (FFT2D), le wavelet e la Trasformata Coseno Discreta (DCT), con l’obiettivo di identificare quale di queste tecniche possa meglio rappresentare le immagini mediche in forma sparsa. 
- Link del codice sviluppato:
- Link del file README.md relativo:

## Campionamento FFT2D per le immagini MRI
Il secondo sotto-progetto si concentra sulla FFT2D, in quanto le immagini MRI sono già rappresentate in questo formato. In particolare, si studiano le modalità di campionamento di questo dominio per ridurre la quantità di dati necessari. Sono state esplorate diverse tipologie di maschere di campionamento: una maschera gaussiana centrata, una maschera randomica e una maschera che seleziona i valori nella FFT2D che superano una determinata soglia di ampiezza. Inoltre, sono state analizzate due modalità di codifica per la trasmissione dei dati campionati: COO (Coordinate List), che rappresenta i dati non nulli in forma sparsa, e RLE (Run-Length Encoding), che codifica le sequenze di valori uguali in modo compatto.
- Link del codice sviluppato:
- Link del file README.md relativo:

## Tecniche di ottimizzazione per la ricostruzione dell'immagine
Il terzo sotto-progetto si concentra sulle tecniche di ottimizzazione utilizzate per ricostruire o avvicinarsi all'immagine originale, partendo dai dati compressi. Sono stati esplorati due algoritmi principali: ADMM (Alternating Direction Method of Multipliers) e ISTA (Iterative Shrinkage Thresholding Algorithm), con l'obiettivo di ridurre l'errore di ricostruzione. Con tali algoritmi sono state utilizzate due tecniche di regularizzazione: la trasformata Wavelet e la Total Variation, che hanno lo scopo di preservare le caratteristiche principali dell'immagine durante il processo di ricostruzione. 
- Link del codice sviluppato:
- Link del file README.md relativo: