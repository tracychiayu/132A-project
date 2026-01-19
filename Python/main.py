import numpy as np
import matplotlib.pyplot as plt
from detector import zf_detector, mmse_detector, maximum_likelihood_detector, simple_sic_detector, sic_detector, mmse_sic_detector
from channel import gen_rayleigh_channel_matrix
from mapping import gen_rand_qam_symbols, qam_symbols_detection

np.random.seed(42)

SNR_dB = np.arange(10, 42, 2)
SNR = 10 ** (0.1 * SNR_dB)
M = 16
decoder = "SIC"  # ZF, LMMSE, ML, simple SIC, SIC, MMSE SIC decoder

# Channel Parameters
Nt = 4
Nr = 4
mu = 0
var = 1

ERROR_THRESHOLD = 1000
SER = np.zeros_like(SNR_dB, dtype=np.float32)
for i, snr in enumerate(SNR):
    print(f'SNR: {SNR_dB[i]} dB')
    
    # Number of symbol vector error
    errors = 0
    total_vecs = 0
    while errors < ERROR_THRESHOLD:
    # Noise varianve
        N0 = Nt / snr

        # Generate N random symbols
        N = Nt
        total_vecs += N // Nt
        symbol, constellation = gen_rand_qam_symbols(N=N, Nt=Nt, M=M)
        x = symbol.reshape(Nt, N // Nt) 

        if N % Nt != 0:
            raise ValueError(f"{N} is not divisible by {Nt}\n")

        # Generate channel matrix
        H = gen_rayleigh_channel_matrix(mu=mu, var=var, Nr=Nr, Nt=Nt)

        # Generate AWGN noise with shape (Nr, 1)
        noise = np.random.normal(loc=0, scale=np.sqrt(N0/2), size=(Nr, N // Nt)) + 1j * np.random.normal(loc=0, scale=np.sqrt(N0/2), size=(Nr, N // Nt)) 

        y = H @ x + noise  

        if (decoder == "ZF"):
            detected_symbols = zf_detector(H=H, y=y, constellation=constellation)
            detected_symbols = detected_symbols.reshape(Nt, N // Nt)

        elif (decoder == "LMMSE"):
            detected_symbols = mmse_detector(H=H, y=y, SNR=snr, constellation=constellation)
            detected_symbols = detected_symbols.reshape(Nt, N // Nt)
        
        elif (decoder == "ML"):
            detected_symbols = maximum_likelihood_detector(H=H, y=y, constellation=constellation)
        
        # (H = Q @ R)
        elif (decoder == "simple SIC"):
            detected_symbols = simple_sic_detector(H=H, y=y, constellation=constellation)
            detected_symbols = detected_symbols.reshape(Nt, N // Nt)

        # sorted SIC (H @ P = Q @ R)
        elif (decoder == "SIC"):
            detected_symbols = sic_detector(H=H, y=y, constellation=constellation)
            detected_symbols = detected_symbols.reshape(Nt, N // Nt)
        
        # regularized sorted SIC ([H^T I^T]^T @ P = Q @ R)
        elif (decoder == "MMSE SIC"):
            detected_symbols = mmse_sic_detector(H=H, y=y, SNR=snr, constellation=constellation)
            detected_symbols = detected_symbols.reshape(Nt, N // Nt)

        # Increment error as long as we detect any errors in a received vector
        # if np.any(x != detected_symbols):
        if not np.allclose(x, detected_symbols):
            errors += 1


    print('errors: ', errors)
    print('total_vecs: ', total_vecs)
    SER[i] = errors / total_vecs
    print('SER: ', SER[i])
    print('---------------------')

# print("SER: ", SER)
# np.save("SER_MMSE_SIC.npy", SER)

plt.figure()
plt.semilogy(SNR_dB, SER, '-o')
plt.xlabel('SNR (dB)')
plt.ylabel('SER')
plt.ylim([1e-4, 1])
plt.grid(which='both', linestyle='--', color='black', alpha=0.3)
plt.show()

    
    




