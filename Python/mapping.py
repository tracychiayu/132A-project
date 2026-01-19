import numpy as np

def gen_rand_qam_symbols(N, Nt, M=4):
    '''
    Generate N random M-QAM symbols with average symbol energy of Es = 1
    For SISO with 1 transmit antenna, the constellation has average symbol energy of 1.
    For MIMO with Nt transmit antennas, the constellation is scaled down by a factor of sqrt(Nt) so that
    the average symbol energy on the transmitter side is 1.
    Input
        - N: number of random symbols to generate
        - Nt: number of transmit antennas
        - M: order of the QAM constellation
    Output
        - symbols: N randomly selected M-QAM symbols
        - constellation: the full M-QAM constellation
    '''

    A = np.sqrt(1.5/(M-1))
    d = int(np.sqrt(M))     
    constellation = np.zeros((d,d), dtype=complex)
    for im in range(0,d):
        for re in range (0,d):
            a = A * (2 * re + 1 - d)
            b = A * (2 * im + 1 - d)
            constellation[d-1-im][re] = a + b * 1j
    
    constellation = constellation.flatten()  # (sqrt(M), sqrt(M)) -> (M,)
    # constellation = constellation / np.sqrt(Nt)
    symbol = np.random.choice(constellation, size=N)

    return symbol, constellation

def qam_symbols_detection(symbols, constellation):
    '''
    Input
        - symbols: array of received symbols (complex)
        - constellation
    Output
        - detected_symbols: array of detected symbols
    '''

    detected_symbols = np.zeros_like(symbols, dtype=complex)

    for i, sym in enumerate(symbols):
        dist = np.abs(sym-constellation)
        nearest_idx = np.argmin(dist)
        detected_symbols[i] = constellation[nearest_idx]
    
    return detected_symbols