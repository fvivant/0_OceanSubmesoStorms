# Coded by Hector Torres (NASA-JPL)
# Adapted by Felix Vivant (Scripps-UCSD and ENS Paris-Saclay)
# Program required for calc_spectral_w.py

def spectra_w(A,B,d1,d2):
    ## Wavenumber spectral analysis
    ### Outputs:
    ## A,B,cs,coh,f1,f2,df1,df2,Atofft,Btofft
    ##
    ## where
    ## A is the power spectrum of A
    ## B is the power spectrum of B
    ## cs is the cospectrum
    ## coh is the spectral coherence
    ## f1,f2 are the spectral coordinates
    ## df1,df2 are the spectral resolution
    ## Atofft,Btofft are the variable in the physical space 
    ## used for the FFt, i.e. after detrending and hanning windowing
    
    import numpy as np
    l1,l2 = A.shape
    df1 = 1./(l1*d1)
    df2 = 1./(l2*d2)
    f1Ny = 1./(2*d1)
    
    f1 = np.arange(-f1Ny,f1Ny,df1)
    f2 = np.arange(0,l2/2+1)*df2
    # Spectral window
    # the spatial window
    w1 = np.hanning(l1)
    wx = np.matrix(w1)
    w2 = np.hanning(l2)
    wy = np.matrix(w2)
    window_s = np.array(wx.T*wy).reshape(l1,l2)

    # ===== Spectral space ======

    corr_fact=d1*d2
    Atofft=window_s*A
    Btofft=window_s*B
    Ahat = np.fft.rfftn(Atofft)*corr_fact
    Bhat = np.fft.rfftn(Btofft)*corr_fact
    
    # Cospectrum of A and B
    cs = (Ahat*Bhat.conjugate()).real
    # cs = (Ahat*Ahat.conjugate()).real
    # Power spectrum of A
    A = (Ahat*Ahat.conjugate()).real
    # Power spectrum of B
    B = (Bhat*Bhat.conjugate()).real
    #:::::
    ### Coherence #####
    coh = cs/(np.sqrt(A)*np.sqrt(B))
  
    # zero-padding
    coh = np.fft.fftshift(coh,axes=(0))
    # cospectrum density
    cs = np.fft.fftshift(cs,axes=(0))

    # power spectrum density
    A = np.fft.fftshift(A,axes=(0))
    B = np.fft.fftshift(B,axes=(0))
    return A,B,cs,coh,f1,f2,df1,df2,Atofft,Btofft

