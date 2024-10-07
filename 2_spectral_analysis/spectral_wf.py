# Coded by Hector Torres (NASA-JPL)
# Adapted by Felix Vivant (Scripps-UCSD and ENS Paris-Saclay)
# Program required for calc_spectral_wf.py

def spectra_wf(A,B,d1,d2,d3):
    ## Wavenumber-frequency spectral analysis
    ### Outputs:
    ## A,B,cs,coh,f1,f2,f3,df1,df2,df3,Atofft,Btofft
    ##
    ## where
    ## A is the power spectrum of A
    ## B is the power spectrum of B
    ## cs is the cospectrum
    ## coh is the spectral coherence
    ## f1,f2,f3 are the spectral coordinates
    ## df1,df2,df3 are the spectral resolution
    ## Atofft,Btofft are the variable in the physical space 
    ## used for the FFt, i.e. after detrending and hanning windowing
    
    import numpy as np
    l1,l2,l3 = A.shape
    df1 = 1./(l1*d1)
    df2 = 1./(l2*d2)
    df3 = 1./(l3*d3)
    f1Ny = 1./(2*d1)
    f2Ny = 1./(2*d2)
    f3Ny = 1./(2*d3)
    f1 = np.arange(-f1Ny,f1Ny,df1)
    f2 = np.arange(-f2Ny,f2Ny,df2)
    f3 = np.arange(0,l3/2+1)*df3
    # Spectral window
    # first, the spatial window
    w1 = np.hanning(l1)
    wx = np.matrix(w1)
    w2 = np.hanning(l2)
    wy = np.matrix(w2)
    window_s = np.repeat(np.array(wx.T*wy),l3).reshape(l1,l2,l3)
    # second, the temporal window
    wt = np.hanning(l3)
    window_t = np.repeat(wt,l1*l2).reshape(l3,l2,l1).T
    # ===== Spectral space ======
    corr_fact=d1*d2*d3 
    Atofft=window_s*window_t*A
    Btofft=window_s*window_t*B
    Ahat = np.fft.rfftn(Atofft)*corr_fact
    Bhat = np.fft.rfftn(Btofft)*corr_fact

    # Cospectrum of A and B
    cs = (Ahat*Bhat.conjugate()).real#/((l1*l2*l3)**2)/(df1*df2*df3)
    #
    # Power spectrum of A
    A = (Ahat*Ahat.conjugate()).real
    # Power spectrum of B
    B = (Bhat*Bhat.conjugate()).real
    #:::::
    ### Coherence #####
    coh = cs/(np.sqrt(A)*np.sqrt(B))
    # zero-padding
    coh = np.fft.fftshift(coh.real,axes=(0,1))
    # cospectrum density
    cs = np.fft.fftshift(cs,axes=(0,1))

    # power spectrum density
    A = np.fft.fftshift(A,axes=(0,1))
    B = np.fft.fftshift(B,axes=(0,1))
    return A,B,cs,coh,f1,f2,f3,df1,df2,df3,Atofft,Btofft
