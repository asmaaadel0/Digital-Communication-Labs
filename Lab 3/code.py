import numpy as np
import matplotlib.pyplot as plt

# later req 4
# documantation
# report

###########################################################
############## Gram-Schmidt Orthogonalization #############
###########################################################


def Gram_Schmidt(S1, S2):

    # later check is it ok or abs(s1)**2
    # calculate the energy of the s1(t)
    Energy1 = np.sum(S1 ** 2) / numOfSamples
    # calculate Gram-Schmidt base function Phi1
    Phi1 = S1 / np.sqrt(Energy1)
    # later divide on numOfSamples or not
    S21 = np.sum(S2 * Phi1) / numOfSamples
    # calculate Gram-Schmidt second orthogonalized
    g2 = S2 - S21 * Phi1
    # later divide on numOfSamples or not
    # calculate the energy of the s2(t)
    Energy2 = np.sum(g2 ** 2) / numOfSamples
    # calculate Gram-Schmidt base function Phi2
    Phi2 = g2 / np.sqrt(Energy2)

    return Phi1, Phi2


###########################################################
############### Signal Space representation ###############
###########################################################
def Signal_Space(S, Phi1, Phi2):

    # later divide on numOfSamples or not
    # just difference on scaling
    # calculate the signal space representation
    V1 = np.sum(S * Phi1)
    V2 = np.sum(S * Phi2)

    return V1, V2


# generate time axis
numOfSamples = 50
timeAxis = np.linspace(0, 1, numOfSamples)

# generate s1(t)
S1 = np.ones(numOfSamples)

# generate s2(t)
S2 = np.ones(numOfSamples)
S2[int(0.75*numOfSamples):numOfSamples] = -1

# plot s1(t) vs time
plt.figure(1)
plt.plot(timeAxis, S1, linewidth=2)
plt.xlabel("time")
plt.ylabel("S1(t)")
# later change them
plt.vlines(x=0, ymin=0, ymax=1)
plt.vlines(x=1, ymin=0, ymax=1)
plt.title('S1(t)')
plt.legend()

# plot s2(t) vs time
plt.figure(2)
plt.plot(timeAxis, S2, linewidth=2)
plt.xlabel("time")
plt.ylabel("S2(t)")
# later change them
plt.vlines(x=0, ymin=0, ymax=1)
plt.vlines(x=1, ymin=-1, ymax=0)
plt.title('S2(t)')
plt.legend()

# get the bases functions of s1(t) & s2(t)
Phi1, Phi2 = Gram_Schmidt(S1, S2)

# plot Phi1 vs time
plt.figure(3)
plt.plot(timeAxis, Phi1, linewidth=2)
plt.vlines(x=0, ymin=0, ymax=1)
plt.vlines(x=1, ymin=0, ymax=1)
plt.title('Gram-Schmidt For S1(t)')
plt.xlabel("time")
plt.ylabel("Phi1(t)")
plt.legend()

# plot Phi2 vs time
plt.figure(4)
plt.plot(timeAxis, Phi2, linewidth=2)
plt.vlines(x=0, ymin=0, ymax=Phi2[0])
plt.vlines(x=1, ymin=Phi2[numOfSamples-1], ymax=0)
plt.title('Gram-Schmidt For S2(t)')
plt.xlabel("time")
plt.ylabel("Phi2(t)")
plt.legend()
plt.show()

# get the signal space representation of s1(t) & s2(t)
V11, V12 = Signal_Space(S1, Phi1, Phi2)
V21, V22 = Signal_Space(S2, Phi1, Phi2)

# Plot the signal space representation
plt.figure(5)
plt.scatter(V11, V12, label='S1', c='g')
plt.scatter(V21, V22, label='S2', c='y')
plt.title('Signal space')
plt.xlabel("Phi1(t)")
plt.ylabel("Phi2(t)")
plt.legend()
plt.show()

# generate E/σ2 array in db
EOverSigma2_db_arr = [-5, 0, 10]
E = 1  # as energy for S1, S2 = 1

# for loop for each element in E/σ2 array
for EOverSigma2_db in EOverSigma2_db_arr:

    # plot noise signal Phit1 vs Phi2
    plt.title('Noise ('+str(EOverSigma2_db)+')dB')
    plt.xlabel("Phi1(t)")
    plt.ylabel("Phi2(t)")
    # for loop 50 times for random noise samples
    for i in range(50):

        # calculate standard deviation
        standardDev = E/(10**(EOverSigma2_db/10))
        # generate random noise samples
        W = np.random.normal(0, np.sqrt(standardDev), numOfSamples)

        # add noise for signals
        r1 = S1 + W
        r2 = S2 + W

        # get the signal space representation of r1(t) & r2(t)
        V11_Req_3, V12_Req_3 = Signal_Space(r1, Phi1, Phi2)
        V21_Req_3, V22_Req_3 = Signal_Space(r2, Phi1, Phi2)

        # plot the signal space representation
        plt.scatter(V11_Req_3, V12_Req_3, c='g')
        plt.scatter(V21_Req_3, V22_Req_3, c='y')

    plt.legend(["r1", "r2"])
    plt.show()
