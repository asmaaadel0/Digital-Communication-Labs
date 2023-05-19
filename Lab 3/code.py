import numpy as np
import matplotlib.pyplot as plt

# later req 4
# documantation
# report

def Gram_Schmidt(S1, S2):

    # later check is it ok or abs(s1)**2
    # Energy1 = np.sum(S1 ** 2)
    Energy1 = abs(S1)**2
    Phi1 = S1 / np.sqrt(Energy1)
    # later divide on numOfSamples or not
    S21 = np.sum(S2 * Phi1) / numOfSamples
    g2 = S2 - S21 * Phi1
    # later divide on numOfSamples or not
    Energy2 = np.sum(g2 ** 2) / numOfSamples
    Phi2 = g2 / np.sqrt(Energy2)

    return Phi1, Phi2


def Signal_Space(S, Phi1, Phi2):

    # later divide on numOfSamples or not
    # just difference on scaling
    V1 = np.sum(S * Phi1)
    V2 = np.sum(S * Phi2)

    return V1, V2


numOfSamples = 50
timeAxis = np.linspace(0, 1, numOfSamples)

S1 = np.ones(numOfSamples)
S2 = np.ones(numOfSamples)

S2[int(0.75*numOfSamples):numOfSamples] = -1

plt.figure(1)
plt.plot(timeAxis, S1, linewidth=2)
plt.xlabel("time")
plt.ylabel("S1")
# later change them
plt.vlines(x=0, ymin=0, ymax=1)
plt.vlines(x=1, ymin=0, ymax=1)

plt.title('S1')
plt.legend()

plt.figure(2)
plt.plot(timeAxis, S2, linewidth=2)
plt.xlabel("time")
plt.ylabel("S2")
# later change them
plt.vlines(x=0, ymin=0, ymax=1)
plt.vlines(x=1, ymin=-1, ymax=0)

plt.title('S2')
plt.legend()

Phi1, Phi2 = Gram_Schmidt(S1, S2)

plt.figure(3)
plt.plot(timeAxis, Phi1, linewidth=2)
plt.vlines(x=0, ymin=0, ymax=1)
plt.vlines(x=1, ymin=0, ymax=1)
plt.title('Gram-Schmidt For S1')
plt.xlabel("time")
plt.ylabel("Phi1")
plt.legend()

plt.figure(4)
plt.plot(timeAxis, Phi2, linewidth=2)
plt.vlines(x=0, ymin=0, ymax=Phi2[0])
plt.vlines(x=1, ymin=Phi2[numOfSamples-1], ymax=0)
plt.title('Gram-Schmidt For S2')
plt.xlabel("time")
plt.ylabel("Phi2")
plt.legend()
plt.show()


V11, V12 = Signal_Space(S1, Phi1, Phi2)
V21, V22 = Signal_Space(S2, Phi1, Phi2)

plt.figure(5)
plt.scatter(V11,V12, label='S1', c='g')
plt.scatter(V21,V22 , label='S2', c='y')
plt.title('Signal space')
plt.xlabel("phi1")
plt.ylabel("phi2")
plt.legend()
plt.show()


EOverSigma2_db_arr = [-5, 0, 10]
E = 1  # as energy for S1, S2 = 1

for EOverSigma2_db in EOverSigma2_db_arr:
    
    plt.title('Noise ('+str(EOverSigma2_db)+')dB')
    plt.xlabel("Phi1")
    plt.ylabel("Phi2")
    for i in range(50):
        
        sigma2 = E/(10**(EOverSigma2_db/10))
        W = np.random.normal(0, np.sqrt(sigma2), numOfSamples)
        
        r1 = S1 + W
        r2 = S2 + W

        V11_Req_3, V12_Req_3 = Signal_Space(r1, Phi1, Phi2)
        V21_Req_3, V22_Req_3 = Signal_Space(r2, Phi1, Phi2)

        plt.scatter(V11_Req_3, V12_Req_3, c='g')
        plt.scatter(V21_Req_3, V22_Req_3, c='y')

    plt.legend(["r1", "r2"])
    plt.show()
