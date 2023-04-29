import random
import numpy as np
import matplotlib.pyplot as plt
import math

numOfBits = 10
numOfSamplesPerBit = 10


############################################################
################## Convolution Function ####################
############################################################
def applyConvolution(noisySamples, receivedFilterValues):

    convolutionResultSampledTp = np.zeros(numOfBits)

    if (receivedFilterValues is not None):
        convolutionResult = np.convolve(
            noisySamples.flatten(), receivedFilterValues)
    else:
        convolutionResult = noisySamples.flatten()
    # later what is this??
    for i in range(numOfBits):
        convolutionResultSampledTp[i] = convolutionResult[(
            numOfSamplesPerBit - 1) + numOfSamplesPerBit * i]
    return convolutionResult, convolutionResultSampledTp



############################################################
######################## BER Simulated #####################
############################################################
def calculateBERSimulated(bits, recievedBits):
    # applying lamda thresholed (lamda optimum = 0)
    # receivedSamples = np.ones(numOfBits)
    # receivedSamples += (-2 * (recievedBits < 0))
    receivedSamples = (np.real(recievedBits) < 0).astype(int)
    # calculate probability of error
    error_probability = np.sum(receivedSamples != bits)
    error_probability /= numOfBits
    # later important always error_probability = 1
    return error_probability


# def BER_sim(bits, recievedBits):
#         decoded=np.zeros(numOfBits)
#         k=0
#         for j in range(numOfSamplesPerBit-1,len(recievedBits),numOfSamplesPerBit):
#             if recievedBits[j]>0 and k<numOfBits:
#                 decoded[k]=1
#             k+=1
#         #------------------------------
#         return np.sum(decoded != bits)/numOfBits


############################################################
######################## main function #####################
############################################################
def func(receivedFilterMatched, ramp):

    experimentalBER = []
    theorticalBER = []

    E = 1
    for EOverN0_db in range(-10, 21):
        EOverN0_dec = 10 ** (EOverN0_db / 10)

        # generating random noise
        generatedNoise = np.random.normal(
            # later waht is this??
            0, E/(2*EOverN0_dec), size=numOfBits*numOfSamplesPerBit).reshape((numOfBits, numOfSamplesPerBit))

        # add noise to samples
        noisySamples = samples + generatedNoise

        # apply convolution for the noisy samples
        filteredSamples, recievedBits = applyConvolution(
            noisySamples, receivedFilterMatched)
        # append BER simulated for plotting
        experimentalBER.append(calculateBERSimulated(numOfBits, recievedBits))
        # append BER theortical for plotting

        if ramp == 1:
            theorticalBER.append(math.erfc((3**0.5/2) * EOverN0_dec ** 0.5))
        #later math.erfc(EOverN0_dec ** 0.5) or 0.5*math.erfc(EOverN0_dec ** 0.5) ??
        else: 
          theorticalBER.append(math.erfc(EOverN0_dec ** 0.5))

    # ploting
    plt.figure()
    plt.plot(range(0, filteredSamples.flatten(
    ).shape[0]), filteredSamples.flatten(), label="bit value")

    plt.xlabel('time')
    plt.ylabel('bit')
    plt.title('output with Delta Filter')

    plt.grid()
    plt.show()
    return experimentalBER, theorticalBER


############################################################
############### Generating bits and samples ################
############################################################

# later who is A in G(t)
# later we want to optain lamda optimal  = 0, We want to make 0, 1 same probability
# later ask elmo3ed np.random.choice not equal probability
# bits = np.asarray([(random.randint(0, 1)*2 - 1) for i in range(numOfBits)])


# generate random bits with equal probability
bits = np.random.choice([-1, 1], size=(numOfBits,), p=[1./2, 1./2])
# generate samples in range [numOfBits, numOfSamplesPerBit]
samples = (np.asarray([[bits[i] for i in range(numOfBits)]
           for _ in range(numOfSamplesPerBit)])).T

# receive with matched filter
filter_1 = np.ones(numOfSamplesPerBit)
experimentalBER_1, theorticalBER_1 = func(filter_1, 0)

# receive with no filter
filter_2 = None
experimentalBER_2, theorticalBER_2 = func(filter_2, 0)

# receive with ramp filter
# filter_2 = np.random.uniform(low=0, high=3**0.5, size=numOfSamplesPerBit)
# later improve 
# later we have to review the filter and the result
filter_3 = np.linspace(0, 10, numOfSamplesPerBit)
for i in range(len(filter_3)):
    filter_3[i] = np.sqrt(3) * filter_3[i]
experimentalBER_3, theorticalBER_3 = func(filter_3, 1)


#ploting
plt.figure()
plt.plot(range(-10, 21), experimentalBER_1, label = "first BER experimental")
plt.plot(range(-10, 21), experimentalBER_1, "--", label = "first BER theortical")

plt.plot(range(-10, 21), experimentalBER_2, label = "second experimental BER")
plt.plot(range(-10, 21), experimentalBER_2, "--", label = "second theortical BER")

plt.plot(range(-10, 21), experimentalBER_3, label = "third experimental BER")
plt.plot(range(-10, 21), experimentalBER_3, "--", label = "third theortical BER")

plt.xlabel('E_Over_N0(DB)')
plt.ylabel('Bit Error Rate')
plt.yscale('log')
plt.ylim(10**(-5))
plt.title('Bit Error Rate')

plt.legend()
plt.grid()
plt.show()
