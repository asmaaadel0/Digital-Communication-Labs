import numpy as np
import matplotlib.pyplot as plt
import math

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
    receivedSamples = np.ones(numOfBits)
    # Done later why * -2
    # receivedSamples += (-2 * (recievedBits < 0))
    receivedSamples = np.sign(recievedBits)
    # receivedSamples = (np.real(recievedBits) < 0).astype(int)
    # calculate probability of error
    error_probability = np.sum(receivedSamples != bits)
    error_probability /= numOfBits
    # DONE later important always error_probability = 1
    return error_probability


############################################################
######################## main function #####################
############################################################
def func(receivedFilterMatched, type, plot):

    experimentalBER = []
    theorticalBER = []

    E = 1
    for EOverN0_db in range(-10, 21):
        EOverN0_dec = 10 ** (EOverN0_db / 10)

        # generating random noise
        generatedNoise = np.random.normal(
            # later waht is this??
            0, E/(2*EOverN0_dec), numOfBits*numOfSamplesPerBit).reshape((numOfBits, numOfSamplesPerBit))

        # add noise to samples
        noisySamples = samples + generatedNoise

        # apply convolution for the noisy samples
        filteredSamples, recievedBits = applyConvolution(
            noisySamples, receivedFilterMatched)
        # append BER simulated for plotting
        experimentalBER.append(calculateBERSimulated(bits, recievedBits))
        # append BER theortical for plotting

        if type == 3:
            theorticalBER.append(0.5*math.erfc((3**0.5/2) * EOverN0_dec ** 0.5))
        #later math.erfc(EOverN0_dec ** 0.5) or 0.5*math.erfc(EOverN0_dec ** 0.5) ??
        else: 
          theorticalBER.append(math.erfc(EOverN0_dec ** 0.5))
    if plot == 1:

        # ploting
        plt.figure()
        plt.plot(range(0, filteredSamples.flatten(
        ).shape[0]), filteredSamples.flatten(), label="bit value")

        plt.xlabel('time')
        plt.ylabel('bit')
        if type == 1:
            plt.title('output with matched Filter')

        elif type == 2:
            plt.title('output with no Filter')    

        else:
            plt.title('output with ramp Filter')        

        plt.grid()
        plt.show()
    return experimentalBER, theorticalBER


#constants
numOfBits = 10
numOfSamplesPerBit = 10

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
experimentalBER_1, theorticalBER_1 = func(filter_1, 1, 1)

# receive with no filter
filter_2 = None
experimentalBER_2, theorticalBER_2 = func(filter_2, 2, 1)

# receive with ramp filter
# filter_3 = np.random.uniform(low=0, high=3**0.5, size=numOfSamplesPerBit)
# later improve 
# later we have to review the filter and the result
filter_3 = np.linspace(0, 10, numOfSamplesPerBit)
for i in range(len(filter_3)):
    filter_3[i] = np.sqrt(3) * filter_3[i]
experimentalBER_3, theorticalBER_3 = func(filter_3, 3, 1)


# constants
numOfBits = 100000
numOfSamplesPerBit = 10

############################################################
############### Generating bits and samples ################
############################################################

# generate random bits with equal probability
bits = np.random.choice([-1, 1], size=(numOfBits,), p=[1./2, 1./2])
# generate samples in range [numOfBits, numOfSamplesPerBit]
samples = (np.asarray([[bits[i] for i in range(numOfBits)]
           for _ in range(numOfSamplesPerBit)])).T


# receive with matched filter
filter_1 = np.ones(numOfSamplesPerBit)
experimentalBER_1, theorticalBER_1 = func(filter_1, 1, 0)

# receive with no filter
filter_2 = None
experimentalBER_2, theorticalBER_2 = func(filter_2, 2, 0)

# receive with ramp filter
filter_3 = np.linspace(0, 10, numOfSamplesPerBit)
for i in range(len(filter_3)):
    filter_3[i] = np.sqrt(3) * filter_3[i]
experimentalBER_3, theorticalBER_3 = func(filter_3, 3, 0)


#ploting
plt.figure()
# DONE later no plotting for theorticalBER_1
plt.plot(range(-10, 21), experimentalBER_1, label = "First BER experimental")
plt.plot(range(-10, 21), theorticalBER_1, "--", label = "First BER theortical")

plt.plot(range(-10, 21), experimentalBER_2, label = "Second experimental BER")
plt.plot(range(-10, 21), theorticalBER_2, "--", label = "Second theortical BER")

plt.plot(range(-10, 21), experimentalBER_3, label = "Third experimental BER")
plt.plot(range(-10, 21), theorticalBER_3, "--", label = "Third theortical BER")

plt.xlabel('E/N0(DB)')
plt.ylabel('Bit Error Rate')
plt.yscale('log')
plt.ylim(10**(-5))
plt.title('Bit Error Rate')

plt.legend()
plt.grid()
plt.show()
