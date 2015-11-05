#!/usr/bin/env python
'''
Created on Mar 31, 2015

@author: thrabe
'''
import numpy
unitStartDictionary = {}

unitStartDictionary['2OMZ'] = [51,73,95,117,139,161,183,205,227,249,271,293] 
unitStartDictionary['1AWC'] = [3,40,73,106]

unitLenghts = {}
unitLenghts['2OMZ'] = 22
unitLenghts['1AWC'] = 33



class NonSolenoidException(Exception):
    pass

def filter(signal,band,doPlot = False):
    from numpy import fft
    
    sfft = fft.fft(signal)
    
    shifted = fft.fftshift(sfft)
    
    filter = numpy.zeros(signal.shape[-1])
    
    filter[signal.shape[-1]/2 - band : signal.shape[-1]/2 + band] = 1
    
#    filter[0:10] = 1
#    filter[-10:-1] = 1

    filtered = shifted * filter * sum(filter)
    
    return numpy.real(fft.ifft(fft.ifftshift(filtered)))
    
    
def frange(start, stop, step):
     r = start
     while not r == stop:
        yield r
        r += step

def peaksFromSignal(signal,sigmalScale,verbose=False):
    import numpy
    
    peaks = numpy.argwhere(signal[:] >= (numpy.mean(signal) + sigmalScale*numpy.std(signal)))
    return peaks
    
def meanLengthFromSignal(signal,verbose=False):
    
    sigmaScale = 2
    length = 0
    
    for sigmaScale in frange(3,0.5,-0.5):
        length = 0
        peaks = peaksFromSignal(signal,sigmaScale,verbose)
        
#         if verbose:
#             print '<!-- ######### meanLengthFromSignal ######### -->'
#             print sigmaScale
#             print peaks
#             print '<!-- ######### END ######### -->'
            
        if len(peaks) == 0:
            continue
        
        counter = 0
        for i,peak in enumerate(peaks[1:]):
            if peaks[i+1][0] - peaks[i][0] <= 18 or peaks[i+1][0] - peaks[i][0] >= 35:
                continue
            #print i,peaks[i+1][0] , peaks[i][0], peaks[i+1][0] - peaks[i][0]

            length += peaks[i+1][0] - peaks[i][0]
            counter += 1
            
        if counter == 0:
            if verbose:
                print '<!-- ######### meanLengthFromSignal ######### -->'
                print '<!-- ######### '+str(counter)+' ######### -->'
            
            continue
        
        length /= counter
        
        if verbose:
            print '<!-- ######### meanLengthFromSignal ######### -->'
            print '<!-- ######### '+str(length)+' ######### -->'
            
        
        if 18 <= length <= 35:
            break #return length if reasonable 
            
    return length,sigmaScale


def argmaxUnit(i,indexes,values,meanLength,maxIrregularity):

    maxV = None
    maxJ = None
    # print indexes
    # for j in xrange(i+meanLength-maxIrregularity,i+meanLength+maxIrregularity):
        
    j = i + 1
    while j < len(indexes) and indexes[j] - indexes[i] <= (meanLength+maxIrregularity):

        unitValue = float(values[indexes[i]] + values[indexes[j]])
        if (not maxV or maxV <= unitValue ) and (15 <= indexes[j] - indexes[i] <= (meanLength+maxIrregularity)):
            maxV = unitValue
            maxJ = int(j)
        #if maxJ:
            #print values[indexes[i]] , values[indexes[j]],values[indexes[maxJ]],maxV,i,j,indexes[j] - indexes[i]
        j += 1
    #print [maxJ,maxV]
    #print ''
    return [maxJ,maxV]

def getUnitsFromSignal(signal,queryOffset,maxIrregularity = 15,sigmaScale = None):

    queryOffset = 0
   #def optimizeUnits(indexes,values,meanLength,maxIrregularity):
    meanLength,newSigmaScale = meanLengthFromSignal(signal[queryOffset:])
    
    if sigmaScale == None:
        indexes = peaksFromSignal(signal,1.5)
    else:
        indexes = peaksFromSignal(signal,sigmaScale)
        
    units = []
    
    i = numpy.argmax(signal[indexes[0:3]])
    isStart = True
    while i < len(indexes):
        if indexes[i] < queryOffset:
            i += 1
            #print 'Skipping!'
            continue
        
        #print indexes[i]
        [j,maxV] = argmaxUnit(i,indexes,signal,meanLength,maxIrregularity)

        if j:
            units.append([indexes[i][0],indexes[j][0]])
            i = j
        else:
            i += 1
    
    if not sigmaScale and len(units) == 0:
        units = getUnitsFromSignal(signal,queryOffset,maxIrregularity,newSigmaScale)
    
    return units

def compressMatrix(matrix,referenceName):
    import numpy
    unitLength = unitLenghts[referenceName]

    compressedMatrix = numpy.zeros(shape=(unitLength,matrix.shape[1]))
    counter = 0

    unitStarts = unitStartDictionary[referenceName]

    for i,unitNumber in enumerate(unitStarts):
        try:
            compressedMatrix = compressedMatrix + matrix[unitStarts[i]:unitStarts[i+1],:]
            counter += 1
        except ValueError:
            continue
        except IndexError:
            continue

    compressedMatrix /= counter


    return compressedMatrix

def diagonalSum(matrix):
    
    sumArray = numpy.zeros(shape=(matrix.shape[1]))

    for i in xrange(matrix.shape[1]):
        sum = 0
        for j in xrange(matrix.shape[0]):
            try:
                sum += matrix[j,i+j]
            except IndexError:
                continue
                
        sumArray[i] = sum / matrix.shape[0]
            
            
    return sumArray

def parseFFASMatrix(ffasMatrixFileName):
    import numpy
    import re
    
    lines = None
    
    with open(ffasMatrixFileName,'r') as f:
        lines = f.readlines()
    
    query1Length = len(lines)
    l = re.sub( '\s+', ' ',lines[0][:-1]).strip()
    query2Length = len(l.split(' '))

    matrix = numpy.zeros(shape=(query1Length,query2Length))
    # print lines[0]
    # print len(lines[1][:-1])
    for i,line in enumerate(lines):

        l = line[1:-1]

        l = re.sub('\s+',' ',l).strip()
        
        for j,value in enumerate(l.split(' ')):
            if value == '':
                continue
            
            #value
            try:
                matrix[i,j] = float(value)
            except ValueError:
#                 print value.split(' ')[0]
                matrix[i,j] = float(value.split(' ')[0])

    return matrix

def filterDiagonal(matrix,kernelLength):

    from scipy import ndimage
    import numpy

    kernel = numpy.zeros(shape=(kernelLength,kernelLength))

    for i in xrange(kernelLength):
        kernel[i,i] = 1

    return ndimage.convolve(matrix,kernel)

def parseFFASAlignment(ffasAlignmentFile):

    lines = None

    with open(ffasAlignmentFile,'r') as f:
        lines = f.readlines()
    if len(lines) == 3 and '>*' in lines[2]:
        raise NonSolenoidException
    
    referenceAlignment = lines[3]
    
    queryAlignment = lines[4]
    
    try:
#         print referenceAlignment.split(' ')
        referenceOffset = int(referenceAlignment.split(' ')[4])
    except ValueError:
        referenceOffset = int(referenceAlignment.split(' ')[5])
      
    try:        
        queryOffset = int(queryAlignment.split(' ')[3])
    except ValueError:
        queryOffset = int(queryAlignment.split(' ')[4])
        
    return referenceOffset,queryOffset,referenceAlignment.split(' ')[4],queryAlignment.split(' ')[4]

def run(fileName,sequence=None,plot=False,verbose=False,referenceName='2OMZ'):
    from scipy.ndimage.filters import laplace
    import numpy
    upperValues = 5
    diagonalLength = 10
    filterBand = 25
    if verbose:
        print '<!--',fileName,'-->'
    try:
        m = parseFFASMatrix(fileName)
    except IOError:
        print 'Error reading file', fileName
        raise Exception
    
    if plot:
        pyplot.imshow(m,interpolation='nearest',cmap = cm.BrBG)
        pyplot.show()
    
#     pyplot.imshow(m,interpolation='nearest',cmap = cm.BrBG)
#     pyplot.show()
    
    m = filterDiagonal(m,diagonalLength)
    
    if plot:
        pyplot.imshow(m,interpolation='nearest',cmap = cm.BrBG)
        pyplot.show()
    m = compressMatrix(m,referenceName)
    if plot:
        pyplot.imshow(m,interpolation='nearest',cmap = cm.BrBG)
        pyplot.show()
        
    m = laplace(m)
    if plot:
        pyplot.imshow(m,interpolation='nearest',cmap = cm.BrBG)
        pyplot.show()
    
    signal = diagonalSum(m[:upperValues,:])
    signal = (signal - numpy.mean(signal)) / numpy.std(signal)
    
    irregularities = {}
    try:
        if verbose:
            print '<!--','Before parseFFASAlignment','-->' 
        referenceOffset,queryOffset,referenceAlignment,queryAlignment = parseFFASAlignment(fileName.replace('.mat','.txt'))
        if verbose:
            print '<!--','After parseFFASAlignment','-->'                
    except :
        queryOffset = 0
    
        
#     print signal
#     queryOffset = 0
    print '<!--','Query Offset', queryOffset,'-->'
    meanLength,sigmaScale = meanLengthFromSignal(signal[queryOffset:],verbose)
    print '<!--',sigmaScale,'-->'
    indexes = peaksFromSignal(signal,sigmaScale)

    units = getUnitsFromSignal(signal,queryOffset,12)

    unitStarts = [u[0] for u in units]
    unitEnds = [u[1] for u in units]
    
    irregularities = [unit[1] - unit[0] for unit in units]
    if sequence:
        subsequences = [sequence.seq[unit[0]:unit[1]] for unit in units] 
    else:
        subsequences = []
    return unitStarts, unitEnds, units,irregularities , subsequences,meanLength,signal

if __name__ == '__main__':
    
    import numpy
    import argparse
    from Bio import SeqIO
    
    parser = argparse.ArgumentParser(description='Extract aperiodicity from sequence.')
    parser.add_argument('-f', action='store',type=str,  help='The ffas matrix file')
    parser.add_argument('-q', action='store',type=str,  help='The query sequence in fasta format')
    parser.add_argument('-p', action='store_true',  help='Plot aperiodicity')
    parser.add_argument('--positionReference', action='store',type=str,  help='PDB ID of profile reference')
    parser.add_argument('--pdb', action='store',type=str,  help='PDB')
    parser.add_argument('--chain', action='store',type=str,  help='Chain') 
    parser.add_argument('--plotFolder', action='store',type=str,  help='Folder for plots')
         
    argumentDictionary = parser.parse_args().__dict__ 

    fileName           = argumentDictionary['f']
    fastaFile          = argumentDictionary['q']
    plot               = argumentDictionary['p']
    pdbID              = argumentDictionary['pdb']
    chainID            = argumentDictionary['chain']
    referenceName      = argumentDictionary['positionReference']
    plotFolder         = argumentDictionary['plotFolder']
    
    if plotFolder:
        import matplotlib
        matplotlib.use('Agg')
    
    from matplotlib import pyplot,cm
    
    try:
        sequence = list(SeqIO.parse(fastaFile, "fasta"))[0]
    except:
        sequence = None
        

    try:
        unitStarts,unitEnds,units,irregularities,subsequences,meanLength,signal = run(fileName,sequence,plot,True,referenceName)
        
        idName = fileName.split('/')[-1].split('.')[0][4:]
        
        #print fileName
        xmlString   = '<Sequences><Sequence ID="'+idName+'" Status="Solenoid">\n'
    #     xmlString  += '<Description>'
    #     xmlString  += sequence.description
    #     xmlString  += '</Description>'
#         print '<OriginalMeanLength Value="'+str(meanLength)+'"/>'
        meanLength = numpy.mean(numpy.array(irregularities)) 
#         print "<Irregularities>",irregularities,"</Irregularities>"
#         print "<MeanLength>",meanLength, "</MeanLength>"
        profileArea = numpy.sum(numpy.abs(numpy.array(irregularities) - meanLength)) / float(len(irregularities))
#         xmlString  += '<Area Value="'+str(profileArea)+'"/>\n'
        
        for i,unit in enumerate(units):
                irregularity = irregularities[i]
                xmlString += '<Unit Boundaries="'+str(unit)+'" UnitLength="'+str(irregularity)+'"/>\n' #'" Sequence="'+sequence.seq[unit[0]:unit[1]]+'"/>\n'
        xmlString +=    '<UnitMeanLength Value="'+str(meanLength)+'"/>\n'
        xmlString += '</Sequence></Sequences>\n'
        
    except NonSolenoidException:
        xmlString = '<Sequences><Sequence ID="'+fileName.split('.')[0][4:]+'" Status="NonSolenoid"/></Sequences>'
        
        
    print xmlString
#     print fileName , meanLength , irregularity
    
    if plot:
        pyplot.plot(signal)
        pyplot.scatter(unitStarts,signal[unitStarts],c = 'b',marker='o')
        pyplot.scatter(unitEnds,signal[unitEnds],c = 'b',marker='o')
        pyplot.show()
    
    if plotFolder:
        signalFileName = plotFolder + '/' + idName.split('_')[0] + '_signal.png' 
        pyplot.plot(signal)
        pyplot.scatter(unitStarts,signal[unitStarts],c = 'b',marker='o')
        pyplot.scatter(unitEnds,signal[unitEnds],c = 'b',marker='o')
        pyplot.savefig(signalFileName)
        pyplot.clf()
        
    if plot:
        pyplot.plot(irregularities)
        pyplot.plot([meanLength] * len(irregularities))
        pyplot.plot([abs(i - meanLength) for i in irregularities])
        pyplot.plot([0] * len(irregularities))
        pyplot.show()
        
    if plotFolder:
        signalFileName = plotFolder + '/' + idName.split('_')[0] + '_irregular.png' 
        pyplot.plot(irregularities)
        pyplot.plot([meanLength] * len(irregularities))
        pyplot.plot([abs(i - meanLength) for i in irregularities])
        pyplot.plot([0] * len(irregularities))
        pyplot.savefig(signalFileName)
        
        
    
    
    
    