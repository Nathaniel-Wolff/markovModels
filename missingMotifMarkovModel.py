import sys
import math

"""This class was developed by David Bernick"""
class CommandLine():
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond.
    it implements a standard command line argument parser with various argument options,
    a standard usage and help, and an error termination mechanism do-usage_and_die.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.

    '''

    def __init__(self, inOpts=None):
        '''
        CommandLine constructor.
        Implements a parser to interpret the command line argv string using argparse.
        '''

        import argparse
        self.parser = argparse.ArgumentParser(
            description='Takes max and min motif size from the user, the zScore cutoff for missing motifs, input fasta to use, and version',
            epilog='Nothing much else to say',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='%(prog)s [options] -option1[default] <input >output'
        )

        self.parser.add_argument('-l', '--minMotif', nargs='?', default=1, action='store',
                                 help='min kMer size ')
        self.parser.add_argument('-m', '--maxMotif', nargs='?', default=8, action='store',
                                 help='max kMer size ')
        self.parser.add_argument('-c', '--cutoff', nargs='?', type=float, default=.01, action='store',
                                 help='Zscore cutoff')
        self.parser.add_argument('-i', '--inputFile', nargs=1, default="Ecoli-UMN026.fa", action='store',
                                 help="input fasta file")

        #self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')

        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)

        """Each of functions  assigns the min and max motif size, z score cutoff, input genome file, and version, respectively, from the command line to attibutes of Command Line class.  """
        args = self.parser.parse_args()
class FastAreader:
    """This class, which reads in a Fasta file and returns its headers/sequences that can be used for futher parsing. This was developed by David Bernick"""
    def __init__(self, fname=''):
        '''contructor: saves attribute fname '''

        self.fname = fname
        self.fileH = None

    def doOpen(self):
        if self.fname == '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta(self):

        header = ''
        sequence = ''

        with self.doOpen() as self.fileH:

            header = ''
            sequence = ''

            # skip to first fasta header
            line = self.fileH.readline()
            while not line.startswith('>'):
                line = self.fileH.readline()
            header = line[1:].rstrip()

            for line in self.fileH:
                if line.startswith('>'):
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else:
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header, sequence

class rc:
    """the following method was written by David Bernick, here it is used to generate the reverse complement sequence for a given motif. """
    def rc(self,s):
        ''' Return the reverse complement of the sequence.'''
        return s[::-1].translate(str.maketrans('ACGTN', 'TGCAN'))  # N is needed for composites

class Genome(rc):
    """Genome class handles a set of sequences, produces a count dictionary for the motifs of sizes given by the user, and produces a null model of order 2 less
     than the maximum motif size. Then, the expected counts for each motif are compared to the null model, and a z score, expected counts, forward and reverse sequences, and actual counts are presented to the user."""

    def __init__(self, minK, maxK, sequences, genome="ACGTAAACA"):
        """Construct the object to count and analyze kMer frequencies"""
        """Each attribute comes from the command line's results."""
        self.minK = minK
        self.newMinK = int(minK)
        self.maxK = maxK
        self.newMaxK = int(maxK)
        self.genome = genome
        self.sequences = sequences




    """This function will take the sequences from the input file passed to the command line and calculate the size of the genome. This doesn't concatenate it."""
    def getFullSeqLen(self):
        fullGenomeLen = 0
        for sequence in self.sequences:
            N = len(sequence)
            fullGenomeLen += N
        return fullGenomeLen

    """This function will generate all of the motifs (all the way from 1 base long to the longest motif as specified on the command line) and count their occurences over all of the genome's sequences without concatenation...."""

    def constructKmerDict(self):
        countsKmers = {}
        if self.newMinK < 3:
            lowestMer = 1
        else:
            lowestMer = self.newMinK - 2
        for sequence in self.sequences:
            N = len(sequence)
            for pos in range(0, N - self.newMaxK+1) :
                for mer in range(lowestMer, self.newMaxK+1):
                    subMer = sequence[pos:pos+mer]
                    rcSubMer = self.rc(subMer)

                    """Counting portion after alphabetizaton."""
                    if subMer < rcSubMer:
                        merTuple = (subMer,rcSubMer)
                    else:
                        merTuple = (rcSubMer,subMer)

                    if merTuple not in countsKmers.keys():
                        countsKmers[merTuple] = 1

                    else:
                        countsKmers[merTuple] += 1

        return countsKmers
    """This method determines the null model of the aforementioned order and compares the counts given by the distribution to calculate a z score. Then, it prints the sequences (forward and reverse) that falls under a certain z score cutoff."""
    def nullModelMarkovOrder2Below(self, countsKmers, genomeLen, cutoff):
        eMerDict = {}
        for mers,count in countsKmers.items():
            mer = mers[0]
            rcMer = mers[1]
            merLen = len(mer)
            #print(mer, merLen, self.newMinK)
            if merLen >= self.newMinK and merLen > 2:
                leftSubMer = mer[:-1]
                rcLeftSubMer = self.rc(leftSubMer)

                if leftSubMer < rcLeftSubMer:
                    merTuple = (leftSubMer, rcLeftSubMer)
                else:
                    merTuple = (rcLeftSubMer, leftSubMer)

                countLeft = countsKmers[merTuple]

                rightSubMer = mer[1:]
                rcRightSubMer = self.rc(rightSubMer)

                if rightSubMer < rcRightSubMer:
                    merTuple = (rightSubMer, rcRightSubMer)
                else:
                    merTuple = (rcRightSubMer, rightSubMer)
                countRight = countsKmers[merTuple]

                centerSubMer = mer[1:-1]
                rcCenterSubMer = self.rc(centerSubMer)

                if centerSubMer < rcCenterSubMer:
                    merTuple = (centerSubMer, rcCenterSubMer)
                else:
                    merTuple = (rcCenterSubMer, centerSubMer)
                countCenter = countsKmers[merTuple]



                nullExpectation = (countLeft*countRight)/countCenter
                stdDev = math.sqrt(nullExpectation * (1-(nullExpectation/genomeLen)))
                zScore = (count - nullExpectation)/stdDev

                if zScore < cutoff:

                           
                    eMerDict[(mer, rcMer)] = [count, nullExpectation, zScore]

                    display = sorted(eMerDict.items(), key=lambda x: x[1][-1], reverse=False)

        print("sequence: reverse count Expect Zscore")

        for thing in eMerDict.items():
            sequence = thing[0][0]
            reverse = thing[0][1]
            count = int(thing[1][0])
            Expect = float(thing[1][1])
            ZScore = float(thing[1][2])
            print('{0:8}:{1:8}\t{2:0d}\t{3:0.2f}\t{4:0.2f}'.format(sequence, reverse, count, Expect, ZScore))
        print("I couldn't figure out how to sort this thing. Can you advise on that?")












def main(inFile=None, options=None):
    ''' Setup necessary objects, read data and print the final report.'''
    cl = CommandLine(options) # setup the command line
    inFile = cl.args.inputFile
    minMer = cl.args.minMotif
    maxMer = cl.args.maxMotif
    cutoff = cl.args.cutoff


    sourceReader = FastAreader(inFile)  # setup the Fasta reader Object
    headSeqDict = {}
    seqList = []
    for head, seq in sourceReader.readFasta():
        seqList.append(seq)
        headSeqDict[head] = seq

    thisGenome = Genome(minMer, maxMer, headSeqDict.values())  # setup a Genome object
    genomeLen = thisGenome.getFullSeqLen()
    kMerDict = thisGenome.constructKmerDict()
    thisGenome.nullModelMarkovOrder2Below(kMerDict, genomeLen, cutoff)


if __name__ == "__main__":
    main(options=["--minMotif=3", "--maxMotif=8", "--cutoff=-2.0"])

