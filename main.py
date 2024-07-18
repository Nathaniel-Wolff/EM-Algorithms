import math
import random
import sys

class cmdLineCRISPR():
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
            description="""Takes motif size (-k), # of iterations of the profile scoring loop (-i), psuedocounts (-p), input Fasta to use (<>), and output file, in which the final entropy score and consensus sequence is written.""",
            epilog='',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='%(prog)s [options] -option1[default] <input > output'
        )

        self.parser.add_argument('-k', '--motifSize', nargs='?', default=13, action='store',
                                 help='motif size')
        self.parser.add_argument('-i', '--iterations', nargs='?', default=10000, action='store',
                                 help='iterations ')
        self.parser.add_argument('-p', '--pseudocounts', nargs='?', type=float, default=1, action='store',
                                 help='psuedocounts')

        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)

        args = self.parser.parse_args()
class FastAreader:
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

class emModelMotifSearch():
    """Instantiates instances of the EMModelMotif with a list of n sequences from a Fasta file, a kmer size (k), and a number of psuedocounts (p)."""

    def __init__(self, sequences,k,p):
        """Takes the user's motif length (kmer size) and number of pseudocounts (p).
         1) Saves the base composition of the sequence
         2) Saves the length of each sequence to the length of the genome.

        Then,
        It counts all the occurences of the motifs, while adding the user specified number of pseudocounts.
        Meanwhile, it generates a random motif from each sequence of size k, which is saved to the instance. """

        self.sequences = sequences
        #base composition of the sequence for the null model
        self.baseComposition = {}
        #n, number of sequences
        self.numberSequences = 0
        #kmer size
        self.k = k
        # number of pseudocounts
        self.p = p
        self.lenGenome = 0
        self.lenSeq = 0
        self.randomMotifs = []
        self.currentMotifs = []
        for sequence in self.sequences:
            sequence = sequence.rstrip()
            #couting the number of sequences in the FASTA
            self.numberSequences += 1
            #finding the total length of the genome
            lengthSeq = len(sequence)
            self.lenSeq = lengthSeq
            self.lenGenome += lengthSeq
            #determining the base composition of the genome and building the null model
            for pos in range(0, lengthSeq):
                    base = sequence[pos:pos+1]
                    if base in self.baseComposition.keys():
                        self.baseComposition[base] += 1
                    else:
                        self.baseComposition[base] = self.p
            #making null model
            self.nullModel = {}
            for base in self.baseComposition:
                self.nullModel[base] = self.baseComposition[base]/self.lenGenome

    def makeRandomKmers(self):
        "Generates a set of random kmers, one from each of the object's sequences."
        randomMotifs = []
        for sequence in self.sequences:
            lenSeq = len(sequence)
            randomMotifStart = random.randint(0, lenSeq-self.k+1)
            randomMotif = sequence[randomMotifStart: randomMotifStart + self.k-1]
            randomMotifs.append(randomMotif)
        return randomMotifs

    def buildProfile(self, motifs):
        """Creates a profile from a set of motifs from the user."""
        """Via the following:
        1) "Aligns" the sequences by iterating through each index
        2) Iterates through all of the sequences and for each, counts the base and updates the column's base composition dictionary, forming the counts matrix.
        3) It then goes through the counts matrix and forms a profile for each row. 
        """

        #the counts matrix consists of a list of dictionaries, in order of sequence index (that is, column)
        countsMatrix = []
        #the profile consists of a list of dictionaries, in the same order as the counts matrix (that is, column)
        profile = []

        #checking to see if the motifs are all the same length, which is a prerequisite for building a functional profile
        motifLength = 0
        for motif in motifs:
            for innerMotif in motifs:
                if len(motif) != len(innerMotif):
                    print("Can't build a functional profile.")
                    return None
                else:
                    motifLength = len(motif)
        #if it passes, the loop goes on
        for sequenceIndex in range(0, motifLength):
            columnDict = {}
            for sequence in motifs:
                base = sequence[sequenceIndex]
                #updating the column dict counts appropriately
                if base not in columnDict.keys():
                    columnDict[base] = self.p
                else:
                    columnDict[base] += 1
            #checking to see if there are any bases that aren't present in this row, and adding blanks to the dictionary
            absentBases = [base for base in ["A", "C", "G", "T"] if base not in columnDict.keys()]
            for absentBase in absentBases:
                #adding pseudocounts to them
                columnDict[absentBase] = self.p
            #alphabetizing the column dict
            sortedColumnDict = dict(sorted(columnDict.items()))
            #appending to counts matrix
            countsMatrix.append(sortedColumnDict)

        for columnDict in countsMatrix:
            profileColumnDict = {}
            totalCounts =  sum(columnDict.values())
            #calculting the prob of each base and saving to the profile column dict
            for base in columnDict:
                probBase = columnDict[base] / totalCounts
                profileColumnDict[base] = probBase
            profile.append(profileColumnDict)
        return profile

    def scoreProfile(self, profile):
        """Creates a null model from the object's base composition dictionary and uses it to score a profile.
        It reads through the columns and computes the point relative entropy for each base, and adds it to the entropy total."""
        relativeEntropy = 0
        for column in profile:
            for base in column.keys():
                nullModel = self.nullModel[base]
                alternativeModel = column[base]
                pointRelativeEntropy = alternativeModel * math.log2(alternativeModel/nullModel)
                relativeEntropy += pointRelativeEntropy
        reScore = relativeEntropy
        return reScore

    def makeNewMotifs (self, profile):
        """This function uses a profile to generate its best scoring motifs.
        It does this by:
        1) Iterating through each of the data DNA sequences
            Building each k-mer for each
            Looking up the probability of obtaining each kmer base in the profile.
            Multiplying them together to get the joing probability of getting said kmer.

            It asks if that probability is higher than the last one built.
            If so, it saves that kmer and probability until the next iteration.
         """
        newMotifs = []
        for sequence in self.sequences:
            currentBestMotif = ""
            currentBestMotifJP = 0
            for baseIndex in range(0, len(sequence)-self.k+1):
                #generating each kmer by iterating through the sequence
                possibleMotif = sequence[baseIndex:baseIndex+self.k]
                jointProbability = 1

                #iterating through each base of the kmer and getting its probability
                for baseIndex in range(0, len(possibleMotif)-1):
                    #finding probability of obtaining each base in kmer at its position
                    base = possibleMotif[baseIndex]

                    column = profile[baseIndex]
                    probBase = column[base]
                    #multiplying the joing probability by said prob
                    jointProbability*=probBase

                #checking to see if the kmer's joint prob is better than the current best
                #if it is, save it and its probability to the bestMotif variables
                if jointProbability > currentBestMotifJP:
                    currentBestMotifJP = jointProbability
                    currentBestMotif = possibleMotif

            #saving the last, bestMotif to the new motifs, and going to the next sequence
            newMotifs.append(currentBestMotif)
        return newMotifs

    def makeConsensus(self, profile):
        """This function builds a consensus sequence from best scoring profile.
        It is called only after the user specified number of iterations of Expectation-Maximization is over (non Gibbs-sampling).
        """
        consensus = ""
        for column in profile:
            maxProbBase = max(column.items(), key=lambda x: x[1])[0]
            maxProb = max(column.values())
            consensus += maxProbBase
        return consensus

def main(inFile=None, options=None):
    ''' Setup necessary objects, read data and print the final report.'''
    #cl contains these attributes: motif size, iterations, psuedocounts, input and output files: use cl.args.attribute to obtain one.
    cl = cmdLineCRISPR(options) # setup the command line
    motifSize = int(cl.args.motifSize)
    pseudocounts = int(cl.args.pseudocounts)
    iterations = int(cl.args.iterations)

    #Fastareader takes standard in.
    thisFastaReader = FastAreader()
    #gets standard in from the Fastareader
    inFile = thisFastaReader.doOpen()

    #Fastareader instance reads its Fasta and saves it internally
    seqList = []
    for head, seq in thisFastaReader.readFasta():
        seqList.append(seq)

    thisEM = emModelMotifSearch(seqList, motifSize, pseudocounts)

    bestScore = 0
    bestMotifs = []
    bestProfile = {}

    for iterationAmount in range(0, iterations):
        #iteration generates a set of random kmers
        currentKmers = thisEM.makeRandomKmers()
        #makes a initial profile from the random kmers and scores it
        currentProfile = thisEM.buildProfile(currentKmers)
        currentScore = thisEM.scoreProfile(currentProfile)

        #optimization loop
        while True:
            #makes motifs from the current profile. The first run of this loop makes motifs from the random kmer's profile
            newMotifs = thisEM.makeNewMotifs(currentProfile)
            #makes a new profile from the motifs that came from the current profile and scores it
            newProfile = thisEM.buildProfile(newMotifs)
            newScore = thisEM.scoreProfile(newProfile)

            #checks to see if the score of the profile from the motifs is better
            #if this is the first run of the loop, it compares the score of the random kmer profile
            #otherwise, it compares it to the last profile
            if newScore > currentScore:
                #if it is better than the score of the last profile
                #it replaces that score
                #the current profile and motifs are saved internally
                currentScore = newScore
                currentProfile = newProfile
                currentMotifs = newMotifs
            else:
                #current score fails to improve and the while loop breaks, ending this iteration
                break
        #after the iteration ends
        #if the iteration's last score is better than the best score over all iterations
        #this score replaces said score
        #profile and motifs are saved to the best of them over all iterations
        if currentScore > bestScore:
            bestScore = currentScore
            bestProfile = currentProfile
    #consensus sequence is made with the best scoring profile over all iterations at the end
    consensus = thisEM.makeConsensus(bestProfile)
    #printing it to the console
    print("Consensus sequence:", consensus)
    print("Best Score:", bestScore)

if __name__ == "__main__":
    main()