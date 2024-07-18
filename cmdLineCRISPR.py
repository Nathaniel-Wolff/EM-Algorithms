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