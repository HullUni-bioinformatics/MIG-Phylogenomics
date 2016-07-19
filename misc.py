### misc ###
def get_two_senses_strs(record):
    from Bio import SeqIO
    """
    record: a SeqRecord object
    returns a list with two strings: the sequence and the revcomp sequence
    """
    return [str(record.seq), str(record.seq.reverse_complement())]

def sed_subset_fastq_gz(full_filename, subset, subset_filename):
    
    values = [full_filename, subset*4, subset_filename]
    
    cline = "zcat {0[0]} | sed -n 1,{0[1]}p | gzip > {0[2]}".format(values)
    
    print "A subset of %i sequences"%(values[1]/4)
    print "from the file %s"%values[0]
    print "will be written to %s"%values[2]
    print "using the command"
    print cline
    
    err, out = execute_cline(cline)
    
    print "stdout:\n%s"%out
    print "stderr:\n%s"%err
    
def parse_nhmmer_tblout(filname):
    lines = [l for l in open(filname,'r').readlines() if not l.startswith('#')]
    matches = []
    for l in lines:
        target_name,accession,query_name,accession,hmmfrom,\
        hmm_to,alifrom,ali_to,envfrom,env_to,sq_len,strand,\
        E_value,score,bias,description_of_target=l.rstrip().split()
        matches.append({
                'target_name':target_name,
                'accession':accession,
                'query_name':query_name,
                'accession':accession,
                'hmmfrom':hmmfrom,
                'hmm_to':hmm_to,
                'alifrom':alifrom,
                'ali_to':ali_to,
                'envfrom':envfrom,
                'env_to':env_to,
                'sq_len':sq_len,
                'strand':strand,
                'E_value':E_value,
                'score':score,
                'bias':bias,
                'description_of_target':description_of_target
                
            })
    return sorted(matches, key=lambda m: float(score), reverse=True)        

def is_overlapping(coords1, coords2, max_overlapp=0):
    
    """
    Check if two ranges are overlapping. The ranges are each
    a list or tupple including of two values, the begining and end.
    The begining and end are included in the range.
    
    # overlapping ranges
    
    >>> coords1 = (10, 100)
    >>> coords2 = (50, 150)
    >>> is_overlapping(coords1, coords2)
    True
    
    >>> is_overlapping(coords2, coords1)
    True
    
    # Second range is revcomp and input is
    # list instead of tuple
    
    >>> coords1 = [10, 100]
    >>> coords2 = [150, 50]
    >>> is_overlapping(coords1, coords2)
    True
    
    # One range is completely nested in the other
    
    >>> coords1 = (10, 100)
    >>> coords2 = (20, 90)
    >>> is_overlapping(coords1, coords2)
    True
    
    >>> is_overlapping(coords2, coords1)
    True
    
    # The ranges overlap by a single position
    
    >>> coords1 = (10, 100)
    >>> coords2 = (100, 200)
    >>> is_overlapping(coords1, coords2)
    True
    
    >>> is_overlapping(coords2, coords1)
    True
    
    # The ranges do not overlap
    
    >>> coords1 = (10, 100)
    >>> coords2 = (200, 300)
    >>> is_overlapping(coords1, coords2)
    False
    """
    
    # inputs need to be lists or tuples
    assert all([isinstance(coords1,(list,tuple)),
                isinstance(coords2,(list,tuple))])
    
    for lst in [coords1, coords2]:
        
        # inputs need to have two values each
        assert len(lst) == 2
        
        # inputs values have to be all int
        for i in lst:
            assert isinstance(i, int)
    
    if any([coords2[0] <= coords1[0] <= coords2[1]-max_overlapp,
            coords2[0]+max_overlapp <= coords1[1] <= coords2[1],
            coords1[0] <= coords2[0] <= coords1[1]-max_overlapp,
            coords1[0]+max_overlapp <= coords2[1] <= coords1[1]]):
        return True
    else:
        return False
    
class Indel:
    
    def __init__(self, pos, length, indel_type, aln):
        
        """
        Characterizes indels in a top sequence alignment
        compared to a bottom sequence alignment
        pos: the position in the top sequence where the indel starts
        length: the number of '-' characters following pos, either in
        the top or in the bottom sequence
        indel_type: insertion: the '-' characters are in the bottom 
        sequence. deletion: the '-' characters are in the top sequence.
        aln: the alignment object.
        """
        
        self.position = pos
        self.length = length
        self.type = indel_type
        self.ref_alignment = aln
        
    def __str__(self):
        return("start: %i,  length: %i, type: %s, sequence: %s, reference: %s"%(
                self.position,
                self.length,
                self.type,
                self.ref_alignment[0].id,
                self.ref_alignment[1].id
               ))

def is_indel(aln, i, indel_type):
    k = None
    if indel_type == 'insertion':
        k = 1
    elif indel_type == 'deletion':
        k = 0
    if aln[k,i] == '-' and (i == 0 or not aln[k,i-1] == '-'):
        # find pos in top sequence
        sub_top_seq_length=len(aln[0,:i])-str(aln[0,:i]).count('-')
        # find insertion length (number of '-' in bottom sequence)
        # or
         # find deletion length (number of '-' in top sequence)
        length = 0
        j = i
        p = aln[k,j]
        while p == '-':
            length += 1
            j += 1
            p = aln[k,j]
        return Indel(sub_top_seq_length-1, length, indel_type, aln)
    else:
        return False
        
def is_insertion(aln, i):
    return is_indel(aln, i, 'insertion')

def is_deletion(aln, i):
    return is_indel(aln, i, 'deletion')

    
def find_indels(aln):
    indels = []
    for i in range(aln.get_alignment_length()):
        indel = is_insertion(aln, i)
        if indel:
            indels.append(indel)
        else:
            indel = is_deletion(aln, i)
            if indel:
                indels.append(indel)
            
    return indels


def __read_kwargs__(suffix='.py', allowed=None, exe_path='',
                    underscore_in_exe='-',
                    underscore_in_keyword='-',
                    **kwargs):
    """
    Given a set of keyword arguments, format as commad line string
    
    # Make up a function with two allowed kwargs: a and bb
    # It will write the function as a command line
    
    >>> def madeup(**kwargs):
    ...     allowed = ['a','bb']
    ...     cline = __read_kwargs__(allowed=allowed, **kwargs)
    ...     return cline
    
    >>> print(madeup(a='a', bb='bb'))
    madeup.py -a a --bb bb
    
    # Now try to use a keyword that is not allowed
    >>> print(madeup(a='a', cc='cc'))
    Traceback (most recent call last):
     ...
    IOError: Keyword cc is not allowed in madeup
    """
    import inspect
    curframe = inspect.currentframe()
    calframe = inspect.getouterframes(curframe, 2)
    caller = calframe[1][3]
    
    cline = ''
    
    if exe_path != '' and not exe_path.endwith('/'):
        exe_path +- '/'
    
    u_in_e = underscore_in_exe
    u_in_k = underscore_in_keyword
    
    for keyword in kwargs:
        if allowed and not keyword in allowed:
            raise IOError("Keyword %s is not allowed in %s"%(keyword, caller))
        prefix = '--'
        if len(keyword) == 1:
            prefix = '-'
        if str(kwargs[keyword]) == 'True':
            cline += "%s%s "%(prefix,
                                 keyword.replace('_',u_in_k))
        else:
            cline += "%s%s %s "%(prefix,
                                 keyword.replace('_',u_in_k),
                                 str(kwargs[keyword]))
    return "%s%s%s %s"%(exe_path,caller.replace('_',u_in_e),suffix,cline[:-1])

def __sufix_kwargs__(cline, arg):
    """
    Add arguments after the keyword arguments
    """
    if type(arg) in [tuple,set,list]:
        for a in arg:
            cline += " %s"%a
    else:
        cline += " %s"%arg
    return cline
    
    
def __prefix_kwargs__(cline, arg):
    """
    Add an arguments between the executable and the next argument/ kwarg
    """
    if type(arg) in [tuple,set,list]:
        raise IOError('Cannot prefix tuple, set ot list to cline')
    cline = cline.partition(' ')
    cline = "%s %s %s"%(cline[0],arg,cline[-1])
    return cline

def execute_cline(cline):
    from subprocess import Popen,PIPE
    p = Popen(cline, shell=True, stdout=PIPE, stderr=PIPE)
    out, err = p.communicate()
    return out, err

def makedir(name, f=False):
    import os
    import warnings
    if os.path.exists(name):
        if f:
            import shutil
            shutil.rmtree(name)
            os.mkdir(name)
        else:
            warnings.warn('keeping existing fpath %s'%name)
    else:
        os.mkdir(name)
        warnings.warn('fpath %s newly created'%name)
        
def extractgz(source, destination):
    
    cline = 'zcat %s > %s'%(source, destination)
    out, err = execute_cline(cline)
    return out, err

def compressgz(source):
    
    cline = 'gzip %s'%source
    out, err = execute_cline(cline)
    return out, err

def pickledump(obj, target):
    import cPickle as pickle
    picout = open(target, 'wb')
    pickle.dump(obj, picout, -1)
    picout.close()

def pickleload(source):
    import cPickle as pickle
    pkl_file = open(source, 'rb')
    obj = pickle.load(pkl_file)
    pkl_file.close()
    return obj

def printoe((out,err)):
    if len(out) > 0:
        print 'STDOUT:'
        print out
    if len(err) > 0:
        print 'STDERR:'
        print err
        
def extract_picload(source, destinationdir='./'):
    
    import os
    
    if not destinationdir.endswith('/'):
        destinationdir += '/'
    pklfile = destinationdir+source.split('/')[-1].replace('.gz','')
    printoe(
        extractgz(source, pklfile)
    )
    
    obj = pickleload(pklfile)
    
    os.remove(pklfile)
    
    return obj

def picdump_compress(obj, destination):
    import gc
    
    pickledump(obj, destination)
    
    printoe(
        compressgz(destination)
    )
    
    del obj
    gc.collect()

# Jellyfish
def jellyfish_count(input_sequence_filenames = None,
                    exe_path='',
                    **kwargs):
    """
    Count k-mers in fasta or fastq files

    Options (default value in (), *required):
     m, mer_len=int                     *Length of mer
     s, size=string                     *Initial hash size (eg 100M)
     t, threads=int                      Number of threads (1)
     F, Files=int                        Number files open simultaneously (1)
     g, generator=path                   File of commands generating fast[aq]
     G, Generators=int                   Number of generators run simultaneously (1)
     S, shell=string                     Shell used to run generator commands ($SHELL or /bin/sh)
     o, output=string                    Output file (mer_counts.jf)
     c, counter_len=Length in bits       Length bits of counting field (7)
        out_counter_len=Length in bytes  Length in bytes of counter field in output (4)
     C, canonical                        Count both strand, canonical representation (false)
        bc=peath                         Bloom counter to filter out singleton mers
        bf_size=string                   Use bloom filter to count high-frequency mers (eg 100M)
        bf_fp=double                     False positive rate of bloom filter (0.01)
        if=path                          Count only k-mers in this files
     Q, min_qual_char=string             Any base with quality below this character is changed to N
     p, reprobes=int                     Maximum number of reprobes (126)
        text                             Dump in text format (false)
        disk                             Disk operation. Do not do size doubling (false)
     L, lower_count=string               Don't output k-mer with count < lower-count (eg 100M)
     U, upper_count=string               Don't output k-mer with count > upper-count (eg 100M)
        timing=Timing file               Print timing information
        usage                            Usage
     h, help                             This message (but for command line, not function)
        full_help                        Detailed help
     V, version                          Version
     
    # Test:
    >>> jellyfish_count(m=21, s='100M', t=10, C=True)
    Traceback (most recent call last):
     ...
    IOError: jellyfish_count: jellyfish count -C -m 21 -s 100M -t 10 : Must provide input_sequence_filenames as first arguments

    """
    allowed = ['m', 'mer_len', 's', 'size', 't', 'threads', 'F', 'Files',
              'g', 'generator', 'G', 'Generators', 'S', 'shell', 'o', 'output', 
              'out_counter_len','C', 'canonical','bc', 'bf_size', 'bf_fp', 'if',
              'Q', 'min_qual_char', 'p', 'reprobes', 'text','disk','L', 'lower_count',
              'U', 'upper_count', 'timing', 'usage', 'h', 'help', 'full_help','V', 'version'] 
    cline = __read_kwargs__(suffix='',
                            exe_path=exe_path,
                            allowed=allowed,
                            underscore_in_exe=' ',
                            underscore_in_keyword='-',
                            **kwargs)
    if input_sequence_filenames:
        cline = __sufix_kwargs__(cline, input_sequence_filenames)
        out, err = execute_cline(cline)
        return out, err
    else:
        raise IOError("jellyfish_count: %s : Must provide input_sequence_filenames "%cline+
                      "as first arguments")

def jellyfish_histo(kmer_hash_bin_filenames = None,
                    exe_path='',
                    **kwargs):
    """
    Create an histogram of k-mer occurrences

    Create an histogram with the number of k-mers having a given
    count. In bucket 'i' are tallied the k-mers which have a count 'c'
    satisfying 'low+i*inc <= c < low+(i+1)*inc'. Buckets in the output are
    labeled by the low end point (low+i*inc).

    The last bucket in the output behaves as a catchall: it tallies all
    k-mers with a count greater or equal to the low end point of this
    bucket.

    Options (default value in (), *required):
     l, low=int                         Low count value of histogram (1)
     h, high=int                        High count value of histogram (10000)
     i, increment=int                   Increment value for buckets (1)
     t, threads=int                     Number of threads (1)
     f, full                            Full histo. Don't skip count 0. (false)
     o, output=string                   Output file
     v, verbose                         Output information (false)
     U, usage                           Usage
        help                            This message
        full_help                       Detailed help
     V, version                         Version
     
    # Test:
    >>> jellyfish_histo(full=True)
    Traceback (most recent call last):
     ...
    IOError: jellyfish_histo: jellyfish histo --full : Must provide kmer_hash_bin_filenames as first arguments

    """
    allowed = ['l', 'low', 'h', 'high', 'i', 'increment', 't',
               'threads', 'f', 'full', 'o', 'output',
               'v', 'verbose', 'U', 'usage', 'help',
               'full_help', 'V', 'version'] 
    
    cline = __read_kwargs__(suffix='',
                            exe_path=exe_path,
                            allowed=allowed,
                            underscore_in_exe=' ',
                            underscore_in_keyword='-',
                            **kwargs)
    if kmer_hash_bin_filenames:
        cline = __sufix_kwargs__(cline, kmer_hash_bin_filenames)
        out, err = execute_cline(cline)
        return out, err
    else:
        raise IOError("jellyfish_histo: %s : Must provide kmer_hash_bin_filenames "%cline+
                      "as first arguments") 
        

# Khmer
def load_into_counting(output_countgraph_filename = None,
                       input_sequence_filenames = None,
                       **kwargs):
    """
    Wrapper for load_into_counting.py
    
    Reuqired arguments:
    output_countgraph_filename - filepath for the output
    input_sequence_filenames   - a list of input fast[q/a] files
    
    optional arguments:
    h, help             show this help message and exit
    version             show program's version number and exit
    ksize,  k
                        k-mer size to use (default: 32)
    n_tables, N
                        number of tables to use in k-mer countgraph (default:
                        4)
    U, unique_kmers
                        approximate number of unique kmers in the input set
                        (default: 0)
    fp_rate             Override the automatic FP rate setting for the current
                        script (default: None)
    max_tablesize, x
                        upper bound on tablesize to use; overrides --max-
                        memory-usage/-M. (default: 1000000.0)
    M, max_memory_usage
                        maximum amount of memory to use for data structure.
                        (default: None)
    threads, T
                        Number of simultaneous threads to execute (default: 1)
    b, no_bigcount      The default behaviour is to count past 255 using
                        bigcount. This flag turns bigcount off, limiting
                        counts to 255. (default: True)
    summary_info, s
                        What format should the machine readable run summary be
                        in? (`json` or `tsv`, disabled by default) (default:
                        None)
    f, force            Overwrite output file if it exists (default: False)

    Note: with `b`/`no_bigcount` the output will be the exact size of the k-mer
    countgraph and this script will use a constant amount of memory. In exchange
    k-mer counts will stop at 255. The memory usage of this script with `b` will
    be about 1.15x the product of the `x` and `N` numbers.

    Example:

        load_into_counting('out', 'data/100k-filtered.fa', k=20 x=5e7)

    Multiple threads can be used to accelerate the process, if you have extra cores
    to spare.

    Example:

        load_into_counting('out', 'data/100k-filtered.fa', k=20, x=5e7, T=4)
   
    # Test:
    >>> load_into_counting(k=20, x=5e7, T=4)
    Traceback (most recent call last):
     ...
    IOError: load_into_counting: load-into-counting.py -x 50000000.0 -k 20 -T 4 : Must provide output_countgraph_filename and input_sequence_filenames as first arguments
    """
    
    allowed = ['ksize',  'k', 'h', 'help', 'version', 'n_tables', 'N',
               'U', 'unique_kmers', 'fp_rate', 'max_tablesize', 'x',
               'M', 'max_memory_usage', 'threads', 'T', 'b', 'no_bigcount',
               'summary_info', 's', 'f', 'force']
    
    cline = __read_kwargs__(suffix='.py', allowed=allowed, **kwargs)
    
    if output_countgraph_filename and input_sequence_filenames:
        if not isinstance(input_sequence_filenames,list):
            raise IOError("input_sequence_filenames must be list")
        cline = __sufix_kwargs__(cline, output_countgraph_filename)
        cline = __sufix_kwargs__(cline, input_sequence_filenames)
        out, err = execute_cline(cline)
        return out, err
        
    else:
        raise IOError("load_into_counting: %s : Must provide output_countgraph_"%cline+
                      "filename and input_sequence_filenames "+
                      "as first arguments")
        
def iupac(pos):
    
    from collections import Counter
    
    if len(pos) > 4:
        raise IOError('up to four pos allowed for IUPAC')
        
    if any([v > 1 for v in Counter(pos).itervalues()]):
        raise IOError('pos cannot repeat in IUPAC')
        
    allowed = 'atgcATGC'
    
    if any([p not in allowed for p in pos]):
        raise IOError('only '+str(allowed)+ ' allowed in IUPAC')
        
    pos = ''.join(sorted(''.join(pos))).upper()
        
    
    IUPAC = { 'AG':'R',
              'GT':'K',
              'AC':'M',
              'AT':'W',
              'CT':'Y',
              'CG':'S',
              'CGT':'B',
              'ACG':'V',
              'AGT':'D',
              'ACT':'H',
              'ACGT':'N',
              'A': 'A',
              'T': 'T',
              'G': 'G',
              'C': 'C'
    }
    
    
    return IUPAC[pos]

def mafft_align_contigs(fpath):
    
    from Bio import AlignIO
    import subprocess
    from Bio.Align.Applications import MafftCommandline
    import warnings
    
    # Align the sequences
    mafft_cline = MafftCommandline(input=fpath, op=3.0)
    child = subprocess.Popen(str(mafft_cline),
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             shell=True)
    
    align = AlignIO.read(child.stdout, "fasta")
    
    if len(align[0].seq) == str(align[0].seq).count('-'):
        warnings.warn('bad alignment pair %s'%input)

    if len(align[1].seq) == str(align[1].seq).count('-'):
        warnings.warn('bad alignment pair %s'%input)
    
    return align

def get_aln_cumulative_entropy(align, char_type = 'dna'):
    # get cumulative entropy for alignment
    from reprophylo import entropy

    entropies =[]
    for i in range(align.get_alignment_length()):
        column = align[:,i]
        entropies.append(entropy(column, char_type))
    cum_entropy = sum(entropies) 
    #print "Entropies:", str(entropies[:10]),'... Total:', cum_entropy
    return cum_entropy

def iter_half_matrix_indices(iterable):
    """
    >>> iterable = [1,2,3]
    >>> for i,j in iter_half_matrix_indices(iterable):
    >>>     print i, j
    0 1
    0 2
    1 2
    """
    length = len(iterable)
    for i in range(len(iterable)):
        for j in range(i+1, len(iterable)):
            yield (i, j)


if __name__ == "__main__":
    import doctest
    doctest.testmod()    