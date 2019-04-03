#! /usr/bin/env python3

import sys

#import and use the find_fird_orf function in find_orf.py to find the first ORF

def find_first_orf(sequence,
        start_codons = ['AUG'],
        stop_codons = ['UAA', 'UAG', 'UGA']):
    """
    Return the first open-reading frame in the DNA or RNA `sequence`.

    An open-reading frame (ORF) is the part of an RNA sequence that is
    translated into a peptide. It must begin with a start codon, followed by
    zero or more codons (triplets of nucleotides), and end with a stop codon.
    If there are no ORFs in the sequence, an empty string is returned.

    Parameters
    ----------
    sequence : str
        A string representing a DNA or RNA sequence (upper or lower-case)
    start_codons : list of strings
        All possible start codons. Each codon must be a string of 3 RNA bases,
        upper or lower-case.
    stop_codons : list of strings
        All possible stop codons. Each codon must be a string of 3 RNA bases,
        upper or lower-case.

    Returns
    -------
    str
       	An uppercase string of the first ORF found in the `sequence` that
        starts with any one of the `start_codons` and ends with any one of the
        `stop_codons`. If no ORF is found an empty string is returned.

    Examples
    --------
    When the whole RNA sequence is an ORF:
    >>> find_first_orf('AUGGUAUAA', ['AUG'], ['UAA'])
    'AUGGUAUAA'

    When the whole DNA sequence is an ORF:
    >>> find_first_orf('ATGGTATAA', ['AUG'], ['UAA'])
    'AUGGUAUAA'

    When there is no ORF:
    >>> find_first_orf('CUGGUAUAA', ['AUG'], ['UAA'])
    ''

    When there is are bases before and after ORF:
    >>> find_first_orf('CCAUGGUAUAACC', ['AUG'], ['UAA'])
    'AUGGUAUAA'
    """
   # Make sure the sequence is valid
    vet_nucleotide_sequence(sequence)

    # Make sure the codons are valid

    # Get copies of everything in uppercase
    seq = sequence.upper()
    starts = [c.upper() for c in start_codons]
    stops = [c.upper() for c in stop_codons]
    # Make sure seq is RNA
    seq = seq.replace('T', 'U')
    orf_pattern_str = r'(' + r'|'.join(starts) + r')([AUGC]{3})*(' + r'|'.join(stops) + r')'

    # Create the regular expression object
    orf_pattern = re.compile(orf_pattern_str)
    # Search the sequence
    match_object = orf_pattern.search(seq)
    if match_object:
        return match_object.group()
    return ''

# import and use the translate_sequence function in translate.py to translate the ORF

def translate_sequence(rna_sequence, genetic_code):
    """Translates a sequence of RNA into a sequence of amino acids.

    Translates `rna_sequence` into string of amino acids, according to the
    `genetic_code` given as a dict. Translation begins at the first position of
    the `rna_sequence` and continues until the first stop codon is encountered
    or the end of `rna_sequence` is reached.

    If `rna_sequence` is less than 3 bases long, or starts with a stop codon,
    an empty string is returned.
    """
    rna_sequence = rna_sequence.upper()
    amino_acid_list = []
    while True:
        if len(rna_sequence) < 3:
            break
        codon, remaining_seq = pop_next_codon(rna_sequence)
        rna_sequence = remaining_seq
        aa = genetic_code[codon]
        if aa == "*":
            break
        amino_acid_list.append(aa)
    return "".join(amino_acid_list)

# import and use the parse_sequence_from_path in find_orf.py to translate the ORF

def parse_sequence_from_path(path):
    # Try to open the path to read from it, and handle exceptions if they
    # arrise
    try:
        file_stream = open(path, 'r')
    except FileNotFoundError as e:
        sys.stderr.write("Sorry, couldn't find path {}".format(path))
        raise e
    except IsADirectoryError as e:
        sys.stderr.write("Sorry, path {} appears to be a directory".format(
                path))
        raise e
    except:
        sys.stderr.write("Sorry, something went wrong when trying to open {}".format(
                path))
        raise
    # If we've reached here, the file is open and ready to read
    sequence = ''
    # A for loop to visit each line in the file
    for line in file_stream:
        # Strip whitespace from the line and concatenate it to the end of the
        # sequence
        sequence += line.strip()
    return sequence

def main():
    import argparse

    # Create a command-line parser object
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)

    # Tell the parser what command-line arguments this script can receive
    parser.add_argument('sequence',
            metavar = 'SEQUENCE',
            type = str,
            help = ('The sequence to search for an open-reading frame. '
                    'If the path flag (\'-p\'/\'--path\') is specified, '
                    'then this should be a path to a file containing the '
                    'sequence to be searched.'))
    parser.add_argument('-p', '--path',
            action = 'store_true',
            help = ('The sequence argument should be treated as a path to a '
                    'containing the sequence to be searched.'))
    parser.add_argument('-s', '--start-codons',
            type = str,
            nargs = '+', # one or more arguments
            default = ['AUG'],
            help = ('One or more possible start codons.'))
    parser.add_argument('-x', '--stop-codons',
            type = str,
            nargs = '+', # one or more arguments
            default = ['UAA', 'UAG', 'UGA'],
            help = ('One or more possible stop codons.'))

    # Parse the command-line arguments into a 'dict'-like container
    args = parser.parse_args()
    # Check to see if the path option was set to True by the caller. If so, parse
    # the sequence from the path
    if args.path:
        sequence = parse_sequence_from_path(args.sequence)
    else:
        sequence = args.sequence

    orf = find_first_orf(sequence = sequence,
            start_codons = args.start_codons,
            stop_codons = args.stop_codons)
    sys.stdout.write('{}\n'.format(orf))

if __name__ == '__main__':
    genetic_code = {'GUC': 'V', 'ACC': 'T', 'GUA': 'V', 'GUG': 'V', 'ACU': 'T', 'AAC': 'N', 'CCU': 'P', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'CAU': 'H', 'AAU': 'N', 'AGU': 'S', 'GUU': 'V', 'CAC': 'H', 'ACG': 'T', 'CCG': 'P', 'CCA': 'P', 'ACA': 'T', 'CCC': 'P', 'UGU': 'C', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'UGC': 'C', 'CAG': 'Q', 'GAU': 'D', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'AGG': 'R', 'GGG': 'G', 'UCC': 'S', 'UCA': 'S', 'UAA': '*', 'GGA': 'G', 'UAC': 'Y', 'GAC': 'D', 'UAG': '*', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'GGC': 'G', 'AUG': 'M', 'CUG': 'L', 'GAG': 'E', 'CUC': 'L', 'AGA': 'R', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'GCU': 'A', 'GAA': 'E', 'AUU': 'I', 'UUG': 'L', 'UUA': 'L', 'UGA': '*', 'UUC': 'F'}
    rna_seq = ("AUG"
            "UAC"
            "UGG"
            "CAC"
            "GCU"
            "ACU"
            "GCU"
            "CCA"
            "UAU"
            "ACU"
            "CAC"
            "CAG"
            "AAU"
            "AUC"
            "AGU"
            "ACA"
            "GCG")
    main()
