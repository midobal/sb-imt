# -*- coding: utf-8 -*-
import sys
import os
from sbimt.sbimt import SBIMT
from sbimt.user import User


def usage():
    """
    This function shows the usage message.
    """
    sys.stderr.write('Usage: ' + sys.argv[0] + ' -s source_file -r '
                     + 'reference_file -m moses_ini -a alignments '
                     + '[options]\n\n')
    sys.stderr.write('Options: \n')
    sys.stderr.write('  -h             show this message.\n')
    sys.stderr.write('  -v             verbose mode on.\n')
    sys.stderr.write('  -XML           show XML Markup.\n')
    sys.stderr.write('  -p theshold    probability threshold. (Default: 0.)\n')
    sys.exit(-1)


def getArguments():
    """
    This function checks the arguments and returns their value.
    """
    src = None
    ref = None
    moses_ini = None
    verbose = False
    XML = False
    alignments = None
    prob_threshold = 0.0

    # Loop through the arguments.
    n = 1
    while n < len(sys.argv):

        if sys.argv[n] == '-s':
            src = sys.argv[n + 1]
            n += 2

        elif sys.argv[n] == '-r':
            ref = sys.argv[n + 1]
            n += 2

        elif sys.argv[n] == '-m':
            moses_ini = sys.argv[n + 1]
            n += 2

        elif sys.argv[n] == '-v':
            verbose = True
            n += 1

        elif sys.argv[n] == '-xml':
            XML = True
            verbose = True
            n += 1

        elif sys.argv[n] == '-a':
            alignments = sys.argv[n + 1]
            n += 2

        elif sys.argv[n] == '-p':
            prob_threshold = float(sys.argv[n + 1])
            n += 2

        else:
            usage()

    # Check that mandatory arguments are present.
    if src is None or ref is None or moses_ini is None or alignments is None:
        usage()

    # Check all paths.
    if not os.path.isfile(src):
        sys.stderr.write('Error opening file ' + src + '\n')
        sys.exit(-1)

    if not os.path.isfile(ref):
        sys.stderr.write('Error opening file ' + ref + '\n')
        sys.exit(-1)

    if not os.path.isfile(moses_ini):
        sys.stderr.write('Error opening file ' + moses_ini + '\n')
        sys.exit(-1)

    if not os.path.isfile(alignments):
        sys.stderr.write('Error opening file ' + alignments + '\n')
        sys.exit(-1)

    # Return arguments.
    return src, ref, moses_ini, verbose, XML, alignments, prob_threshold


if __name__ == "__main__":
    """
    Start of the simulation.
    """

    # Check arguments.
    src, ref, moses_ini, verbose, XML, alignments_path, prob_threshold = \
        getArguments()

    # Session set-up.
    sys.stderr.write("\x1b[2J\x1b[H")
    sys.stderr.write('Preparing systems  (it may take a while)...\n')
    session = SBIMT(moses_ini, alignments_path, prob_threshold)
    reference_file = open(ref, 'r')
    total_sentences = 0
    for s in open(src):
        total_sentences += 1
    current_sentence = 1

    # IMT session (sentence by sentence).
    for source in open(src):

        # Data initialization.
        reference = reference_file.readline().strip().split()
        sim = User(reference)

        # Show progess.
        sys.stderr.write("\x1b[2J\x1b[H")
        sys.stderr.write('Progress: ' + str(current_sentence) + '/'
                         + str(total_sentences)
                         + ' [' + "{0:.2f}".format(current_sentence
                                                   / float(total_sentences)
                                                   * 100) + ' %]\n')
        current_sentence += 1

        # Session initialization.
        session.newSentence(source.split())
        validated_translation = False

        if verbose:
            print('SOURCE: ' + source.strip())
            print('REFERENCE: ' + ' '.join(reference))

        # Iterative process.
        while not validated_translation:

            # Load new hypothesis.
            session.newHypothesis()

            if verbose:
                print('TRANSLATION:', session.getTranslation())

            # Check if new hypothesis is the desired translation.
            if session.getTranslation() == ' '.join(reference):
                session.validateTranslation()
                validated_translation = True
                continue

            # Segment Validation.
            sim.validateSegments(session)

            # Segment Merging.
            sim.mergeSegments(session)

            # Word Correction / Sentence Validation.
            new_word = sim.wordCorrection(session)

            if new_word == '':
                session.validateTranslation()
                if verbose:
                    print('')
                    print('')
                    print('CORRECTED WORD: ')
                    print('WORD SEGMENTS:', session.getWordSegments())
                    print('DELETED WORDS:', session.getDeletedWords())
                    print('')
                break

            # XML Generation.
            session.generateXML()

            if verbose:
                print('')
                print('')
                print('CORRECTED WORD:',  new_word)
                print('WORD SEGMENTS:', session.getWordSegments())
                print('DELETED WORDS:', session.getDeletedWords())
                print('')
                if XML:
                    print('XML:', session.getXML())

        if verbose:
            print('')
            print('Word Strokes: ', session.getWordStrokes())
            print('Mouse Actions:', session.getMouseActions())

            print("-----------------------------------------")
            print('')

    # Show metrics.
    print('WSR:', "{0:.1f}".format(session.getWSR()))
    print('MAR:', "{0:.1f}".format(session.getMAR()))
    print('WDR:', "{0:.1f}".format(session.getWDR()))
