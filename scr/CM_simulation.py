#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import argparse
import SBIMT as SB

###############################################################################
###############################   CLASS Simulation   ##########################
###############################################################################


class Simulation:
    """
    This class simulates a user in a segment-based IMT framework.
    """

    def __init__(self, ref):
        """
        This method initializes the simulation of a new sentence. The method
        receives the reference of the sentence.
        """
        self.reference = ref
        self.reference_in_segment = [None for r in ref]
        self.reference_size = len(ref)

    def validateSegments(self, session):
        """
        This method simulates a user validating segments. The method receives
        an object of the SBIMT class that contains the current session.
        """

        # Current hypothesis and validated segments.
        hyp = session.getTranslation().split()
        segments = session.getSegments()

        # New segments are obtained by appliying the longest common
        # subsequence algorithm to hypothesis and reference.
        target_index, reference_index = self.commonWords(hyp)

        # Due to matching errors with the LCS algorithm, we have to make
        # sure that all previous validated segments are part of the
        # algorithm output (if they're not, that means the algorithm
        # has validated a segment inconsistently).
        pending_references = [False for ref in self.reference]
        for n in range(self.reference_size):
            if self.reference_in_segment[n] is not None:
                pending_references[n] = True

        # Check new validated segments word by word.
        validated_segments = []
        validated_segments_reference = []
        current_segments = []
        current_segments_reference = []
        current_segment_index = -1
        old_reference_in_segment = [n for n in self.reference_in_segment]
        for n in range(len(reference_index)):

            # If current word belongs to a previous validated segment:
            if self.reference_in_segment[reference_index[n]] is not None:
                # If the corresponding target is not the one that is
                # supposed to be:
                if target_index[n] not in segments[
                        old_reference_in_segment[reference_index[n]]]:
                    # Discard the word as an error.
                    continue
                # Else:
                if current_segments != []:
                    # Store the current saved words as a new validated segment.
                    validated_segments.append(current_segments)
                    validated_segments_reference.append(
                        current_segments_reference)
                    current_segments = []
                    current_segments_reference = []
                # Unmark the word as pending.
                current_segment_index = self.reference_in_segment[
                    reference_index[n]]
                pending_references[reference_index[n]] = False
                continue

            # Otherwise, check if there has been a matching error with
            # LCS algorithm. (For a given word, it is considered to be
            # errouneus either if there exists a previous pending word--a
            # word that should have been validated before that one--or
            # a future pending word that is not on the list of new
            # validated words.)
            erroneus_segment = False
            for index in range(reference_index[n] - 1, -1, -1):
                if pending_references[index]:
                    erroneus_segment = True
                    break
            if not erroneus_segment:
                for index in range(reference_index[n] + 1,
                                   self.reference_size):
                    if (pending_references[index]
                            and index not in reference_index[n + 1:]):
                        erroneus_segment = True
            if erroneus_segment:  # Discard the word if an error is detected.
                continue

            # If no error is detected, add the word to the new
            # validated segments:

            if current_segments == []:  # If the list of current new segment
                # is empty:
                # Add the word as a beggining of a new segment.
                current_segments.append(target_index[n])
                current_segments_reference.append(reference_index[n])
                self.reference_in_segment[reference_index[n]] = \
                    current_segment_index + 1
                for index in range(reference_index[n] + 1,
                                   self.reference_size):
                    if self.reference_in_segment[index] is not None:
                        self.reference_in_segment[index] += 1

            elif (target_index[n] == current_segments[-1] + 1
                  and reference_index[n] == current_segments_reference[
                      -1] + 1):  # Else,
                # if the word is following one of the current new segment:
                # Add the word to the current new segment.
                current_segments.append(target_index[n])
                current_segments_reference.append(reference_index[n])
                self.reference_in_segment[reference_index[n]] = \
                    current_segment_index + 1

            else:  # Otherwise:
                # Add the current new segment to the new segment lists and
                # create a current new segment with the word as a beggining
                # of the segment.
                validated_segments.append(current_segments)
                validated_segments_reference.append(current_segments_reference)
                current_segments = [target_index[n]]
                current_segments_reference = [reference_index[n]]
                current_segment_index += 1
                self.reference_in_segment[reference_index[n]] = \
                    current_segment_index + 1
                for index in range(reference_index[n] + 1,
                                   self.reference_size):
                    if self.reference_in_segment[index] is not None:
                        self.reference_in_segment[index] += 1

        # After checking all words, make sure there isn't a current new
        # segment pending to be added to the list.
        if current_segments != []:
            validated_segments.append(current_segments)
            validated_segments_reference.append(current_segments_reference)

        # Finally, make sure that there are no pending words.
        for n in range(self.reference_size):
            if pending_references[n]:
                # Discard all new segments if an error is detected.
                self.reference_in_segment = old_reference_in_segment
                return False

        # Otherwise, if there are segments to validate,
        if validated_segments == []:
            return False
        # validate them
        session.segmentValidation(validated_segments)
        return True

    def validateSegmentsCM(self, session):
        """
        This method simulates a user validating segments. The user trustes
        the system blindly and validates those segments for which the
        system has a confidence greater than a threshold. The method
        receives an object of the SBIMT class that contains the current
        session.
        """

        # Current hypothesis, validated segments and new segments.
        hyp = session.getTranslation().split()
        new_segments = session.mostLikelySegments()

        # See which words from the reference have been validated.
        # (Use longest common subsequence to relate reference and
        # target--the target is used as an intermediate to know
        # in which segment each reference word is located.)
        target_in_new_segments = [(hyp[trg], n)
                                  for n in range(len(new_segments))
                                  for trg in new_segments[n]]
        target_index, reference_index = self.commonWords(
            [trg[0] for trg in target_in_new_segments])
        segment_to_reference = [[] for segment in new_segments]
        for n in range(len(reference_index)):
            self.reference_in_segment[int(reference_index[n])] = \
                target_in_new_segments[int(target_index[n])][1]
            segment_to_reference[target_in_new_segments[
                int(target_index[n])][1]].append(int(reference_index[n]))

        # Missing words in a segment are considered as if it were in
        # the segment.
        for segment in segment_to_reference:
            if segment != []:
                for n in range(segment[0] + 1, segment[-1]):
                    self.reference_in_segment[n] = \
                        self.reference_in_segment[n - 1]

        # If there are new segments,
        if new_segments == []:
            return False
        # validate them.
        session.segmentValidation(new_segments)
        return True

    def commonWords(self, hyp):
        """
        This function computes the longest common subsequence between
        two sequences and returns a vector with the indexes of the
        position of the subsequence in the first sequence, and a vector
        with the indexes of the position of the subsequence in the second
        sequence. The first sequence corresponds to the current hypothesis
        and the second corresponds to the reference.
        """

        m = len(hyp)
        n = self.reference_size

        x_list = hyp
        x_path = [str(tr) for tr in range(m)]
        y_list = self.reference
        y_path = [str(ref) for ref in range(n)]

        LCS = [['' for i in range(n + 1)] for j in range(m + 1)]
        LCS_xpath = [['' for i in range(n + 1)] for j in range(m + 1)]
        LCS_ypath = [['' for i in range(n + 1)] for j in range(m + 1)]

        for i in range(m - 1, -1, -1):
            for j in range(n - 1, -1, -1):
                LCS[i][j] = LCS[i + 1][j + 1]
                LCS_xpath[i][j] = LCS_xpath[i + 1][j + 1]
                LCS_ypath[i][j] = LCS_ypath[i + 1][j + 1]

                if x_list[i] == y_list[j]:
                    LCS[i][j] += ' ' + x_list[i]
                    LCS_xpath[i][j] += ' ' + x_path[i]
                    LCS_ypath[i][j] += ' ' + y_path[j]

                if len(LCS[i][j + 1]) > len(LCS[i][j]):
                    LCS[i][j] = LCS[i][j + 1]
                    LCS_xpath[i][j] = LCS_xpath[i][j + 1]
                    LCS_ypath[i][j] = LCS_ypath[i][j + 1]

                if len(LCS[i + 1][j]) > len(LCS[i][j]):
                    LCS[i][j] = LCS[i + 1][j]
                    LCS_xpath[i][j] = LCS_xpath[i + 1][j]
                    LCS_ypath[i][j] = LCS_ypath[i + 1][j]

        return ([int(index) for index in list(
            reversed(LCS_xpath[0][0].split()))],
                [int(index) for index in list(
                    reversed(LCS_ypath[0][0].split()))])

    def mergeSegments(self, session):
        """
        This method simulates a user merging two consecutive segments.
        The method receives an object of the SBIMT class that contains
        the current session.
        """

        # Current validated segments.
        segments = session.getSegments()

        # SPECIAL CASE: if there are incorrect words before the first
        # segment (the first segment containing the begin of the reference),
        # those words must be deleted so that the next hypothesis starts with
        # the first segment.
        if self.reference_in_segment[0] is not None and segments[0][0] != 0:
            session.mergeSegments(-1, 0)

        # Word by word, check if a word and its previous word belong to two
        # different segments.
        for n in range(1, self.reference_size):
            # If they do:
            if (self.reference_in_segment[n] is not None
                and self.reference_in_segment[n - 1] is not None
                and self.reference_in_segment[
                    n - 1] == self.reference_in_segment[n] - 1):
                # Merge those segments.
                session.mergeSegments(self.reference_in_segment[n - 1],
                                      self.reference_in_segment[n])
                # Update reference_in_segment list.
                old_index = self.reference_in_segment[n]
                for index in range(n, self.reference_size):
                    if (self.reference_in_segment[index] is not None
                            and self.reference_in_segment[index] >= old_index):
                        self.reference_in_segment[index] -= 1

    def wordCorrection(self, session):
        """
        This method simulates a user correcting a word. Without loss of
        generality and for simplicity's sake, the first word corrected is
        the leftmost wrong word of the hypothesis. The method receives an
        object of the SBIMT class that contains the current session.
        """

        # Current hypothesis and validated segments.
        hyp = session.getTranslation().split()
        segments = session.getSegments()

        # Target in segments
        target_in_segments = [None for trg in hyp]
        for n in range(len(segments)):
            for trg in segments[n]:
                target_in_segments[trg] = n

        # Search for the first wrong word:
        pending_words = []
        if segments == []:
            trg = 0
        else:
            trg = segments[0][-1] + 1
        for n in range(self.reference_size):
            if self.reference_in_segment[n] is None:
                # If the first wrong word is the reference's first word:
                if n == 0:
                    # Ensure the hypothesis doesn't start as it should.
                    # (The one validating is the system, not the user!)
                    if hyp[0] == self.reference[0]:
                        # If it is, mark it as correct and look for the next
                        # wrong word.
                        self.reference_in_segment[0] = 0
                        # pending_words.append([0])
                        session.segmentValidation([[0]])
                        trg = 1
                        continue
                    # Correct the hypothesis' first word,
                    self.reference_in_segment[0] = 0
                    session.wordCorrection(0, self.reference[0])
                    # and update reference_in_segment list (a new segment
                    # is going to be added at the beggining).
                    for index in range(1, self.reference_size):
                        if self.reference_in_segment[index] is not None:
                            self.reference_in_segment[index] += 1
                    return self.reference[0]

                # Else, check if the word to correct is actually right.
                # (The one validating is the system, not the user!)
                if trg < len(hyp) and hyp[trg] == self.reference[n]:
                    # If it is, mark it as correct and look for the next
                    # wrong word.
                    self.reference_in_segment[n] = 0
                    # pending_words.append([trg])
                    session.segmentValidation([[trg]])
                    session.mergeSegments(0, 1)
                    trg += 1
                    # Make sure that the next word is not in a segment
                    if trg < len(hyp) and target_in_segments[trg] is not None:
                        session.mergeSegments(0, 1)
                        for r in range(self.reference_size):
                            if (self.reference_in_segment[r] is not None
                                    and self.reference_in_segment[r] > 0):
                                self.reference_in_segment[r] -= 1
                        first_segment = False
                        for r in range(self.reference_size - 1, -1, -1):
                            if self.reference_in_segment[r] == 0:
                                first_segment = True
                            elif first_segment:
                                self.reference_in_segment[r] = 0
                        trg = segments[0][-1] + 1
                    continue

                # Otherwise, correct the word,
                self.reference_in_segment[n] = max(1,
                                                   min(1, len(
                                                       pending_words) + 1))
                # self.reference_in_segment[n] = 0
                session.wordCorrection(trg, self.reference[n])
                # update reference_in_segment list (a new segment is
                # going to be added at the beginning)
                for index in range(n + 1, self.reference_size):
                    if self.reference_in_segment[index] is not None:
                        self.reference_in_segment[index] += 1
                # and validate the pending words.
                return self.reference[n]

        # Finally, if all words in the reference have already been validated
        # and hypothesis and reference differ:
        if segments[-1][-1] != len(hyp) - 1:
            # Input an end of sentence.
            session.endOfSentence()
            return '#'

        # Otherwise, there's no correction to be made.
        return ''

    def wordCorrectionCM(self, session):
        """
        This method simulates a user correcting a word. The user corrects
        the word which the system indicates that it has the lest confidence
        in that word being correct. The method receives an object of the
        SBIMT class that contains the current session.
        """

        # Current hypothesis, validated segments and word to correct.
        hyp = session.getTranslation().split()
        segments = session.getSegments()
        trg = session.mostUnlikelyTarget()

        # There is no word to correct
        if trg is None:
            for n in range(self.reference_size):
                if self.reference_in_segment[n] is None:
                    if n == 0:
                        session.wordCorrection(0, self.reference[0])
                        self.reference_in_segment[n] = 0
                        for r in range(n + 1, self.reference_size):
                            if self.reference_in_segment[r] is not None:
                                self.reference_in_segment[r] += 1
                        return self.reference[0]
                    else:
                        session.wordCorrection(
                            segments[self.reference_in_segment[
                                n - 1]][-1] + 1, self.reference[n])
                        self.reference_in_segment[n] = \
                            self.reference_in_segment[n - 1] + 1
                        for r in range(n + 1, self.reference_size):
                            if self.reference_in_segment[r] is not None:
                                self.reference_in_segment[r] += 1
                        return self.reference[n]
            return ''

        print('Word to correct:', hyp[trg])

        if (segments != [] and trg < segments[0][0]
                and self.reference_in_segment[0] == 0):
            session.mergeSegments(-1, 0)
            return 'Ò.Ó'

        # Word to correct is the first word of the hypothesis.
        if trg == 0:
            self.reference_in_segment[0] = 0
            for n in range(1, self.reference_size):
                if self.reference_in_segment[n] is not None:
                    self.reference_in_segment[n] += 1
            session.wordCorrection(trg, self.reference[0])
            return self.reference[0]

        # Word to correct is the last word of the hypothesis.
        if trg == len(hyp) - 1:
            # Check if the reference's last word is already in the hypothesis.
            if self.reference_in_segment[-1] is None:
                session.wordCorrection(trg, self.reference[-1])
                self.reference_in_segment[-1] = len(segments) - 1
                return self.reference[-1]
            # If it is, input an end of sentence.
            else:
                session.endOfSentence()
                return '#'
                session.wordCorrection(trg, '')
                return ''

        # Target position in segments.
        target_in_segments = [None for tr in hyp]
        for n in range(len(segments)):
            for tr in segments[n]:
                target_in_segments[tr] = n

        # Word to correct is at the beginning of a segment.
        if target_in_segments[trg + 1] is not None:
            segment = target_in_segments[trg + 1]
            ref = -1
            for n in range(1, self.reference_size):
                if self.reference_in_segment[n] == segment:
                    ref = n - 1
                    break
            if ref != -1 and self.reference_in_segment[ref] is not None:
                session.mergeSegments(self.reference_in_segment[ref],
                                      self.reference_in_segment[ref] + 1)
                for rf in range(self.reference_size):
                    if (self.reference_in_segment[rf] is not None
                        and self.reference_in_segment[
                            rf] > self.reference_in_segment[ref]):
                        self.reference_in_segment[rf] -= 1
                return 'Ò.Ó'
            # There may not be any reference in the segment,
            # in which case we don't know the correction. The
            # same happens if the reference is in another segment.
            if ref == -1 or self.reference_in_segment[ref] is not None:
                session.endOfSentence()
                return '#'
            self.reference_in_segment[ref] = segment
            for n in range(ref + 1, self.reference_size):
                if self.reference_in_segment[n] is not None:
                    self.reference_in_segment[n] += 1
            session.wordCorrection(trg, self.reference[ref])
            return self.reference[ref]

        # Word to correct is at the end of a segment.
        if target_in_segments[trg - 1] is not None:
            segment = target_in_segments[trg - 1]
            ref = -1
            for n in range(self.reference_size - 1):
                if self.reference_in_segment[n] == segment:
                    ref = n + 1
            if ref != -1 and self.reference_in_segment[ref] is not None:
                session.mergeSegments(self.reference_in_segment[ref] - 1,
                                      self.reference_in_segment[ref])
                for rf in range(self.reference_size):
                    if (self.reference_in_segment[rf] is not None
                        and self.reference_in_segment[
                            rf] >= self.reference_in_segment[ref]):
                        self.reference_in_segment[rf] -= 1
                return 'Ò.Ó'
            # There may not be any reference in the segment,
            # in which case we don't know the correction. The
            # same happens if the reference is in another segment.
            if ref == -1 or self.reference_in_segment[ref] is not None:
                session.endOfSentence()
                return '#'
            self.reference_in_segment[ref] = segment + 1
            for n in range(ref + 1, self.reference_size):
                if self.reference_in_segment[n] is not None:
                    self.reference_in_segment[n] += 1
            session.wordCorrection(trg, self.reference[ref])
            return self.reference[ref]

        # Otherwise, there's no correction to be made.
        return ''

##############################################################################
###############################   END Simulation   ###########################
##############################################################################


def parse_args():
    parser = argparse.ArgumentParser(description='Taking a test to translate, \
                                     its reference file, the Moses init file \
                                     of a trained system and an HMM alignment \
                                     model; this software simulates a user \
                                     working on a segment-based IMT \
                                     framework.')
    parser.add_argument('-s', '--sources', metavar='source_file',
                        required=True, help='file containing the source \
                        segments.')
    parser.add_argument('-r', '--references', metavar='reference_file',
                        required=True, help='file containing the reference\
                        segments.')
    parser.add_argument('-c', '--config', metavar='moses_ini',
                        required=True, help='file containing moses\
                        configuration.')
    parser.add_argument('-a', '--alignments', metavar='alignments_file',
                        required=True, help='file containing the alignments\
                        generated by alignments.sh')
    parser.add_argument('-v', '--verbose', default=False,
                        action='store_true', required=False, help='Activate \
                        verbose mode.')
    parser.add_argument('-x', '--xml', default=False,
                        action='store_true', required=False, help='Show \
                        XML markup.')
    parser.add_argument('-p', '--probability', metavar='threshold',
                        required=False, type=float, help='probability \
                        threshold. (Default 0.)', default=0)
    parser.add_argument('-t', '--confidence', metavar='threshold',
                        required=False, type=float, help='condidence \
                        threshold. (Default 0.)', default=0)
    parser.add_argument('-l', '--level', metavar='[1-3]',
                        required=False, type=int, help='what to use CM for: \
                        1- use CM to validate segments. 2- use CM to \
                        correct words. 3- use CM to both validate and \
                        correct.\n(Default 3.)', default=3,
                        choices=range(1, 4))
    parser.add_argument('-m', '--moses', metavar='moses_bin',
                        required=False, default='/opt/moses/bin/moses',
                        help='Path to moses bin. (Default: \
                        /opt/moses/bin/moses.)')

    return parser.parse_args()


if __name__ == "__main__":
    """
    Start of the simulation.
    """

    # Check arguments.
    args = parse_args()

    # Session set-up.
    sys.stderr.write("\x1b[2J\x1b[H")
    sys.stderr.write('Preparing systems  (it may take a while)...\n')
    session = SB.SBIMT(args.config, args.alignments, args.probability,
                       args.confidence)
    reference_file = open(args.references, 'r')
    translation_file = open('test.hyp', 'w')
    total_sentences = 0
    for s in open(args.sources):
        total_sentences += 1
    current_sentence = 1

    # IMT session (sentence by sentence).
    for source in open(args.sources):

        # Data initialization.
        reference = reference_file.readline().strip().split()
        sim = Simulation(reference)

        # Show progess.
        sys.stderr.write("\x1b[2J\x1b[H")
        sys.stderr.write('Progress: ' + str(current_sentence) + '/'
                         + str(total_sentences) + ' [' + "{0:.2f}".format(
                             current_sentence / float(total_sentences) * 100)
                         + ' %]\n')
        current_sentence += 1

        # Session initialization.
        session.newSentence(source.split())
        validated_translation = False

        if args.verbose:
            print('SOURCE: ' + source.strip())
            print('REFERENCE: ' + ' '.join(reference))

        # Iterative process.
        while not validated_translation:

            # Load new hypothesis.
            session.newHypothesis()

            if args.verbose:
                print('TRANSLATION:', session.getTranslation())

            # Check if new hypothesis is the desired translation.
            if session.getTranslation() == ' '.join(reference):
                session.validateTranslation()
                validated_translation = True
                continue

            # Segment Validation.
            if args.level == 1 or args.level == 3:
                segment_validated = sim.validateSegmentsCM(session)
            else:
                segment_validated = sim.validateSegments(session)

            # Word Correction / Sentence Validation.
            if args.level == 2 or args.level == 3:
                new_word = 'Ò.Ó'
                while new_word == 'Ò.Ó':
                    new_word = sim.wordCorrectionCM(session)
            else:
                new_word = sim.wordCorrection(session)

            if new_word == '':  # and not segment_validated:
                session.validateTranslation()
                if args.verbose:
                    print('')
                    print('')
                    print('CORRECTED WORD: ')
                    print('WORD SEGMENTS:', session.getWordSegments())
                    print('DELETED WORDS:', session.getDeletedWords())
                    print('')
                validated_translation = True
                continue

            # XML Generation.
            session.generateXML()

            if args.verbose:
                print('')
                print('')
                print('CORRECTED WORD:',  new_word)
                print('WORD SEGMENTS:', session.getWordSegments())
                print('DELETED WORDS:', session.getDeletedWords())
                print('')
                if args.xml:
                    print('XML:', session.getXML())

        if args.verbose:
            print('')
            print('Word Strokes: ', session.getWordStrokes())
            print('Mouse Actions:', session.getMouseActions())

            print("-----------------------------------------")
            print('')

        translation_file.write(session.getValidatedTranslation() + '\n')

    # Show metrics.
    print('WSR:', "{0:.1f}".format(session.getWSR()))
    print('MAR:', "{0:.1f}".format(session.getMAR()))
    print('WDR:', "{0:.1f}".format(session.getWDR()))
    translation_file.close()
