#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''Taking a test to translate, its reference file, the Moses init file of
a trained system and an HMM alignment model; this software simulates a user
working on a segment-based IMT framework.'''

__author__ = "Miguel Domingo"
__copyright__ = "Copyright (C) 2016 PRHLT"
__license__ = "Apache 2.0"
__email__ = "midobal@prhlt.upv.es"
__version__ = "1.0"

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import sys
import os
import SBIMT as SB

##############################################################################
###############################   CLASS Simulation   #########################
##############################################################################


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

        # New segments are obtained by applied the longest common subsequence
        # algorithm to hypothesis and reference.
        target_index, reference_index = self.commonWords(hyp)

        # Due to matching errors with the LCS algorithm, we have to make sure
        # that all previous validated segments are part of the algorithm output
        # (if they're not, that means the algorithm has validated a segment
        # inconsistently).
        pending_references = [False for ref in self.reference]
        for n in range(self.reference_size):
            if self.reference_in_segment[n] is not None:
                pending_references[n] = True

        # Check new validated segments word by word.
        validated_segments = []
        validated_segments_reference = []
        current_segments = []
        current_segments_ref = []
        current_segment_index = -1
        old_reference_in_segment = [n for n in self.reference_in_segment]
        for n in range(len(reference_index)):

            # If current word belongs to a previous validated segment:
            if self.reference_in_segment[reference_index[n]] is not None:
                # If the corresponding target is not the one that is supposed
                # to be:
                if target_index[n] not in segments[old_reference_in_segment[
                        reference_index[n]]]:
                    # Discard the word as an error.
                    continue
                if current_segments != []:
                    # Store the current saved words as a new validated segment.
                    validated_segments.append(current_segments)
                    validated_segments_reference.append(current_segments_ref)
                    current_segments = []
                    current_segments_ref = []
                # Unmark the word as pending.
                current_segment_index = self.reference_in_segment[
                    reference_index[n]]
                pending_references[reference_index[n]] = False
                continue

            # Otherwise, check if there has been a matching error with LCS
            # algorithm. (For a given word, it is considered to be errouneus
            # either if there exists a previous pending word--a word that
            # should have been validated before that one--or or a future
            # pending word that is not on the list of new validated words.)
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

            # If no error is detected, add the word to the new validated
            # segments:
            if current_segments == []:  # If the list of current new segment
                # is empty:
                # Add the word as a beggining of a new segment.
                current_segments.append(target_index[n])
                current_segments_ref.append(reference_index[n])
                self.reference_in_segment[reference_index[n]] = \
                    current_segment_index + 1
                for index in range(reference_index[n] + 1,
                                   self.reference_size):
                    if self.reference_in_segment[index] is not None:
                        self.reference_in_segment[index] += 1

            elif (target_index[n] == current_segments[-1] + 1
                  and reference_index[n] == current_segments_ref[-1] + 1):
                # Else, if the word is following one of the current new
                # segment:
                # Add the word to the current new segment.
                current_segments.append(target_index[n])
                current_segments_ref.append(reference_index[n])
                self.reference_in_segment[reference_index[n]] = \
                    current_segment_index + 1

            else:  # Otherwise:
                # Add the current new segment to the new segment lists and
                # create a current new segment with the word as a beggining
                # of the segment.
                validated_segments.append(current_segments)
                validated_segments_reference.append(current_segments_ref)
                current_segments = [target_index[n]]
                current_segments_ref = [reference_index[n]]
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
            validated_segments_reference.append(current_segments_ref)

        # Finally, make sure that there are no pending words.
        for n in range(self.reference_size):
            if pending_references[n]:
                # Discard all new segments if an error is detected.
                self.reference_in_segment = old_reference_in_segment
                return []

        # Validate them otherwise.
        session.segmentValidation(validated_segments)

    def commonWords(self, hyp):
        """
        This function computes the longest common subsequence between two
        sequences and returns a vector with the indexes of the position of
        the subsequence in the first sequence, and a vector with the indexes
        of the position of the subsequence in the second sequence. The first
        sequence corresponds to the current hypothesis and the second
        corresponds to the reference.
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

        return [int(index) for index in list(reversed(
            LCS_xpath[0][0].split()))], [int(index) for index in list(
                reversed(LCS_ypath[0][0].split()))]

    def mergeSegments(self, session):
        """
        This method simulates a user merging two consecutive segments. The
        method receives an object of the SBIMT class that contains the
        current session.
        """

        # Current hypothesis and validated segments.
        segments = session.getSegments()
        merged = False

        # SPECIAL CASE: if there are incorrect words before the first segment
        # (the first segment containing the begin of the reference), those
        # words must be deleted so that the next hypothesis starts with
        # the first segment.
        if self.reference_in_segment[0] != None and segments[0][0] != 0:
            session.mergeSegments(-1, 0)
            return True

        # Word by word, check if a word and its previous word belong to two
        # different segments.
        for n in range(1, self.reference_size):

            if self.reference_in_segment[n] == None:
                break
                continue

            # If they do:
            if (self.reference_in_segment[n - 1] is not None
                and self.reference_in_segment[n - 1]
                    != self.reference_in_segment[n]):
                merged = True
                # Merge those segments.
                session.mergeSegments(self.reference_in_segment[n - 1],
                                      self.reference_in_segment[n])
                # Update reference_in_segment list.
                old_index = self.reference_in_segment[n]
                for index in range(n, self.reference_size):
                    if (self.reference_in_segment[index] is not None
                            and self.reference_in_segment[index] >= old_index):
                        self.reference_in_segment[index] -= 1
        return merged

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

        # If the reference's first word hasn't been validated:
        if self.reference_in_segment[0] is None:
            # Correct the hypothesis' first word.
            self.reference_in_segment[0] = 0
            # Update reference_in_segment list (a new segment is going to
            # be added at the beggining).
            for index in range(1, self.reference_size):
                if self.reference_in_segment[index] is not None:
                    self.reference_in_segment[index] += 1
            session.wordCorrection(0, self.reference[0])
            return self.reference[0]

        # Otherwise, look for the reference's first unvalidated word.
        for n in range(1, self.reference_size):
            if self.reference_in_segment[n] is None:
                # Once located, find the target which should be corrected.
                # (The target in the hypothesis that follows the segment in
                # which the previous reference--n - 1--was located.)
                self.reference_in_segment[n] = self.reference_in_segment[
                    n - 1] + 1
                # Update reference_in_segment (a new segment is going to be
                # added after the segment in which the previous reference was
                # located).
                for index in range(n + 1, self.reference_size):
                    if self.reference_in_segment[index] is not None:
                        self.reference_in_segment[index] += 1
                session.wordCorrection(segments[
                    self.reference_in_segment[n - 1]][-1] + 1,
                                        self.reference[n])
                return self.reference[n]

        # If all word in the reference have already been validated:

        if hyp[-1] != self.reference[-1]:  # If hypothesis and reference are
            # different: Input an end of sentence.
            session.endOfSentence()
            return '#'

        # Otherwise, there's no correction to be made.
        return ''

##############################################################################
###############################   END Simulation   ###########################
##############################################################################


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
    session = SB.SBIMT(moses_ini, alignments_path, prob_threshold)
    reference_file = open(ref, 'r')
    total_sentences = 0
    for s in open(src):
        total_sentences += 1
    current_sentence = 1

    # IMT session (sentence by sentence).
    for source in open(src):

        # Data initialization.
        reference = reference_file.readline().strip().split()
        sim = Simulation(reference)

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
