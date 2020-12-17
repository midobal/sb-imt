#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''This software implements the segment-based methodology described by Domingo et al. (2016).

On using this software, please cite the following paper:

  Domingo, M., Peris, Á., and Casacuberta, F. (2016). Interactive-predictive translation
based on multiple word-segments. In Proceedings of the Annual Meeting of
the European Association for Machine Translation, pages 282–291.'''

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
# limitations under the license.

import sys, bisect, subprocess, pty, os, operator
from random import shuffle

########################################################################################
###############################   CLASS Phrase   #######################################
########################################################################################

class Phrase:
    """
    This class stores phrases information (mainly, sources and translations).
    It also contains a flag to know if a phrase is contained in a segment (validated) and
    a pointer to the phrase in which the segment is stored. Phrases data structures
    are ordered according to source.
    """

    def __init__(self, src):
        """
        This method initializes a new Phrase. The method receives the source to which
        is originally associated.
        """
        self.sources = [src]
        self.translation = ''
        self.segment_position = None
        self.validated = False

    def addTranslation(self, trans):
        """
        This method ads more translation to the current translation of the phrase.
        The method receives a string with the new translation.
        """
        if self.translation != '':
            self.translation += ' '
        self.translation += trans

    def addSource(self, src):
        """
        This method ads a source to the list of sources contained in the phrase.
        The method receives the position of the new source.
        """
        if src not in self.sources:
            bisect.insort(self.sources, src)

########################################################################################
###############################   END Phrase   #########################################
########################################################################################



########################################################################################
###############################   CLASS SBIMT   ########################################
########################################################################################

class SBIMT:
    """
    This class implements the segment-based IMT methodology.
    """

    def __init__(self, moses_ini, HMM_alignments, HMM_probability = 0.0, CM_threshold = 0.0, max_phrase_size = 20):
        """
        This method initializes a new segment-based IMT session. The method receives the path
        containing the Moses initialization file, the path to a file containing HMM alignments
        between source and target, a threshold of the probability with which a source and a new
        target should be aligned (default: 0.0), a threshold for the confidence the system have
        for a certaing translation to be correct (default: 0.0) and the maximum size of a phrase
        (default: 20).
        """

        self.ensureMosesVariableDefined()
        self.moses = subprocess.Popen('$MOSES/bin/moses -xml-input '
                                      + 'exclusive -f ' + moses_ini
                                      + ' -print-alignment-info 2> /dev/null',
                                      shell=True, stdin=subprocess.PIPE,
                                      stdout=subprocess.PIPE, bufsize=1,
                                      universal_newlines=True)

        self.loadAlignments(HMM_alignments)
        self.alignment_probability = HMM_probability

        self.CM_threshold = CM_threshold

        self.max_sources = max_phrase_size

        self.word_strokes = 0
        self.mouse_actions = 0
        self.word_deletions = 0
        self.words = 0
        self.characters = 0

    def getSegments(self):
        """
        This method returns a list with all the validated segments.
        """
        return self.segments

    def getWordSegments(self):
        """
        This method returns a string with all the validated segments. It is intended for
        logging purposes.
        """

        return '    '.join([' '.join(segment) for segment in self.segments_string])

    def getDeletedWords(self):
        """
        This method returns a string with the words deleted when merging two segments.
        It is intended for logging purposes.
        """

        #return '    '.join([' '.join(words) for words in self.deleted_words])
        deleted_words = [[] for words in self.deleted_words]
        for n in range(len(self.deleted_words)):
            for w in self.deleted_words[n]:
                deleted_words[n].append(self.target[w])
        return '    '.join([' '.join(words) for words in deleted_words])

    def getXML(self):
        """
        This method returns the current xml version of the source sentence.
        """
        return self.xml

    def getTranslation(self):
        """
        This method returns the current hypothesis.
        """
        return ' '.join(self.target)

    def getValidatedTranslation(self):
        """
        This method returns the translation validated by the user (the current
        translation without the words the user has removed).
        """

        target_deleted = [False for trg in self.target]
        for words in self.deleted_words:
            for word in words:
                target_deleted[word] = True
        return ' '.join([self.target[n] for n in range(self.target_size) if not target_deleted[n]])

    def getWordStrokes(self):
        """
        This method returns the word strokes for the current sentence.
        """

        return self.current_word_strokes

    def getMouseActions(self):
        """
        This method returns the mouse action for the current sentence.
        """

        return self.current_mouse_actions

    def getWSR(self):
        """
        This method computes the WSR for the current session.
        """
        return (self.word_strokes / float(self.words) * 100)

    def getMAR(self):
        """
        This method computes the MAR for the current session.
        """
        return (self.mouse_actions / float(self.characters) * 100)

    def getWDR(self):
        """
        This method computes the Word Deletion Ratio for the current session.
        """
        return (self.word_deletions / float(self.words) * 100)

    def segmentsInTarget(self):
        """
        This method returns the segments position in the new hypothesis, using
        a variation of the longest common subsequence algorithm.
        """
        sgmnts = [w for s in self.segments_string for w in s]
        sgmnts_size = len(sgmnts)

        trgs = [self.target[n] for n in range(self.target_size) if self.target_to_phrase[n] == []]
        trgs_index = [str(n) for n in range(self.target_size) if self.target_to_phrase[n] == []]
        trgs_size = len(trgs)

        LCS = [['' for i in range(sgmnts_size + 1)] for j in range(trgs_size + 1)]
        LCS_path = [['' for i in range(sgmnts_size + 1)] for j in range(trgs_size + 1)]

        for i in range(trgs_size - 1, -1, -1):
            for j in range(sgmnts_size - 1, -1, -1):
                LCS[i][j] = LCS[i + 1][j + 1]
                LCS_path[i][j] = LCS_path[i + 1][j + 1]

                if trgs[i] == sgmnts[j]:
                    LCS[i][j] += ' ' + trgs[i]
                    LCS_path[i][j] += ' ' + trgs_index[i]

                if len(LCS[i][j + 1]) >= len(LCS[i][j]):
                    LCS[i][j] = LCS[i][j + 1]
                    LCS_path[i][j] = LCS_path[i][j + 1]

                if len(LCS[i + 1][j]) >= len(LCS[i][j]):
                    LCS[i][j] = LCS[i + 1][j]
                    LCS_path[i][j] = LCS_path[i + 1][j]

        return [int(n) for n in list(reversed(LCS_path[0][0].split()))]

    def ensureMosesVariableDefined(self):
        """
        This method checks if there is a variable pointing to Moses path and
        aborts the execution if there is not.
        """
        if os.getenv('MOSES') == '':
            sys.stderr.write('Error loading Moses. Please define a variable "$MOSES" pointing to Moses path (e.g., export MOSES=/opt/moses)\n')
            sys.exit(-1)

    def loadAlignments(self, alignments_path):
        """
        This method creates a dictionary with the HMM word alignments between source
        and target. The method receives the path to a file containing the alignments.
        """
        self.alignments = {}

        for lines in open(alignments_path):
            try:
                self.alignments[lines.split()[0]]
                self.alignments[lines.split()[0]][lines.split()[1]] = float(lines.split()[2])
            except:
                self.alignments[lines.split()[0]] = {lines.split()[1] : float(lines.split()[2])}

    def newSentence(self, src):
        """
        This method sets the data structures to translate a new sentence. The method
        receives a list containing the new sentence.
        """

        #Source and phrases structures are initialized according to the new sentence.
        self.source = src
        self.source_size = len(src)
        self.artificial_sources = 0
        self.artificial_source_at_beggining = False
        self.phrases = [Phrase(n) for n in range(len(src))]

        #The xml is initially initialized as the new sentence.
        self.xml = ' '.join(src)

        #The rest of the structures are initialized empty.
        self.target = []
        self.target_size = 0
        self.target_to_phrase = []

        self.segments = []
        self.segments_string = []
        self.segments_size = 0

        self.current_word_strokes = 0
        self.current_mouse_actions = 0

    def newHypothesis(self):
        """
        This method updates the data structures according to the new hyp. This is done
        using the alignment information provided by Moses when using XML Markup.
        """

        #New hypothesis.
        self.moses.stdin.write((self.xml + '\n'))
        hyp = self.moses.stdout.readline().strip()

        #Update data structures.
        self.deleted_words = []
        self.target = hyp.split('|||')[0].split()
        self.target_size = len(self.target)
        alignments = hyp.strip().split('|||')[1].split()
        old_target_to_phrase = self.target_to_phrase
        self.target_to_phrase = [[] for n in range(self.target_size)]

        #Source-target alignments given by Moses.
        for w in alignments:
            src = int(w.split('-')[0])
            trg = int(w.split('-')[1])

            if self.artificial_source_at_beggining and not self.phrases[0].validated:
                #Exception for using prefix-based as a particula case of segment-based.
                src -= 1

            if not self.phrases[src].validated:
                self.target_to_phrase[trg].append(src)

        #Moses doesn't give alignment information about the sources containted in an xml tag.
        #Therefore, information related to segments must be updated independently.
        segments_in_target = self.segmentsInTarget()
        index = 0
        for n in range(self.segments_size):
            for m in range(len(self.segments[n])):
                self.target_to_phrase[segments_in_target[index]] = old_target_to_phrase[self.segments[n][m]]
                self.segments[n][m] = segments_in_target[index]
                index += 1

        #Ensure there is no error.
        for n in range(self.segments_size):
            for trg in range(1, len(self.segments[n])):
                if self.segments[n][trg] != self.segments[n][trg - 1] + 1:
                    self.segments[n][trg] = self.segments[n][trg - 1] + 1

    def segmentValidation(self, segmnts):
        """
        This method updates the data structures according to the list of segments validated
        by the user. The method receives a list of list with the position of the targets
        that conform the segments.
        """

        for n in range(len(segmnts)):

            #Search for the phrase in which the segment is going to be stored.
            phrase_index = self.source_size + self.artificial_sources
            for trg in segmnts[n]:
                for src in self.target_to_phrase[trg]:
                    if src < phrase_index and not self.phrases[src].validated:
                        phrase_index = src

            #If the segment isn't aligned to any source, an artificial one is created.
            #if phrase_index == None:
            if phrase_index == self.source_size + self.artificial_sources:
                self.phrases.append(Phrase(-1))
                phrase_index = self.source_size + self.artificial_sources
                self.artificial_sources += 1

            #Phrase is validated. (The segment is stored in it.)
            self.phrases[phrase_index].validated = True
            self.phrases[phrase_index].segment_position = phrase_index
            self.phrases[phrase_index].translation = ' '.join([self.target[tr] for tr in segmnts[n]])

            #Phrases aligned with any word of the segment is "merged" into the phrase that represents the segment.
            for trg in segmnts[n]:
                for src in self.target_to_phrase[trg]:
                    if src != phrase_index and not self.phrases[src].validated:
                        self.phrases[src].validated = True
                        self.phrases[src].translation = ''
                        self.phrases[src].segment_position = phrase_index
                        self.phrases[phrase_index].addSource(src)
                self.target_to_phrase[trg].insert(0, phrase_index)

            #The segment is added into the segment list:
            if self.segments_size == 0: #At the end of the list, if the list is empty.
                self.segments.append(segmnts[n])
                self.segments_string.append(self.phrases[phrase_index].translation.split())
                self.segments_size += 1

            else:#At the corresponding position, otherwise.
                inserted = False
                for index in range(self.segments_size):
                    if segmnts[n][-1] < self.segments[index][0]:
                        self.segments.insert(index, segmnts[n])
                        self.segments_string.insert(index, self.phrases[phrase_index].translation.split())
                        self.segments_size += 1
                        inserted = True
                        break
                if not inserted:
                    self.segments.append(segmnts[n])
                    self.segments_string.append(self.phrases[phrase_index].translation.split())
                    self.segments_size += 1

            #Finally, mouse actions are computed (1 action for single-word segments, 2 otherwise).
            if len(segmnts[n]) > 1:
                self.current_mouse_actions += 2
            else:
                self.current_mouse_actions += 1

    def mergeSegments(self, left_segment, right_segment):
        """
        This method updates the data structures to merge two consecutives segments
        into one. Words between those segments are stored in a list for logging purposes.
        Them method receives the position in the segments list of both segments to merge.
        """

        #SPECIAL CASE: If left_segment equal to '-1', the user is deleting words at the beggining so that the
        #hypothesis starts with the first segment.
        if left_segment == -1:
            self.artificial_source_at_beggining = True
            #self.deleted_words.append(self.target[:self.segments[0][0]])
            self.deleted_words.append([n for n in range(self.segments[0][0])])
            if self.segments[0][0] > 1:
                self.current_mouse_actions += 2
            else:
                self.current_mouse_actions += 1
            return

        #First of all, deleted words are stored and mouse actions are computed (1 action when deleting
        #a single word and 2 when deleting more than one).
        deleted_words = []
        for n in range(self.segments[left_segment][-1] + 1, self.segments[right_segment][0]):
            #deleted_words.append(self.target[n])
            deleted_words.append(n)
        if deleted_words != []:
            self.deleted_words.append(deleted_words)
            if len(deleted_words) > 1:
                self.current_mouse_actions += 2
            else:
                self.current_mouse_actions += 1
            self.word_deletions += len(deleted_words)

        #Location of the phrases where the segments are stored.
        right_phrase = self.phrases[self.target_to_phrase[self.segments[right_segment][0]][0]].segment_position
        left_phrase = self.phrases[self.target_to_phrase[self.segments[left_segment][0]][0]].segment_position

        #If one of the segments contains an artificial source, it is treated differently (to remove the artificial source):
        if self.artificial_sources > 0 and (right_phrase >= self.source_size or left_phrase >= self.source_size):
            #If the articial source is in the right segment, the segments are merged into the left phrase and
            #the right phrase is deleted (target_to_phrase list must be updated).
            if right_phrase >= self.source_size:
                self.phrases[left_phrase].addTranslation(self.phrases[right_phrase].translation)
                del self.phrases[right_phrase]
                self.artificial_sources -= 1
                for trg in range(self.target_size):
                    if self.target_to_phrase[trg] != [] and self.target_to_phrase[trg][0] == right_phrase:
                        self.target_to_phrase[trg][0] = left_phrase
                    elif self.target_to_phrase[trg] != [] and self.target_to_phrase[trg][0] > right_phrase:
                        self.target_to_phrase[trg][0] -= 1
                for source in range(self.source_size, self.source_size + self.artificial_sources):
                    self.phrases[source].segment_position = source

            #If the articial source is in the left segment, the segments are merged into the right phrase and
            #the left phrase is deleted (target_to_phrase list must be updated).
            else:
                self.phrases[right_phrase].translation = self.phrases[left_phrase].translation + ' ' + self.phrases[right_phrase].translation
                del self.phrases[left_phrase]
                self.artificial_sources -= 1
                for trg in range(self.target_size):
                    if self.target_to_phrase[trg] != [] and self.target_to_phrase[trg][0] == left_phrase:
                        self.target_to_phrase[trg][0] = right_phrase
                    elif self.target_to_phrase[trg] != [] and self.target_to_phrase[trg][0] > left_phrase:
                        self.target_to_phrase[trg][0] -= 1
                for source in range(self.source_size, self.source_size + self.artificial_sources):
                    self.phrases[source].segment_position = source

        #Otherwise, both segments are merged into the leftmost phrase.
        elif right_phrase < left_phrase:
            for src in self.phrases[left_phrase].sources:
                self.phrases[right_phrase].addSource(src)
            self.phrases[right_phrase].translation = self.phrases[left_phrase].translation + ' ' + self.phrases[right_phrase].translation
            self.phrases[left_phrase].segment_position = self.phrases[right_phrase].segment_position
            self.phrases[left_phrase].translation = ''
            for trg in range(self.target_size): #target_to_phrase is updated so that all targets pointing to left phrase point now to right phrase.
                if self.target_to_phrase[trg] != [] and self.target_to_phrase[trg][0] == left_phrase:
                    self.target_to_phrase[trg][0] = right_phrase

        else:
            for src in self.phrases[right_phrase].sources:
                self.phrases[left_phrase].addSource(src)
            self.phrases[left_phrase].addTranslation(self.phrases[right_phrase].translation)
            self.phrases[right_phrase].segment_position = self.phrases[left_phrase].segment_position
            self.phrases[right_phrase].translation = ''
            for trg in range(self.target_size): #target_to_phrase is updated so that all targets pointing to right phrase point now to left phrase.
                if self.target_to_phrase[trg] != [] and self.target_to_phrase[trg][0] == right_phrase:
                    self.target_to_phrase[trg][0] = left_phrase

        #Finally, segments are merged into the left segment.
        self.segments[left_segment] += self.segments[right_segment]
        self.segments_string[left_segment] += self.segments_string[right_segment]
        del self.segments[right_segment]
        del self.segments_string[right_segment]
        self.segments_size -= 1

    def getSources(self, word):
        """
        This method returns a list of sources aligned with a target word with a probability
        higher than the threshold. The method receives a string containing the target word.
        """

        probabilities = [-1 for n in range(self.source_size)]

        #Probabilities are computed only by non-validated phrases. Computations is made
        #through the HMM aligments. Both source word and target word are lowercase before
        #computing alignment probability.
        for n in range(self.source_size):
            if not self.phrases[n].validated:
                try:
                    self.alignments[self.source[n].lower()][word.lower()]
                    probabilities[n] = self.alignments[self.source[n].lower()][word.lower()]
                except:
                    continue

        sources = []

        for n in range(self.source_size):
            if probabilities[n] > self.alignment_probability:
                sources.append(n)

        return sources

    def wordCorrection(self, wrong_word, corrected_word):
        """
        This method corrects a wrong word. The method receives the index of the target word
        to correct and a string containing the corrected word.
        """

        #Corrections are made accordinly to the possible scenarios:
        if wrong_word == self.target_size: #The user is inserting a word at the end.
            #Target data structures are increased in one element at the end.
            trg = wrong_word
            self.target.append(corrected_word)
            self.target_to_phrase.append([])
            self.target_size += 1

        #elif self.target_to_phrase[wrong_word] == [] or not self.phrases[self.target_to_phrase[wrong_word][0]].validated: #The user is replacing a word.
        elif self.target_to_phrase[wrong_word] == [] or sum([1 for src in self.target_to_phrase[wrong_word] if not self.phrases[src].validated]) > 0: #The user is replacing a word.
            #(Nothing special to consider.)
            trg = wrong_word

        else: #The user is inserting a word between segments.
            #Target data structures are increased in one element at the position of the target word corrected.
            trg = wrong_word
            self.target = self.target[:trg] + [corrected_word] + self.target[trg:]
            self.target_to_phrase = self.target_to_phrase[:trg] + [[]] + self.target_to_phrase[trg:]
            for sgmnts in range(self.segments_size):
                for w in range(len(self.segments[sgmnts])):
                    if self.segments[sgmnts][w] >= trg:
                       self.segments[sgmnts][w]  += 1
            self.target_size += 1

        #The validated phrase into which the new segment (composed of the  word correction) belongs is
        #conformed by all the source words aligned with the segment.
        srcs = self.getSources(corrected_word)
        src = None
        for source in srcs: #The phrase of one of those words will represent the segment.
            if not self.phrases[source].validated:
                src = source
                break

        #If no source is aligned with the segment, an artificial source is created.
        if src == None:
            self.phrases.append(Phrase(-1))
            src = self.source_size + self. artificial_sources
            self.artificial_sources += 1

        #Otherwise, phrases corresponding to the sources aligned with the new word are "merged" into the phrase
        #that represents the segment.
        else:
            for source in srcs:
                if source != src:
                    self.phrases[source].validated = True
                    self.phrases[source].translation = ''
                    self.phrases[source].segment_position = src

        #The segment is stored into the phrase which represents it.
        self.phrases[src].validated = True
        self.phrases[src].translation = corrected_word
        self.phrases[src].segment_position = src
        self.target_to_phrase[trg] = [src]

        #The segment is added into the segment list:
        if self.segments_size == 0 or trg > self.segments[-1][-1]: #At the end of the list if the list is empty or the word belongs there.
            self.segments.append([trg])
            self.segments_string.append([corrected_word])
            self.segments_size += 1

        else:#At the corresponding position, otherwise.
            for index in range(self.segments_size + 1):
                if index == self.segments_size:
                    self.segments.insert(index, [trg])
                    self.segments_string.insert(index, [corrected_word])
                    self.segments_size += 1
                    break
                if self.segments[index][0] > trg:
                    self.segments.insert(index, [trg])
                    self.segments_string.insert(index, [corrected_word])
                    self.segments_size += 1
                    break
        #Finally, word strokes and mouse actions are computed (1 stroke and 1 action per correction).
        self.current_word_strokes += 1
        self.current_mouse_actions += 1

    def xmlTag(self, target, source):
        """
        This method generates an xml tag. The method receives a string containing the target
        part of the tag, and a string containing the source part of the tag.
        """
        if target == '':
            return '<x translation=" ">' + source + '</x><wall/>'

        return '<x translation="' + target + '">' + source + '</x><wall/>'

    def generateXML(self):
        """
        This method generates the xml source sentence for obtaining the new hypothesis.
        """

        self.xml = ''

        #No matter the current phrase, segments must be included in the xml in the desired order.
        #(Moses doesn't reorder translations belonging to a tag.) Sourcer must be added to the xml
        #in it's original order. Therefore, they must be split out of the phrase if there is some
        #more word in between (in that case, a new tag with an "empty" translation would be created).
        #No matter into which phrase they belong, phrases with an "empty" translation must be grouped
        #into as few tags as possible.

        segment_index = 0
        pending_phrases = []
        pending_phrases_size = 0
        segment_translation = ''
        empty_translation = True
        phrase_to_avoid = -1

        #Take into account the special merge case in which the first segment must be located at the
        #beggining of the new hypothesis.
        if self.artificial_source_at_beggining and not self.phrases[0].validated:
            segment_index = 1
            phrase_to_avoid = self.target_to_phrase[self.segments[0][0]][0]
            self.xml += self.xmlTag(' '.join(self.segments_string[0]), '.')

        #Phrases are visited in the source order.
        for n in range(self.source_size):

            '''if segment_index >= self.segments_size:
                print segment_index, self.segments_size, self.segments
                break'''

            #If the maximum number of sources is reached, the tag is created (no matter its content).
            if pending_phrases_size == self.max_sources:
                self.xml += self.xmlTag(segment_translation, ' '.join(pending_phrases)) + ' '
                pending_phrases = []
                pending_phrases_size = 0
                segment_translation = ''
                empty_translation = True

            #If the current phrase is a non-validated phrase:
            if not self.phrases[n].validated:
                if pending_phrases_size > 0: #Add to a tag phrases that are pending.
                    self.xml += self.xmlTag(segment_translation, ' '.join(pending_phrases)) + ' '
                #Add the sources belonging to the phrase into the xml.
                self.xml += ' '.join([self.source[w] for w in self.phrases[n].sources if not self.phrases[w].validated]) + ' '
                #Reset the list of pending phrases and segment translation.
                pending_phrases = []
                pending_phrases_size = 0
                segment_translation = ''
                empty_translation = True

            #If the current phrases is validated and represents a segment:
            elif self.phrases[n].segment_position == n and n != phrase_to_avoid:
                if not empty_translation: #Add to a tag phrases that are pending.
                    self.xml += self.xmlTag(segment_translation, ' '.join(pending_phrases)) + ' '
                    pending_phrases = []
                    pending_phrases_size = 0
                #Add current phrase to the list of pending phrases and current segment to the segment translation.
                pending_phrases.append(self.source[n])
                pending_phrases_size += 1
                segment_translation = ' '.join(self.segments_string[segment_index])
                empty_translation = False
                segment_index += 1

            #If the current phrase is a validated phrase included in other phrase:
            else:
                #Add the phrase to the list of pending phrases.
                pending_phrases.append(self.source[n])
                pending_phrases_size += 1

        if pending_phrases_size > 0: #Add to a tag phrases that are pending.
            self.xml += self.xmlTag(segment_translation, ' '.join(pending_phrases)) + ' '

        #Finally, add all phrases with artificial sources into xml tags.
        for n in range(self.source_size, self.source_size + self.artificial_sources):
            if n != phrase_to_avoid and segment_index < self.segments_size:
                self.xml += self.xmlTag(' '.join(self.segments_string[segment_index]), '.')
                segment_index += 1

    def generateXMLPB(self):
        """
        This method generates the xml source sentence for obtaining the new hypothesis working in a
        prefix-based environment. (Prefix-based is a particular case of segment-based--where there's
        only one validated segment--with the additional restriction that the segment always has to
        be at the beggining of the hypothesis.)
        """

        self.xml = ''

        #No matter the current phrase, segments must be included in the xml in the desired order.
        #(Moses doesn't reorder translations belonging to a tag.) Sourcer must be added to the xml
        #in it's original order. Therefore, they must be split out of the phrase if there is some
        #more word in between (in that case, a new tag with an "empty" translation would be created).
        #No matter into which phrase they belong, phrases with an "empty" translation must be grouped
        #into as few tags as possible.

        if not self.phrases[0].validated:
            self.xml += self.xmlTag(' '.join(self.segments_string[0]), '.')
            segment_translation = ''
        else:
            segment_translation = ' '.join(self.segments_string[0])

        pending_phrases = []
        pending_phrases_size = 0
        empty_translation = True

        #Phrases are visited in the source order.
        for n in range(self.source_size):

            #If the maximum number of sources is reached, the tag is created (no matter its content).
            if pending_phrases_size == self.max_sources:
                self.xml += self.xmlTag(segment_translation, ' '.join(pending_phrases)) + ' '
                pending_phrases = []
                pending_phrases_size = 0
                segment_translation = ''
                empty_translation = True

            #If the current phrase is a non-validated phrase:
            if not self.phrases[n].validated:
                if pending_phrases_size > 0: #Add to a tag phrases that are pending.
                    self.xml += self.xmlTag(segment_translation, ' '.join(pending_phrases)) + ' '
                #Add the sources belonging to the phrase into the xml.
                self.xml += ' '.join([self.source[w] for w in self.phrases[n].sources if not self.phrases[w].validated]) + ' '
                #Reset the list of pending phrases and segment translation.
                pending_phrases = []
                pending_phrases_size = 0
                segment_translation = ''
                empty_translation = True

            #If the current phrase is a validated phrase:
            else:
                #Add the phrase to the list of pending phrases.
                pending_phrases.append(self.source[n])
                pending_phrases_size += 1

        if pending_phrases_size > 0: #Add to a tag phrases that are pending.
            self.xml += self.xmlTag(segment_translation, ' '.join(pending_phrases)) + ' '

        self.artificial_source_at_beggining = True

    def endOfSentence(self):
        """
        This method accounts for the user indicating that the current validated segments
        conform their desired translation.
        """

        #Non-validated phrases become validated. (With an "empty" translation.)
        for src in range(self.source_size):
            if not self.phrases[src].validated:
                self.phrases[src].validated = True
                self.phrases[src].translation = ''

        #Word stroke and mouse actions are computed (1 stroke and 1 action--as if it were
        #a correction).
        self.current_word_strokes += 1
        self.current_mouse_actions += 1

    def validateTranslation(self):
        """
        This method validates the current hypothesis as the desired translation.
        """

        #Mouse actions are increased (by 1).
        self.current_mouse_actions += 1

        #Total word strokes, mouse actions, number of characters and number of words
        #are increased.
        self.word_strokes += self.current_word_strokes
        self.mouse_actions += self.current_mouse_actions
        self.characters += sum([len(w) for w in ' '.join(self.target).decode('utf-8').split()])
        self.words += self.target_size

    def mostUnlikelyTarget(self):
        """
        This method returns the word which the system has the least confidence on having
        translated correctly. For simplicity's sake, only the words next to segments
        (either before or after a segment) are considered. The first and last words
        from the hypothesis are also considered. If all words are validated, the method
        returns a 'None' value.
        """

        #Initialize probabilities to 1.
        probabilities = [9000 for trg in self.target]

        #Look for words next to segments.
        target_in_segments = [None for trg in self.target]
        for n in range(self.segments_size):
            for trg in self.segments[n]:
                target_in_segments[trg] = n

        #Deleted words must be ignored.
        target_deleted = [False for trg in self.target]
        for words in self.deleted_words:
            for word in words:
                target_deleted[word] = True

        for n in range(self.target_size):
            #Ignore deleted words.
            if target_deleted[n]:
                continue

            #Beginning and ending of the hypothesis.
            if (n == 0 or n == self.target_size - 1) and target_in_segments[n] == None:
                p = 0
                #for src in self.target_to_phrase[n]:
                for src in range(self.source_size):
                    try:
                        if self.alignments[self.source[src].lower()][self.target[n].lower()] > p:
                            p = self.alignments[self.source[src].lower()][self.target[n].lower()]
                    except:
                        continue
                probabilities[n] = p

            #Word before a segment.
            #if n > 0 and target_in_segments[n] and not target_in_segments[n - 1]:
            elif target_in_segments[n] == None and target_in_segments[n + 1] != None and target_in_segments[n - 1] != target_in_segments[n + 1]:
                p = 0
                #for src in self.target_to_phrase[n]:
                for src in range(self.source_size):
                    try:
                        if self.alignments[self.source[src].lower()][self.target[n].lower()] > p:
                            p = self.alignments[self.source[src].lower()][self.target[n].lower()]
                    except:
                        continue
                probabilities[n] = p

            #Word after a segment.
            #if n < self.target_size - 1 and target_in_segments[n] and not target_in_segments[n + 1]:
            elif target_in_segments[n] == None and target_in_segments[n - 1] != None and target_in_segments[n - 1] != target_in_segments[n + 1]:
                p = 0
                #for src in self.target_to_phrase[n]:
                for src in range(self.source_size):
                    try:
                        if self.alignments[self.source[src].lower()][self.target[n].lower()] > p:
                            p = self.alignments[self.source[src].lower()][self.target[n].lower()]
                    except:
                        continue
                probabilities[n] = p

        #Find the word with lowest probability.
        p = 9000
        trg = None
        for n in range(self.target_size):
            if probabilities[n] != 9000 and probabilities[n] < p:
                p = probabilities[n]
                trg = n
        new_targets = []
        for n in range(self.target_size):
            if probabilities[n] != 9000:
                new_targets.append(n)
        if new_targets == []:
            return None
        shuffle(new_targets)
        return new_targets[0]
        return trg

    def mostLikelySegments(self):
        """
        This method returns a list with those segments which the system considers that have
        a translation confidence greater than a threshold.
        """

        #Compute confidence of each word of the hypothesis.
        confidence = [-1 for trg in self.target]
        target_in_segments = [None for trg in self.target]
        for n in range(self.segments_size):
            for trg in self.segments[n]:
                target_in_segments[trg] = n
        for n in range(self.target_size):
            #for src in self.target_to_phrase[n]:
            if target_in_segments[n] == None:
                for src in range(self.source_size):
                    #if src < len(self.phrases) and not self.phrases[src].validated:
                    try:
                        self.alignments[self.source[src].lower()][self.target[n].lower()]
                        confidence[n] = self.alignments[self.source[src].lower()][self.target[n].lower()]
                    except:
                        continue

        #Divide the words with a confidence greater than the threshold into segments.
        current_segment = []
        segments = []
        for n in range(self.target_size):
            if confidence[n] > self.CM_threshold:
                if current_segment == [] or n == current_segment[-1] + 1:
                    current_segment.append(n)
                else:
                    segments.append(current_segment)
                    current_segment = [n]
        if current_segment != []:
            segments.append(current_segment)

        #Return the list of segments.
        return segments



########################################################################################
###############################   END SBIMT   ##########################################
########################################################################################
