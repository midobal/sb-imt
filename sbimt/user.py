# -*- coding: utf-8 -*-


class User:
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
        if self.reference_in_segment[0] is not None and segments[0][0] != 0:
            session.mergeSegments(-1, 0)
            return True

        # Word by word, check if a word and its previous word belong to two
        # different segments.
        for n in range(1, self.reference_size):

            if self.reference_in_segment[n] is None:
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
