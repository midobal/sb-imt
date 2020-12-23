# -*- coding: utf-8 -*-
import bisect


class Phrase:
    """
    This class stores phrases information (mainly, sources and translations).
    It also contains a flag to know if a phrase is contained in a segment
    (validated) and a pointer to the phrase in which the segment is stored.
    Phrases data structures are ordered according to source.
    """

    def __init__(self, src):
        """
        This method initializes a new Phrase. The method receives the source
        to which is originally associated.
        """
        self.sources = [src]
        self.translation = ''
        self.segment_position = None
        self.validated = False

    def addTranslation(self, trans):
        """
        This method ads more translation to the current translation of the
        phrase. The method receives a string with the new translation.
        """
        if self.translation != '':
            self.translation += ' '
        self.translation += trans

    def addSource(self, src):
        """
        This method ads a source to the list of sources contained in the
        phrase. The method receives the position of the new source.
        """
        if src not in self.sources:
            bisect.insort(self.sources, src)
