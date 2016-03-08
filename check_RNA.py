import sys 


class RNASequence(object):

	alphabet = set('GAUC')
	
	def __init__(self, identifier, sequence):

		self.identifier = identifier
		self.sequence = sequence

		for letter in self.sequence:
			if letter not in self.alphabet:
				e = IncorrectSequenceLetter(letter, self.__class__.__name__)
				raise e
	def get_seq(self):

		return dict(identifier => sequence)

class IncorrectSequenceLetter(ValueError):
	
	def __init__(self, letter, class_name):

		self.letter = letter
		self.class_name = class_name
	def __str__(self):
		return "The sequence item %s is not found in the alphabet of class %s" %(self.letter, self.class_name)