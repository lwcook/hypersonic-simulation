import numpy as np
import math
import pdb

class obj_function:
	"""A class for objective functions, mapping x to y"""
	def __init__(self, x = 1):
		
		self.x = x
		self.evaluate(x)

	def __call__(self,x):
		return self.evaluate(x)

	def evaluate(self,x):
		# The actual objective function

		x = np.array(x)

		if x.size == 1:
			a = 2 + (9.0/20.0)*x
			b = 1
		elif x.size == 2:
			a = 2 + (9.0/20.0)*x[0]
			b = x[1]


		y = (a + b)**a

		self.y = y
		return y

	eval = evaluate



if __name__ == "__main__":

	fobj = obj_function()

	print fobj.eval(2)



