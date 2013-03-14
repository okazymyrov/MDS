#  Created on: June 12, 2012
#		Author: Oleksandr Kazymyrov <oleksandr.kazymyrov@uib.no>

class MDS(SageObject):
	def __init__(self, **kwargs):
		r'''
			Parameters
				generator	- a generator polynomial of MDS code
				modulus		- a primitive polynomial over GF(2)
				L			- linear representation of MDS matrix by matrix in GF(2)
				G			- MDS matrix
		'''
		self._g			= kwargs.get('generator',None)
		self._G			= kwargs.get('G',None)
		self._L			= kwargs.get('L',None)
		self._f			= kwargs.get('modulus',None)
		
		if self._f is None:
			raise TypeError("'modulus' is not defined")

		if self._G is None and self._L is None and self._g is None:
			raise TypeError("None of 'G', 'L' or 'generator' are not defined")

		R.<x>=ZZ[]
		self._n 		= R(self._f).degree()
		self._length 	= 1 << self._n		
		self._K 		= GF(1 << R(self._f).degree(),'a',modulus=R(self._f))
		self._P 		= PolynomialRing(self._K,'x')

		if self._g is not None:
			self._nbr 		= R(self._g).degree() + 1 # number of bytes in row or column
			self._nbc 		= self._nbr*self._n		  # number of bits in row or column
		
		if self._L is not None:
			self._nbc = self._L.nrows()
			self._nbr = self._nbc/self._n

		if self._G is not None:
			self._nbr = self._G.nrows()
			self._nbc = self._nbr*self._n

	def __convert_g_P(self):
		r'''
			Apply PolinomialRing to generator
		'''

		R.<x> = ZZ[]
		self._g 		= R(self._g).coeffs()
		self._g 		= [self._K.fetch_int(i) for i in self._g]
	
	def __convert_g_G(self):
		r'''
			Convert generator to circulant matrix
		'''

		self._G = matrix(self._K,self._nbr)

		for i in xrange(self._nbr):
			for j in xrange(self._nbr):
				self._G[j,i] = self._g[(j-i)%self._nbr]

	def __convert_g_L(self):
		r'''
			Convert generator to the linear representation of MDS code

			B = self._L * A
		'''
		self._L = matrix(GF(2),self._nbc)

		if self._G is None:
			self.__convert_g_G()

		self.__convert_G_L()

		for i in xrange(self._nbr):
			for j in xrange(self._nbr):
				self._L.set_block(i*self._n,j*self._n,self._G[i][j]._matrix_())	

	def __convert_G_g(self):
		r'''
			Convert circulant matrix to the generator, otherwise g is None.
		'''
		G = self._G[:]
		self._g = self._G.column(0).list()
		self.__convert_g_G()

		if G != self._G:
			self._g = None
			self._G = G

	def __convert_G_K(self):
		r'''
			Convert integer elements of G to the field K
		'''
		
		G = matrix(self._K,self._nbr)

		for i in xrange(self._nbr):
			for j in xrange(self._nbr):
				G[i,j] = self._K.fetch_int(self._G[i,j])
		
		self._G = G

	def __convert_G_L(self):
		r'''
			Convert MDS matrix to the linear representation
			
			B = self._L * A
		'''
		self._L = matrix(GF(2),self._nbc)

		for i in xrange(self._nbr):
			for j in xrange(self._nbr):
				self._L.set_block(i*self._n,j*self._n,self._G[i][j]._matrix_())	

	def __convert_L_G(self):
		r'''
			Convert the linear representation to the MDS matrix
		'''

		self._G = matrix(self._K,self._nbr)

		for i in xrange(self._nbr):
			for j in xrange(self._nbr):
				self._G[i,j] = self._K((self._L[i*self._n:(i+1)*self._n,j*self._n:(j+1)*self._n]).column(0))

	def __convert_L_system(self):
		r'''
			Convert the linear representation to the system of equations
		'''
		P = BooleanPolynomialRing(self._n+self._nbc, ["c%d"%i for i in range(self._n)] + ["d%d"%i for i in range(self._nbc)])
		gens = list(P.gens())
		
		self._system = range(self._nbc)

		for i in xrange(self._nbc):
			self._system[i] = sum([self._L[i][g]*gens[self._n + g] for g in xrange(self._nbc) ])

	def __convert_P_g(self):
		r'''
			Apply PolinomialRing to the generator
		'''
		self._g = P(self._g)

	def is_MDS(self):
		r'''
			Check MDS property of matrix G
		'''

		for m in xrange(1,self._nbr+1):
			if (0 in self._G.minors(m)) == True:
				#print "{0}".format(self._G.minors(m))
				return False

		return True

	def __MDS_mul(self,state):
		r'''
			Perform MDS transformation as in AES using self._g
		'''
		g = self._P(self._g)
		state = ZZ(state).digits(2^8,padto=self._nbr^2)

		for i in xrange(self._nbr):
			t=self._P([self._K.fetch_int(h) for h in state[i*self._nbr:i*self._nbr+self._nbr]])
			t=g*t
			t=t.mod(self._P("x^{0}+1".format(self._nbr)))
			for h in xrange(self._nbr):
				state[i*self._nbr+h]=t[h]

		state = [i.integer_representation() for i in state]
		state = ZZ(state,2^self._n)
		
		return state

	def __MDS_matrix(self,state):
		r'''
			Perform MDS transformation as in AES using self._G
		'''
		S = matrix(self._K,self._nbr,[self._K.fetch_int(i) for i in ZZ(state).digits(2^8,padto=self._nbr^2)])
		S = S.transpose()

		S = self._G*S

		S = S.transpose()
		state = ZZ([i.integer_representation() for i in S.list()],2^8)
		
		return state

	def __MDS_matrix_2(self,state):
		r'''
			Perform MDS transformation as in AES using self._L
		'''

		state = ZZ(state).digits(2^self._nbc,padto=self._nbr)
		
		for i in xrange(self._nbr):
			state[i] = vector(ZZ(state[i]).digits(2,padto=self._nbc))
			state[i] = self._L*state[i]
			state[i]=ZZ(state[i].list(),2)

		state = sum([state[i]*(2^self._nbc)^i for i in xrange(self._nbr)])

		return state

	def __MDS_system(self,state):
		r'''
			Perform MDS transformation as in AES using self._system
		'''

		P = BooleanPolynomialRing(self._n+self._nbc, ["c%d"%i for i in range(self._n)] + ["d%d"%i for i in range(self._nbc)])
		gens = list(P.gens())

		state = ZZ(state).digits(2^self._nbc,padto=self._nbr)
		for i in xrange(self._nbr):
			state[i] = state[i].digits(2,padto=self._nbc)
			state[i] = [self._system[g].subs(dict(zip(gens[self._n:self._n+self._nbc],state[i]))) for g in xrange(self._nbc)]
			state[i] = ZZ(state[i],2)

		state = sum([state[i]*(2^self._nbc)^i for i in xrange(self._nbr)])

		return state

	def convert(self):
		r'''
			Convert given L, G or generator(g) to other representations
		'''

		if self._g is not None:
			self.__convert_g_P()
			self.__convert_g_G()
			self.__convert_g_L()
			self.__convert_L_system()
			return

		if self._G is not None:
			self.__convert_G_K()
			self.__convert_G_g()
			self.__convert_G_L()
			self.__convert_L_system()
			return

		if self._L is not None:
			self.__convert_L_system()
			self.__convert_L_G()
			self.__convert_G_g()
			return

	def get_g(self):
		r'''
			Return generator self._g
		'''
		
		if self._g is None:
			self.convert()
		
		return self._g

	def get_G(self):
		r'''
			Return MDS matrix self._G
		'''
		if self._G is None:
			self.convert()
		
		return self._G

	def get_L(self):
		r'''
			Return linear matirx over GF(2) self._L
		'''
		if self._L is None:
			self.convert()
		
		return self._L

	def print_g(self, form='integer'):
		r'''
			Print self._g in numerical notation
		'''
		if self._g is None:
			print "g = None"
			return 
		
		s = "g = "
		if form == 'integer':
			for i in xrange(self._nbr-1,1,-1):
				t = self._g[i].integer_representation()
				if t >= 2:
					s += "{0}*x^{1}+".format(t,i)
				elif t == 1:
					s += "x^{0}+".format(i)

			t = self._g[1].integer_representation()
			if t >= 2:
				s += "{0}*x+".format(t,i)
			elif t == 1:
				s += "x+"
			t = self._g[0].integer_representation()
			if t >= 1:
				s += "{0}".format(t,i)
			print s

	def print_G(self):
		r'''
			Print self._G in numerical notation
		'''
		print "G:"
		for i in xrange(self._nbr):
			for j in xrange(self._nbr):
				print "%.2X"%self._G[i][j].integer_representation(),
			print ""

	def print_L(self):
		r'''
			Print self._L
		'''
		print "L:"
		for i in xrange(self._nbc):
			for j in xrange(self._nbc):
				print "%d"%self._L[i][j],
			print ""

	def print_system(self):
		r'''
			Print self._system
		'''
		print "system:"
		for i,s in enumerate(self._system):
			print "system[{0}] = {1}".format(i,s)

	def selftesting(self,**kwargs):
		r'''
			Check equality of self._G, self._L, self._system and/or self._g
			Parameters
				verbose	- print additional information 
				g_test	- check g on equality
		'''
		verbose = kwargs.get('verbose',False)
		g_test = kwargs.get('g_test',False)

		state = randint(0,2^(self._nbr*self._nbc)-1)

		s1 = self.__MDS_matrix(state)
		s2 = self.__MDS_matrix_2(state)
		s3 = self.__MDS_system(state)
		if g_test == True:
			s4 = self.__MDS_mul(state)

		if verbose == True:
			print "MDS		= {0:x}".format(s1)
			print "L		= {0:x}".format(s2)
			print "system		= {0:x}".format(s3)
			if g_test == True:
				print "mul		= {0:x}".format(s4)

		if g_test == True:
			if s1 != s2 or s1 != s3 or s1 != s4 or s2 != s3 or s2 != s4 or s3 != s4:
				return False
		else:
			if s1 != s2 or s1 != s3 or s2 != s3:
				return False

		return True
