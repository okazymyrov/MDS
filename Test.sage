#  Created on: June 12, 2012
#		Author: Oleksandr Kazymyrov <oleksandr.kazymyrov@uib.no>

load Data.sage 

def test_g():
	m = MDS(generator=AES['g'],modulus=AES['modulus'])
	#m = MDS(generator=Kalyna['g'],modulus=Kalyna['modulus'])

	t1=cputime()

	m.convert()

	m.print_G()
	m.print_g()
	m.print_L()
	m.print_system()

	print "Selftesting	= {0}".format(m.selftesting(g_test=True))
	print "is_MDS		= {0}".format(m.is_MDS())

	t2=cputime()

	print "Time = {0}".format(t2-t1)

def test_L():
	m = MDS(L=AES['L'],modulus=AES['modulus'])
	#m = MDS(L=Kalyna['L'],modulus=Kalyna['modulus'])
	#m = MDS(L=Stribog['L'],modulus=Stribog['modulus'])

	t1=cputime()

	m.convert()

	m.print_G()
	m.print_g()
	m.print_L()
	m.print_system()

	print "Selftesting	= {0}".format(m.selftesting())
	print "is_MDS		= {0}".format(m.is_MDS())

	t2=cputime()

	print "Time = {0}".format(t2-t1)

def test_G():
	m = MDS(G=AES['G'],modulus=AES['modulus'])
	#m = MDS(G=Kalyna['G'],modulus=Kalyna['modulus'])
	#m = MDS(G=Stribog['G'],modulus=Stribog['modulus'])

	t1=cputime()

	m.convert()

	m.print_G()
	m.print_g()
	m.print_L()
	m.print_system()

	print "Selftesting	= {0}".format(m.selftesting())
	print "is_MDS		= {0}".format(m.is_MDS())

	t2=cputime()

	print "Time = {0}".format(t2-t1)

def test_L_Stribog():
	pol=["x^8+x^4+x^3+x^2+1","x^8+x^6+x^5+x^4+x^2+x^1+1","x^8+x^7+x^6+x^5+x^4+x^1+1","x^8+x^6+x^5+x^3+1","x^8+x^7+x^5+x^4+x^3+x^2+1",
		 "x^8+x^7+x^6+x^5+x^2+x^1+1","x^8+x^5+x^3+x^1+1","x^8+x^7+x^6+x^4+x^2+x^1+1","x^8+x^6+x^5+x^2+1","x^8+x^7+x^3+x^1+1",
		 "x^8+x^6+x^5+x^1+1","x^8+x^4+x^3+x^1+1","x^8+x^5+x^4+x^3+x^2+x^1+1","x^8+x^7+x^3+x^2+1","x^8+x^5+x^3+x^2+1",
		 "x^8+x^6+x^4+x^3+x^2+x^1+1","x^8+x^7+x^6+x^5+x^4+x^3+1","x^8+x^7+x^6+x^1+1","x^8+x^5+x^4+x^3+1","x^8+x^7+x^5+x^3+1",
		 "x^8+x^7+x^2+x^1+1","x^8+x^7+x^5+x^4+1","x^8+x^6+x^3+x^2+1","x^8+x^7+x^6+x^3+x^2+x^1+1","x^8+x^7+x^6+x^4+x^3+x^2+1",
		 "x^8+x^7+x^5+x^1+1","x^8+x^7+x^6+x^5+x^4+x^2+1","x^8+x^7+x^4+x^3+x^2+x^1+1","x^8+x^6+x^5+x^4+x^3+x^1+1","x^8+x^6+x^5+x^4+1"
		]

	t1=cputime()

	for i,p in enumerate(pol):
		m = MDS(L=Stribog['L'],modulus=p)
		m.convert()
		if m.get_G() is not None:
			print "polynomial = {0}".format(p)
			break

	t2=cputime()

	print "Time = {0}".format(t2-t1)
	