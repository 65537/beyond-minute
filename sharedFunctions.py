import sage.all
import argparse
import sys
import functools
import itertools

def argumentParser():
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('CartanType', default='A2', help='Cartan type of the finite root system, such as A2, B3, F4 etc.')
	parser.add_argument('--frobenius', default='split',help="Frobenius action. Should be 'split', 'quasi-split' or a permutation in cycle notation. The word 'quasi-split' is only available if there is a unique non-trivial automorphism of the finite Dynkin diagram")
	return parser

def initializeFromArgparse(args=None):
	if args is None:
		parser = argumentParser()
		args = parser.parse_args()
	global dynkin_diagram, cartan_matrix, root_system, coweight_lattice, weyl_group, extended_affine_weyl_group, sigma, sigma_orbits
	dynkin_diagram = sage.all.DynkinDiagram(args.CartanType)
	root_system = sage.all.RootSystem(dynkin_diagram.cartan_type())
	cartan_matrix = root_system.cartan_matrix()
	coweight_lattice = root_system.coweight_lattice()
	weyl_group = coweight_lattice.weyl_group(prefix="s")
	extended_affine_weyl_group = sage.all.ExtendedAffineWeylGroup(cartan_matrix.cartan_type())
	sigma = args.frobenius
	if sigma is None or sigma in ['s', 'split']:
		sigma = sage.all.PermutationGroupElement(())
	elif sigma in ['qs', 'quasi-split']:
		automorphisms = dynkin_diagram.automorphism_group()
		assert(automorphisms.order() == 2)
		sigma = next(iter(s for s in automorphisms if not s.is_one()))
	else:
		sigma = sage.all.PermutationGroupElement(str(sigma))
	sigma_orbits = set(tuple(sorted(sigma.orbit(i))) for i in root_system.index_set())


def enumerateAdmissibleLocus(mu,min_length_set = (), record_bruhat_lower_covers = False):
	candidates = []
	worklist = [extended_affine_weyl_group.W0P().from_translation(w.action(mu)) for w in weyl_group]
	worklist = [c for c in worklist if not any(c.has_descent(d, side='right') for d in min_length_set)]
	found = set(worklist)
	worklist = list(found)
	fundamental_element = None
	output = {}
	while worklist:
		candidate = worklist.pop()
		if record_bruhat_lower_covers:
			c_covers = []
			output[candidate] = c_covers
		candidate = extended_affine_weyl_group.FW()(candidate)
		candidate_aff = candidate.to_affine_weyl_right()
		if fundamental_element is None:
			fundamental_element = candidate * candidate.parent().from_affine_weyl(candidate_aff).inverse()
		for cand2_aff in candidate_aff.bruhat_lower_covers():
			cand2 = extended_affine_weyl_group.W0P()(fundamental_element * candidate.parent().from_affine_weyl(cand2_aff))
			if any(cand2.has_descent(d, side='right') for d in min_length_set):
				continue
			if record_bruhat_lower_covers:
				c_covers.append(cand2)
			if cand2 in found:
				continue
			found.add(cand2)
			worklist.append(cand2)
	return output if record_bruhat_lower_covers else found
def newtonPoint(x):
	x_power = x
	sig_power = sigma
	exponent = 1
	while not x_power.to_classical_weyl().is_one() or not sig_power.is_one():
		sig_power *= sigma
		if sigma.is_one():
			x_power *= x
		else:
			sigma_x_power_classical = x_power.parent().from_reduced_word([sigma(s) for s in x_power.to_classical_weyl().reduced_word()])
			sigma_x_power_translation = sum(coef * coweight_lattice.basis()[sigma(i+1)] for (i, coef) in enumerate(x_power.to_translation_right().dense_coefficient_list()))
			x_power = x * sigma_x_power_classical * x_power.parent().from_translation(sigma_x_power_translation)
		exponent += 1
	return root_system.coweight_space()(x_power.to_translation_right()).to_dominant_chamber() / exponent
def newtonPointAndDefect(x):
	newton_point = newtonPoint(x)
	defect_vector = sage.all.vector((newton_point - newtonPoint(extended_affine_weyl_group.W0P().from_translation(x.to_translation_right()))).dense_coefficient_list())
	defect_vector = cartan_matrix.transpose().inverse() * defect_vector
	defect = sum(int(not sum(defect_vector[i-1] for i in orbit).is_integer()) for orbit in sigma_orbits)
	return (newton_point, defect)
@functools.lru_cache(maxsize=None)
def tupleLength(xs):
	if isinstance(xs, tuple) or isinstance(xs, list):
		return sum(tupleLength(x) for x in xs)
	return xs.length()


if __name__ == '__main__':
	initializeFromArgparse()
	for el in enumerateAdmissibleLocus(root_system.coweight_lattice().an_element()):
		print(el)
