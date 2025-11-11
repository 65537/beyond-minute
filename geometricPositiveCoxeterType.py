from sharedFunctions import *

parser = argumentParser()
parser.add_argument("mu", help="The cocharacter defining an admissible set. Pass a tuple to study Weil restrictions of scalars")
parser.add_argument("--level", default="hyperspecial", help="The parahoric level used for the admissible locus. Can be 'hyperspecial' (default), 'iwahori' or a list of affine indices")
args = parser.parse_args()
initializeFromArgparse(args)
from sharedFunctions import *

omega = coweight_lattice.basis()
tworho_vee = 2*sum(omega)
alpha = coweight_lattice.alpha()
mu = eval("("+args.mu+")")
if isinstance(mu, tuple):
	d = len(mu)
	mu = tuple(coweight_lattice(m) for m in mu)
else:
	d=1
	mu = (coweight_lattice(mu),)
index_set = None
if args.level == 'hyperspecial':
	index_set = root_system.index_set()
elif args.level.lower() == 'iwahori':
	index_set=()
else:
	index_set = tuple(sorted(eval('('+args.level+')')))

@functools.lru_cache(maxsize=None)
def lengthPositive(x):
	translation = x.to_translation_right()
	classical = weyl_group(x.to_classical_weyl())
	classical_coweight = classical.inverse().action(tworho_vee)
	translation_dom,rw = root_system.coweight_space()(translation).to_dominant_chamber(reduced_word=True)
	translation_dom = translation_dom.dense_coefficient_list()
	vmin = weyl_group.from_reduced_word(rw).coset_representative(tuple(i for i in root_system.index_set() if translation_dom[i-1] == 0),side='right')
	worklist = [vmin]
	result = set(worklist)
	while worklist:
		v = worklist.pop()
		vd = root_system.weight_lattice().weyl_group()(v)
		for (i,alpha) in enumerate(root_system.root_lattice().simple_roots()):
			valpha = alpha.weyl_action(vd)
			lfun = translation.scalar(valpha) +int(tworho_vee.scalar(valpha)>0) - int(classical_coweight.scalar(valpha)>0)
			#print("%s, %s, %s has valpha=%s, lfun = %d" % (x, v, alpha, valpha, lfun))
			assert(lfun >= 0)
			if lfun > 0:
				continue
			v_new = v.apply_simple_reflection(i+1,side='right')
			if v_new in result:
				continue
			result.add(v_new)
			worklist.append(v_new)
	return result
@functools.lru_cache(maxsize=None)
def descentSets(x):
	return {
		'left' : [i for i in [0]+list(root_system.index_set()) if x.has_descent(i, side='left' )],
		'right': [i for i in [0]+list(root_system.index_set()) if x.has_descent(i, side='right')]
	}

@functools.lru_cache(maxsize=None)
def tupleSupport(xs):
	if isinstance(xs, tuple) or isinstance(xs, list):
		result = set()
		for x in xs:
			result = result.union(tupleSupport(x))
		return result
	return xs.support()


@functools.lru_cache(maxsize=None)
def geometricCoxeterType(xs):
	degree = len(xs)
	simple_reflections = extended_affine_weyl_group.W0P().simple_reflections()
	one = extended_affine_weyl_group.W0P().one()
	worklist = [xs]
	cyclic_shift_class = set(worklist)
	while worklist:
		xs = worklist.pop()
		descents = [descentSets(x) for x in xs]
		for i in range(len(xs)):
			for j in descents[i]['right']:
				conj = tuple( (simple_reflections[sigma.inverse()(j)] if (k+1-i) % degree == 0 else one) * xs[k] * (simple_reflections[j] if k == i else one) for k in range(degree))
				if conj in cyclic_shift_class:
					continue
				cyclic_shift_class.add(conj)
				if j in descents[(i-1) % degree]['left']:
					assert(tupleLength(conj) < tupleLength(xs))
					gct_2 = geometricCoxeterType(conj)
					if gct_2 is False:
						return False
					right_mul = tuple(xs[k] * (simple_reflections[j] if k == i else one) for k in range(degree))
					gct_1 = geometricCoxeterType(right_mul)
					if gct_1 is False:
						return False
					if any(b in gct_2 for b in gct_1):
						return False
					return gct_1.union(gct_2)
				worklist.append(conj)
	xlength = tupleLength(xs)
	x = extended_affine_weyl_group.W0P().one()
	for xel in reversed(xs):
		x *= xel
	newton,defect = newtonPointAndDefect(x)
	newton2rho = sum(newton.scalar(alpha) for alpha in root_system.root_space().positive_roots())
	clsigma_matrix = []
	x_cl = weyl_group(x.to_classical_weyl())
	for i in root_system.index_set():
		omega = coweight_lattice.basis()[i]
		sigma_omega = coweight_lattice.basis()[sigma(i)]
		clsigma_omega = x_cl.action(sigma_omega)
		clsigma_matrix.append( (clsigma_omega - omega).dense_coefficient_list())
	clsigma_matrix = sage.all.matrix(clsigma_matrix).transpose()
	sigma_reflection_length = clsigma_matrix.rank()
	xlength_gct_bound = newton2rho + sigma_reflection_length - defect
	assert(xlength >= xlength_gct_bound)
	if xlength == xlength_gct_bound:
		return set([newton])
	else:
		return False

total_size = 0
pct_size = 0
gct_size = 0
for xs in itertools.product(*[enumerateAdmissibleLocus(m, index_set) for m in mu]):
	print("%s has length %d" % (xs, tupleLength(xs)))
	total_size += 1
	for vs in itertools.product(*[lengthPositive(x) for x in xs]):
		conj_classical_part = tuple(vs[i-1].inverse() * xs[i].to_classical_weyl() * vs[i] for i in range(len(xs)))
		#print("%s -> %s" % (vs, conj_classical_part))
		support = tupleSupport(conj_classical_part)
		if tupleLength(conj_classical_part) == sum(int(any(i in support for i in orbit)) for orbit in sigma_orbits):
			print("positive Coxeter with v = %s" % repr(vs))
			pct_size += 1
			break
	if geometricCoxeterType(xs) is not False:
		print("This element is of geometric Coxeter type")
		gct_size += 1

print("Admissible set has %d elements. Positive Coxeter: %d, Geometric Coxeter: %d" % (total_size, pct_size, gct_size))

