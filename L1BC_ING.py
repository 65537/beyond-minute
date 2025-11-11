from sharedFunctions import *

parser = argumentParser()
parser.add_argument("mu", help="The cocharacter defining an admissible set. Pass a tuple to study Weil restrictions of scalars")
parser.add_argument("--level", default="iwahori", help="The parahoric level used for the admissible locus. Can be 'hyperspecial', 'iwahori' or a list of affine indices")
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


Adm_covers = {}
for xs_covers in itertools.product(*[enumerateAdmissibleLocus(m, index_set, True).items() for m in mu]):
	xs = tuple(t[0] for t in xs_covers)
	covers = []
	Adm_covers[xs] = covers
	for i in range(len(xs)):
		covers += [tuple(xs[:i]) + (y,) + tuple(xs[i+1:]) for y in xs_covers[i][1]]
print(Adm_covers)

xs_products = {}
xs_nd = {}
for xs in Adm_covers:
	product = extended_affine_weyl_group.W0P().one()
	for x in reversed(xs):
		product *= x
	xs_products[xs] = product
	xs_nd[xs] = newtonPointAndDefect(product)

xs_newton2rho = {xs: 2*sum(newton.dense_coefficient_list()) for (xs, (newton,d)) in xs_nd.items()}
xs_fundamental = set([xs for (xs, n2r) in xs_newton2rho.items() if n2r == tupleLength(xs)])
print("Found %d fundamental elements" % len(xs_fundamental))
print(xs_fundamental)
newton_hnindec = set()
for (newton,_) in set(xs_nd.values()):
	vector = cartan_matrix.transpose().inverse() * sage.all.vector((sum(mu) - newton).dense_coefficient_list())
	print(vector)
	if all(sum(vector[i-1] for i in orbit)>0 for orbit in sigma_orbits):
		newton_hnindec.add(newton)
print("%d / %d sigma conjugacy classes are HN indec" % (len(newton_hnindec), len(set(xs_nd.values()))))

xs_hnindec = set(xs for (xs, (newton,d)) in xs_nd.items() if newton in newton_hnindec)
print("%d / %d admissible set elements are HN indec" % (len(xs_hnindec), len(xs_products)))

xs_hnindec_max = set(xs for xs in xs_hnindec if all(not xs in Adm_covers[ys] for ys in xs_hnindec))
print("%d maximal HN indec elements" % len(xs_hnindec_max))

xs_gnp = dict()

for xs in sorted(xs_products.keys(), key=tupleLength):
	best_nd = xs_nd[xs]
	best_n2r = 2*sum(best_nd[0].dense_coefficient_list())
	for ys in Adm_covers[xs]:
		yn, yd = xs_gnp[ys]
		y_n2r = 2*sum(yn.dense_coefficient_list())
		if best_n2r is None or best_n2r < y_n2r:
			best_n2r = y_n2r
			best_nd = (yn, yd)
	xs_gnp[xs] = best_nd

max_gnp = None
for xs in xs_hnindec_max:
	n,d = xs_gnp[xs]
	print("%s is Bruhat maximal in the HN-indec locus and has generic Newton point %s" % (xs, n))
	if max_gnp is None:
		max_gnp = n
	elif max_gnp != n:
		print("Surprising difference")
		if len(index_set) > 0:
			print("(Might not be a counterexample to (ING), make sure to check the actual order)")

l1bc_failures = set()
for xs in xs_fundamental:
	xs_n, xs_d = xs_nd[xs]
	for ys in Adm_covers[xs]:
		ys_n, ys_d = xs_gnp[ys]
		l = 2*sum((xs_n-ys_n).dense_coefficient_list()) + ys_d - xs_d
		assert(l>0)
		if l > 2:
			l1bc_failures.add(xs_n)
if not l1bc_failures:
	print("As far as this admissible set is concerned, the property (L1BC) is always satisfied.")
else:
	print("In this admissible set, we found the following %d failures of (L1BC):" % len(l1bc_failures))
	for n in l1bc_failures:
		print("%s (HN %sdecomposable)" % (n, "in" if n in newton_hnindec else ""))
