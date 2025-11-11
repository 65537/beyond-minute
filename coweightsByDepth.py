from sharedFunctions import *

parser = argumentParser()
parser.add_argument("Depth", type=float, help="Maximum Depth allowed")
args = parser.parse_args()
initializeFromArgparse(args)
from sharedFunctions import *

queue = [coweight_lattice.zero()]
fundamental_weight_coweight_pairings = cartan_matrix.inverse().columns()
depth_pairings = [sum(fundamental_weight_coweight_pairings[i-1] for i in orbit) for orbit in sigma_orbits]


while queue:
	mu = queue.pop()
	coefficients = mu.dense_coefficient_list()
	depth = max(sum(c * p for (c,p) in zip(coefficients, pairings)) for pairings in depth_pairings)
	if depth >= args.Depth:
		continue
	print("%s has depth %s" % (mu,depth))
	for fundamental_coweight in coweight_lattice.basis():
		new_mu = mu + fundamental_coweight
		if new_mu not in queue:
			queue = [new_mu] + queue


