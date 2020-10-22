
from surface_dynamics import AbelianStratum
from sage.functions.other import sqrt
from sage.arith.misc import xgcd, gcd
from sage.rings.infinity import Infinity
from sage.modular.arithgroup.arithgroup_perm import ArithmeticSubgroup, ArithmeticSubgroup_Permutation
from multiprocessing import Pool as ThreadPool




_SL2Z = ArithmeticSubgroup_Permutation()
_SL2Z_Farey = _SL2Z.farey_symbol()




def max_systole_h11( n, lower_bound=0, n_threads=2 ):
    r"""
    This function computes the maximal length of a shortest systole on Origamis
    in the stratum H(1,1) with exactly `n` squares.

    INPUT:

    - ``n`` -- a positive integer. The number of squares of the Origamis in
      H(1,1) to consider.

    - ``lower_bound`` -- a positive real number (default: `0`). A known lower
      bound for the maximal length of a shortest systole on Origamis in H(1,1)
      with `n` squares. If the algorithm encounters an Origami with a systole
      of smaller length than `lower_bound`, it discards the Origami immediately
      without determining the exact systolic ratio.

    - ``n_threads`` -- a positive integer (default: `2`). The number of threads
      that should be used for the computation. The algorithm decomposes the set
      of Origamis with `n` squares into a disjoint union of orbits under the
      action of SL(2,ZZ) on Origamis and deals with each orbit in a separate
      thread.


    OUTPUT: A tuple `(o,m)` consisting of the maximal length `m` of a shortest
    systole and an Origami `o` that achieves this bound.

    """

    component = AbelianStratum(1,1).hyperelliptic_component()
    curves = component.arithmetic_teichmueller_curves(n)

    pool = ThreadPool(n_threads)
    curves_max = pool.map( _max_systole_h11_packed, [(c,lower_bound) for c in curves])
    pool.close()
    pool.join()

    o, m = max( curves_max, key=lambda l: l[1][0] )
    return o, m




def _max_systole_h11_packed( x ):
    r""" 
    This function is a facade to the actual algorithm implemented in
    `shortest_systles_on_h11_orbit` that takes a single tuple as argument and
    computes the maximum length of a systole.

    INPUT:

    - ``x`` -- a tuple `(c,lower_bound)` of an arithmetic Teichmüller
      curve/SL(2,ZZ)-orbit of Origamis in H(1,1) and a real number
      `lower_bound` that is a known lower bound for the maximal length of a
      systole on Origamis in `c`.

    OUTPUT: A tuple `(o,m)` consisting of the maximal systole length `m` and an
    Origami `o` in `c` that achieves this bound.

    """

    curve, lower_bound = x 
    return max(shortest_systoles_on_h11_orbit( curve, lower_bound ).iteritems(), key=lambda l: l[1][0])




def shortest_systoles_on_h11_orbit( curve, lower_bound ):
    # R = S * L^-1
    origami = curve.origami()
    l_action, _, s_action = origami.sl2z_edges()
    r_action = dict((lo, s_action[o]) for o, lo in l_action.iteritems())
    action = dict((o,dict()) for o in s_action.iterkeys())
    for o, so in s_action.iteritems():
        ro = r_action[o]
        action[o][ 1] = so
        action[so][-1] = o
        action[o][ 2] = ro
        action[ro][-2] = o
    horizontal_saddles = dict((o, _shortest_horizontal_saddles(o)) for o in action.iterkeys())

    # we compute systoles only up to action of S as S preserves all lengths.
    i, s_orbits = 0, s_action.keys()
    while i < len(s_orbits):
        o = s_orbits[i]
        so = set( [o, s_action[o], s_action[s_action[o]], action[o][-1] ] )
        so.remove(o)
        for o in so:
            s_orbits.remove(o)
        i += 1

    systoles = dict()
    for o in s_orbits:
        l, e0, e1 = horizontal_saddles[o]
        best_edges = [ (l, (1,0)), (e0, (1,0)), (e1, (1,0)) ]

        # treat the case (x,y) = (0,1) separately.
        _update(best_edges, horizontal_saddles[s_action[o]], (0,1))

        mincycle = min(best_edges[0][0], best_edges[1][0]+best_edges[2][0])

        # TODO: iterate over coprime pairs directly and optimize for not doing work twice! Need to solve word problem only half of the time with
        # A S (y,-x) = (1,0) for A (x,y) = (1,0)

        y = 1
        while y < mincycle and mincycle >= lower_bound:
            x = 1
            while x**2 + y**2 < mincycle**2 and mincycle >= lower_bound:
                d, a, b = xgcd(x, y)
                if d == 1:
                    dlen = sqrt( x**2 + y**2 )

                    xy_o = o
                    for p in reversed(_SL2Z_Farey.word_problem( _SL2Z( [[a,b], [-y,x]]))):
                        xy_o = action[xy_o][p]
                    _update(best_edges,\
                        (dlen * l for l in horizontal_saddles[xy_o]),\
                        (x,y))

                    xy_o = o
                    for p in reversed(_SL2Z_Farey.word_problem( _SL2Z( [[-a,b], [-y,-x]]))):
                        xy_o = action[xy_o][p]
                    _update(best_edges,\
                        (dlen * l for l in horizontal_saddles[xy_o]),\
                        (x,y))

                    mincycle = min(best_edges[0][0], best_edges[1][0]+best_edges[2][0])
                x += 1
            y+=1

        if mincycle == best_edges[0][0]:
            systoles[o] = (mincycle, 'loop', best_edges[0][1] )
        else:
            systoles[o] = (mincycle, 'cycle', best_edges[1], best_edges[2] )

    return systoles




def _shortest_horizontal_saddles( origami ):
    r"""
    This helper function computes the shortest horizontal saddle connections on
    an Origami with exactly two singularities.

    INPUT:

    - ``origami`` -- an Origami with exactly two singularities.

    OUTPUT: A tuple `(l,e0,e1)`, where `l` is the length of the shortest
    horizontal loop on the given origami and `e0` and `e1` are the lengths of
    the two shortest horizontal edges between the two singularities.
    """

    # The singularities correspond to cycles in the cycle decomposition of the
    # commutator of the "right" and "up" permutations of the origami.
    r, u = origami.r(), origami.u()
    comm = r * u * r.inverse() * u.inverse()
    singularities = comm.cycle_tuples()

    # This code snippet works only for Origamis with two singularities.
    assert len(singularities) == 2

    # If (i1 i2 ... ir) is such a cycle then the squares that contain the
    # corresponding singularity at their lower left corner are i1, ..., ir.
    singularity_squares = [ i for c in singularities for i in c  ]

    # Mapping from singularity_squares to the singularity (either 0 or 1) that
    # they contain.
    square_to_singularity = dict(\
            (s,0) if s in singularities[0] else (s,1)\
            for s in singularity_squares\
            )

    # We compute the shortest horizontal saddle connection by going to the
    # right from each square that contains one of the two singularities until
    # we hit another singularity.
    #
    # As there are only two singularities, we keep track of the shortest loop
    # and the two shortest edges between the two singularities.
    l, e0, e1 = Infinity, Infinity, Infinity
    for square in singularity_squares:
        cur_len = 1
        cur_square = r(square)
        while not cur_square in singularity_squares:
            cur_len, cur_square = cur_len+1, r(cur_square)

        s1, s2 = square_to_singularity[square], square_to_singularity[cur_square]
        if s1 == s2:
            l = min(l, cur_len)
        elif cur_len <= e0:
            e0, e1 = cur_len, e0
        elif cur_len < e1:
            e1 = cur_len
    return (l, e0, e1)




# code snippet for update of short edges and loops.
def _update(old, new, direction):
    # pass old as list as [ (loop length, loop direction), ... ]
    l, e0, e1 = old
    new_l, new_e0, new_e1 = new

    if new_l < old[0][0]:
        old[0] = (new_l, direction)
    if new_e0 <= old[1][0]:
        if new_e1 < old[1][0]:
            old[2] = (new_e1, direction)
        else:
            old[2] = old[1]
        old[1] = (new_e0, direction)
    elif new_e0 < old[2][0]:
        old[2] = (new_e0, direction)


# vim:ft=python
