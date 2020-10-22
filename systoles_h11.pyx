
from multiprocessing import Pool as ThreadPool
from sage.functions.other import sqrt
from sage.arith.misc import xgcd
from sage.rings.infinity import Infinity
from sage.modular.arithgroup.arithgroup_perm import ArithmeticSubgroup_Permutation
from surface_dynamics import AbelianStratum




_SL2Z = ArithmeticSubgroup_Permutation()
_SL2Z_Farey = _SL2Z.farey_symbol()




def max_systole_h11(n_squares, lower_bound=0, n_threads=2):
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

    component = AbelianStratum(1, 1).hyperelliptic_component()
    curves = component.arithmetic_teichmueller_curves(n_squares)

    pool = ThreadPool(n_threads)
    curves_max = pool.map(_max_systole_h11_packed, [(c, lower_bound) for c in curves])
    pool.close()
    pool.join()

    return max(curves_max, key=lambda l: l[1][0])




def _max_systole_h11_packed(data):
    r"""
    This function is a facade to the actual algorithm implemented in
    `shortest_systles_on_h11_orbit` that takes a single tuple as argument and
    computes the maximum length of a systole.

    INPUT:

    - ``data`` -- a tuple `(c,lower_bound)` of an arithmetic Teichmueller
      curve/SL(2,ZZ)-orbit of Origamis in H(1,1) and a real number
      `lower_bound` that is a known lower bound for the maximal length of a
      systole on Origamis in `c`.

    OUTPUT: A tuple `(o,m)` consisting of the maximal systole length `m` and an
    Origami `o` in `c` that achieves this bound.

    """

    curve, lower_bound = data
    return max(\
            shortest_systoles_on_h11_curve(curve, lower_bound).iteritems(),\
            key=lambda l: l[1][0]\
            )




def shortest_systoles_on_h11_curve(curve, lower_bound):
    r"""
    This function computes the maximal length of a shortest systoles on the
    Origamis in a given SL(2,ZZ) orbit in the stratum H(1,1).

    INPUT:

    - ``curve`` -- an arithmetic Teichmueller curve in H(1,1), i.e. a SL(2,ZZ)
      orbit of Origamis in H(1,1).
    - ``lower_bound`` -- a real number that is a known to be a lower bound for
      the maximal length of a systole on Origamis in `curve`.

    OUTPUT: A dictionary mapping each origami `o` in `curve` to a description
    `(l, t, d)` of the shortest (relevant) systole on `o`. If only a single
    systole of length less than `lower_bound` or the current maximal length of
    a shortest systole is found on `o`, then this systole is taken regardless
    of the fact that there may exist shorter systoles.

    The constituents `l`, `t` and `d` of the values of the dictionary have the
    following semantics: The value `l` is the length of the shortest relevant
    systole on `o` and `t` is either the string `'loop'` or `'cycle'` depending
    on whether the shortest systole is a loop connecting a singularity with
    itself or a cycle running through both singularites of `o`, respectively.

    In the case that the shortest systole is a loop, the final entry `d` is a
    vector indicating the direction of the loop.

    If the systole is a cycle, then `d` is a tuple `(edge_0, edge_1)` of
    tuples, where the tuples `edge_0` and `edge_1` contain length and direction
    of the edges in the cycle that make up the shortest systole.

    """


    origami = curve.origami()
    _, _, s_action = origami.sl2z_edges()

    # Precompute information on shortest horizontal saddle connections for all
    # origamis in this orbit.
    horizontal_saddles = dict((o, _shortest_horizontal_saddles(o)) for o in s_action.iterkeys())

    # We compute systoles only up to action of the matrix S as all lengths are
    # preserved by this action. To this end, we compute representatives of the
    # orbits under S.
    i, s_orbit_reps = 0, s_action.keys()
    while i < len(s_orbit_reps):
        o = s_orbit_reps[i]
        s_o = s_action[o]
        s_orbit_reps.remove(s_o)
        s_o = s_action[s_o]
        if s_o != o:
            s_orbit_reps.remove(s_o)
            s_o = s_action[s_o]
            s_orbit_reps.remove(s_o)
        i += 1

    systoles = dict()
    for o in s_orbit_reps:
        systole = shortest_systole_on_origami(o, horizontal_saddles, lower_bound)
        lower_bound = max(lower_bound, systole[0])
        systoles[o] = systole

    return systoles




def shortest_systole_on_origami(origami, horizontal_saddles, lower_bound):
    r"""
    This function checks whether the shortest systole on the given origami in
    H(1,1) has length greater or equal to `lower_bound` and returns such a
    systole if possible.

    INPUT:

    - ``origami`` -- the origami under consideration.
    - ``horizontal_saddles`` -- a dictionary mapping each origami in the orbit
      of `origami` to its shortest horizontal saddles. See the documentation of
      `shortest_horizontal_saddles` for more information.
    - ``lower_bound`` -- a real number. If a sysole of length less than
      `lower_bound` is found, the search is aborted and this systole is
      returned.

    OUTPUT: The function returns a tuple `(l,t,...)`, where `l` is the length
    of the shortest relevant systole and `t` is either the string `'loop'` or
    `'cycle'` depending on whether the shortest relevant systole connects a
    singularity on the given origami to itself or runs through both
    singularities.

    If the shortest relevant systole connects a singularity to itself, then the
    return value is a tuple `(l,t,d)` consisting of the length `l`, the type
    `t` being `'loop'` and the direction of the loop on the origami.

    If the shortest relevant systole is a cycle passing through both
    singularities, then the return value is a tuple `(l,t,edge_0,edge_1)`
    consisting of the length `l` of the sysole, the type `t` being `'cycle'`
    and the directions of the two saddle connections that constitute the
    systole.

    """

    _, _, s_action = origami.sl2z_edges()

    # We keep track of the shortest loop and the two shortest edges between the
    # two singularities (ordered by their length) in the list shortest_edges.
    # The auxiliary function _update_shortest_edges_ is used to update that
    # data as the implementation is a bit ugly.
    loop, edge_0, edge_1 = horizontal_saddles[origami]
    shortest_edges = [(loop, (1, 0)), (edge_0, (1, 0)), (edge_1, (1, 0))]

    # For each primitive direction (+/-dir_x, dir_y) in ZZ^2 with positive
    # dir_x and dir_y whose length is less than the current length of a minimal
    # cycle of saddle connections, we choose a matrix SL(2,ZZ) making the
    # direction (dir_x, dir_y) a horizontal direction and transport back the
    # precomputed information about horizontal saddle connections. The case
    # (dir_x,dir_y) = (1,0) are just the horizontal saddle connections on
    # origami and we treat the case (dir_x,dir_y) = (0,1) separately:
    min_cycle = _update_shortest_edges_(\
            shortest_edges,\
            horizontal_saddles[s_action[origami]], (0, 1))

    dir_y = 1
    while dir_y < min_cycle and min_cycle >= lower_bound:
        dir_x = 1
        while dir_x**2 + dir_y**2 < min_cycle**2 and min_cycle >= lower_bound:
            dir_origami, primitive = _origami_with_horizontal_saddle_direction_(\
                    origami, dir_x, dir_y, 1)
            if primitive:
                dir_len = sqrt(dir_x**2 + dir_y**2)
                _update_shortest_edges_(shortest_edges,\
                    (dir_len * l for l in horizontal_saddles[dir_origami]),\
                    (dir_x, dir_y))
                dir_origami, _ = _origami_with_horizontal_saddle_direction_(\
                        origami, dir_x, dir_y, -1)
                min_cycle = _update_shortest_edges_(shortest_edges,\
                    (dir_len * l for l in horizontal_saddles[dir_origami]),\
                    (-dir_x, dir_y))

            dir_x += 1
        dir_y += 1

    if min_cycle == shortest_edges[0][0]:
        return (min_cycle, 'loop', shortest_edges[0][1])
    return (min_cycle, 'cycle', shortest_edges[1], shortest_edges[2])


def _origami_with_horizontal_saddle_direction_(origami, dir_x, dir_y, sign):
    r"""
    This function computes a further origami in the SL(2,ZZ)-orbit of origami
    such that the direction (+/- x,y) on origami becomes the horizontal direction
    (1,0) on the new origami.

    INPUT:

    - ``origami`` -- the origami on which the direction dir_x, dir_y is under
      consideration.
    - ``dir_x`` -- an integer.
    - ``dir_y`` -- an integer.
    - ``sign`` -- an integer. If `sign` is equal to 1, then this function
      returns an origami with horizontal direction corresponding to (x,y) and
      otherwise an origami with horizontal direction corresponding to (-x,y).

    OUTPUT: A tuple `(o,p)` consisting of a boolean `p` and possibly an origami
    `o`. If the greatest common divisor of the input `dir_x` and `dir_y` is not
    equal to one, then no origami is returned and `p` is set to `False`.
    Otherwise, `p` is `True` and the sought-after origami is `o`.
    """

    d, a, b = xgcd(dir_x, dir_y)
    if d != 1:
        return (None, False)

    # Consider the matrices
    #
    # L = [1  1] [0 1]
    # R = [1  0] [1 1]
    # S = [0 -1] [1 0]
    #
    # The action of SL(2,ZZ) on the curve is encoded as three dictionaries
    #   l_action, r_action, s_action
    # mapping origamis to its image under the respective matrix.
    #
    # The object _SL2Z_Farey represents SL(2,ZZ) as being generated by the
    # matrices S and Q = [0 -1] [1 -1] though. We write Q = RS and Q^-1 =
    # LS^{-1} and may thus use the dictionaries as provided by the
    # surface_dynamics package.

    l_action, r_action, s_action = origami.sl2z_edges()
    result = origami
    sq_word = _SL2Z_Farey.word_problem(_SL2Z([[a, b], [-dir_y, dir_x]])) if sign == 1 \
            else _SL2Z_Farey.word_problem(_SL2Z([[-a, b], [-dir_y, -dir_x]]))
    for p in reversed(sq_word):
        if p == 1:
            result = s_action[result]
        elif p == 2:
            result = r_action[s_action[result]]
        elif p == -1:
            result = s_action[s_action[s_action[result]]]
        elif p == -2:
            result = l_action[s_action[s_action[s_action[result]]]]
    return (result, True)


def _shortest_horizontal_saddles(origami):
    r"""
    This auxiliary function computes the shortest horizontal saddle connections
    on an Origami with exactly two singularities.

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
    singularity_squares = [i for c in singularities for i in c]

    # Mapping from singularity_squares to the singularity (either 0 or 1) that
    # they contain.
    square_to_singularity = dict(\
            (s, 0) if s in singularities[0] else (s, 1)\
            for s in singularity_squares\
            )

    # We compute the shortest horizontal saddle connection by going to the
    # right from each square that contains one of the two singularities until
    # we hit another singularity.
    #
    # As there are only two singularities, we keep track of the shortest loop
    # and the two shortest edges between the two singularities.
    loop, edge_0, edge_1 = Infinity, Infinity, Infinity
    for square in singularity_squares:
        cur_len = 1
        cur_square = r(square)
        while not cur_square in singularity_squares:
            cur_len, cur_square = cur_len+1, r(cur_square)

        source, target = square_to_singularity[square], square_to_singularity[cur_square]
        if source == target:
            loop = min(loop, cur_len)
        elif cur_len <= edge_0:
            edge_0, edge_1 = cur_len, edge_0
        elif cur_len < edge_1:
            edge_1 = cur_len
    return (loop, edge_0, edge_1)


def _update_shortest_edges_(shortest_edges, update, direction):
    r"""
    This function is an auxiliary helper function that deals with some ugly
    implementation details concerning the way shortest systoles on origamis in
    H(1,1) are represented in the function `shortest_systole_on_origami`.
    """

    new_loop, new_edge_0, new_edge_1 = update

    if new_loop < shortest_edges[0][0]:
        shortest_edges[0] = (new_loop, direction)
    if new_edge_0 <= shortest_edges[1][0]:
        if new_edge_1 < shortest_edges[1][0]:
            shortest_edges[2] = (new_edge_1, direction)
        else:
            shortest_edges[2] = shortest_edges[1]
        shortest_edges[1] = (new_edge_0, direction)
    elif new_edge_0 < shortest_edges[2][0]:
        shortest_edges[2] = (new_edge_0, direction)

    return min(shortest_edges[0][0], shortest_edges[1][0]+shortest_edges[2][0])


# vim:ft=python
