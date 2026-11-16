"""Fusion ring automorphisms: all permutations pi of the simple objects, fixing the unit, with
N_{pi X, pi Y}^{pi Z} = N_{XY}^Z (backtracking).

As a library:  from fusion_aut import fusion_automorphisms;  fusion_automorphisms(obj, unit, N)
               -> list of dicts pi (N: dict (X, Y, Z) -> multiplicity, zero entries omitted).
"""


def fusion_automorphisms(obj, unit, N):
    n = lambda X, Y, Z: N.get((X, Y, Z), 0)
    others = [x for x in obj if x != unit]
    out = []

    def ok(pi):
        dom = list(pi)
        for X in dom:
            for Y in dom:
                for Z in dom:
                    if n(X, Y, Z) != n(pi[X], pi[Y], pi[Z]):
                        return False
        return True

    def rec(pi, rest):
        if not ok(pi):
            return
        if len(pi) == len(obj):
            out.append(dict(pi))
            return
        X = others[len(pi) - 1]
        for Y in rest:
            pi[X] = Y
            rec(pi, rest - {Y})
            del pi[X]

    rec({unit: unit}, set(others))
    return out
