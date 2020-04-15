"""
FILE: discretisation.py
LAST MODIFIED: 24-12-2015
DESCRIPTION:
Modules for discretising meshes.
General Inputs: discretisation scheme
                discretisation density
                geometric field mesh
                
General outputs: xi coordinates per mesh element
    
===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

import numpy


def discretiseRegularXi(GF, D):
    """
    regular discretisation in Xi space
    """

    pass


def discretiseRegularGeoN(GF, D):
    """
    regular discretisation in geometric space given fixed number of 
    points per xi direction
    """
    pass


def discretiseRegularGeoD(d, elemEval):
    """
    regular discretisation in geometric space given a fixed spacing
    distance
    """

    d2 = d * d

    if 'quad' in elemEval.element.type:
        if elemEval.element.dimensions == 2:
            # initial points = 4 corners
            pointsXiInit = [(0.0, 0.0),
                            (1.0, 0.0),
                            (0.0, 1.0),
                            (1.0, 1.0),
                            ]

            # subdivide
            pointsXi = _subdivideQuad2D(pointsXiInit, elemEval, d2)
        elif elemEval.element.dimensions == 3:
            # initial points = 4 corners
            pointsXiInit = numpy.array([[0.0, 0.0, 0.0],
                                        [1.0, 0.0, 0.0],
                                        [0.0, 1.0, 0.0],
                                        [1.0, 1.0, 0.0],
                                        [0.0, 0.0, 1.0],
                                        [1.0, 0.0, 1.0],
                                        [0.0, 1.0, 1.0],
                                        [1.0, 1.0, 1.0],
                                        ])

            # subdivide
            pointsXi = _subdivideQuad3D(pointsXiInit, elemEval, d2)
        else:
            raise NotImplementedError('only 2 and 3 dimensional quad elements supported')

    elif 'tri' in elemEval.element.type:
        if elemEval.element.dimensions == 2:
            # initial points = 4 corners
            pointsXiInit = [(0.0, 0.0),
                            (1.0, 0.0),
                            (0.0, 1.0),
                            ]

            # subdivide
            pointsXi = _subdivideTri2D(pointsXiInit, elemEval, d2)
        else:
            raise NotImplementedError('only 2 dimensional tri elements supported')
    else:
        raise TypeError('unrecognised element type: ' + elemEval.elemType)

    return numpy.array(list(pointsXi))


def _subdivideQuad2D(cornerXi, elemEval, maxDistanceSquared):
    cornerCoord = elemEval.eval(cornerXi)
    divideMode = _checkQuad2D(cornerCoord, maxDistanceSquared)

    if divideMode == 0:
        # ~ print 0
        return set(cornerXi)
    elif (divideMode == 1) or (divideMode == 4):
        # subdivide
        # G---h---I
        # |   |   |
        # d---e---f
        # |   |   |
        # A---b---C
        # caps are the original 
        # ~ print 1
        A, C, G, I = cornerXi

        b = (0.5 * (A[0] + C[0]), 0.5 * (A[1] + C[1]))
        d = (0.5 * (A[0] + G[0]), 0.5 * (A[1] + G[1]))
        f = (0.5 * (C[0] + I[0]), 0.5 * (C[1] + I[1]))
        h = (0.5 * (G[0] + I[0]), 0.5 * (G[1] + I[1]))
        e = (0.5 * (d[0] + f[0]), 0.5 * (d[1] + f[1]))

        pointsQ1 = _subdivideQuad2D([A, b, d, e], elemEval, maxDistanceSquared)
        pointsQ2 = _subdivideQuad2D([b, C, e, f], elemEval, maxDistanceSquared)
        pointsQ3 = _subdivideQuad2D([d, e, G, h], elemEval, maxDistanceSquared)
        pointsQ4 = _subdivideQuad2D([e, f, h, I], elemEval, maxDistanceSquared)
        return set(list(pointsQ1) + list(pointsQ2) + list(pointsQ3) + list(pointsQ4))

    elif divideMode == 2:
        # subdivide
        # D---e---F
        # |   |   |
        # A---b---C
        # caps are the original 
        # ~ print 2
        A, C, D, F = cornerXi

        b = (0.5 * (A[0] + C[0]), 0.5 * (A[1] + C[1]))
        e = (0.5 * (D[0] + F[0]), 0.5 * (D[1] + F[1]))

        pointsQ1 = _subdivideQuad2D([A, b, D, e], elemEval, maxDistanceSquared)
        pointsQ2 = _subdivideQuad2D([b, C, e, F], elemEval, maxDistanceSquared)
        return set(list(pointsQ1) + list(pointsQ2))

    elif divideMode == 3:
        # subdivide
        # E---F
        # |   |
        # c---d
        # |   |
        # A---B
        # caps are the original 
        # ~ print 3
        A, B, E, F = cornerXi

        c = (0.5 * (A[0] + E[0]), 0.5 * (A[1] + E[1]))
        d = (0.5 * (B[0] + F[0]), 0.5 * (B[1] + F[1]))

        pointsQ1 = _subdivideQuad2D([A, B, c, d], elemEval, maxDistanceSquared)
        pointsQ2 = _subdivideQuad2D([c, d, E, F], elemEval, maxDistanceSquared)
        return set(list(pointsQ1) + list(pointsQ2))


def _checkQuad2D(corners, maxDistanceSquared):
    # c---d
    # |   |
    # a---b

    # a-b, b-d, d-c, c-a
    a, b, c, d = corners
    X1 = numpy.array([a, b, d, c])
    X2 = numpy.array([b, d, c, a])

    edgeDistsSquared = ((X1 - X2) ** 2.0).sum(1)

    lengthCheck = (edgeDistsSquared > maxDistanceSquared).astype(int)

    if lengthCheck.sum() > 2:
        # ~ print 1, lengthCheck
        return 1
    elif (lengthCheck[0] + lengthCheck[2]) == 2:
        # ~ print 2, lengthCheck
        return 2
    elif (lengthCheck[1] + lengthCheck[3]) == 2:
        # ~ print 3, lengthCheck
        return 3
    # ~ elif (lengthCheck[0] + lengthCheck[1]) == 2:
    # ~ print 4, lengthCheck
    # ~ return 4
    # ~ elif (lengthCheck[2] + lengthCheck[3]) == 2:
    # ~ print 4, lengthCheck
    # ~ return 4
    else:
        # ~ print 0, lengthCheck
        return 0


def _subdivideQuad3D(cornerXi, elemEval, maxDistanceSquared):
    def nodiv(cornerXi, elemEval, maxDistanceSquared):
        return set([tuple(p) for p in cornerXi])

    def xdiv(cornerXi, elemEval, maxDistanceSquared):
        # subdivide
        # D---e---F J---k---L
        # |   |   | |   |   |
        # |   |   | |   |   |
        # |   |   | |   |   |
        # A---b---C G---h---I
        # caps are the original 

        A, C, D, F, G, I, J, L = cornerXi
        b = 0.5 * (A + C)
        e = 0.5 * (D + F)
        h = 0.5 * (G + I)
        e = 0.5 * (J + L)

        pointsQ1 = _subdivideQuad3D([A, b, D, e, G, h, J, k], elemEval, maxDistanceSquared)
        pointsQ2 = _subdivideQuad3D([b, C, e, F, h, I, k, L], elemEval, maxDistanceSquared)
        return set([tuple(p) for p in pointsQ1] + [tuple(p) for p in pointsQ2])

    def ydiv(cornerXi, elemEval, maxDistanceSquared):
        # subdivide
        # E-------F K-------L
        # |       | |       |
        # c-------d i-------j
        # |       | |       |
        # A-------B G-------H
        # caps are the original 

        A, B, E, F, G, H, K, L = cornerXi
        c = 0.5 * (A + E)
        d = 0.5 * (B + F)
        i = 0.5 * (G + K)
        j = 0.5 * (H + L)

        pointsQ1 = _subdivideQuad3D([A, B, c, d, G, H, i, j], elemEval, maxDistanceSquared)
        pointsQ2 = _subdivideQuad3D([c, d, E, F, i, j, K, L], elemEval, maxDistanceSquared)
        return set([tuple(p) for p in pointsQ1] + [tuple(p) for p in pointsQ2])

    def zdiv(cornerXi, elemEval, maxDistanceSquared):
        # subdivide
        # C-------D g-------h K-------L
        # |       | |       | |       |
        # |       | |       | |       |
        # |       | |       | |       |
        # A-------B e-------f I-------J
        # caps are the original 

        A, B, C, D, I, J, K, L = cornerXi
        e = 0.5 * (A + I)
        f = 0.5 * (B + J)
        g = 0.5 * (C + K)
        h = 0.5 * (D + L)

        pointsQ1 = _subdivideQuad3D([A, B, C, D, e, f, g, h], elemEval, maxDistanceSquared)
        pointsQ2 = _subdivideQuad3D([e, f, g, h, I, J, K, L], elemEval, maxDistanceSquared)
        return set([tuple(p) for p in pointsQ1] + [tuple(p) for p in pointsQ2])

    def alldiv(cornerXi, elemEval, maxDistanceSquared):
        # subdivide
        # G---h---I p---q---r Y---z---A1
        # |   |   | |   |   | |   |   |
        # d---e---f m---n---o v---w---x
        # |   |   | |   |   | |   |   |
        # A---b---C j---k---l S---t---U
        # caps are the original 

        A, C, G, I, S, U, Y, A1 = cornerXi
        b = 0.5 * (A + C)
        d = 0.5 * (A + G)
        f = 0.5 * (C + I)
        e = 0.5 * (G + C)
        h = 0.5 * (G + I)

        t = 0.5 * (S + U)
        v = 0.5 * (S + Y)
        x = 0.5 * (U + A1)
        w = 0.5 * (Y + U)
        z = 0.5 * (Y + A1)

        j = 0.5 * (A + S)
        k = 0.5 * (A + U)
        l = 0.5 * (C + U)
        m = 0.5 * (G + S)
        n = 0.5 * (A + A1)
        o = 0.5 * (I + U)
        p = 0.5 * (G + Y)
        q = 0.5 * (G + A1)
        r = 0.5 * (I + A1)

        pointsQ1 = _subdivideQuad3D([A, b, d, e, j, k, m, n], elemEval, maxDistanceSquared)
        pointsQ2 = _subdivideQuad3D([b, C, e, f, k, l, n, o], elemEval, maxDistanceSquared)
        pointsQ3 = _subdivideQuad3D([d, e, G, h, m, n, p, q], elemEval, maxDistanceSquared)
        pointsQ4 = _subdivideQuad3D([e, f, h, I, n, o, q, r], elemEval, maxDistanceSquared)

        pointsQ5 = _subdivideQuad3D([j, k, m, n, S, t, v, w], elemEval, maxDistanceSquared)
        pointsQ6 = _subdivideQuad3D([k, l, n, o, t, U, w, x], elemEval, maxDistanceSquared)
        pointsQ7 = _subdivideQuad3D([m, n, p, q, v, w, Y, z], elemEval, maxDistanceSquared)
        pointsQ8 = _subdivideQuad3D([n, o, q, r, w, x, z, A1], elemEval, maxDistanceSquared)

        return set([tuple(p) for p in pointsQ1] + [tuple(p) for p in pointsQ2] + \
                   [tuple(p) for p in pointsQ3] + [tuple(p) for p in pointsQ4] + \
                   [tuple(p) for p in pointsQ5] + [tuple(p) for p in pointsQ6] + \
                   [tuple(p) for p in pointsQ7] + [tuple(p) for p in pointsQ8])

    divideModes = {0: nodiv,
                   1: xdiv,
                   3: ydiv,
                   4: alldiv,  # xydiv,
                   5: zdiv,
                   6: alldiv,  # xzdiv,
                   8: alldiv,  # yzdiv,
                   9: alldiv
                   }

    cornerCoord = elemEval.eval(cornerXi)
    divideMode = _checkQuad3D(cornerCoord, maxDistanceSquared)
    return divideModes[divideMode](cornerXi, elemEval, maxDistanceSquared)


def _checkQuad3D(corners, maxDistanceSquared):
    # 2---3 6---7
    # | 1 | | 2 |
    # 0---1 4---5
    #
    # returnCode:
    # 0: no div
    # 1: x div
    # 3: y div
    # 4: x,y div
    # 5: z div
    # 6: x,z div
    # 8: y,z div
    # 9: all div

    a, b, c, d, e, f, g, h = corners
    # 0-1, 2-3, 4-5, 6-7, x direction lengths
    x1 = corners[[0, 2, 3, 6]]
    x2 = corners[[1, 3, 5, 7]]
    xLengths = ((x1 - x2) ** 2.0).sum(1)

    # 0-2, 1-3, 4-6, 5-7, y direction lengths
    y1 = corners[[0, 1, 4, 5]]
    y2 = corners[[2, 3, 6, 7]]
    yLengths = ((y1 - y2) ** 2.0).sum(1)

    # 0-4, 1-5, 2-6, 3-7, z direction lengths
    z1 = corners[[0, 1, 2, 3]]
    z2 = corners[[4, 5, 6, 7]]
    zLengths = ((z1 - z2) ** 2.0).sum(1)

    xLengthCheck = (xLengths > maxDistanceSquared).astype(int)
    yLengthCheck = (yLengths > maxDistanceSquared).astype(int)
    zLengthCheck = (zLengths > maxDistanceSquared).astype(int)

    returnCode = 0

    if xLengthCheck.sum() > 2:
        returnCode += 1
    if yLengthCheck.sum() > 2:
        returnCode += 3
    if zLengthCheck.sum() > 2:
        returnCode += 5

    return returnCode


def _subdivideTri2D(cornerXi, elemEval, maxDistanceSquared):
    cornerCoord = elemEval.eval(cornerXi)
    if _checkTri2D(cornerCoord, maxDistanceSquared):
        # subdivide
        # F
        # | \
        # d--e 
        # | \| \  
        # A--b--C
        # caps are the original 

        A, C, F = cornerXi
        # ~ b = (A+C)/2.0
        # ~ d = (A+F)/2.0
        # ~ e = (F+C)/2.0

        b = (0.5 * (A[0] + C[0]), 0.5 * (A[1] + C[1]))
        d = (0.5 * (A[0] + F[0]), 0.5 * (A[1] + F[1]))
        e = (0.5 * (F[0] + C[0]), 0.5 * (F[1] + C[1]))

        # ~ b = 0.5*(A+C)
        # ~ d = 0.5*(A+F)
        # ~ e = 0.5*(F+C)

        pointsT1 = _subdivideTri2D([A, b, d], elemEval, maxDistanceSquared)
        pointsT2 = _subdivideTri2D([e, d, b], elemEval, maxDistanceSquared)
        pointsT3 = _subdivideTri2D([b, C, e], elemEval, maxDistanceSquared)
        pointsT4 = _subdivideTri2D([d, e, F], elemEval, maxDistanceSquared)

        return set(list(pointsT1) + list(pointsT2) + list(pointsT3) + list(pointsT4))
    else:
        return set(cornerXi)


def _checkTri2D(corners, maxDistanceSquared):
    # c
    # |\ 
    # a-b

    # a-b, b-c, c-a
    a, b, c = corners
    X1 = numpy.array([a, b, c])
    X2 = numpy.array([b, c, a])

    edgeDistsSquared = ((X1 - X2) ** 2.0).sum(1)

    # ~ if any(edgeDistsSquared > maxDistanceSquared):
    if sum(edgeDistsSquared > maxDistanceSquared) > 1:
        return 1
    else:
        return 0
