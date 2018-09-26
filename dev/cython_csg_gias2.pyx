
# cython: profile=True

# import cython
# import math
import operator
from functools import reduce

# from cpython cimport array
# import array

import cython
cimport cython
# cimport numpy as np

# increase the max number of recursive calls
# sys.setrecursionlimit(10000)  # my default is 1000, increasing too much may cause a seg fault

# cdef array.array int_array_template = array.array('i', [])

from libc.math cimport sqrt, sin, cos, atan2, abs

cdef double PI = 3.14159265358979323846264338327

cdef double dot_vectors(Vector a, Vector b):
    return a.x * b.x + a.y * b.y + a.z * b.z

@cython.nonecheck(False)
cdef class Vector(object):
    """
    class Vector

    Represents a 3D vector.

    Example usage:
         Vector(1, 2, 3);
         Vector([1, 2, 3]);
         Vector({ 'x': 1, 'y': 2, 'z': 3 });
    """

    cdef public double x
    cdef public double y
    cdef public double z

    # cdef tuple __slots__ = ('x','y','z')
    # def __cinit__(self, double x, double y, double z):
    #     self.x = x
    #     self.y = y
    #     self.z = z

    def __init__(self, double x, double y, double z):
        self.x = x
        self.y = y
        self.z = z

        # self.x, self.y, self.z = 0., 0., 0.
        # if len(args) == 3:
        #     self.x = args[0]
        #     self.y = args[1]
        #     self.z = args[2]
        # elif len(args) == 1:
        #     a = args[0]
        #     if isinstance(a, dict):
        #         self.x = a.get('x', 0.0)
        #         self.y = a.get('y', 0.0)
        #         self.z = a.get('z', 0.0)
        #     elif a is not None and len(a) == 3:
        #         self.x = a[0]
        #         self.y = a[1]
        #         self.z = a[2]

    def __repr__(self):
        return '({0}, {1}, {2})'.format(self.x, self.y, self.z)

    cpdef Vector clone(self):
        """ Clone. """
        return Vector(self.x, self.y, self.z)

    cpdef Vector negated(self):
        """ Negated. """
        return Vector(-self.x, -self.y, -self.z)

    def __neg__(self):
        return self.negated()

    cpdef Vector plus(self, Vector a):
        """ Add. """
        return Vector(self.x + a.x, self.y + a.y, self.z + a.z)

    def __add__(self, Vector a):
        return self.plus(a)

    cpdef Vector minus(self, Vector a):
        """ Subtract. """
        return Vector(self.x - a.x, self.y - a.y, self.z - a.z)

    def __sub__(self, Vector a):
        return self.minus(a)

    cpdef Vector times(self, double a):
        """ Multiply. """
        return Vector(self.x * a, self.y * a, self.z * a)

    def __mul__(self, double a):
        return self.times(a)

    cpdef Vector dividedBy(self, double a):
        """ Divide. """
        return Vector(self.x / a, self.y / a, self.z / a)

    def __truediv__(self, double a):
        return self.dividedBy(float(a))

    def __div__(self, double a):
        return self.dividedBy(a)

    cpdef double dot(self, Vector a):
        """ Dot. """
        return self.x * a.x + self.y * a.y + self.z * a.z

    cpdef Vector lerp(self, Vector a, double t):
        """ Lerp. Linear interpolation from self to a"""
        return self.plus(a.minus(self).times(t));

    cpdef double length(self):
        """ Length. """
        # return math.sqrt(self.dot(self))
        return sqrt(dot_vectors(self, self))

    cpdef Vector unit(self):
        """ Normalize. """
        return self.dividedBy(self.length())

    cpdef Vector cross(self, Vector a):
        """ Cross. """
        return Vector(
            self.y * a.z - self.z * a.y,
            self.z * a.x - self.x * a.z,
            self.x * a.y - self.y * a.x)

    def __getitem__(self, int key):
        return (self.x, self.y, self.z)[key]

    def __setitem__(self, int key, double value):
        l = [self.x, self.y, self.z]
        l[key] = value
        self.x, self.y, self.z = l

    def __len__(self):
        return 3

    def __iter__(self):
        return iter((self.x, self.y, self.z))

    def __repr__(self):
        return 'Vector(%.2f, %.2f, %0.2f)' % (self.x, self.y, self.z)


cdef class Vertex(object):
    """
    Class Vertex

    Represents a vertex of a polygon. Use your own vertex class instead of this
    one to provide additional features like texture coordinates and vertex
    colors. Custom vertex classes need to provide a `pos` property and `clone()`,
    `flip()`, and `interpolate()` methods that behave analogous to the ones
    defined by `Vertex`. This class provides `normal` so convenience
    functions like `CSG.sphere()` can return a smooth vertex normal, but `normal`
    is not used anywhere else.
    """

    cdef public Vector pos
    cdef public Vector normal

    def __init__(self, Vector pos, Vector normal):
        self.pos = pos
        self.normal = normal

    cpdef Vertex clone(self):
        return Vertex(self.pos.clone(), self.normal.clone())

    cpdef void flip(self):
        """
        Invert all orientation-specific data (e.g. vertex normal). Called when the
        orientation of a polygon is flipped.
        """
        self.normal = self.normal.negated()

    cpdef Vertex interpolate(self, Vertex other, double t):
        """
        Create a new vertex between this vertex and `other` by linearly
        interpolating all properties using a parameter of `t`. Subclasses should
        override this to interpolate additional properties.
        """
        return Vertex(self.pos.lerp(other.pos, t),
                      self.normal.lerp(other.normal, t))

    def __repr__(self):
        return repr(self.pos)

cpdef Plane planeFromPoints(Vector a, Vector b, Vector c):
    n = b.minus(a).cross(c.minus(a)).unit()
    # return Plane(n, n.dot(a))
    return Plane(n, dot_vectors(n, a))

cdef class Plane(object):
    """
    class Plane

    Represents a plane in 3D space.
    """

    """
    `Plane.EPSILON` is the tolerance used by `splitPolygon()` to decide if a
    point is on the plane.
    """
    # cdef tuple __slots__ = ('normal', 'w')


    cdef public double EPSILON
    cdef public Vector normal
    cdef public double w

    def __init__(self, Vector normal, double w):
        self.normal = normal
        # w is the (perpendicular) distance of the plane from (0, 0, 0)
        self.w = w
        self.EPSILON = 1e-5

    cpdef Plane clone(self):
        return Plane(self.normal.clone(), self.w)

    cpdef void flip(self):
        self.normal = self.normal.negated()
        self.w = -self.w

    def __repr__(self):
        return 'normal: {0} w: {1}'.format(self.normal, self.w)

    @cython.nonecheck(False)
    cpdef void splitPolygon(self, Polygon polygon, list coplanarFront,
        list coplanarBack, list front, list back):
        """
        Split `polygon` by this plane if needed, then put the polygon or polygon
        fragments in the appropriate lists. Coplanar polygons go into either
        `coplanarFront` or `coplanarBack` depending on their orientation with
        respect to this plane. Polygons in front or in back of this plane go into
        either `front` or `back`
        """
        cdef int COPLANAR = 0  # all the vertices are within EPSILON distance from plane
        cdef int FRONT = 1  # all the vertices are in front of the plane
        cdef int BACK = 2  # all the vertices are at the back of the plane
        cdef int SPANNING = 3  # some vertices are in front, some in the back
        cdef int numVertices = len(polygon.vertices)
        cdef Py_ssize_t i, j, ti, tj
        cdef int loc
        cdef double t, normalDotPlaneNormal, t1, t2
        cdef double EPS = self.EPSILON
        cdef double NEPS = -1.0*self.EPSILON
        cdef Vertex vi, vj, v
        cdef Polygon poly
        cdef Vector pos, pos2

        # Classify each point as well as the entire polygon into one of the above
        # four classes.
        cdef int polygonType = 0
        cdef list vertexLocs = []  # list of ints
        cdef list f = []  # list of vertices
        cdef list b = []  # list of vertices
        
        for i in range(numVertices):
            # t = self.normal.dot(polygon.vertices[i].pos) - self.w
            v = polygon.vertices[i]
            t = dot_vectors(self.normal, v.pos) - self.w
            loc = -1
            if t < NEPS:
                loc = BACK
            elif t > EPS:
                loc = FRONT
            else:
                loc = COPLANAR
            polygonType |= loc
            vertexLocs.append(loc)  # list op

        # Put the polygon in the correct list, splitting it when necessary.
        if polygonType == COPLANAR:
            # normalDotPlaneNormal = self.normal.dot(polygon.plane.normal)
            normalDotPlaneNormal = dot_vectors(self.normal, polygon.plane.normal)
            if normalDotPlaneNormal > 0.0:
                coplanarFront.append(polygon)  # list op
            else:
                coplanarBack.append(polygon)  # list op
        elif polygonType == FRONT:
            front.append(polygon)  # list op
        elif polygonType == BACK:
            back.append(polygon)  # list op
        elif polygonType == SPANNING:
            f = []
            b = []
            for i in range(numVertices):
                j = (i + 1) % numVertices
                ti = vertexLocs[i]
                tj = vertexLocs[j]
                vi = polygon.vertices[i]
                vj = polygon.vertices[j]
                if ti != BACK:
                    f.append(vi)  # list op
                if ti != FRONT:
                    if ti != BACK:
                        b.append(vi.clone())  # list op
                    else:
                        b.append(vi)  # list op
                if (ti | tj) == SPANNING:
                    # interpolation weight at the intersection point
                    # t = (self.w - self.normal.dot(vi.pos)) / self.normal.dot(vj.pos.minus(vi.pos))
                    pos = vi.pos
                    t1 = dot_vectors(self.normal, pos)
                    pos2 = vj.pos.minus(pos)
                    t2 = dot_vectors(self.normal, pos2)
                    t = (self.w - t1) / t2
                    # intersection point on the plane
                    v = vi.interpolate(vj, t)
                    f.append(v)  # list op
                    b.append(v.clone())  # list op
            if len(f) >= 3:
                poly = Polygon(f, polygon.shared)
                front.append(poly)  # list op
            if len(b) >= 3:
                poly = Polygon(b, polygon.shared)
                back.append(poly)  # list op


cdef class Polygon(object):
    """
    class Polygon

    Represents a convex polygon. The vertices used to initialize a polygon must
    be coplanar and form a convex loop. They do not have to be `Vertex`
    instances but they must behave similarly (duck typing can be used for
    customization).

    Each convex polygon has a `shared` property, which is shared between all
    polygons that are clones of each other or were split from the same polygon.
    This can be used to define per-polygon properties (such as surface color).
    """

    cdef public list vertices
    cdef public bint shared
    cdef public Plane plane

    # cdef tuple __slots__ = (
    #     'vertices',
    #     'shared',
    #     'plane'
    #     )

    def __init__(self, list vertices, bint shared):
        self.vertices = vertices
        self.shared = shared
        self.plane = planeFromPoints(vertices[0].pos, vertices[1].pos, vertices[2].pos)

    cpdef Polygon clone(self):
        # vertices = list(map(lambda v: v.clone(), self.vertices))
        cdef list vertices
        cdef Py_ssize_t vi
        cdef int nverts
        cdef Vertex v

        vertices = []
        nverts = len(self.vertices)
        for vi in range(nverts):
            v = self.vertices[vi]
            vertices.append(v.clone())

        return Polygon(vertices, self.shared)

    cpdef void flip(self):
        cdef list vertices
        cdef Py_ssize_t vi
        cdef int nverts = len(self.vertices)
        cdef Vertex v

        self.vertices.reverse()
        
        # map(lambda v: v.flip(), self.vertices)
        for vi in range(nverts):
            v = self.vertices[vi]
            v.flip()

        self.plane.flip()

    def __repr__(self):
        return reduce(lambda x, y: x + y,
                      ['Polygon(['] + [repr(v) + ', ' \
                                       for v in self.vertices] + ['])'], '')


cdef class BSPNode(object):
    """
    class BSPNode

    Holds a node in a BSP tree. A BSP tree is built from a collection of polygons
    by picking a polygon to split along. That polygon (and all other coplanar
    polygons) are added directly to that node and the other polygons are added to
    the front and/or back subtrees. This is not a leafy BSP tree since there is
    no distinction between internal and leaf nodes.
    """

    # __slots__ = ('plane',
    #              'front',
    #              'back',
    #              'polygons')

    cdef public Plane plane
    cdef public BSPNode front
    cdef public BSPNode back
    cdef public list polygons

    def __init__(self, list polygons=None):
        self.plane = None  # Plane instance
        self.front = None  # BSPNode
        self.back = None  # BSPNode
        self.polygons = []
        if polygons:
            self.build(polygons)

    cpdef BSPNode clone(self):
        cdef int npolys
        cdef Py_ssize_t pi
        cdef BSPNode node = BSPNode()
        cdef Polygon poly
        
        if self.plane:
            node.plane = self.plane.clone()
        if self.front:
            node.front = self.front.clone()
        if self.back:
            node.back = self.back.clone()

        # node.polygons = list(map(lambda p: p.clone(), self.polygons))
        node.polygons = []
        npolys = len(self.polygons)
        for pi in range(npolys):
            poly = self.polygons[pi]
            node.polygons.append(poly.clone())

        return node

    cpdef void invert(self):
        """
        Convert solid space to empty space and empty space to solid space.
        """
        # Polygon([Vector(-14.00, 0.00, 37.00), Vector(-14.00, 796.87, 37.00), Vector(-16.00, 0.00, 37.00), ])
        cdef Polygon poly
        cdef BSPNode temp
        cdef int npolys = len(self.polygons)
        cdef Py_ssize_t pi

        for pi in range(npolys):
            poly = self.polygons[pi]
            poly.flip()

        # for poly in self.polygons:
        #     poly.flip()
            # x = 100
        self.plane.flip()
        if self.front:
            self.front.invert()
        if self.back:
            self.back.invert()
        temp = self.front
        self.front = self.back
        self.back = temp

    cpdef list clipPolygons(self, list polygons):
        """
        Recursively remove all polygons in `polygons` that are inside this BSP
        tree.
        """
        cdef list front = []
        cdef list back = []
        cdef Polygon poly
        cdef int npolys
        cdef Py_ssize_t pi

        if not self.plane:
            return polygons[:]

        npolys = len(polygons)
        for pi in range(npolys):
            poly = polygons[pi]
            self.plane.splitPolygon(poly, front, back, front, back)
        # for poly in polygons:
        #     self.plane.splitPolygon(poly, front, back, front, back)

        if self.front:
            front = self.front.clipPolygons(front)

        if self.back:
            back = self.back.clipPolygons(back)
        else:
            back = []

        front.extend(back)
        return front

    cpdef void clipTo(self, BSPNode bsp):
        """
        Remove all polygons in this BSP tree that are inside the other BSP tree
        `bsp`.
        """
        self.polygons = bsp.clipPolygons(self.polygons)
        if self.front:
            self.front.clipTo(bsp)
        if self.back:
            self.back.clipTo(bsp)

    cpdef list allPolygons(self):
        """
        Return a list of all polygons in this BSP tree.
        """
        cdef list polygons = self.polygons[:]
        if self.front:
            polygons.extend(self.front.allPolygons())
        if self.back:
            polygons.extend(self.back.allPolygons())
        return polygons

    cpdef build(self, list polygons):
        """
        Build a BSP tree out of `polygons`. When called on an existing tree, the
        new polygons are filtered down to the bottom of the tree and become new
        nodes there. Each set of polygons is partitioned using the first polygon
        (no heuristic is used to pick a good split).
        """
        cdef list front, back
        cdef int npolys
        cdef Polygon poly
        cdef Py_ssize_t pi

        if len(polygons) == 0:
            return
        if not self.plane:
            poly = polygons[0]
            self.plane = poly.plane.clone()

        # add polygon to this node
        self.polygons.append(polygons[0])
        front = []
        back = []
        # split all other polygons using the first polygon's plane
        npolys = len(polygons)
        for pi in range(1, npolys):
            # coplanar front and back polygons go into self.polygons
            poly = polygons[pi]
            self.plane.splitPolygon(
                poly, self.polygons, self.polygons,
                front, back
                )
        # for poly in polygons[1:]:
        #     # coplanar front and back polygons go into self.polygons
        #     self.plane.splitPolygon(poly, self.polygons, self.polygons,
        #                             front, back)
        # recursively build the BSP tree
        if len(front) > 0:
            if not self.front:
                self.front = BSPNode()
            self.front.build(front)
        if len(back) > 0:
            if not self.back:
                self.back = BSPNode()
            self.back.build(back)


cpdef CSG csgFromPolygons(list polygons):
    cdef CSG csg

    csg = CSG()
    csg.polygons = polygons
    return csg

cdef Vector rotateVector(Vector v, Vector ax, double cosAngle, double sinAngle):
    cdef double vA, vPerpLen, vCosA, vSinA
    cdef Vector vPerp, u1, u2

    # vA = v.dot(ax)
    vA = dot_vectors(v, ax)
    vPerp = v.minus(ax.times(vA))
    vPerpLen = vPerp.length()
    if vPerpLen == 0:
        # vector is parallel to axis, no need to rotate
        return v
    u1 = vPerp.unit()
    u2 = u1.cross(ax)
    vCosA = vPerpLen * cosAngle
    vSinA = vPerpLen * sinAngle
    return ax.times(vA).plus(u1.times(vCosA).plus(u2.times(vSinA)))

cdef class CSG(object):
    """
    Constructive Solid Geometry (CSG) is a modeling technique that uses Boolean
    operations like union and intersection to combine 3D solids. This library
    implements CSG operations on meshes elegantly and concisely using BSP trees,
    and is meant to serve as an easily understandable implementation of the
    algorithm. All edge cases involving overlapping coplanar polygons in both
    solids are correctly handled.

    Example usage::

        from csg.core import CSG

        cube = CSG.cube();
        sphere = CSG.sphere({'radius': 1.3});
        polygons = cube.subtract(sphere).toPolygons();

    ## Implementation Details

    All CSG operations are implemented in terms of two functions, `clipTo()` and
    `invert()`, which remove parts of a BSP tree inside another BSP tree and swap
    solid and empty space, respectively. To find the union of `a` and `b`, we
    want to remove everything in `a` inside `b` and everything in `b` inside `a`,
    then combine polygons from `a` and `b` into one solid::

        a.clipTo(b);
        b.clipTo(a);
        a.build(b.allPolygons());

    The only tricky part is handling overlapping coplanar polygons in both trees.
    The code above keeps both copies, but we need to keep them in one tree and
    remove them in the other tree. To remove them from `b` we can clip the
    inverse of `b` against `a`. The code for union now looks like this::

        a.clipTo(b);
        b.clipTo(a);
        b.invert();
        b.clipTo(a);
        b.invert();
        a.build(b.allPolygons());

    Subtraction and intersection naturally follow from set operations. If
    union is `A | B`, subtraction is `A - B = ~(~A | B)` and intersection is
    `A & B = ~(~A | ~B)` where `~` is the complement operator.

    ## License

    Copyright (c) 2011 Evan Wallace (http://madebyevan.com/), under the MIT license.

    Python port Copyright (c) 2012 Tim Knip (http://www.floorplanner.com), under the MIT license.
    Additions by Alex Pletzer (Pennsylvania State University)
    """

    # cdef tuple __slots__ = ('polygons')
    cdef public list polygons

    def __init__(self):
        self.polygons = []

    cpdef CSG clone(self):
        cdef CSG csg
        cdef int npolys
        cdef Polygon poly
        cdef Py_ssize_t pi

        csg = CSG()
        # csg.polygons = list(map(lambda p: p.clone(), self.polygons))
        csg.polygons = []
        npolys = len(self.polygons)
        for pi in range(npolys):
            poly = self.polygons[pi]
            csg.polygons.append(poly.clone())
        return csg

    cpdef list toPolygons(self):
        return self.polygons

    cpdef CSG refine(self):
        """
        Return a refined CSG. To each polygon, a middle point is added to each edge and to the center
        of the polygon
        """
        cdef CSG newCSG
        cdef Polygon poly, newPoly
        cdef list verts, newVerts, vs
        cdef int numVerts
        cdef int npolys = len(self.polygons)
        cdef Py_ssize_t pi, i, vi
        cdef Vector midPos, midNormal, midPosAccum
        cdef Vertex midVert

        newCSG = CSG()
        for pi in range(npolys):
            poly = self.polygon[pi]
            verts = poly.vertices
            numVerts = len(verts)

            if numVerts == 0:
                continue

            midPos = reduce(operator.add, [v.pos for v in verts]) / float(numVerts)
            # midPosAccum = verts[0].clone()
            # for vi in range(1, numVerts):
            #     with cython.nonecheck(False):
            #         midPosAccum = midPosAccum + verts[vi].pos
            # midPos = midPosAccum/float(numVerts)

            midNormal = None
            if verts[0].normal is not None:
                midNormal = poly.plane.normal
            midVert = Vertex(midPos, midNormal)

            newVerts = verts + \
                       [verts[i].interpolate(verts[(i + 1) % numVerts], 0.5) for i in range(numVerts)] + \
                       [midVert]

            i = 0
            vs = [newVerts[i], newVerts[i + numVerts], newVerts[2 * numVerts], newVerts[2 * numVerts - 1]]
            newPoly = Polygon(vs, poly.shared)
            newPoly.shared = poly.shared
            newPoly.plane = poly.plane
            newCSG.polygons.append(newPoly)

            for i in range(1, numVerts):
                vs = [newVerts[i], newVerts[numVerts + i], newVerts[2 * numVerts], newVerts[numVerts + i - 1]]
                newPoly = Polygon(vs, poly.shared)
                newCSG.polygons.append(newPoly)

        return newCSG

    cpdef void translate(self, list disp):
        """
        Translate Geometry.
           disp: displacement (array of floats)
        """
        cdef Vector d
        cdef Vertex v
        cdef int nverts
        cdef int npolys = len(self.polygons)
        cdef Py_ssize_t pi, vi

        d = Vector(disp[0], disp[1], disp[2])
        for pi in range(npolys):
            poly = self.polygon[pi]
            nverts = len(poly.vertices)
            for vi in range(nverts):
                v = poly.vertices[vi]
                v.pos = v.pos.plus(d)
                # no change to the normals

    cpdef void rotate(self, list axis, double angleDeg):
        """
        Rotate geometry.
           axis: axis of rotation (array of floats)
           angleDeg: rotation angle in degrees
        """
        cdef Vector ax
        cdef Vertex vert
        cdef Polygon poly
        cdef double cosAngle, sinAngle
        cdef int npolys = len(self.polygons)
        cdef Py_ssize_t pi, vi

        ax = Vector(axis[0], axis[1], axis[2]).unit()
        # cosAngle = math.cos(math.pi * angleDeg / 180.)
        # sinAngle = math.sin(math.pi * angleDeg / 180.)
        cosAngle = cos(PI * angleDeg / 180.)
        sinAngle = sin(PI * angleDeg / 180.)

        for pi in range(npolys):
            poly = self.polygon[pi]
            nverts = len(poly.vertices)
            for vi in range(nverts):
                vert = poly.vertices[vi]
                vert.pos = rotateVector(vert.pos, ax, cosAngle, sinAngle)
                if vert.normal.length() > 0:
                    vert.normal = rotateVector(vert.normal, ax, cosAngle, sinAngle)

    cpdef int toVerticesAndPolygons(self, list verts, list polys):
        """
        Return list of vertices, polygons (cells), and the total
        number of vertex indices in the polygon connectivity list
        (count).
        """
        
        cdef dict vertexIndexMap
        cdef double offset = 1.234567890
        cdef int npolys = len(self.polygons)
        cdef int count, index, nSortedVertexIndices, nWords
        cdef Polygon poly
        cdef Vertex v
        cdef Vector p
        cdef Py_ssize_t pi, vi, si, wi
        cdef list cell, sortedVertexIndex, _p, words
        cdef str vKey
        
        verts = []
        polys = []
        vertexIndexMap = {}
        count = 0
        for pi in range(npolys):
            poly = self.polygon[pi]
            verts = poly.vertices
            cell = []
            nverts = len(poly.vertices)
            for vi in range(nverts):
                v = poly.vertices[vi]
                p = v.pos
                # use string key to remove degeneracy associated
                # very close points. The format %.10e ensures that
                # points differing in the 11 digits and higher are
                # treated as the same. For instance 1.2e-10 and
                # 1.3e-10 are essentially the same.
                vKey = '%.10e,%.10e,%.10e' % (p.x + offset,
                                              p.y + offset,
                                              p.z + offset)
                if not vKey in vertexIndexMap:
                    vertexIndexMap[vKey] = len(vertexIndexMap)
                index = vertexIndexMap[vKey]
                cell.append(index)
                count += 1
            polys.append(cell)
        # sort by index
        sortedVertexIndex = sorted(vertexIndexMap.items(),
                                   key=operator.itemgetter(1))
        verts = []
        nSortedVertexIndices = len(sortedVertexIndex)
        for si in range(nSortedVertexIndices):
            v = sortedVertexIndex[si][0]
            i = sortedVertexIndex[si][1]
            _p = []
            words = v.split(',')
            nWords = len(words)
            for wi in range(nWords):
                c = words[wi]
                _p.append(float(c) - offset)
            verts.append(tuple(p))
        # return verts, polys, count
        return count

    def saveVTK(self, filename):
        """
        Save polygons in VTK file.
        """
        with open(filename, 'w') as f:
            f.write('# vtk DataFile Version 3.0\n')
            f.write('pycsg output\n')
            f.write('ASCII\n')
            f.write('DATASET POLYDATA\n')
            verts = []
            cells = []

            count = self.toVerticesAndPolygons(verts, cells)

            f.write('POINTS {0} float\n'.format(len(verts)))
            for v in verts:
                f.write('{0} {1} {2}\n'.format(v[0], v[1], v[2]))
            numCells = len(cells)
            f.write('POLYGONS {0} {1}\n'.format(numCells, count + numCells))
            for cell in cells:
                f.write('{0} '.format(len(cell)))
                for index in cell:
                    f.write('{0} '.format(index))
                f.write('\n')

    cpdef CSG union(self, CSG csg):
        """
        Return a new CSG solid representing space in either this solid or in the
        solid `csg`. Neither this solid nor the solid `csg` are modified.::

            A.union(B)

            +-------+            +-------+
            |       |            |       |
            |   A   |            |       |
            |    +--+----+   =   |       +----+
            +----+--+    |       +----+       |
                 |   B   |            |       |
                 |       |            |       |
                 +-------+            +-------+
        """
        cdef BSPNode a, b

        a = BSPNode(self.clone().polygons)
        b = BSPNode(csg.clone().polygons)
        a.clipTo(b)
        b.clipTo(a)
        b.invert()
        b.clipTo(a)
        b.invert()
        a.build(b.allPolygons());
        return csgFromPolygons(a.allPolygons())

    def __add__(self, CSG csg):
        return self.union(csg)

    cpdef CSG subtract(self, CSG csg):
        """
        Return a new CSG solid representing space in this solid but not in the
        solid `csg`. Neither this solid nor the solid `csg` are modified.::

            A.subtract(B)

            +-------+            +-------+
            |       |            |       |
            |   A   |            |       |
            |    +--+----+   =   |    +--+
            +----+--+    |       +----+
                 |   B   |
                 |       |
                 +-------+
        """
        cdef BSPNode a, b

        a = BSPNode(self.clone().polygons)
        b = BSPNode(csg.clone().polygons)
        a.invert()
        a.clipTo(b)
        b.clipTo(a)
        b.invert()
        b.clipTo(a)
        b.invert()
        a.build(b.allPolygons())
        a.invert()
        return csgFromPolygons(a.allPolygons())

    def __sub__(self, CSG csg):
        return self.subtract(csg)

    cpdef CSG intersect(self, CSG csg):
        """
        Return a new CSG solid representing space both this solid and in the
        solid `csg`. Neither this solid nor the solid `csg` are modified.::

            A.intersect(B)

            +-------+
            |       |
            |   A   |
            |    +--+----+   =   +--+
            +----+--+    |       +--+
                 |   B   |
                 |       |
                 +-------+
        """
        cdef BSPNode a, b

        a = BSPNode(self.clone().polygons)
        b = BSPNode(csg.clone().polygons)
        a.invert()
        b.clipTo(a)
        b.invert()
        a.clipTo(b)
        b.clipTo(a)
        a.build(b.allPolygons())
        a.invert()
        return csgFromPolygons(a.allPolygons())

    def __mul__(self, CSG csg):
        return self.intersect(csg)

    cpdef CSG inverse(self):
        """
        Return a new CSG solid with solid and empty space switched. This solid is
        not modified.
        """
        cdef CSG csg
        cdef int npolys
        cdef Py_ssize_t pi

        csg = self.clone()
        # map(lambda p: p.flip(), csg.polygons)
        npolys = len(csg.polygons)
        for pi in range(npolys):
            csg.polygons[pi].flip()

        return csg


def cube(center=[0, 0, 0], radius=[1, 1, 1]):
    """
    Construct an axis-aligned solid cuboid. Optional parameters are `center` and
    `radius`, which default to `[0, 0, 0]` and `[1, 1, 1]`. The radius can be
    specified using a single number or a list of three numbers, one for each axis.

    Example code::

        cube = CSG.cube(
          center=[0, 0, 0],
          radius=1
        )
    """
    cdef Vector c
    cdef list r, polygons

    c = Vector(0, 0, 0)
    r = [1, 1, 1]
    if isinstance(center, list): c = Vector(center)
    if isinstance(radius, list):
        r = radius
    else:
        r = [radius, radius, radius]

    polygons = list(map(
        lambda v: Polygon(
            list(map(lambda i:
                     Vertex(
                         Vector(
                             c.x + r[0] * (2 * bool(i & 1) - 1),
                             c.y + r[1] * (2 * bool(i & 2) - 1),
                             c.z + r[2] * (2 * bool(i & 4) - 1)
                         ),
                         None
                     ), v[0]))),
        [
            [[0, 4, 6, 2], [-1, 0, 0]],
            [[1, 3, 7, 5], [+1, 0, 0]],
            [[0, 1, 5, 4], [0, -1, 0]],
            [[2, 6, 7, 3], [0, +1, 0]],
            [[0, 2, 3, 1], [0, 0, -1]],
            [[4, 5, 7, 6], [0, 0, +1]]
        ]))
    return csgFromPolygons(polygons)


def sphere(**kwargs):
    """ Returns a sphere.

        Kwargs:
            center (list): Center of sphere, default [0, 0, 0].

            radius (float): Radius of sphere, default 1.0.

            slices (int): Number of slices, default 16.

            stacks (int): Number of stacks, default 8.
    """
    cdef list polygon, vertices, verticesN, verticesS, verticesE, verticesW
    cdef Vector c, d
    cdef int slices, stacks, j0, i0
    cdef float r, j1, j2, i1, dTheta, dPhi

    center = kwargs.get('center', [0.0, 0.0, 0.0])
    if isinstance(center, float):
        center = [center, center, center]
    c = Vector(*center)
    r = kwargs.get('radius', 1.0)
    # if isinstance(r, list) and len(r) > 2:
    #     r = r[0]
    slices = kwargs.get('slices', 16)
    stacks = kwargs.get('stacks', 8)
    polygons = []

    def appendVertex(vertices, theta, phi):
        # d = Vector(
        #     math.cos(theta) * math.sin(phi),
        #     math.cos(phi),
        #     math.sin(theta) * math.sin(phi))
        d = Vector(
            cos(theta) * sin(phi),
            cos(phi),
            sin(theta) * sin(phi))
        vertices.append(Vertex(c.plus(d.times(r)), d))

    dTheta = PI * 2.0 / float(slices)
    dPhi = PI / float(stacks)

    j0 = 0
    j1 = j0 + 1
    for i0 in range(0, slices):
        i1 = i0 + 1
        #  +--+
        #  | /
        #  |/
        #  +
        vertices = []
        appendVertex(vertices, i0 * dTheta, j0 * dPhi)
        appendVertex(vertices, i1 * dTheta, j1 * dPhi)
        appendVertex(vertices, i0 * dTheta, j1 * dPhi)
        polygons.append(Polygon(vertices, False))

    j0 = stacks - 1
    j1 = j0 + 1
    for i0 in range(0, slices):
        i1 = i0 + 1
        #  +
        #  |\
        #  | \
        #  +--+
        vertices = []
        appendVertex(vertices, i0 * dTheta, j0 * dPhi)
        appendVertex(vertices, i1 * dTheta, j0 * dPhi)
        appendVertex(vertices, i0 * dTheta, j1 * dPhi)
        polygons.append(Polygon(vertices, False))

    for j0 in range(1, stacks - 1):
        j1 = j0 + 0.5
        j2 = j0 + 1
        for i0 in range(0, slices):
            i1 = i0 + 0.5
            i2 = i0 + 1
            #  +---+
            #  |\ /|
            #  | x |
            #  |/ \|
            #  +---+
            verticesN = []
            appendVertex(verticesN, i1 * dTheta, j1 * dPhi)
            appendVertex(verticesN, i2 * dTheta, j2 * dPhi)
            appendVertex(verticesN, i0 * dTheta, j2 * dPhi)
            polygons.append(Polygon(verticesN, False))
            verticesS = []
            appendVertex(verticesS, i1 * dTheta, j1 * dPhi)
            appendVertex(verticesS, i0 * dTheta, j0 * dPhi)
            appendVertex(verticesS, i2 * dTheta, j0 * dPhi)
            polygons.append(Polygon(verticesS, False))
            verticesW = []
            appendVertex(verticesW, i1 * dTheta, j1 * dPhi)
            appendVertex(verticesW, i0 * dTheta, j2 * dPhi)
            appendVertex(verticesW, i0 * dTheta, j0 * dPhi)
            polygons.append(Polygon(verticesW, False))
            verticesE = []
            appendVertex(verticesE, i1 * dTheta, j1 * dPhi)
            appendVertex(verticesE, i2 * dTheta, j0 * dPhi)
            appendVertex(verticesE, i2 * dTheta, j2 * dPhi)
            polygons.append(Polygon(verticesE, False))

    return csgFromPolygons(polygons)


def cylinder(**kwargs):
    """ Returns a cylinder.

        Kwargs:
            start (list): Start of cylinder, default [0, -1, 0].

            end (list): End of cylinder, default [0, 1, 0].

            radius (float): Radius of cylinder, default 1.0.

            slices (int): Number of slices, default 16.
    """
    cdef int slices, i
    cdef double t0, i1, t1, dt
    cdef list polygons
    cdef Vector s, e, ray, axisZ, axisX, axisY, out, pos, normal
    cdef Vertex end

    _s = kwargs.get('start', Vector(0.0, -1.0, 0.0))
    _e = kwargs.get('end', Vector(0.0, 1.0, 0.0))
    if isinstance(_s, list):
        s = Vector(*_s)
    if isinstance(_e, list):
        e = Vector(*_e)
    r = kwargs.get('radius', 1.0)
    slices = kwargs.get('slices', 16)
    ray = e.minus(s)

    axisZ = ray.unit()
    # isY = (math.fabs(axisZ.y) > 0.5)
    isY = (abs(axisZ.y) > 0.5)
    axisX = Vector(float(isY), float(not isY), 0).cross(axisZ).unit()
    axisY = axisX.cross(axisZ).unit()
    start = Vertex(s, axisZ.negated())
    end = Vertex(e, axisZ.unit())
    polygons = []

    def point(stack, angle, normalBlend):
        out = axisX.times(cos(angle)).plus(
            axisY.times(sin(angle)))
        pos = s.plus(ray.times(stack)).plus(out.times(r))
        # normal = out.times(1.0 - math.fabs(normalBlend)).plus(
        normal = out.times(1.0 - abs(normalBlend)).plus(
            axisZ.times(normalBlend))
        return Vertex(pos, normal)

    dt = PI * 2.0 / float(slices)
    for i in range(0, slices):
        t0 = i * dt
        i1 = (i + 1) % slices
        t1 = i1 * dt
        polygons.append(Polygon([start.clone(),
                                 point(0., t0, -1.),
                                 point(0., t1, -1.)], False))
        polygons.append(Polygon([point(0., t1, 0.),
                                 point(0., t0, 0.),
                                 point(1., t0, 0.),
                                 point(1., t1, 0.)], False))
        polygons.append(Polygon([end.clone(),
                                 point(1., t1, 1.),
                                 point(1., t0, 1.)], False))

    return csgFromPolygons(polygons)


def cone(**kwargs):
    """ Returns a cone.

        Kwargs:
            start (list): Start of cone, default [0, -1, 0].

            end (list): End of cone, default [0, 1, 0].

            radius (float): Maximum radius of cone at start, default 1.0.

            slices (int): Number of slices, default 16.
    """
    
    cdef list polygons
    cdef double r, taperAngle, sinTaperAngle, cosTaperAngle, dt, t0, i1, t1
    cdef int slices, i
    cdef Vector s, e, ray, axisZ, axisX, axisY, startNormal, out, pos, normal, p0, n0, p1, n1, nAvg
    cdef Vertex start
    cdef Polygon polyStart, polySide

    _s = kwargs.get('start', Vector(0.0, -1.0, 0.0))
    _e = kwargs.get('end', Vector(0.0, 1.0, 0.0))
    if isinstance(_s, list):
        s = Vector(*_s)
    if isinstance(_e, list):
        e = Vector(*_e)
    r = kwargs.get('radius', 1.0)
    slices = kwargs.get('slices', 16)
    ray = e.minus(s)

    axisZ = ray.unit()
    # isY = (math.fabs(axisZ.y) > 0.5)
    isY = (abs(axisZ.y) > 0.5)
    axisX = Vector(float(isY), float(not isY), 0).cross(axisZ).unit()
    axisY = axisX.cross(axisZ).unit()
    startNormal = axisZ.negated()
    start = Vertex(s, startNormal)
    polygons = []

    # taperAngle = math.atan2(r, ray.length())
    taperAngle = atan2(r, ray.length())
    sinTaperAngle = sin(taperAngle)
    cosTaperAngle = cos(taperAngle)

    def point(angle):
        # radial direction pointing out
        out = axisX.times(cos(angle)).plus(
            axisY.times(sin(angle)))
        pos = s.plus(out.times(r))
        # normal taking into account the tapering of the cone
        normal = out.times(cosTaperAngle).plus(axisZ.times(sinTaperAngle))
        return pos, normal

    dt = PI * 2.0 / float(slices)
    for i in range(0, slices):
        t0 = i * dt
        i1 = (i + 1) % slices
        t1 = i1 * dt
        # coordinates and associated normal pointing outwards of the cone's
        # side
        p0, n0 = point(t0)
        p1, n1 = point(t1)
        # average normal for the tip
        nAvg = n0.plus(n1).times(0.5)
        # polygon on the low side (disk sector)
        polyStart = Polygon([start.clone(),
                             Vertex(p0, startNormal),
                             Vertex(p1, startNormal)], False)
        polygons.append(polyStart)
        # polygon extending from the low side to the tip
        polySide = Polygon([Vertex(p0, n0), Vertex(e, nAvg), Vertex(p1, n1)], False)
        polygons.append(polySide)

    return csgFromPolygons(polygons)
