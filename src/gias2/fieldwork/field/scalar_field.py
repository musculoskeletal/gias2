"""
FILE: scalar_field.py
LAST MODIFIED: 24-12-2015 
DESCRIPTION: class for a ensemble_field_function representing a scalar field.

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""

import copy
import logging
import shelve

import itertools
import numpy
import sys
from scipy import sparse
from scipy.optimize import leastsq
from scipy.spatial import cKDTree

from gias2.fieldwork.field import ensemble_field_function as EFF

log = logging.getLogger(__name__)


def loadScalarField(filename, ensFilename=None, meshFilename=None, path=None):
    if path is not None:
        filename = path + filename

    try:
        S = shelve.open(filename, 'r')
    except:
        raise IOError(filename + ' not found')
    else:
        if ensFilename:
            F = EFF.load_ensemble(ensFilename, meshFilename=meshFilename, path=path)
        else:
            try:
                F = EFF.load_ensemble(S['ensemble_field'], meshFilename=meshFilename, path=path)
            except:
                raise ValueError('no ensemble field function defined')

        sf = scalarField(S['name'], F)
        sf.setFieldParameters(S['field_parameters'])
        return sf


class scalarField(object):
    def __init__(self, name, ensembleFieldFunction, parameters=None):
        self.name = name
        self.ensembleFieldFunction = ensembleFieldFunction
        if parameters is not None:
            self.setFieldParameters(parameters)
        else:
            self.parameters = None

    # ==================================================================#
    def saveScalarField(self, filename, fieldFilename=None, meshFilename=None):
        """ only saves the field and field parameters
        """
        S = shelve.open(filename + '.scf')

        S['name'] = self.name
        if fieldFilename:
            S['ensemble_field'] = self.ensembleFieldFunction.save_ensemble(fieldFilename, mesh_filename=meshFilename)
        # ~ else:
        # ~ S['ensemble_field'] = self.ensemble_field_function.save_ensemble( self.ensemble_field_function.name, mesh_filename=mesh_filename )
        else:
            S['ensemble_field'] = None

        S['field_parameters'] = self.fieldParameters

        S.close()
        return filename + '.scf'

    # ==================================================================#
    def setFieldParameters(self, p):
        nEnsemblePoints = self.ensembleFieldFunction.get_number_of_ensemble_points()
        if len(p) != nEnsemblePoints:
            raise ValueError('mismatch in number of points {} and parameters {}'.format(
                self.ensembleFieldFunction.get_number_of_ensemble_points(), len(p)))
        else:
            # if nothings wrong so far
            self.fieldParameters = numpy.array(p, dtype=float)

    # ==================================================================#
    def getFieldParameters(self):
        return self.fieldParameters.copy()

    # ==================================================================#
    def evaluateField(self, density, derivs=None):
        """ evaluates the field for all parameter components.
        Returns a list of self.dimension lists
        """
        if not derivs:
            V = self.ensembleFieldFunction.evaluate_field_in_mesh(density, parameters=self.fieldParameters)
            return V
        else:
            V, D = self.ensembleFieldFunction.evaluate_field_in_mesh(density, parameters=self.fieldParameters,
                                                                     derivs=derivs)
            return V, D

    # ==================================================================#
    def evaluateFieldAtElementPoints(self, element, XI):
        V = self.ensembleFieldFunction.evaluate_field_at_element_point(element, XI, parameters=self.fieldParameters)
        return V

    def evaluateFieldAtMaterialPoints(self, matpoints):
        """
        evaluate in multiple elements. matpoints is list of element points expressed
        as a tuple (element_number, xi_coordinates). Example of input EP:
        [(0, [0.5, 0.1]), (1, [0.1, 0.1]),]
        """
        evaluator = makeScalarFieldEvaluatorSparse(self, [1, 1], matPoints=matpoints)
        X = evaluator(self.fieldParameters)
        return X

    def getFieldParameters(self):
        return self.fieldParameters.copy()

    def flattenEnsembleFieldFunction(self):
        self.ensembleFieldFunctionOld = copy.deepcopy(self.ensembleFieldFunction)
        flatEFF = self.ensembleFieldFunction.flatten()[0]
        self.ensembleFieldFunction = flatEFF


# ======================================================================#
def makeScalarFieldEvaluatorSparse(S, evalD, matPoints=None):
    """ create a function for evaluation the scalar field values,
    taking advantage of a precomputed sparse matrix of basis function
    values at fixed element coordinates
    
    This is about 10x faster than the same setup with a dense matrix
    implementation. csc is a little bit faster than csr sparse matrix.
    """

    f = S.ensembleFieldFunction
    if not f.is_flat():
        f = f.flatten()[0]
    ep = S.evaluateField(evalD)

    if matPoints is not None:
        nEPs = len(matPoints)
    else:
        ep = S.evaluate_geometric_field(evalD)
        nEPs = ep.shape[1]

    A = numpy.zeros((nEPs, f.get_number_of_ensemble_points()), dtype=float)

    if matPoints is not None:
        # ~ pdb.set_trace()
        elemEnsNodes = {}
        for mpI, (elem, xi) in enumerate(matPoints):
            element = f.mesh.elements[elem]
            b = f.basis[element.type].eval(xi)
            ensNodes = elemEnsNodes.get(elem)
            if ensNodes is None:
                emap = f.mapper._element_to_ensemble_map[elem]
                ensNodes = [emap[k][0][0] for k in list(emap.keys())]
                elemEnsNodes[elem] = ensNodes

            A[mpI, ensNodes] = b
    else:
        # calculate static basis values for the required evalD and assemble
        # matrix
        basisValues = {}
        A = numpy.zeros((ep.shape[0], f.get_number_of_ensemble_points()), dtype=float)
        row = 0
        for elementNumber in numpy.sort(list(f.mesh.elements.keys())):
            # get element
            element = f.mesh.elements[elementNumber]
            # calculate basis values
            if basisValues.get(element.type) == None:
                evalGrid = element.generate_eval_grid(evalD).squeeze()
                basisValues[element.type] = f.basis[element.type].eval(evalGrid.T).T

            # fill in A matrix
            b = basisValues[element.type]  # basis values
            emap = f.mapper._element_to_ensemble_map[elementNumber]  # element to ensemble map
            for n in range(b.shape[1]):
                A[row:row + b.shape[0], emap[n][0][0]] = b[:, n]

            row += b.shape[0]

    # sparsify A
    As = sparse.csc_matrix(A)

    def evaluator(P):
        E = As * P
        return E

    return evaluator


def makeScalarFieldDerivativesEvaluatorSparse(S, evalD):
    """ create a function for evaluating the scalar field derivatives,
    taking advantage of a precomputed sparse matrix of basis function
    values at fixed element coordinates
    
    This is about 10x faster than the same setup with a dense matrix
    implementation. csc is a little bit faster than csr sparse matrix.
    """

    f = S.ensembleFieldFunction.flatten()[0]
    ep, epd = S.evaluateField(evalD, derivs=-1)

    # calculate static basis values for the required evalD and assemble
    # matrices
    basisValues = {}
    A = [numpy.zeros((ep.shape[0], f.get_number_of_ensemble_points()), dtype=float) for i in range(epd.shape[0])]
    row = 0
    for elementNumber in numpy.sort(list(f.mesh.elements.keys())):
        # get element
        element = f.mesh.elements[elementNumber]
        # calculate basis values
        if basisValues.get(element.type) == None:
            evalGrid = element.generate_eval_grid(evalD)
            basisValues[element.type] = f.basis[element.type].eval_derivatives(evalGrid.T, None)

        # fill in A matrix
        b = basisValues[element.type].copy()  # basis values dongdong
        emap = f.mapper._element_to_ensemble_map[elementNumber]  # element to ensemble map
        for d in range(b.shape[0]):
            for n in range(b.shape[1]):
                A[d][row:row + b.shape[2], emap[n][0][0]] = b[d, n, :]

        row += b.shape[2]

    # sparsify matrices in A
    # ~ As = [ sparse.csc_matrix(a) for a in A ]
    # ~ dx1, dx2, ddx1, ddx2, dx1x2 = As

    # ~ def evaluatorOld( P ):
    # ~ Pd = P.reshape((dim,-1)).T
    # ~
    # ~ Dx1 = (dx1*Pd).T
    # ~ Dx2 = (dx2*Pd).T
    # ~ DDx1 = (ddx1*Pd).T
    # ~ DDx2 = (ddx2*Pd).T
    # ~ Dx1x2 = (dx1x2*Pd).T
    # ~
    # ~ return numpy.array([Dx1,Dx2,DDx1,DDx2,Dx1x2]).swapaxes(0,1)

    # stack all derivative A matrices
    AStackedSparse = sparse.csc_matrix(numpy.vstack(A))
    nDerivs = len(A)

    def evaluator(P):
        """ uses a A matrix that is the vstack of all derivative A matrices
        """

        D = AStackedSparse * P
        return D.reshape((nDerivs, -1))

    return evaluator


# ======================================================================#
def norm(x):
    return x / numpy.sqrt((x ** 2.0).sum(1))[:, numpy.newaxis]


class normalSmoother(object):
    """ for 2D elements
    """

    def __init__(self, F):

        self.F = F
        self.en2el, self.el2en = self.F.get_mapping()  # ens2elem, elem2ens maps
        self._procMap()

    def _procMap(self):
        """ find element points of shared nodes
        """
        self.edgePoints = []
        self.vertexPoints = []  # not implemented

        # for each ensemble point
        for p in list(self.en2el.keys()):
            # edge nodes
            if len(self.en2el[p]) == 2:
                E = list(self.en2el[p].items())  # element points mapped to by ensemble point p
                self.edgePoints.append(((E[0][0], list(E[0][1].keys())[0]), (
                    E[1][0], list(E[1][1].keys())[0])))  # ((element1, elementnode1), (element2, elementnode2))
            # vertex nodes
            elif len(self.en2el[p]) > 2:
                pass

        return

    def _procEdge(self, D):
        """ for each pair in edgePoints, find the element and edge shared
        and generate eval points along the shared edges
        """
        self.edgeEvalBasis = []
        self.edgeEvalPoints = []
        self.nPairs = 0  # number of pairs of edge points
        doneEdges = []  # list of element edges that have been done
        # for each pair
        for ep1, ep2 in self.edgePoints:
            # get elements
            e1 = self.F.mesh.elements[ep1[0]]
            e2 = self.F.mesh.elements[ep2[0]]
            # get edges for element points
            edgeList1 = e1.get_point_edge(ep1[1])
            edgeList2 = e2.get_point_edge(ep2[1])

            # ignore corners
            if len(edgeList1) > 1 or len(edgeList2) > 2:
                pass
            else:
                (ei1, pi1, edge1) = edgeList1[0]
                (ei2, pi2, edge2) = edgeList2[0]

                # ignore repeat element edges
                if ((ep1[0], ei1) in doneEdges) and ((ep2[0], ei2)) in doneEdges:
                    pass
                else:
                    # check if element edges are reversed
                    ### assumes no hanging nodes and same number of nodes per edge!!! ###
                    if pi1 != pi2:
                        eval1 = edge1.get_elem_coord(numpy.linspace(0.0, 1.0, D))
                        eval2 = edge2.get_elem_coord(numpy.linspace(1.0, 0.0, D))
                    else:
                        eval1 = edge1.get_elem_coord(numpy.linspace(0.0, 1.0, D))
                        eval2 = edge2.get_elem_coord(numpy.linspace(0.0, 1.0, D))

                    # get basis values for these
                    basis1 = [self.F.basis[e1.type].eval_derivatives(eval1.T, d).T for d in ((1, 0), (0, 1))]
                    basis2 = [self.F.basis[e2.type].eval_derivatives(eval2.T, d).T for d in ((1, 0), (0, 1))]

                    self.edgeEvalPoints.append((ep1[0], eval1, ep2[0], eval2))  # ( element1, ep1, element2, ep2 )
                    self.edgeEvalBasis.append(
                        (ep1[0], basis1, ep2[0], basis2))  # ( element1, basis1, element2, basis2 )
                    doneEdges.append((ep1[0], ei1))
                    doneEdges.append((ep2[0], ei2))
                    self.nPairs += eval1.shape[0]

        return

    def makeObj(self, D):
        """ make a lagrange multiplier element edge smoothing objective 
        function with each edge discretised at D
        """

        # calculate edge point basis values
        self._procEdge(D)

        # make derivatives evaluation matrices
        A1dxi1 = numpy.zeros((self.nPairs, self.F.get_number_of_ensemble_points()), dtype=float)
        A1dxi2 = numpy.zeros((self.nPairs, self.F.get_number_of_ensemble_points()), dtype=float)
        A2dxi1 = numpy.zeros((self.nPairs, self.F.get_number_of_ensemble_points()), dtype=float)
        A2dxi2 = numpy.zeros((self.nPairs, self.F.get_number_of_ensemble_points()), dtype=float)
        row1 = 0
        row2 = 0

        for e1, b1, e2, b2 in self.edgeEvalBasis:
            emap1 = self.el2en[e1]
            emap2 = self.el2en[e2]

            # make dxi1 and dxi2 matrices for points on one side
            for n in range(b1[0].shape[1]):
                A1dxi1[row1:row1 + b1[0].shape[0], emap1[n][0][0]] = b1[0][:, n]
                A1dxi2[row1:row1 + b1[1].shape[0], emap1[n][0][0]] = b1[1][:, n]

            row1 += b1[0].shape[0]

            # make dxi1 and dxi2 matrices for points on the other side
            for n in range(b2[0].shape[1]):
                A2dxi1[row2:row2 + b2[0].shape[0], emap2[n][0][0]] = b2[0][:, n]
                A2dxi2[row2:row2 + b2[1].shape[0], emap2[n][0][0]] = b2[1][:, n]

            row2 += b2[0].shape[0]

        # sparsify matrices
        sA1dxi1 = sparse.csc_matrix(A1dxi1)
        sA1dxi2 = sparse.csc_matrix(A1dxi2)
        sA2dxi1 = sparse.csc_matrix(A2dxi1)
        sA2dxi2 = sparse.csc_matrix(A2dxi2)

        V1 = numpy.ones([self.nPairs, 3], dtype=float)
        V2 = numpy.ones([self.nPairs, 3], dtype=float)

        def obj(x):

            # evaluate normal on one side
            # evaluate  dxi1
            d1dxi1 = sA1dxi1 * x
            # evaluate  dxi2
            d1dxi2 = sA1dxi2 * x
            # cross product and normalise
            V1[:, 0] = -d1dxi1
            V1[:, 1] = -d1dxi2
            n1 = norm(V1)

            # evaluation normal of the other side
            # evaluate  dxi1
            d2dxi1 = sA2dxi1 * x
            # evaluate  dxi2
            d2dxi2 = sA2dxi2 * x
            # cross product and normalise
            V2[:, 0] = -d2dxi1
            V2[:, 1] = -d2dxi2
            n2 = norm(V2)

            # dot product normals
            err = 1.0 - (n1[:, 0] * n2[:, 0] + n1[:, 1] * n2[:, 1] + n1[:, 2] * n2[:, 2])

            return err

        return obj


# ======================================================================#
def fitScalarFieldEPDP(SF, data, evalD, ftol=None, xtol=None, epsfcn=None, maxit=None):
    if maxit:
        maxfev = maxit * SF.ensembleFieldFunction.get_number_of_ensemble_points()
    else:
        maxfev = 0

    evaluator = makeScalarFieldEvaluatorSparse(SF, evalD)

    # ~ it = 0
    def obj(x):
        v = evaluator(x)
        d = data - v
        d *= d
        # ~ print 'it:', it, ' rms:', numpy.sqrt( d.mean() )
        log.debug('rms:', numpy.sqrt(d.mean()))
        # ~ it+=1
        return d * d

    pOpt = leastsq(obj, SF.fieldParameters.copy().ravel(), ftol=ftol, xtol=xtol, epsfcn=epsfcn, maxfev=maxfev)[0][:,
           numpy.newaxis]
    SF.setFieldParameters(pOpt)
    return pOpt, SF


# ======================================================================#
def fitScalarFieldDPEP(SF, data, evalD, dpepi, ftol=None, xtol=None, epsfcn=None, maxit=None):
    if maxit:
        maxfev = maxit * SF.ensembleFieldFunction.get_number_of_ensemble_points()
    else:
        maxfev = 0

    evaluator = makeScalarFieldEvaluatorSparse(SF, evalD)

    # ~ it = 0
    def obj(x):
        v = evaluator(x)[dpepi]
        d = data - v
        d *= d
        # ~ print 'it:', it, ' rms:', numpy.sqrt( d.mean() )
        log.debug('rms:', numpy.sqrt(d.mean()))
        # ~ it+=1
        return d

    pOpt = leastsq(obj, SF.fieldParameters.copy().ravel(), ftol=ftol, xtol=xtol, epsfcn=epsfcn, maxfev=maxfev)[0][:,
           numpy.newaxis]
    SF.setFieldParameters(pOpt)
    return pOpt, SF


# ======================================================================#

class fitGeomFieldToScalar(object):
    def __init__(self, G, dataGeom, dataScalar):
        self.G = G
        self.dataG = dataGeom
        self.dataS = dataScalar
        self.xtol = 1e-8
        self.ftol = 1e-8
        self.epsfcn = 0
        self.maxIt = 5

    def initialiseScalarField(self, fieldName=None, p0=None):
        # get initialise values for curvField. At each node, find curvature of closest datapoint
        if p0 == None:
            geoNodeCoord = self.G.field_parameters.squeeze().T
            self.dataGTree = cKDTree(self.dataG)
            closestDI = self.dataGTree.query(list(geoNodeCoord))[1]
            p0 = self.dataS[closestDI][:, numpy.newaxis]

        # create a new field curvField to fit on the mesh ensemble
        self.SF = scalarField(fieldName, self.G.ensemble_field_function, parameters=p0)
        return self.SF

    def fitEPDP(self, evalD):
        # = EPDP FIT ===========================================================#
        # for each ep find curvature of closest datapoint
        ep = self.G.evaluate_geometric_field(evalD).T
        epdpi = self.dataGTree.query(list(ep))[1]
        data = self.dataS[epdpi]

        # visualise initial curvature field
        # ~ SFValue = self.SF.evaluate_field( evalD )
        # ~ G.display_geometric_field( dispD, scalar=curvFieldValue, data=dispData.T )
        # ~ initialErrorField =  curvField.evaluate_field( evalD ) - data
        # ~ G.display_geometric_field( evalD, scalar=initialErrorField )

        # fit fC to data curvature
        pOpt, SFFit = fitScalarFieldEPDP(self.SF, self.dataS, evalD, ftol=self.ftol, xtol=self.xtol, epsfcn=self.epsfcn,
                                         maxit=self.maxIt)
        return SFFit

    def fitDPEP(self, evalD):
        # = DPEP FIT ===========================================================#
        # for each dp find curvature of closest ep
        ep = self.G.evaluate_geometric_field(evalD).T
        epTree = cKDTree(ep)
        dpepi = epTree.query(list(self.dataG))[1]

        # visualise initial curvature field
        # ~ curvFieldValue = curvField.evaluate_field( dispD )
        # ~ G.display_geometric_field( dispD, scalar=curvFieldValue, data=dispData.T )
        # ~ initialErrorField =  curvField.evaluate_field( evalD )[dpepi] - data
        # ~ G.display_geometric_field( evalD, scalar=initialErrorField )

        # fit fC to data curvature
        pOpt, SFFit = fitScalarFieldDPEP(self.SF, self.dataS, evalD, dpepi, ftol=self.ftol, xtol=self.xtol,
                                         epsfcn=self.epsfcn, maxit=self.maxIt)
        return SFFit


# ======================================================================#
def makeSobelovPenalty2D(S, evalD, w):
    SDEval = makeScalarFieldDerivativesEvaluatorSparse(S, evalD)

    def obj(p):
        D1, D2, D11, D22, D12 = SDEval(p)

        S = numpy.sqrt(w[0] * D1 * D1 + w[1] * D11 * D11 + w[2] * D2 * D2 + w[3] * D22 * D22 + w[4] * D12 * D12)

        return S

    return obj


# ======================================================================#
def makeMaskedDataObj(S, evalD, maskedData, sobW=None, sobD=None, nW=None, nD=None):
    evaluator = makeScalarFieldEvaluatorSparse(S, evalD)
    c = itertools.count(0)

    if (sobW == None) and (nW == None):
        # no sobelov or normal smoothing
        def obj(x):
            v = evaluator(x)
            d = maskedData - v
            sys.stdout.write('\rit: %(it)05i rms: %(rms)8.6f' % {'it': next(c), 'rms': numpy.sqrt((d * d).mean())})
            sys.stdout.flush()
            return d

        return obj
    else:
        sobObj = makeSobelovPenalty2D(S, sobD, sobW)
        if nW == None:
            # only sobelov smoothing
            def obj(x):
                v = evaluator(x)
                dData = abs(maskedData - v)
                dSob = sobObj(x)
                d = numpy.hstack([dData, dSob])
                sys.stdout.write(
                    '\rit: %(it)05i rms: %(rms)8.6f' % {'it': next(c), 'rms': numpy.sqrt((dData * dData).mean())})
                sys.stdout.flush()
                return d

            return obj
        else:
            # sobelov and normal smoothing
            nObjMaker = normalSmoother(S.ensembleFieldFunction.flatten()[0])
            nObj = nObjMaker.makeObj(nD)

            def obj(x):
                v = evaluator(x)
                dData = abs(maskedData - v)
                # ~ d = numpy.hstack( [d*d, sobObj(x)] )
                dSob = sobObj(x)
                dN = nObj(x)
                d = numpy.hstack([dData, dSob, dN])
                sys.stdout.write(
                    '\rit: %(it)05i rms: %(rms)8.6f' % {'it': next(c), 'rms': numpy.sqrt((dData * dData).mean())})
                sys.stdout.flush()
                return d

            return obj

    return obj


def fitGeomFieldToMaskedScalar(G, GD, maskedData, fieldName=None, sobW=None, sobD=None, nW=None, nD=None, xtol=1e-9):
    """ given a geometricField, creates from it a scalar field, and fits
    it to the data provided. Data is masked by dataMask (true for valid
    values, false for values to be masked
    """
    fieldFitter = fitGeomFieldToScalar(G, G.evaluate_geometric_field(GD).T, maskedData)
    field = fieldFitter.initialiseScalarField(fieldName=fieldName)
    fieldObj = makeMaskedDataObj(field, GD, maskedData, sobW=sobW, sobD=sobD, nW=nW, nD=nD)
    fieldPOpt = leastsq(fieldObj, field.getFieldParameters().ravel(), xtol=xtol)
    field.setFieldParameters(fieldPOpt[0][:, numpy.newaxis])
    fittedRMS = numpy.sqrt(fieldObj(fieldPOpt[0]).mean())
    return field, fieldPOpt, fittedRMS
