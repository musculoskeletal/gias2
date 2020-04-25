import logging

from numpy import ndarray, zeros, array

log = logging.getLogger(__name__)


class PolydataReader:
    def __init__(self, filename: str):
        self.fileName = filename
        self.file = None
        self.fileEnd = None
        self.numPoints = None
        self.numFaces = None
        self.points = None
        self.faces = None

    def loadData(self) -> None:
        log.debug("opening %s", self.fileName)
        self.file = open(self.fileName, "r")
        self.file.seek(0, 2)
        self.fileEnd = self.file.tell()

        # determine number of points in the data
        self.file.seek(0)
        found_num_points = False
        while (not found_num_points) and (self.file.tell() < self.fileEnd):
            line = self.file.readline()
            if line.find("POINTS") > -1:
                found_num_points = True
                line = line.split()
                self.numPoints = int(line[1])
                log.debug(str(self.numPoints) + " points in file")

        if not found_num_points:
            log.debug("Unable to find number of points!")

            # determine number of faces in the data
        self.file.seek(0)
        found_num_faces = False
        while (not found_num_faces) and (self.file.tell() < self.fileEnd):
            line = self.file.readline()
            if line.find("POLYGONS") > -1:
                found_num_faces = True
                line = line.split()
                self.numFaces = int(line[1])
                log.debug("%s faces in file", self.numFaces)

        if not found_num_faces:
            log.debug("Unable to find number of Faces!")

    def getPoints(self) -> ndarray:

        try:
            self.numPoints
        except AttributeError:
            log.debug("number of points unknown. Run .loadData() first")
        else:
            self.file.seek(0)
            found_points = False
            while (not found_points) and (self.file.tell() < self.fileEnd):

                # find points section header
                line = self.file.readline()
                if line.find("POINTS") > -1:
                    found_points = True
                    log.debug("getting points...")
                    self.points = zeros([self.numPoints, 3], dtype=float)
                    point_counter = 0

                    while point_counter < self.numPoints:
                        point_line = array(self.file.readline().split(), dtype=float)

                        for a in range(0, len(point_line) / 3):
                            self.points[point_counter, :] = point_line[0 + 3 * a:3 + 3 * a]
                            point_counter += 1

                    log.debug("Got %s points", self.points.shape[0])
                    return self.points

            if not found_points:
                log.debug("No points found!")

    def getPointNormals(self) -> ndarray:

        try:
            self.numPoints
        except AttributeError:
            log.debug("number of points unknown. Run .loadData() first")
        else:
            self.file.seek(0)
            found_normals = False
            while (not found_normals) and (self.file.tell() < self.fileEnd):

                # find points section header
                line = self.file.readline()
                if line.find("NORMALS") > -1:
                    found_normals = True
                    log.debug("getting point normals...")
                    self.pointNormals = zeros([self.numPoints, 3], dtype=float)
                    normal_counter = 0

                    while normal_counter < self.numPoints:
                        normal_line = array(self.file.readline().split(), dtype=float)

                        for a in range(0, len(normal_line) / 3):
                            self.pointNormals[normal_counter, :] = normal_line[0 + 3 * a:3 + 3 * a]
                            normal_counter += 1

                    log.debug("Got %s point normals", self.pointNormals.shape[0])
                    return self.pointNormals

            if not found_normals:
                log.debug("No point Normals found!")

    def getFaces(self):
        try:
            self.numFaces
        except AttributeError:
            log.debug("number of faces unknown. Run .loadData() first")
        else:
            self.file.seek(0)
            foundFaces = False
            while (foundFaces == False) and (self.file.tell() < self.fileEnd):

                # find points section header
                line = self.file.readline()
                if line.find("POLYGONS") > -1:
                    foundFaces = True
                    log.debug("getting faces...")
                    self.faces = zeros([self.numFaces, 3], dtype=int)
                    faceCounter = 0

                    while faceCounter < self.numFaces:
                        faceLine = array(self.file.readline().split(), dtype=int)

                        self.faces[faceCounter, :] = faceLine[1:4]
                        faceCounter += 1

                    log.debug("Got " + str(self.faces.shape[0]) + " faces")
                    return self.faces

            if foundFaces == False:
                log.debug("No faces found!")

    def getCurvature(self):

        try:
            self.numPoints
        except AttributeError:
            log.debug("number of points unknown. Run .loadData() first")
        else:
            self.file.seek(0)
            foundCurv = False
            while (foundCurv == False) and (self.file.tell() < self.fileEnd):

                # find curvature section header
                line = self.file.readline()
                if line.find("Curvature") > -1:
                    foundCurv = True
                    log.debug("Getting curvature values...")

                    self.curv = zeros([self.numPoints], dtype=float)
                    curvCounter = 0
                    self.file.readline()

                    while curvCounter < self.numPoints:
                        curvLine = array(self.file.readline().split(), dtype=float)
                        for i in curvLine:
                            self.curv[curvCounter] = i
                            curvCounter += 1

                    self.curvMean = self.curv.mean()
                    self.curvSD = self.curv.std()
                    log.debug("Got " + str(self.curv.shape[0]) + " curvature values")
                    return self.curvMean, self.curvSD

            if foundCurv == False:
                log.debug("No curvature values found!")

    def getEdgePoints(self, sd):

        # gets datapoints with local curvature greater than sd standard
        # deviations away from the mean

        try:
            self.points
        except AttributeError:
            self.getPoints()
            self.getCurvature()
            self.getEdgePoints(sd)
        else:
            self.edgePoints = []
            limit = sd * self.curvSD
            counter = 0

            for i in range(0, self.numPoints):
                if abs(self.curv[i] - self.curvMean) > limit:
                    self.edgePoints.append(self.points[i])
                    counter += 1

            log.debug(str(counter) + " edge points found at " + str(sd) + " SD")
            return self.edgePoints
