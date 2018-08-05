"""
PLY writer that supports point normals
"""

HEADER = """ply
format ascii 1.0
comment GIAS2 PLYWriter generated
element vertex {nvertices}
property float x
property float y
property float z
{normal_block}{vcolour_block}element face {nfaces}
property list uchar int vertex_indices
end_header
"""

HEADER_NORMAL = """property float nx
property float ny
property float nz
"""

HEADER_VCOLOUR = """property uchar red
property uchar green
property uchar blue
"""

VERTEX_ROW_PAT = "{vcoord}{vnormal}{vcolour}\n"
VCOORD_PAT = "{:10.8f} {:10.8f} {:10.8f}"
VNORMAL_PAT = " {:10.8f} {:10.8f} {:10.8f}"
VCOLOUR_PAT = " {:d} {:d} {:d}"
FACE_PAT = "{} {} {} {}\n"

class PLYWriter(object):

	def __init__(self, v, f, filename=None, vn=None, vcolours=None):
		self.v = v
		self.f = f
		self.vnormals = vn
		self.vcolours = vcolours
		self.filename = filename

	def write(self, filename=None):
		if filename is None:
			filename = self.filename

		with open(filename, 'w') as file_:

			# write header block
			if self.vnormals is not None:
				header_normal_block = HEADER_NORMAL
			else:
				header_normal_block = ''

			if self.vcolours is not None:
				header_vcolour_block = HEADER_VCOLOUR
			else:
				header_vcolour_block = ''	

			header_block = HEADER.format(
				nvertices=len(self.v),
				nfaces=len(self.f),
				normal_block=header_normal_block,
				vcolour_block=header_vcolour_block,
				)

			file_.write(header_block)

			# write vertex coord, normal, and colour
			if (self.vnormals is not None) and (self.vcolours is not None):
				vrow_pat = "{vcoord}{vnormal}{vcolour}\n"
				for vi, v in enumerate(self.v):
					vcoord = VCOORD_PAT.format(*v)
					vnormal = VNORMAL_PAT.format(*self.vnormals[vi])
					vcolour = VCOLOUR_PAT.format(*self.vcolours[vi])
					file_.write(vrow_pat.format(
						vcoord=vcoord, vnormal=vnormal, vcolour=vcolour
						))
			# write vertex coord, normal
			elif (self.vnormals is not None) and (self.vcolours is None):
				vrow_pat = "{vcoord}{vnormal}\n"
				for vi, v in enumerate(self.v):
					vcoord = VCOORD_PAT.format(*v)
					vnormal = VNORMAL_PAT.format(*self.vnormals[vi])
					file_.write(vrow_pat.format(vcoord=vcoord, vnormal=vnormal))
			# write vertex coord, colour
			elif (self.vnormals is None) and (self.vcolours is not None):
				vrow_pat = "{vcoord}{vcolour}\n"
				for vi, v in enumerate(self.v):
					vcoord = VCOORD_PAT.format(*v)
					vcolour = VCOLOUR_PAT.format(*self.vcolours[vi])
					file_.write(vrow_pat.format(vcoord=vcoord, vcolour=vcolour))
			# write vertex coord
			elif (self.vnormals is None) and (self.vcolours is None):
				vrow_pat = "{vcoord}\n"
				for vi, v in enumerate(self.v):
					vcoord = VCOORD_PAT.format(*v)
					file_.write(vrow_pat.format(vcoord=vcoord))

			# write faces
			for fi, f in enumerate(self.f):
				file_.write(FACE_PAT.format(3, *f))