"""
FILE: plywriter.py
LAST MODIFIED: 05-08-2018 
DESCRIPTION: Class for writing PLY files

===============================================================================
This file is part of GIAS2. (https://bitbucket.org/jangle/gias2)

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
===============================================================================
"""
from typing import List, Optional

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

# VERTEX_ROW_PAT = "{vcoord}{vnormal}{vcolour}\n"
VCOORD_PAT = "{:10.8f} {:10.8f} {:10.8f}"
VNORMAL_PAT = " {:10.8f} {:10.8f} {:10.8f}"
VCOLOUR_PAT = " {:d} {:d} {:d}"
FACE_PAT = "{} {} {} {}\n"

VCOORD_PAT_BIN = b"%10.8f %10.8f %10.8f"
VNORMAL_PAT_BIN = b" %10.8f %10.8f %10.8f"
VCOLOUR_PAT_BIN = b" %d %d %d"
FACE_PAT_BIN = b"%d %d %d %d\n"


def _write_ascii(f, s):
    f.write(s)


def _write_bin(f, s):
    f.write(s.encode('utf-8'))


class PLYWriter(object):

    def __init__(
            self,
            v: List[List[float]],
            f: List[List[int]],
            filename: Optional[str] = None,
            vn: Optional[List[List[float]]] = None,
            vcolours: Optional[List[List[float]]] = None,
            ascenc: bool = True):
        self.v = v
        self.f = f
        self.vnormals = vn
        self.vcolours = vcolours
        self.filename = filename
        self.ascenc = ascenc

    def _header_block(self) -> str:
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
        return header_block

    def write(self, filename: Optional[str] = None) -> None:
        if filename is None:
            filename = self.filename

        # define ascii or binary writer
        if self.ascenc:
            write_func = _write_ascii
            file_mode = 'w'
        else:
            raise RuntimeError('Binary write not supported')
        # write_func = _write_bin
        # file_mode = 'wb'

        # write file
        with open(filename, file_mode) as file_:
            # writer header
            header_block = self._header_block()
            write_func(file_, header_block)

            # write vertex coord, normal, and colour
            if (self.vnormals is not None) and (self.vcolours is not None):
                vrow_pat = "{vcoord}{vnormal}{vcolour}\n"
                for vi, v in enumerate(self.v):
                    vcoord = VCOORD_PAT.format(*v)
                    vnormal = VNORMAL_PAT.format(*self.vnormals[vi])
                    vcolour = VCOLOUR_PAT.format(*self.vcolours[vi])
                    row_str = vrow_pat.format(
                        vcoord=vcoord, vnormal=vnormal, vcolour=vcolour
                    )
                    write_func(file_, row_str)
            # write vertex coord, normal
            elif (self.vnormals is not None) and (self.vcolours is None):
                vrow_pat = "{vcoord}{vnormal}\n"
                for vi, v in enumerate(self.v):
                    vcoord = VCOORD_PAT.format(*v)
                    vnormal = VNORMAL_PAT.format(*self.vnormals[vi])
                    row_str = vrow_pat.format(vcoord=vcoord, vnormal=vnormal)
                    write_func(file_, row_str)
            # write vertex coord, colour
            elif (self.vnormals is None) and (self.vcolours is not None):
                vrow_pat = "{vcoord}{vcolour}\n"
                for vi, v in enumerate(self.v):
                    vcoord = VCOORD_PAT.format(*v)
                    vcolour = VCOLOUR_PAT.format(*self.vcolours[vi])
                    row_str = vrow_pat.format(vcoord=vcoord, vcolour=vcolour)
                    write_func(file_, row_str)
            # write vertex coord
            elif (self.vnormals is None) and (self.vcolours is None):
                vrow_pat = "{vcoord}\n"
                for vi, v in enumerate(self.v):
                    vcoord = VCOORD_PAT.format(*v)
                    row_str = vrow_pat.format(vcoord=vcoord)
                    write_func(file_, row_str)

            # write faces
            for fi, f in enumerate(self.f):
                # file_.write(FACE_PAT.format(3, *f))
                write_func(file_, FACE_PAT.format(3, *f))

# def write_bin(self, filename=None):
# 	if filename is None:
# 		filename = self.filename

# 	# define ascii or binary writer
# 	write_func = _write_bin
# 	file_mode = 'wb'

# 	# write file
# 	with open(filename, file_mode) as file_:
# 		# writer header
# 		header_block = self._header_block()
# 		file_.write(header_block.encode('utf-8'))

# 		# write vertex coord, normal, and colour
# 		if (self.vnormals is not None) and (self.vcolours is not None):
# 			vrow_pat = b"%b%b%b\n"
# 			for vi, v in enumerate(self.v):
# 				vcoord = VCOORD_PAT_BIN%(v[0], v[1], v[2])
# 				vnormal = VNORMAL_PAT_BIN%(self.vnormals[vi][0], self.vnormals[vi][1], self.vnormals[vi][2])
# 				vcolour = VCOLOUR_PAT_BIN%(self.vcolours[vi][0], self.vcolours[vi][1], self.vcolours[vi][2])
# 				row_str = vrow_pat%(vcoord, vnormal, vcolour)
# 				file_.write(row_str)
# 		# write vertex coord, normal
# 		elif (self.vnormals is not None) and (self.vcolours is None):
# 			vrow_pat = b"%b%b\n"
# 			for vi, v in enumerate(self.v):
# 				vcoord = VCOORD_PAT_BIN%(v[0], v[1], v[2])
# 				vnormal = VNORMAL_PAT_BIN%(self.vnormals[vi][0], self.vnormals[vi][1], self.vnormals[vi][2])
# 				row_str = vrow_pat%(vcoord, vnormal)
# 				file_.write(row_str)
# 		# write vertex coord, colour
# 		elif (self.vnormals is None) and (self.vcolours is not None):
# 			vrow_pat = b"%b%b\n"
# 			for vi, v in enumerate(self.v):
# 				vcoord = VCOORD_PAT_BIN%(v[0], v[1], v[2])
# 				vcolour = VCOLOUR_PAT_BIN%(self.vcolours[vi][0], self.vcolours[vi][1], self.vcolours[vi][2])
# 				row_str = vrow_pat%(vcoord, vcolour)
# 				file_.write(row_str)
# 		# write vertex coord
# 		elif (self.vnormals is None) and (self.vcolours is None):
# 			vrow_pat = b"%b\n"
# 			for vi, v in enumerate(self.v):
# 				vcoord = VCOORD_PAT_BIN%(v[0], v[1], v[2])
# 				row_str = vrow_pat%(vcoord)
# 				file_.write(row_str)

# 		# write faces
# 		for fi, f in enumerate(self.f):
# 			file_.write(FACE_PAT_BIN%(3, f[0], f[1], f[2]))
