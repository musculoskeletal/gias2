"""
Script to update binary fieldwork files to json text format
"""

import argparse
import logging
import os
import sys
from glob import glob

from gias2.fieldwork.field import ensemble_field_function as eff
from gias2.fieldwork.field import geometric_field as gf
from gias2.fieldwork.field.topology import mesh

log = logging.getLogger(__name__)


# =============================================================================#
def is_json(inpath):
    with open(inpath, 'r') as f:
        head = f.read(1)
        if head == '{':
            return True
        else:
            return False


def _make_outpath(inpath, keepold):
    outpath = os.path.splitext(inpath)[0]
    if keepold:
        outpath = outpath + '_json'
    return outpath


def convert_geof(inpath, outpath, keepold, verbose=False):
    obj = gf.load_geometric_field(inpath, force=True)
    if outpath is None:
        outpath = _make_outpath(inpath, keepold)

    if verbose:
        log.debug('output filename: {}'.format(outpath))
    obj.save_geometric_field(outpath)


def convert_ens(inpath, outpath, keepold, verbose=False):
    obj = eff.load_ensemble(inpath, force=True)
    if outpath is None:
        outpath = _make_outpath(inpath, keepold)

    if verbose:
        log.debug('output filename: {}'.format(outpath))
    obj.save_ensemble(outpath)


def convert_mesh(inpath, outpath, keepold, verbose=False):
    obj = mesh.load_mesh(inpath)
    if outpath is None:
        outpath = _make_outpath(inpath, keepold)

    if verbose:
        log.debug('output filename: {}'.format(outpath))
    obj.save_mesh(outpath)


converters = {
    '.geof': convert_geof,
    '.ens': convert_ens,
    '.mesh': convert_mesh,
}


def convert_file(inpath, outpath=None, keepold=False, verbose=False):
    if verbose:
        log.debug('converting {}'.format(inpath))

    # check is not json already
    if is_json(inpath):
        if verbose:
            log.debug('Input file is already in json text format. Aborting.')
        return

    # identify file type
    ext = os.path.splitext(inpath)[1].lower()
    converter = converters.get(ext)
    if converter is None:
        raise ValueError('Unknown file extension {}'.format(ext))
        sys.exit(1)

    converter(inpath, outpath, keepold, verbose=verbose)


def convert_dir(inpath, keepold=False, verbose=False):
    # get all .geof files
    geof_files = glob(os.path.join(inpath, '*.geof'))
    if verbose:
        log.debug('Converting {} .geof files'.format(len(geof_files)))
    for gff in geof_files:
        convert_file(gff, keepold=keepold, verbose=verbose)

    # get all .ens files
    ens_files = glob(os.path.join(inpath, '*.ens'))
    if verbose:
        log.debug('Converting {} .ens files'.format(len(ens_files)))
    for enf in ens_files:
        convert_file(enf, keepold=keepold, verbose=verbose)

    # get all .ens files
    mesh_files = glob(os.path.join(inpath, '*.mesh'))
    if verbose:
        log.debug('Converting {} .mesh files'.format(len(mesh_files)))
    for mf in mesh_files:
        convert_file(mf, keepold=keepold, verbose=verbose)


# =============================================================================#
def make_parser():
    parser = argparse.ArgumentParser(
        description='Convert binary fieldwork files (.geof, .ens, .mesh) to json text format'
    )
    parser.add_argument(
        'inpath',
        help='Path to file or directory. If directory, all fieldwork files will be converted'
    )
    parser.add_argument(
        '-o', '--outfile',
        help='Output filename. Default is to overwrite original file. If inpath is a directory, this is ignored.'
    )

    parser.add_argument(
        '-k', '--keepold',
        action='store_true',
        help='Keep old binary file. New file will be named with a _json suffix'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Extra info.'
    )
    return parser


def main():
    parser = make_parser()
    args = parser.parse_args()

    # =============================================================================#
    # check path is valid
    if not os.path.exists(args.inpath):
        raise ValueError('Invalid input path.')

    # check output path is valid
    if args.outfile is not None:
        out_dir = os.path.split(args.outfile)
        if not os.path.isdir(out_dir):
            raise ValueError('Invalid output file directory.')

    # check if inpath is a file or a directory
    if os.path.isdir(args.inpath):
        convert_dir(args.inpath, keepold=args.keepold, verbose=args.verbose)
    else:
        convert_file(args.inpath, outpath=args.outfile, keepold=args.keepold, verbose=args.verbose)


if __name__ == '__main__':
    try:
        main()
        exit(0)
    except Exception as e:
        log.exception('Unexpected exception', e)
        exit(99)
