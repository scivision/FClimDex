#!/usr/bin/env python
"""
Runs Fortran FClimDex on each text data file, putting output in separate directories
Michael Hirsch, Ph.D.

First you need to compile the fclimdex.f code (once) by:
    gfortran fclimdex.f -o FClimDex
"""
import shutil
import subprocess
from pathlib import Path

EXE='FClimDex'

def main(path:Path, pat:str):
    path = Path(path).expanduser()

    flist = sorted(path.glob(pat))

    if not flist:
        raise FileNotFoundError(f'No {pat} files found in {path}')

    for f in flist:
        f = f.resolve()

        cdir = f.parent/f.stem
        if not f.suffix:
            cdir = cdir.with_suffix('.dir')

        sitefn = f.with_suffix('.site')

        cdir.mkdir(exist_ok=True) # create output directory per input file
        try:
            (cdir/EXE).resolve().symlink_to(f.parents[1]/EXE)
        except FileExistsError:
            pass

        shutil.copy2(f,cdir) # copy input file to output dir
        shutil.copy2(sitefn,cdir) # copy site def. file to output dir

        subprocess.check_call(['./'+EXE, sitefn.name, f.name],cwd=cdir)



if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser()
    p.add_argument('datadir',help='directory where text files to process in Fortran are')
    p.add_argument('glob',help='pattern to glob files by',nargs='?',default='*.dat')
    p = p.parse_args()

    main(p.datadir, p.glob)
