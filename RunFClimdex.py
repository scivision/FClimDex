#!/usr/bin/env python
"""
Runs Fortran FClimDex on each text data file, putting output in separate directories
Michael Hirsch, Ph.D.
"""
import shutil
import subprocess
from pathlib import Path

EXE='FClimdex'

def main(path:Path, pat:str):
    flist = finddata(path,pat)

    for f in flist:
        if not f.is_file():
            continue

        f = f.resolve()

        cdir = f.parent/f.stem
        if not f.suffix:
            cdir = cdir.with_suffix('.dir')

        sitefn = f.with_suffix('.site')
        if not sitefn.is_file(): # fallbcak to one-line para.txt
            sitefn = f.parent/'para.txt'


        cdir.mkdir(exist_ok=True) # create output directory per input file
        shutil.copy2(EXE,cdir)
#        try:
#            (cdir/EXE).resolve().symlink_to(f.parents[1]/EXE)
#        except FileExistsError:
#            pass

        shutil.copy2(f,cdir) # copy input file to output dir
        shutil.copy2(sitefn,cdir) # copy site def. file to output dir

        subprocess.check_call(['./'+EXE, sitefn.name, f.name],cwd=cdir)


def finddata(path:Path, pat:str) -> list:
    path = Path(path).expanduser()
# %% single file specified
    if path.is_file():
        return [path]
# %%
    if not path.is_dir():
        raise FileNotFoundError(f'{path} is not a directory')

    flist = path.glob(pat)
    flist = [f for f in flist if f.is_file()] # get rid of directories

    if not flist:
        raise FileNotFoundError(f'No {pat} files found in {path}')

    return flist


if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser()
    p.add_argument('datadir',help='directory where text files to process in Fortran are')
    p.add_argument('glob',help='pattern to glob files by',nargs='?',default='*.dat')
    p = p.parse_args()

    main(p.datadir, p.glob)
