#!/usr/bin/env python
"""
Runs Fortran FClimDex on each text data file, putting output in separate directories
Michael Hirsch, Ph.D.

First you need to compile the fclimdex.f code (once) by:
    gfortran fclimdex.f -o FClimDex
"""
import subprocess
from pathlib import Path


