#!/usr/bin/env python
from ez_setup import use_setuptools
use_setuptools()

from setuptools import setup, find_packages, Extension

from Cython.Distutils import build_ext
from Cython.Build import cythonize

from os import path
from glob import glob
import subprocess

#Cython extensions can be built in the standard fashion via distribute
#but files including numpy headers must have numpy include path specified.
import numpy

#Versioning logic
def get_git_version():
    """Get current version tag via call to git describe.

    Version tags are of the format:
        v<major>.<minor>.<patch>
    """

    version_string = subprocess.check_output(["git", "describe", "--match", "v*"])
    return version_string.strip().lstrip("v")

#Cython build logic
cython_subdirs = [
    "interface_fragment_matching/structure_database",
    "interface_fragment_matching/transform_types",
    "interface_fragment_matching/fragment_fitting/rmsd_calc",
    "interface_fragment_matching/fragment_fitting/clustering",
    "interface_fragment_matching/fragment_fitting/lookup",
    "interface_fragment_matching/fragment_profiling",
    "interface_fragment_matching/parallel",
    "interface_fragment_matching/residue_stubs/energy_evaluation",
]

cython_modules= cythonize([pyxfile for d in cython_subdirs for pyxfile in glob(path.join(d, "*.pyx"))])

for r in cython_modules:
    r.extra_compile_args=["-fopenmp"]
    r.extra_link_args=["-fopenmp"]

def pkgconfig(*packages, **kw):
    flag_map = {'-I': 'include_dirs', '-L': 'library_dirs', '-l': 'libraries'}
    for token in subprocess.check_output(["pkg-config", "--libs", "--cflags"] + list(packages)).strip().split(" "):
        kw.setdefault(flag_map.get(token[:2]), []).append(token[2:])
    return kw

include_dirs = [numpy.get_include(), path.abspath("."), path.abspath("interface_fragment_matching")] + pkgconfig("eigen3")["include_dirs"]

setup(
        name="interface_fragment_matching",
        description="PyRosetta modules for interface fragment matching.",
        author="Alex Ford",
        author_email="fordas@uw.edu",
        url="https://github.com/fordas/interface_fragment_matching/",
        version=get_git_version(),
        provides=["interface_fragment_matching"],
        packages=find_packages(where="interface_fragment_matching"),
        install_requires=[
            "numpy>=1.7",
            "numexpr>=2.0",
            "pandas>=0.12",
            "scipy>=0.13",
            "biopython",
            "tables>=3.0.0",
            "posix_ipc>=0.9.8",
            "fastcluster",
            "lockfile",
            "futures",
            "ipython>=1.0.0",
            "ipython_glmol",
            "pyrosetta>=3.4.20140216004611",
            "Jinja2",
            "SQLAlchemy>=0.9.4",
            "MySQL-python>=1.2.5",
            "bz2file",
            "toolz",
            "decorator"
        ],
        setup_requires=[
            "Cython>=0.18",
        ],
        tests_require=[
            "nose",
            "nose-html",
            "coverage",
            "nose-progressive"
        ],
        test_suite = "nose.collector",
        cmdclass = {'build_ext': build_ext},
        #Extra options to cythonize are *not* the same 'Extension'
        #See Cython.Compiler.Main.CompilationOptions for details
        #
        #In particular, include_dirs must be specified in setup
        #   as opposed to Extension if globs are used.
        include_dirs=include_dirs,
        ext_modules = cython_modules + [Extension(
                        'interface_fragment_matching/external/_transformations',
                        ['interface_fragment_matching/external/transformations.c'],
                       include_dirs=[numpy.get_include()])],
        entry_points = {
                "console_scripts" : [
                    "close_structure_chainbreaks=interface_fragment_matching.apps.close_structure_chainbreaks:main"
                    ],
                },
        )
