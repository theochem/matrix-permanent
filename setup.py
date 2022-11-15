#!/usr/bin/env python3

from os.path import dirname

from setuptools import Extension, setup

import numpy as np


name = "permanent"


version = "0.1.0"


licence = "GPLv3"


author = "AUTHOR"


author_email = "AUTHOR_EMAIL"


url = "URL"


description = "DESCRIPTION"


long_description = open("README.md", "r", encoding="utf-8").read()


classifiers = [
    "Environment :: Console",
    "Intended Audience :: Science/Research",
    "Programming Language :: Python :: 3",
    # "Topic :: Science/Engineering :: Molecular Science",
    # Figure out what to put here before publishing the package
]


install_requires = [
    "numpy>=1.13",
]


packages = [
    "permanent",
    "permanent.test",
]


ext_modules = [
    Extension(
        "permanent.permanent",
        ["permanent/permanent.c", "permanent/py_permanent.c"],
        include_dirs=[np.get_include()],
    ),
]


if __name__ == "__main__":

    setup(
        name=name,
        version=version,
        license=licence,
        author=author,
        author_email=author_email,
        url=url,
        description=description,
        long_description=long_description,
        classifiers=classifiers,
        install_requires=install_requires,
        packages=packages,
        ext_modules=ext_modules,
    )
