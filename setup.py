#!/usr/bin/env python3

from setuptools import setup

import numpy as np


name = "permanent"


version = "0.1.0"


licence = "GPLv3"


author = "QC-Devs"


author_email = "email@address.com"


url = "https://permanent.qcdevs.org/"


description = "Functions to compute the permanent of arbitrary matrices."


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


extras_require = {
    "test": ["pytest"],
    "docs": ["sphinx", "sphinx_rtd_theme"],
}


packages = [
    "permanent",
    "permanent.test",
]


package_data = {
    "permanent": ["permanent.so", "tuning.h"],
}


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
        extras_require=extras_require,
        packages=packages,
        package_data=package_data,
        include_package_data=True,
    )
