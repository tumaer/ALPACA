#!/usr/bin/env python3
import io
import os
import re

from setuptools import find_packages, setup

# Get version


def read(*names, **kwargs):
    with io.open(os.path.join(os.path.dirname(__file__), *names), encoding=kwargs.get("encoding", "utf8")) as fp:
        return fp.read()


def get_scripts():
    folder = os.path.join(os.getcwd(), "scripts")
    return [os.path.join(folder, name) for name in os.listdir(folder) if os.path.isfile(os.path.join(folder, name)) and name.endswith(".py")]


def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")


package_name = "alpacapy"
readme = open("README.md").read()
version = find_version(package_name, "__init__.py")

install_requires = [
    "h5py",
    "mpmath",
    "numpy",
    "pandas",
    "python-dateutil",
    "pytz",
    "scipy",
    "six",
    "sympy"
]

# Run the setup
setup(
    name=package_name,
    version=version,
    description="A module to use the Alpaca framework in python",
    long_description=readme,
    long_description_content_type="text/markdown",
    author="Chair of Aerodynamics and Fluid Mechanics",
    author_email="nanoshock@aer.mw.tum.de",
    url="git@gitlab.lrz.de:nanoshock/alpaca_aer.git",
    project_urls={
        "Git repository": "git@gitlab.lrz.de:nanoshock/alpaca_aer.git",
    },
    license="LGPLv3",
    classifiers=["Development Status :: 4 - Beta", "Programming Language :: Python :: 3"],
    packages=find_packages(),
    python_requires=">=3.6",
    install_requires=install_requires,
    extras_require={
    },
    scripts=get_scripts()
)
