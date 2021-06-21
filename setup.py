from setuptools import find_packages, setup

setup(
    name="mizlab_tools",
    version="0.0.1",
    description="some codes",
    long_description="README.md",
    author="Omochice",
    author_email="h20ms419@hirosaki-u.ac.jp",
    url="",
    licence="LICENCE",
    install_requires=["biopython"],
    packages=find_packages(exclude=("tests", "docs")),
    tests_require=["pytest"],
)
